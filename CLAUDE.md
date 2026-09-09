# specex: C++ -> Python/JAX PSF-fitting port

This repo is DESI's `specex` PSF-fitting code. Branch `python-gpu-port` (the
one you're almost certainly on) adds a from-scratch Python/JAX reimplementation
of the fitting pipeline (`py/specex/`) alongside the original C++ (`src/`),
targeting numerical parity with C++ plus a large GPU speedup. This file exists
so a new Claude Code session (or a new human collaborator) can get oriented
fast without re-deriving months of investigation from scratch.

## Read these, in this order

`how-to-run.md` and `env_setup.sh` live at the repo root; everything else
about the port's history, status, and known edge cases lives under
`docs/python-port/` -- moved there 2026-09-09 to keep the repo root down to
just this file, `how-to-run.md`, and `env_setup.sh`.

1. **`how-to-run.md`** (repo root) -- how to actually run the Python port
   (single bundle / full CCD / full night, GPU and CPU), environment setup,
   full CLI reference. Start here for "how do I invoke this."
2. **`docs/python-port/current-status.txt`** -- a dated bottom-line
   snapshot. Good for "where does this project currently stand," but always
   dated -- treat `porting-notes.md`'s most recent entries as more current
   than anything here if they conflict.
3. **`docs/python-port/edge-cases.md`** -- quick-reference index of known
   non-standard inputs/divergence modes (bad amps, missing preproc files,
   fiber-crossing regions, etc.), each with a one-paragraph summary and a
   pointer into `porting-notes.md` for the full story.
4. **`docs/python-port/porting-notes.md`** -- the real source of truth: an
   append-only, chronological session log of every investigation, bug, fix,
   and measurement across the whole port (6000+ lines). Don't read it front
   to back -- `grep` it by topic/filename/date. This is where the actual
   evidence for any claim in the other docs here lives.
5. **`docs/python-port/guide_fits_output.md`** -- the output PSF FITS file
   format (extensions, columns, what STATUS values mean).
6. **`docs/python-port/python-vs-cpp-diff.txt`** -- a standing reference
   diff/comparison summary (timing + correctness + the concrete algorithmic
   differences between the two backends), used by some of the `testing/`
   parity scripts and kept current as new campaigns land.
7. **`docs/python-port/algorithm-paper-map.md`** -- maps the algorithms in
   Julien Guy's specex paper (`2209.14482v2.pdf`, repo root) to where they're
   implemented in both `src/` (C++) and `py/specex/` (Python), including
   where/why the port's behavior deviates from a literal reading of either.

**Do not go looking for local Claude Code session transcripts as a history
source.** They're per-user, not portable across machines/accounts, and are
raw noisy tool-call logs rather than curated findings -- everything worth
keeping from every session has already been distilled into the docs above as
part of this project's own working discipline. If something feels missing
from the docs, that's a real gap worth flagging, not a sign to go hunting for
transcripts.

## Repo layout

- **`py/specex/`** -- the Python/JAX port (the actual subject of this branch).
  - `specex.py` -- CLI entry point (`python -m specex.specex`) and the
    process-pool driver (`fit_ccd_native`, `fit_bundle_task`, `main`). Also
    contains `run_specex()`, a wrapper around the compiled C++ pybind11
    extension (`_libspecex`). **Not legacy/one-off** -- confirmed
    2026-09-06 that this is the actual real-C++-production codepath:
    desispec's `desi_compute_psf` entry point (`desispec/scripts/specex.py`)
    does `from specex.specex import run_specex` and calls it directly, so
    every `--backend cpp`/`cpp-direct` run this project has ever done (via
    `desi_proc --mpi` or direct invocation) goes through this exact
    function. It's also still used by the old one-off comparison scripts
    (`testing/example_specex.py`, `testing/full_analysis.py`), but that's
    not its only or primary role. `read_preproc_cpp()`/`read_preproc()`
    (`io.py`) -- the C++ and Python image-loading paths this feeds --
    apply *identical* masking: `ivar[mask != 0] = 0.0` before either
    fitter ever runs, so cosmic-ray/bad-pixel exclusion (preproc MASK bit
    4 etc.) is handled the same way in both backends by construction (see
    docs/python-port/porting-notes.md 2026-09-06).
  - `fitter.py` -- the actual fit engine: spot selection and the staged
    bundle fit (`PSF_Fitter.fit()`). Despite the name, this is **not** a
    single joint solve: `fit()` advances through fixed stages -- `'flux'`
    (flux + continuum only) -> `'trace'` (adds trace-correction
    coefficients) -> `'sigma'` (adds GHSIGX/GHSIGY) -> `'full'` (every
    remaining GH-shape coefficient) -- with trace and GHSIGX/GHSIGY
    **frozen, never revisited**, once their own stage ends. This
    structurally mirrors real C++ (`specex_psf_fitter.cc`'s `FitEverything`
    never solves trace and PSF shape together -- confirmed the one place a
    combined `fit_trace=true; fit_psf=true` call exists in the C++ source
    is permanently commented-out dead code) and was merged into this
    branch from the since-retired `experiment/cpp-alternating-solve` branch
    (`docs/python-port/porting-notes.md`, 2026-07-29 through 2026-08-05) -- there is no
    separate branch to check out for this anymore, it's simply how `fit()`
    behaves by default. Production defaults as of that merge:
    `trace_per_fiber_deg=6` (each fiber gets its own independent
    7-coefficient trace basis, not a basis shared across the bundle) and
    `trace_prior_deg=1` (a soft prior pulling a *dead-column-flagged*
    fiber's degree>=1 trace coefficients toward the bundle's cross-fiber
    mean -- C++ has the identical mechanism coded but its own CLI default
    leaves it off in real production, `src/specex_pyoptions.h`; see
    `docs/python-port/porting-notes.md`'s 2026-08-25/26 fiber-0-investigation entries for
    what this asymmetry does and doesn't explain). The default
    `--line-search grid` path is what's actually used in production; the
    alternate `'brent'`/`'cpp'` line-search modes and their helper
    `_cpp_brent()` are experimental, kept only for reference (tested
    correctness-neutral, never the default). Several `SPECEX_*` env vars
    gate genuinely optional one-off diagnostic code paths (e.g.
    `SPECEX_DEBUG_MEM`, `SPECEX_FREEZE_GH10`, `SPECEX_MATCH_CPP_DEAD_COLUMN`)
    -- distinct from `SPECEX_MIXED_PRECISION`, `SPECEX_TRACE_PRIOR_WEIGHT`/
    `_NDEAD_THRESHOLD`, which are just the internal plumbing for real,
    on-by-default CLI flags. `_fit_one_spot_jax` is dead code (superseded
    by the batched `_fit_all_spots_batch`/`_get_spot_stats_jax`), left over
    from an earlier version.
  - `psf.py` -- the Gauss-Hermite PSF model (`GaussHermitePSF`) and the
    per-fiber `PSF`/`PSF_Params` containers. `single_pix_value_np` is a
    dead reference implementation next to the actually-used
    `single_pix_value_jnp`.
  - `io.py` -- FITS I/O: reading preproc/input-PSF files, writing the
    merged output PSF (`write_python_psf`, including its own inline QA
    pass -- see below).
  - `math.py` -- Legendre/Hermite polynomial bases shared by the trace and
    PSF-shape fits.
  - `qa.py` -- **fully dead code.** Adapted from upstream specex#91;
    `specex.py` imports it but the one call site is commented out.
    `io.py`'s `write_python_psf` has its own independently-adapted inline
    reimplementation of the same trace-crossing QA logic, which is what
    actually runs.
  - `fitter_old.py`, `fitter_old_gpt.py` -- abandoned early drafts, not
    imported anywhere. Kept on disk (not deleted) but untracked from git;
    ignore them.
- **`src/`** -- the original C++ implementation (plus vendored `pybind11`).
  Still built and used as the correctness baseline for comparisons.
- **`testing/`** -- validation/comparison tooling, not unit tests in the
  pytest sense (except `test_math.py`, `test_math_psf.py`,
  `test_vectorization.py`, which are real pytest suites CI actually runs).
  `how-to-run.md` documents which scripts are the current, maintained
  entry points (`run_night.py`, `select_test_case.py`,
  `instrumentation_analysis.py`, `validate_all_modes.py`,
  `full_ccd_campaign.py`, `compare_correctness.py`, `per_fiber_breakdown.py`,
  `stage_preproc.py`) -- treat any `testing/*.py` not mentioned there as a
  one-off debugging script, not a supported tool.
- **`2209.14482v2.pdf`** -- Julien Guy's original specex paper (arXiv). The
  algorithms it describes are what both the C++ code and this Python port
  implement.

## Who's involved

Stephen and Julien Guy (specex's original C++ author) are reviewing this
branch. Julien's two standing asks for the port going forward are (1) a guide
to reading the code so he can follow the port's structure, and (2) a very
complete test suite -- both open, ongoing work, not one-shot deliverables.
