# specex: C++ -> Python/JAX PSF-fitting port

This repo is DESI's `specex` PSF-fitting code. Branch `python-gpu-port` (the
one you're almost certainly on) adds a from-scratch Python/JAX reimplementation
of the fitting pipeline (`py/specex/`) alongside the original C++ (`src/`),
targeting numerical parity with C++ plus a large GPU speedup. This file exists
so a new Claude Code session (or a new human collaborator) can get oriented
fast without re-deriving months of investigation from scratch.

## Read these, in this order

1. **`how-to-run.md`** -- how to actually run the Python port (single bundle /
   full CCD / full night, GPU and CPU), environment setup, full CLI reference.
   Start here for "how do I invoke this."
2. **`current-status.txt`** -- a dated bottom-line snapshot (written
   2026-07-21, with a later pointer at the top to its own section 14/15 as the
   most current finding). Good for "where does this project currently stand,"
   but always dated -- treat `porting-notes.md`'s most recent entries as more
   current than anything here if they conflict.
3. **`edge-cases.md`** -- quick-reference index of known non-standard
   inputs/divergence modes (bad amps, missing preproc files, fiber-crossing
   regions, etc.), each with a one-paragraph summary and a pointer into
   `porting-notes.md` for the full story.
4. **`porting-notes.md`** -- the real source of truth: an append-only,
   chronological session log of every investigation, bug, fix, and
   measurement across the whole port (3500+ lines). Don't read it front to
   back -- `grep` it by topic/filename/date. This is where the actual
   evidence for any claim in the three docs above lives.
5. **`guide_fits_output.md`** -- the output PSF FITS file format (extensions,
   columns, what STATUS values mean).
6. **`python-vs-cpp-diff.txt`** -- a standing reference diff/comparison
   artifact used by some of the `testing/` parity scripts.

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
    contains `run_specex()`, a **legacy wrapper around the compiled C++
    pybind11 extension** (`_libspecex`) used only for old one-off
    comparison scripts (`testing/example_specex.py`,
    `testing/full_analysis.py`) -- not part of the production pipeline.
  - `fitter.py` -- the actual fit engine: spot selection, the joint
    bundle fit (`PSF_Fitter.fit()`), trace-prior/dead-column/masked-amp
    handling. The default `--line-search grid` path is what's actually
    used in production; the alternate `'brent'`/`'cpp'` line-search modes
    and their helper `_cpp_brent()` are experimental, kept only for
    reference (tested correctness-neutral, never the default). Several
    `SPECEX_*` env vars gate genuinely optional one-off diagnostic code
    paths (e.g. `SPECEX_DEBUG_MEM`, `SPECEX_FREEZE_GH10`,
    `SPECEX_MATCH_CPP_DEAD_COLUMN`) -- distinct from `SPECEX_MIXED_PRECISION`,
    `SPECEX_TRACE_PRIOR_WEIGHT`/`_NDEAD_THRESHOLD`, which are just the
    internal plumbing for real, on-by-default CLI flags. `_fit_one_spot_jax`
    is dead code (superseded by the batched `_fit_all_spots_batch`/
    `_get_spot_stats_jax`), left over from an earlier version.
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
