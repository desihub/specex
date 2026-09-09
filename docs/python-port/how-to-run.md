# Specex Python/GPU Usage Guide

This guide documents how to run the ported Python/JAX version of Specex and the various validation scripts used to verify numerical parity with the C++ baseline.

## Quick Start

One-time environment setup, then the three run scopes most people need. Full documentation -- every flag, timing table, and operational gotcha -- is below this section.

```bash
# Once per session (sets PYTHONPATH/LD_LIBRARY_PATH, loads cudatoolkit)
source env_setup.sh
```

### Full night/expid (production scale, all 30 cameras)

```bash
python testing/run_night.py --night 20260401 --expid 00344649 \
    --backend python --worker-mode persistent
```

`--worker-mode persistent` is the current fastest, fully-validated config: **~483-503s/night** (1 node/4 GPUs, warm cache) vs. C++'s **~592-648s** -- see Section 2.3 for the full breakdown and why `--worker-mode` still defaults to the older, slower `subprocess` mode instead of this one. Swap `--backend cpp` for the real production C++/MPI driver, or `--backend cpp-direct` for a `desi_proc`-free C++ baseline (Section 2.3 explains the difference between the three backends).

### One camera (full CCD, all 20 bundles)

```bash
python -m specex.specex \
    -a /path/to/preproc-z8-00344649.fits.gz \
    --in-psf /path/to/shifted-input-psf-z8-00344649.fits \
    --out-psf $SCRATCH/pyfit-psf-z8-00344649.fits \
    --broken-fibers 473,474 \
    --gpu 4 --workers-per-gpu 5
```

~33-45s warm (Section 2.2). Drop `--broken-fibers` if the camera has none; find the right value with `testing/select_test_case.py` (Section 3).

### One bundle (25 fibers, fast iteration/debugging)

```bash
python -m specex.specex \
    -a /path/to/preproc-z8-00344649.fits.gz \
    --in-psf /path/to/shifted-input-psf-z8-00344649.fits \
    --out-psf $SCRATCH/pyfit-psf-z8-00344649_05.fits \
    --first-bundle 5 --last-bundle 5 \
    --first-fiber 125 --last-fiber 149 \
    --broken-fibers 473,474 \
    --gpu 1
```

A few seconds warm (Section 2.1). `--first-fiber`/`--last-fiber` is optional (25 fibers/bundle = `bundle_id*25` to `bundle_id*25+24`) but keeps the run scoped to exactly the fibers you're inspecting.

No GPU available? Every command above also runs with `--backend cpu` (Section 2.4) -- correctness-equivalent, ~11x slower, no `run_night.py` support yet so loop it per-camera for a full night.

---

## 0. Environment Creation (One-time setup)

To recreate the environment used for development (`specex_env`):

```bash
# Create venv -- use the desiconda Python below, or any Python 3.11+
python -m venv /path/to/your/specex_env
source /path/to/your/specex_env/bin/activate

# Install core dependencies
pip install --upgrade pip
pip install numpy fitsio astropy scipy

# Install JAX with CUDA support (for Perlmutter A100s)
pip install --upgrade "jax[cuda13]"

# Note: no "-f ..." index URL and no "_pip" suffix on the extra -- both are
# obsolete. The jax-cuda13-plugin/jax-cuda13-pjrt wheels and all needed
# nvidia-* CUDA runtime libs now publish straight to PyPI via the plain
# "cuda13" extra. "cuda13_pip" is NOT a valid extra as of jax 0.10.x/0.11.x
# -- pip only warns (does not error) and silently falls back to installing
# a CPU-only jaxlib, which then fails at runtime with "Unknown backend:
# 'gpu' requested... Platforms present are: cpu" (as of 2026-08-18, that
# specific crash is now a clear actionable error instead -- see the box at
# the end of this section). Watch for the pip warning if you ever see it
# again after a jax upgrade.
```

### Known-good versions (this environment, confirmed working 2026-08-18)

For reference/reproducibility -- these are the exact versions `specex_env` currently uses on a Perlmutter GPU node. Nothing here is strictly pinned by the code; if you land on nearby versions of any of these, that's fine.

| Component | Version |
|---|---|
| Base Python (via `python -m venv`) | 3.13.12, from `desiconda` (see below) |
| DESI environment | `source $CFS/desi/software/desi_environment.sh 26.3` -- optional, see note below |
| jax / jaxlib | 0.10.1 |
| jax-cuda13-plugin / jax-cuda13-pjrt | 0.10.1 |
| numpy | 2.3.5 |
| fitsio | 1.3.0 |
| astropy | 7.2.0 |
| scipy | 1.16.3 |
| GPU | NVIDIA A100-PCIE-40GB |
| NVIDIA driver | 580.159.04 (reports CUDA 13.0) |
| `module load cudatoolkit` (loaded by `env_setup.sh`) | 12.9 |

Notes:
*   **`desi_environment.sh` is optional, not a prerequisite.** The core fit driver (`python -m specex.specex`, everything in Sections 2-4 below) never imports `desiutil`/`desispec` and doesn't need it. Only `py/specex/qa.py` (a separate, optional QA-plotting helper, not part of any run mode documented here) imports `desispec.io.xytraceset` and `desiutil.log` -- source `desi_environment.sh` (or otherwise have those two packages on `PYTHONPATH`) only if you plan to use `qa.py`.
*   **Base Python matters less than you'd think.** `specex_env` was built via `python -m venv` from `desiconda`'s Python 3.13.12 (`/global/common/software/desi/perlmutter/desiconda/.../conda/bin/python3.13`) with `include-system-site-packages = false` -- a clean, non-inheriting venv. Any reasonably modern Python 3 (3.11+) that can install the pip packages above should work equally well; there's nothing desiconda-specific baked into the pip-installed dependency set itself.
*   **Driver vs. loaded module CUDA version mismatch is expected and harmless.** The node's NVIDIA driver reports CUDA 13.0 (via `nvidia-smi`), while `env_setup.sh`'s `module load cudatoolkit` loads 12.9 -- these don't need to match. JAX's CUDA support comes entirely from the pip-installed `jax-cuda13-*`/`nvidia-*` wheels (self-contained, see `env_setup.sh`'s `LD_LIBRARY_PATH` derivation below); the driver only needs to be new enough to run CUDA 13 code (it is), and the loaded `cudatoolkit` module isn't actually load-bearing for the Python/JAX path at all.

**If `--gpu`/`--backend gpu` (the default) is requested but no CUDA-enabled jaxlib is installed**, `python -m specex.specex` now fails immediately with a clear, actionable `RuntimeError` pointing back at this section, rather than either JAX's own opaque `Unknown backend: 'gpu' requested... Platforms present are: cpu` or (worse) silently falling back to CPU and running ~10x slower with no indication anything is wrong. This is deliberate fail-fast behavior, not a bug: since `--gpu` is the default and this is a performance-critical batch pipeline, a silent CPU fallback would be a much nastier trap than a loud failure. Pass `--backend cpu` explicitly if you ever want to run on CPU on purpose.

### Alternative to `specex_env`: a shared DESI environment with JAX (tested 2026-09-13)

Stephen is assembling a shared DESI module environment that bundles JAX,
as a path toward not needing a personal venv at all for production. Tested
`source /global/cfs/cdirs/desi/software/desi_environment.sh test-26.9`
(note: must `deactivate` any active personal venv and `unset PYTHONPATH`
*before* sourcing it, or the venv's own `python`/`PYTHONPATH` silently wins
and you're not actually testing the new environment):

*   **The GPU-native Python/JAX path (`main()`/`fit_ccd_native()`, everything
    in Sections 2-4 below) works with zero code changes.** This environment
    provides its own Python 3.14.7 (desiconda 20260908-3.0.0), jax/jaxlib
    0.10.2 (vs. `specex_env`'s 0.10.1), numpy 2.5.3, scipy 1.18.0, astropy
    8.0.1, fitsio 1.4.2 -- all close to or newer than `specex_env`'s
    versions -- and `jax.devices()` correctly reports all 4 A100s with no
    extra setup (no `env_setup.sh`-style `LD_LIBRARY_PATH`/`NVLIBS` dance
    needed; this environment's own CUDA plumbing already works). A real
    single-bundle fit (`z8/00344649` bundle 5) ran clean end-to-end,
    `SPECEX_RESULT: OK 1/1 bundles`. **You still need to prepend this
    repo's own `py/` to `PYTHONPATH`**
    (`export PYTHONPATH=/path/to/specex/py:$PYTHONPATH`,
    same idea as `env_setup.sh`, just without the venv-specific
    `NVLIBS`/`LD_LIBRARY_PATH` piece, which isn't needed here) --
    **this environment already bundles its own separate `specex` install**
    (`.../desiconda/.../code/specex/main/py`, the old pre-port C++-only
    version -- no `fitter.py`/`psf.py`/`math.py` at all), which silently
    shadows this branch's code if you don't put this repo's `py/` first.
*   **The C++-wrapper path (`run_specex()`, `--backend cpp`/`cpp-direct`)
    does NOT work as-is.** The compiled pybind11 extension
    (`py/specex/_libspecex.cpython-313-x86_64-linux-gnu.so`) is built
    against `specex_env`'s CPython 3.13 ABI; this environment's Python 3.14
    can't load it (`ModuleNotFoundError: No module named
    'specex._libspecex'` the moment `run_specex()` actually tries the
    lazy `from ._libspecex import ...`, even though `from specex.specex
    import run_specex` itself succeeds -- the import is lazy, inside the
    function body). **Needs a rebuild against this environment's Python/
    toolchain before `--backend cpp` can run under it** -- `cmake`
    (4.4.3) and `g++` (via `PrgEnv-gnu/8.7.0`) are both present, so a
    rebuild looks straightforward, just not attempted yet (not needed for
    the GPU-native path, which is the actual subject of this port).
*   **Bottom line**: this environment is a viable `specex_env` replacement
    for the Python/JAX path today, once `PYTHONPATH` is set correctly. The
    C++ path is the one piece of "readying an environment for production"
    work still open, and it's independent of anything about this branch's
    own Python code.

---

## 1. Environment Setup (Every session)

Before running any scripts, ensure your environment is set up correctly on Perlmutter.

```bash
# From the project root
source env_setup.sh
```

This sets `PYTHONPATH` to include the local `py` directory, points `LD_LIBRARY_PATH` at the pip-installed NVIDIA CUDA libraries JAX needs, and loads the `cudatoolkit` module if on a compute node.

---

## 2. Running the Python CCD Fit -- Three Modes

All modes go through the same entry point, `python -m specex.specex` (or the `fit_ccd_native()` function directly). What changes is scope: a single bundle, one camera's full CCD, or all 30 cameras of an exposure.

**Before you start timing anything, read this:** the *first* `python -m specex.specex` invocation in a fresh environment (or after a `~/.cache/specex/jax_compilation_cache` wipe) pays JAX's JIT-compilation cost on top of the real fit -- expect it to be several times slower than every run after it. `docs/python-port/porting-notes.md`'s persistent compilation cache means that cost is paid once, not once per run: once a given (function, array-shape) pair has been compiled anywhere on this filesystem, every later run reuses it. So "cold" vs "warm" below isn't about caching your specific inputs, it's about whether *any* prior run has already compiled the shapes this run needs. Don't judge real throughput off a first run.

### 2.1 Mode 1: Single Bundle (25 fibers)

Useful for fast iteration/debugging. Restrict fitting to one bundle with `--first-bundle`/`--last-bundle` (0-19) and, optionally, a matching `--first-fiber`/`--last-fiber` range (25 fibers/bundle, `bundle_id * 25` to `bundle_id * 25 + 24`). A single bundle doesn't benefit from more than one GPU/worker, so `--gpu 1` (no `--workers-per-gpu` needed) is the right call here:

```bash
python -m specex.specex \
    -a /dvs_ro/cfs/cdirs/desi/spectro/redux/matterhorn/preproc/20260401/00344649/preproc-z8-00344649.fits.gz \
    --in-psf /dvs_ro/cfs/cdirs/desi/spectro/redux/matterhorn/exposures/20260401/00344649/shifted-input-psf-z8-00344649.fits \
    --out-psf $SCRATCH/pyfit-psf-z8-00344649_05.fits \
    --first-bundle 5 --last-bundle 5 \
    --first-fiber 125 --last-fiber 149 \
    --broken-fibers 473,474 \
    --gpu 1
```

### 2.2 Mode 2: Full CCD (one camera, all 20 bundles)

Drop `--first-bundle`/`--last-bundle` (they default to the full 0-19 range) to fit an entire camera. Unlike Mode 1, this scope *does* benefit from spreading bundles across every GPU on the node -- on a standard 4-GPU Perlmutter node, `--gpu 4` (the CLI default) is the right starting point, not `--gpu 1`:

```bash
python -m specex.specex \
    -a /dvs_ro/cfs/cdirs/desi/spectro/redux/matterhorn/preproc/20260401/00344649/preproc-z8-00344649.fits.gz \
    --in-psf /dvs_ro/cfs/cdirs/desi/spectro/redux/matterhorn/exposures/20260401/00344649/shifted-input-psf-z8-00344649.fits \
    --out-psf $SCRATCH/pyfit-psf-z8-00344649.fits \
    --broken-fibers 473,474 \
    --gpu 4 --workers-per-gpu 5
```

`--gpu N` spreads the 20 bundles across `N` GPUs on the *current node*; `--workers-per-gpu` controls how many bundle-fit worker processes are packed onto each GPU concurrently. `5` is a good starting point for one-off/ad-hoc single-camera runs on any band (`4 GPUs x 5 workers = 20` -- exactly one wave for a 20-bundle CCD) -- if omitted it's auto-detected per band instead, defaulting lower (3) for z-band specifically once `--trace-per-fiber-deg` (on by default) is active, since the z-band per-fiber design matrix is the most memory-hungry combination and the auto-default errs conservative. See the callout below for when to use the lower, more conservative values instead.

**Cold vs. warm, and what "optimal" actually means here** -- measured on this exact command (z8/00344649, one node/4 A100s):
| Run | Flags | Time |
|---|---|---|
| 1st ever run (cold JIT cache) | `--gpu 1 --workers-per-gpu 4` | 270s |
| 2nd run (warm cache) | `--gpu 4` (workers-per-gpu auto) | 42s |
| 3rd run (warm cache) | `--gpu 4 --workers-per-gpu 5` | **33s** |

The 270s->42s jump is dominated by the JIT cache going from cold to warm, not the GPU-count change alone (don't read that as "`--gpu 1` is 6x slower than `--gpu 4`" -- a cold `--gpu 4` run would still be slow). The 42s->33s jump is the real, repeatable effect of `--workers-per-gpu 5` over the lower auto-default for this camera/band. **33s with a warm cache is the best full-CCD single-camera number seen so far and a reasonable target to expect.**

> **`--workers-per-gpu 5` vs. the lower production table in 2.3 below -- these are not in conflict.** `5` (this section) is safe and fast for testing *one specific camera at a time*, where you can watch the log for `WARNING: Bundle`/`RESOURCE_EXHAUSTED` and just re-run at a lower value if you ever see one. The lower per-band values in Section 2.3 (`10 b / 7 r / 4 z`, and the code's own even more conservative z auto-default of `3`) exist because an *unattended* full-night batch has no one watching logs camera-by-camera, and z-band's OOM risk at `workers-per-gpu 5` is data-dependent -- confirmed to hit real, silent `RESOURCE_EXHAUSTED` failures on roughly 4 of 10 z-band cameras tested (higher-than-typical spot density combined with the per-fiber trace's larger design matrix), while the other 6/10 ran clean at `5` just like this section's example. So: for a single ad-hoc camera, start at `5` and drop to `4` only if you see an OOM warning; for a full unattended night across cameras/exposures you haven't pre-screened, use Section 2.3's validated-safe table instead.

For a batch of full-CCD runs (many cameras, one at a time, each compared against C++), see `testing/full_ccd_campaign.py`:
```bash
python testing/full_ccd_campaign.py --night 20260401 --expid 00344649 \
    --cameras b5,b4,b2,r3,r5,r1,z1,z6,z9 --outdir $SCRATCH/specex/full_ccd
```
This launches the real C++ wrapper (`srun -n 20 desi_compute_psf --mpi`) and the Python port concurrently per camera (they don't contend for the same resources -- C++ is CPU/MPI, Python is GPU), and reports wall time plus X/Y trace RMS and wavelength-residual RMS for both.

### 2.3 Mode 3: Full Night/Expid (all 30 cameras, production scale)

**`testing/run_night.py` is the single entry point for this mode**, switchable between the C++ and Python/JAX backends via one flag (or the `SPECEX_BACKEND` env var), so the same command works for either pipeline:

```bash
python testing/run_night.py --night 20260401 --expid 00344649 --backend python
python testing/run_night.py --night 20260401 --expid 00344649 --backend cpp

# or set it once for the session:
export SPECEX_BACKEND=python
python testing/run_night.py --night 20260401 --expid 00344649
```

**Scope:** this fits an *already-preprocessed* exposure (`preproc-*.fits.gz` + `shifted-input-psf-*.fits` must already exist -- true for any real matterhorn production night/expid). It is not a `desi_proc` replacement:

*   `--backend cpp` runs the real production driver, `desi_proc --mpi`, which does its own preprocessing (idempotent -- skips it if outputs already exist) then calls the C++ `desi_psf_fit` binary per camera via MPI ranks. One `srun` call; scales via `--nodes` using the validated rank formula (`100*nodes + 1` -- `-N1`/`-n101` measured ~11 min, `-N3`/`-n301` measured ~7 min, see `docs/python-port/porting-notes.md`). Output goes to a **private** `$DESI_SPECTRO_REDUX/$SPECPROD` tree (`--redux-dir`/`--specprod`, default under `--outdir`), never the real production `matterhorn` tree.
*   `--backend python` runs `python -m specex.specex` once per camera, each **pinned to a dedicated GPU** via `CUDA_VISIBLE_DEVICES` (no GPU sharing across cameras), using a **per-node dynamic work queue** so however many GPUs you have stay busy. Auto-detects node count (from the SLURM allocation) and GPUs/node (`nvidia-smi -L`); one node runs locally, multiple nodes launch via `srun -N1 -n1 -w <hostname>` per node (same pattern validated this session). Does **not** call `desi_proc` at all -- it goes straight to `specex.specex` on the existing preprocessed files.
*   `--backend cpp-direct` runs the real C++ `desi_compute_psf --mpi` binary once per camera (`--cpp-ranks`, default 20 -- same invocation `testing/full_ccd_campaign.py` already validates single-camera), reading the exact same preprocessed inputs `--backend python` does and writing the same `fit-psf-<cam>-<expid>.{fits,log}` naming, so the two are directly comparable file-for-file. Like `--backend python`, does **not** call `desi_proc` -- no idempotent-preprocessing pass, and no single 101+-rank MPI collective for one bad camera to hang (see 7.2 below). CPU-only; sequential by default (`--cpp-concurrency 1`) for clean per-camera timing. Added 2026-08-11 specifically to get a `desi_proc`-free C++ timing/correctness baseline without its MPI-hang or missing-calib-state failure modes.

**`--worker-mode {subprocess, persistent}`** (default: `subprocess`) controls how `--backend python` executes across cameras:
*   **`subprocess`** (default): a fresh `python -m specex.specex` process per camera -- the original, most-tested methodology, and still the default for backward compatibility with every prior campaign's numbers.
*   **`persistent`**: one long-lived worker process per GPU, pulling cameras off a shared queue and calling `fit_ccd_native()` in-process instead of paying interpreter/JAX-import startup cost per camera. Also reuses the per-bundle `multiprocessing.Pool` across a GPU worker's whole stream of cameras (`--no-pool-reuse` to disable, A/B-testing only, no correctness effect). **This is the fastest validated configuration and the recommended choice for any real timing-sensitive run** (see the campaign numbers below) -- it isn't the CLI default only because it's newer and less battle-tested than `subprocess`, not because of any known downside.

**Validated `--workers-per-gpu-{b,r,z}` settings per band** (concurrent bundle-fit workers packed onto one GPU for one camera -- override if needed):
| Band | `subprocess` default | `persistent` default | Why |
|------|---------|---------|-----|
| b    | 10      | 12      | Validated sweet spot, no OOM |
| r    | 7       | 8       | Naive value silently OOMs a handful of bundles (`RESOURCE_EXHAUSTED`) on some r-band cameras |
| z    | 4       | 5       | Larger per-fiber design matrix (per-fiber trace default) OOMs on higher values even in isolation |

`persistent` mode's higher per-band values are safe specifically because pool reuse frees up extra GPU memory headroom that per-camera Pool teardown/recreation didn't leave available; don't reuse the `persistent` column's values with `--worker-mode subprocess`.

**Full 10-night campaign, `subprocess` -> `persistent` (30 cameras/night, 1 node/4 GPUs, same 10 nights, correctness held constant to 4 decimal places throughout)**:
| Configuration | 10-night mean wall time | vs. C++ (592.3s) |
|---|---|---|
| C++ (production MPI, single node) | 592.3s | -- |
| Python, `subprocess` (pre-persistent) | 673.8s | +13.8% (slower) |
| Python, `persistent`, naive camera-to-GPU split | 579.5s | -2.2% |
| + bundle-pool reuse | 514.0s | -13.2% |
| **+ `workers-per-gpu` 12/8/5 (current `persistent`-mode default)** | **483.3s** | **-18.4%** |

A second, independent 10-night set (2026-09-05, bringing the cumulative validated sample to 20 nights) reproduced this: Python `persistent` mean **502.5s warm / 591.2s cold** (first-touch-on-a-fresh-node cost, ~15% higher, consistent across all 10 nights) vs. C++ mean **648.3s**. Correctness across the full 20-night sample: **xrms=0.0120px, yrms≈0.0129px** vs. C++ (see `docs/python-port/porting-notes.md`'s 2026-09-02/05 entries for full per-night/per-camera tables). **Before trusting any timing number on this project**, confirm `$HOME` isn't at its NERSC disk quota -- a full quota produces silent per-bundle failures at `rc=0` and inflated wall time, not an obvious error (see `docs/python-port/porting-notes.md`, 2026-09-04).

**Multi-node camera splitting:** the default is a naive alternating split (band-diverse but not load-balanced -- there's no timing prior for an arbitrary fresh night/expid). Pass `--lpt-profile cameras.json` (a `{"b0": 69.2, ...}` map of measured per-camera wall times, e.g. parsed from a prior run's own logs) to get an **LPT (longest-processing-time-first) balanced split** instead -- sorts cameras descending by known duration and greedily assigns each to whichever GPU-slot currently has the least total load, closing most of the gap a naive split leaves on the table (the slowest, most variable band, z, otherwise gets queued last with nothing to fill the tail).

**Output naming and location:** `--backend python` writes `fit-psf-<cam>-<expid>.fits` per camera (matches `desi_proc`'s/`desi_compute_psf`'s own real output naming -- was bare `<cam>.fits` before 2026-08-11, which didn't line up with anything C++ produces), plus a same-stem `fit-psf-<cam>-<expid>.log` per camera (was `py-<cam>.log` -- no `expid`, so it collided across different expids of the same camera sharing one outdir, and didn't visually pair with its own fits file). Rerunning does **not** skip existing output -- there's no idempotency check, files are simply overwritten in place; no need to clear the output directory between runs. Where those files land is resolved in three tiers:
1. `--outdir <dir>`, if given -- used exactly as-is (a flat directory, for ad-hoc/scratch runs).
2. `$DESI_SPECTRO_REDUX` (+ `$SPECPROD`, default `$USER`), if set in the environment and `--outdir` is not given -- resolves to `$DESI_SPECTRO_REDUX/$SPECPROD/exposures/<night>/<expid>/`, mirroring `desi_proc`'s/`--backend cpp`'s own layout, so a bare `export DESI_SPECTRO_REDUX=...` before running either backend now puts both in the same, directly comparable location.
3. Otherwise, `$SCRATCH/specex/run_night_<night>_<expid>/` (unchanged default).

Refuses to resolve inside the real production `matterhorn` tree under any of the three tiers -- fails fast rather than writing fit-psf output there. Note `--backend cpp` does **not** participate in tier 2 -- it always manages its own private redux tree under `--outdir`/`--redux-dir`, deliberately ignoring any pre-set `$DESI_SPECTRO_REDUX` so a real production value left in the shell can never get written to.

**Other flags:** `--cameras` (restrict to a subset, default all 30), `--nodes` (default: full SLURM allocation), `--gpus-per-node` (default: auto-detect), `--worker-mode {subprocess,persistent}` (default `subprocess`, see above), `--no-pool-reuse` (persistent mode only, A/B-testing), `--footprint-margin` (passthrough to `specex.specex --footprint-margin`, default 7), `--dry-run` (print planned commands without executing).

**The multi-node/LPT numbers below predate `--worker-mode persistent` and were measured under `subprocess` mode only** -- they haven't been separately re-validated with `persistent` mode's per-camera-startup savings stacked on top, so treat "N nodes + LPT" as still using `subprocess` mode until that combination is tested. For single-node runs, `persistent` mode (table above) is faster than any of the multi-node `subprocess` numbers below and doesn't need multiple nodes at all.

Measured with this exact tool's predecessor scripts (2026-08-10, before consolidation into `run_night.py`, `subprocess`-equivalent methodology): **13.2 min / 30 cameras** on 1 node/4 GPUs (0 bundle failures); **6.7 min** on 2 nodes/8 GPUs with a naive split, **5.92 min** LPT-rebalanced -- both correctness-verified (xrms/yrms match a rebuilt/official C++ reference to <0.03px mean), and the result generalizes across independent nights/exposures (confirmed on a second night; absolute timing varies by exposure since real exposures differ in total compute needed, but the technique transfers). Always grep worker logs for `WARNING: Bundle` to catch silent per-bundle GPU-OOM failures -- a nonzero process exit code is *not* a reliable failure signal, `fit_ccd_native` logs a warning and keeps merging on a per-bundle failure. If `workers-per-gpu`, the camera set, or the node count change, an `--lpt-profile` needs to be recomputed from a fresh timing profile -- it isn't portable across settings changes, but the profile-then-rebalance *technique* is.

Re-measured through the consolidated `run_night.py` itself on 1 node/4 GPUs across 4 independent nights (2026-08-10, `subprocess` mode): 13.2, 11.4, 11.8, 12.3 min, all 30/30 cameras, 0 failures -- confirms the tool's own overhead is negligible and the timing is stable across different exposures/nights at this scale.

### Clean failure reporting (`--backend python` only)

Because each camera is an independent subprocess (not one MPI collective), a bad camera never blocks the rest, and failures are reported immediately with a reason instead of a bare nonzero exit code:
*   **Missing input file(s)** (e.g. a preproc file production never generated for that camera/expid): reported as `SKIPPED` with the missing path(s) -- confirmed on `z7@20250822/00307722`, which resolved in **1.9s**. The equivalent `--backend cpp` run on the same case took **~12 minutes to hang** (see the MPI-hang gotcha below) before it had to be killed manually.
*   **Per-camera subprocess failure** (nonzero rc): the run summary prints the last few lines of that camera's log (`tail_error()`) inline under a `PROBLEM:` entry, so the cause is visible without opening individual log files.

### 2.4 Running on CPU (no GPU available)

Every mode above also runs on CPU -- useful on a GPU-less machine (a laptop, a login node, a non-Perlmutter dev box), or just to sanity-check a result without touching a GPU at all. **Pass `--backend cpu` explicitly** (the CLI defaults to `--backend gpu`, and per Section 0's fail-fast check, `--backend gpu` on a machine with no CUDA-enabled jaxlib now errors immediately rather than silently doing the wrong thing):

```bash
python -m specex.specex \
    -a .../preproc-z8-00344649.fits.gz --in-psf .../shifted-input-psf-z8-00344649.fits \
    --out-psf $SCRATCH/pyfit-psf-z8-00344649.fits \
    --broken-fibers 473,474 \
    --backend cpu --cpu-workers 20
```

`--cpu-workers` is the CPU-backend analog of `--gpu`/`--workers-per-gpu` combined -- it's the total number of concurrent bundle-fit worker processes (defaults to the `--gpu` value, 4, if unset, which is usually too low; set it to roughly your machine's real core count instead, e.g. `20` on the interactive-node cases below). Each worker's own thread budget (`OMP_NUM_THREADS` etc.) is auto-computed as `available_cores // cpu_workers`, so you don't need to hand-tune per-worker threading on top of this -- just pick a sensible `--cpu-workers` for the box you're on. Correctness is validated equivalent to the GPU backend (bit-for-bit-consistent to within the standing mixed-precision float32/float64 difference documented in Section 6).

**GPU vs. CPU timing** (both warm, i.e. past the first-run JIT cost from Section 2's intro):
| Scope | GPU | CPU | Slowdown |
|---|---|---|---|
| Single bundle, solo | ~23s | ~80s (`--cpu-workers` >= bundle count, no queueing) | ~3.5x |
| Full night, 30 cameras | 13.2 min (1 node/4 GPUs, `subprocess` mode, pinned `wpg=10 b/7 r/4 z`, Section 2.3) -- or ~8-8.4 min with `--worker-mode persistent` | 147.9 min (1 node, `--cpu-workers 20`, one `python -m specex.specex --backend cpu` call per camera) | ~11.2x (subprocess) / ~17-18x (persistent) |

**`testing/run_night.py` does not have a CPU-only mode** -- its `--backend` choices are `cpp`/`cpp-direct`/`python`, and `python` always assigns each camera a GPU. The 147.9 min full-night CPU number above was measured with a predecessor one-off script (pre-`run_night.py` consolidation) looping `python -m specex.specex --backend cpu --cpu-workers 20` over all 30 cameras sequentially; there's no single documented command for it today -- for a full CPU-only night, write the same kind of loop over Section 3's `select_test_case.py` cases.

CPU is a genuine, correctness-equivalent fallback, not a performance option -- expect roughly an order of magnitude slower at production scale. **Don't try to combine CPU and GPU workers on the same node to "help" a GPU run go faster**: this was tested extensively (six distinct concurrency/pinning designs, `docs/python-port/porting-notes.md`'s CPU+GPU hybrid investigation) and every design either left the GPU run unaffected at best or measurably slowed it down (up to ~3.6x on the cameras that overlapped) -- CPU-side host compute contends with the GPU workers' own host-side work for memory bandwidth. Pure-GPU-pinned is the standing, validated production recommendation; use CPU-only when there's truly no GPU, not alongside one.

### Runtime guidelines (warm JIT cache, all modes)

| Mode | Recommended flags | Typical warm time |
|---|---|---|
| 2.1 Single bundle | `--gpu 1` | a few seconds |
| 2.2 Full CCD, one camera | `--gpu 4 --workers-per-gpu 5` | ~33-45s (this section's measurement; band-dependent, see `docs/python-port/porting-notes.md`) |
| 2.3 Full night, 30 cameras, 1 node/4 GPUs, `--worker-mode persistent` (recommended) | `run_night.py --backend python --worker-mode persistent` | ~483-503s (8-8.4 min) |
| 2.3 Full night, 30 cameras, 1 node/4 GPUs, `--worker-mode subprocess` (default) | `run_night.py --backend python` | ~11-13 min |
| 2.3 Full night, 30 cameras, 2 nodes/8 GPUs, LPT-balanced (`subprocess` mode) | `run_night.py --backend python --lpt-profile ...` | ~6 min |
| 2.4 Full night, 30 cameras, CPU-only | `specex.specex --backend cpu --cpu-workers 20`, looped per camera (no `run_night.py` support yet) | ~148 min (~11x slower than GPU) |

The very first run in a fresh environment (empty `~/.cache/specex/jax_compilation_cache`) will be several times slower than this table for whichever mode you run first -- that cost only has to be paid once per machine/environment, not once per run.

---

## 3. Finding Test Data

If you want to run on a specific night or exposure, use `select_test_case.py` to find the correct file paths and parameters (like `--broken-fibers`).

```bash
# List all cases for a specific night
python testing/select_test_case.py --night 20260401 --list

# Select a random case for testing
python testing/select_test_case.py --night 20260401 --random
```

For picking several *distinct* random nights (e.g. to build an independent test set that avoids reusing the same night twice), see `testing/random_case_picker.py`.

---

## 4. Running Parity Comparisons (Recommended)

The primary tool for verifying Python vs C++ metrics is `testing/instrumentation_analysis.py`. It runs both the production C++ code and the new JAX-GPU code on a single bundle and compares results.

### Basic Usage:
```bash
python testing/instrumentation_analysis.py --night 20260401 --expid 00344649 --cameras b0,r3,z8 --bundle 5
```

### Metrics Produced:
The script generates a table (default output: `instrumentation_analysis.txt`) containing:
*   **Time(s):** Wall-clock time for the fit.
*   **Chi2:** Final chi-squared value (lower is generally better).
*   **Spots:** Number of spots identified and used in the fit.
*   **XT RMS / YT RMS:** Root-Mean-Square difference in pixels between the Python-fitted traces and the C++-fitted traces.

To compare **CPP** vs **JAX-CPU** vs **JAX-GPU** simultaneously, use `testing/validate_all_modes.py`:
```bash
python testing/validate_all_modes.py --cameras b0 --bundle 5 --output comparison_results.txt
```
This is useful for verifying that GPU acceleration doesn't introduce numerical divergence from the CPU version of the same code.

For a full night/expid, whole-CCD C++-vs-Python parity+timing sweep across many cameras at once, see `testing/full_ccd_campaign.py` (Section 2.2 above covers it as a production-scale driver too).

---

## 5. Full CLI Reference

```
python -m specex.specex -h
```

**Required:**
| Flag | Description |
|------|-------------|
| `-a`, `--arc`, `--input-image` | Input preproc arc image |
| `--in-psf`, `--input-psf` | Input (shifted) PSF file |
| `--out-psf`, `--output-psf` | Output PSF file path |

**Scope (bundle/fiber range):**
| Flag | Default | Description |
|------|---------|-------------|
| `--first-bundle` | 0 | First bundle to fit (0-19) |
| `--last-bundle` | 19 | Last bundle to fit (0-19) |
| `--first-fiber` | (unset) | First fiber to fit |
| `--last-fiber` | (unset) | Last fiber to fit |
| `--broken-fibers` | (unset) | Comma-separated fiber IDs to exclude from the fit |

**Fit basis / degrees:**
| Flag | Default | Description |
|------|---------|-------------|
| `--legendre-deg-wave` | auto (3 for z-band, 1 otherwise, from the CAMERA header) | Legendre degree for the joint fit's PSF-shape wavelength basis |
| `--trace-legendre-deg-wave` | auto per axis | Legendre degree for the trace-position wavelength basis, both axes at once; overridden per-axis by the two flags below if given |
| `--trace-legendre-deg-wave-x` | auto (same as `--legendre-deg-wave`) | Trace-position X basis degree only |
| `--trace-legendre-deg-wave-y` | auto (2 for b/r, same as `--legendre-deg-wave` for z) | Trace-position Y basis degree only |
| `--trace-per-fiber-deg` | **6** (production default since 2026-08-05) | Block-diagonal-by-fiber trace basis instead of a shared one, paired with the ndead-gated trace prior below. `0` reverts to the old shared basis. Validated: 30/30 cases improved on xrms (-37.5% mean) and yrms (-60.7% mean) vs. shared-basis, for a ~9% timing cost. |
| `--trace-prior-deg` | 1 | Degree at/above which per-fiber trace coefficients are pulled toward the bundle's cross-fiber consensus. Only active with `--trace-per-fiber-deg` on. Negative disables the prior while keeping per-fiber trace on. |
| `--trace-prior-weight` | 1e5 | Trace-prior penalty weight (C++'s own 1e8 measurably harms healthy bundles applied blanket-style; this is the corrected value) |
| `--trace-prior-ndead-threshold` | 500 | A fiber's dead-pixel count (ndead) above this triggers the trace prior for that fiber only; normal fibers (ndead ~20-120) are unaffected |
| `--fit-continuum` / `--no-fit-continuum` | auto (on for z-band, off otherwise, matching C++) | Fit a per-bundle continuum background |

**Compute / concurrency:**
| Flag | Default | Description |
|------|---------|-------------|
| `--backend` | `gpu` | `cpu` or `gpu` |
| `--gpu` | 4 | Number of GPUs to use (spreads bundles across them on the current node) |
| `--workers-per-gpu` | auto (5, or 3 for z-band when per-fiber trace is active) | Concurrent bundle-fit worker processes per GPU. See Section 2.3's table for the validated per-band values (10 b / 7 r / 4 z) used in production-scale multi-camera runs. |
| `--cpu-workers` | `--gpu` count | Concurrent worker processes for `--backend cpu` |
| `--gpu-worker-threads` | unconstrained | Diagnostic: force an OMP/BLAS/XLA thread cap on each GPU-backend worker's host-side computation. Confirmed *not* load-bearing for the CPU+GPU hybrid-scheduling investigation (see `docs/python-port/porting-notes.md`) -- left in as a diagnostic knob, no effect on a normal run. |
| `--double-precision` | off (mixed float32/float64) | Force full float64 for the joint-fit Jacobian. Validated equivalent accuracy; mixed precision uses ~71% less GPU memory/worker. |

**Spot selection:**
| Flag | Default | Description |
|------|---------|-------------|
| `--sn-threshold` | 3.0 | S/N threshold for spot selection |
| `--max-lines` | 200 | Maximum number of lines to keep per bundle |
| `--h-size-y` | 5 | Override PSF stamp half-size in Y |
| `--force-spots` | (unset) | Path to a file of spots to fit (`fiber,wave,xc,yc`), bypassing spot selection |

**Other:**
| Flag | Default | Description |
|------|---------|-------------|
| `--lamp-lines` | (unset) | Lamp lines file path |
| `--line-search` | `grid` | EXPERIMENTAL: final joint fit's per-iteration step-size search. `grid` (default, coarse 3-point), `brent` (continuous, not C++-faithful), `cpp` (faithful replica of C++'s Numerical Recipes brent() + mode-dependent skip logic). Both alternates tested correctness-neutral vs `grid` -- kept for reference only. |
| `--debug-spots` | off | Write per-pass spot-selection debug dump files (`.pyrawspots.txt`, `.pyspots_pass*.txt`, etc. -- the Python analog of C++'s own `--debug-spots`). Off by default: adds per-bundle-worker I/O overhead with no effect on the fitted output. |

Alternatively, call it from within another Python script via `fit_ccd_native()`:

```python
from specex.specex import fit_ccd_native

fit_ccd_native(
    arc_file='path/to/preproc.fits',
    in_psf_file='path/to/input-psf.fits',
    out_psf_file='output-psf.fits',
    lamp_lines_file='py/specex/data/specex_linelist_desi.txt',
    broken_fibers="473,474",
    gpu=1, workers_per_gpu=10,
)
```

---

## 6. Understanding the Metrics

*   **Spots:** If the Python version identifies significantly fewer spots than C++, check the `--sn-threshold` parameter.
*   **Chi2:** A 2-4% difference is currently expected due to "Dead Column Masking" differences in the pre-processor.
*   **Trace Deltas (XT/YT RMS):** Typical whole-CCD results are ~0.02-0.03 pixels (well under the 0.05px target), using the `--trace-per-fiber-deg 6` production default.
*   **Chi2: -1.0:** Usually indicates a crash or a failure to find any spots.
*   **`WARNING: Bundle N failed` in a run's log:** A per-bundle failure (commonly GPU `RESOURCE_EXHAUSTED`/OOM) that `fit_ccd_native` logs and *continues past*, merging the rest of the camera anyway -- a `0` process exit code alone does **not** mean every bundle succeeded. Always grep multi-camera run logs for this string (or `RESOURCE_EXHAUSTED` directly) before trusting a "0 failures" summary.

---

## 7. Known Operational Gotchas (private-SPECPROD reruns)

These apply when using `run_night.py` (either backend) against a night/expid you haven't run through real production yourself -- i.e. almost any from-scratch rerun via a private `DESI_SPECTRO_REDUX`/`SPECPROD`.

### 7.1 `--backend cpp`: private reruns are missing calibration state real production already has

A from-scratch `desi_proc` rerun only processes the single requested exposure -- it never runs the full calibration-night pipeline that real production uses to pre-generate certain per-night calibration products. Two distinct manifestations confirmed this session:
*   **Missing CTE (charge-transfer-efficiency) correction file**: `RuntimeError: Missing .../calibnight/<night>/ctecorr-<night>.yaml`. **Free, zero-compute pre-screen**: read that yaml directly on the real `matterhorn` production tree -- `[]` (empty list) means no camera on that night needs it (clean rerun); a populated list names the exact `CAMERA`s that will fail. Confirmed accurate across 3 tested nights; ~76% of scanned 2026 nights are clean.
*   **Missing calibration darks**: `Didn't find matching <cam> calibration darks in $DESI_SPECTRO_DARK` -- same root cause, no cheap pre-screen found yet. Even a fully production-validated arc exposure hit this on one camera (b0) in testing -- **29/30 or 30/30 looks like a realistic ceiling for a "clean" private rerun**, not a bug worth chasing further.

`--backend python` is immune to both: it reads production's already-generated `preproc-*`/`shifted-input-psf-*` files directly and never invokes `desi_proc`, so it doesn't care whether the private redux tree has calibration state or not.

### 7.2 `--backend cpp`: any single MPI rank failure hangs the whole job

Confirmed repeatedly, across 3 distinct trigger types (CTE gap, missing preproc input, missing calibration darks): when any one of the ~100 MPI ranks in a `desi_proc --mpi` job fails, the **entire job hangs indefinitely** (log stops updating, zero further output) instead of exiting cleanly -- even though the processes stay in state `R` burning ~98%+ CPU, not zombie/D-state. Waiting does not resolve it; it must be killed manually:
1. `ps aux | grep srun` -- find the `srun` frontend PID(s) for the job.
2. `kill -9 <pid...>` on those. **This is sometimes not enough** -- it can fail to cascade to the actual worker processes.
3. Verify: `ps aux | grep desi_proc | grep -v grep | wc -l`. If nonzero, `pkill -9 -f "desi_proc -n <night> -e <expid>"`, then re-check until it's 0.
4. Do **not** `scancel` the whole SLURM allocation unless you intend to end the whole interactive session -- the hang is a job-level problem, not a node-level one.

This is exactly what motivated `--backend python`'s per-camera clean-failure reporting in Section 2.3 above -- it's structurally immune to this failure mode since cameras are independent subprocesses, not MPI ranks in one collective.

### 7.3 `find_cases()` (used by `--backend python`) only sees exposures processed through the normal nightly pipeline

`run_night.py --backend python` locates each camera's input file paths and `--broken-fibers` list by scraping `arc*.log` files under matterhorn's `run/scripts/night/<night>/` directory -- the logs written by production's own SLURM job-script workflow. An expid that was never run through that workflow (e.g. one flagged in the exposure table as failing, and consequently skipped by production) has **no arc log entry at all**, so `find_cases()` silently skips it (`WARNING: no arc log entry found ... skipping`), even if the raw data and a preprocessed version exist somewhere.

Confirmed on `20211028/00106396`: no arc log exists anywhere in production for that expid (consistent with its exposure-table comment, "fails psf fitting" -- production apparently never attempted it). Worked around by building the camera-to-file-path map manually, pointing at a private-redux preprocessing pass (from a `--backend cpp` run of the same expid, which does its own preprocessing regardless of the scripts/night logs), and borrowing the `--broken-fibers` list from the *next* arc exposure taken the same night (fiber breakage is a persistent hardware property, not a per-exposure one -- confirmed identical camera-by-camera across the two nearby expids where compared). There's no CLI flag for this yet; it requires a short one-off script (see `run_night.py`'s `find_cases()`/`run_node_python()` for the pieces to reuse).

## 8. Comparing C++ vs Python: Common Scripts

Three small scripts, used together to get a true apples-to-apples C++ vs Python comparison for a night/expid. This is the same methodology behind the cumulative 20-night campaign in Section 2.3 above (current figures: xrms=0.0120px/yrms≈0.0129px, Python ~18% faster with `--worker-mode persistent`) -- see `docs/python-port/porting-notes.md`'s 2026-08-11/12 entry for the original 10-night campaign this workflow was built for (superseded numbers, from before the 2026-09-01 footprint-margin fix; kept for history, not current figures) and its 2026-09-02/05 entries for the current ones.

### 8.1 `testing/stage_preproc.py` -- make a real `--backend cpp` run skip its own (buggy/incomplete) preprocessing

A from-scratch `--backend cpp` run's own idempotent preprocessing pass isn't reliable (see 7.1) -- it can silently produce fewer than 30 cameras' `preproc-*.fits.gz` even on a clean, CTE-free night. Since `--backend python` already succeeds against the exact same night/expid, the real production preproc files demonstrably exist; the fix is to hand `desi_proc` copies of them at the exact paths its own idempotency checks look for, so it skips regenerating them entirely and goes straight to the real `-n101` fit stage:

```bash
python testing/stage_preproc.py --night 20260316 --expid 00342128 \
    --redux-dir /pscratch/.../desiproc2/20260316_00342128/redux --specprod cdwarner
python testing/run_night.py --night 20260316 --expid 00342128 --backend cpp \
    --outdir /pscratch/.../desiproc2/20260316_00342128 \
    --redux-dir /pscratch/.../desiproc2/20260316_00342128/redux --specprod cdwarner
```

Idempotent (safe to re-run); errors out up front if any camera has no arc-log entry (see 7.3) rather than silently staging a partial set.

### 8.2 `testing/compare_correctness.py` -- per-camera and per-night xrms/yrms across a list of nights

Same trace-RMS methodology as `full_ccd_campaign.py` (Legendre trace polynomials evaluated on a 100-point wavelength grid, `--broken-fibers` excluded), extended to loop over as many night/expid pairs as you give it:

```bash
python testing/compare_correctness.py \
    --py-dir /pscratch/.../cpptest/redux/cdwarner \
    --cpp-base /pscratch/.../cpptest/desiproc2 --specprod cdwarner \
    --nights 20260316:00342128,20260401:00344649,20220120:00119496
```

`--py-dir` is wherever `--backend python` wrote `fit-psf-<cam>-<expid>.fits` (its flat `--outdir`); `--cpp-base` is the parent of each night's `<night>_<expid>/redux/<specprod>/exposures/<night>/<expid>/fit-psf-<cam>-<expid>.fits` (i.e. each night/expid's own `--outdir` from the `stage_preproc.py`+`--backend cpp` step above, one level up from `redux`).

### 8.3 `testing/per_fiber_breakdown.py` -- drill into one camera's elevated RMS

When `compare_correctness.py` flags one camera's mean xrms/yrms as elevated, this breaks it down per-fiber to distinguish "one or a few genuinely bad fibers" (a data-quality issue, exclude and move on) from "a systematic offset across the whole camera" (worth investigating further):

```bash
python testing/per_fiber_breakdown.py --night 20241021 --expid 00259030 --camera z6 \
    --py-dir /pscratch/.../cpptest/redux/cdwarner \
    --cpp-base /pscratch/.../cpptest/desiproc2 --specprod cdwarner --top 20
```

Prints the top-N fibers by yrms and by xrms, with their bundle number (`fiber // 25`). A cluster of top offenders at literal bundle-boundary positions (first/last of a 25-fiber bundle) was the signature of the bundle-boundary trace divergence -- root-caused and closed by the `--footprint-margin` fix (`docs/python-port/porting-notes.md`, 2026-09-01; see `docs/python-port/python-vs-cpp-diff.txt` section 1.3(a) for the before/after), so seeing that pattern again on current code points to a regression, not the old known issue.
