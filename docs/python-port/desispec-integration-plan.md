# Linking the GPU/JAX port into desispec: how it's called today, and a plan

Stephen's ask (2026-09-13): desispec currently hands PSF-fitting work off to
a bunch of MPI ranks running the C++ code. How would the GPU-native
Python/JAX port play with that, and what's the integration plan? This
investigates `../desispec` (checked out locally, `main` branch) and proposes
a concrete, phased plan. Nothing here has been implemented in desispec
itself -- this is a plan, not a patch.

## How desispec calls specex today

`desispec/py/desispec/scripts/specex.py` is the whole of it. Two entry
points:

- **`main(args, comm)`** -- fits one camera. Splits the camera's 500 fibers
  into 20 bundles of 25, assigns bundles round-robin across `comm`'s ranks
  (`nproc`), and **each rank calls `run_specex(com)` once per bundle it
  owns, in a plain Python loop** -- `com` is the same CLI-argument-list
  convention (`-a`, `--in-psf`, `--out-psf`, `--first-bundle`,
  `--last-bundle`, `--first-fiber`, `--last-fiber`, `--legendre-deg-wave`,
  `--broken-fibers`, ...) our own `main()` accepts, since our Python CLI was
  built to mirror the original C++ one. Each bundle writes its own file
  (`{outroot}_{b:02d}.fits`); rank 0 then merges them into one final PSF
  file via `merge_psf()` (overlays each bundle's fitted `XTRACE`/`YTRACE`/
  `PARAM` coefficients onto a copy of the input PSF, keyed by `STATUS==0`
  fibers).
- **`run(comm, cmds, cameras)`** -- the real per-exposure driver, called
  from `desispec.scripts.proc` (`desi_proc`). Uses
  `desispec.workflow.schedule.Schedule` to split the whole MPI allocation
  into **one 20-rank group per camera** (`group_size=20`, matching the
  20-bundle structure above 1:1) and runs `main()` once per camera inside
  its own group.

So today's architecture is: **20 MPI ranks per camera, each rank fitting
~1 bundle by calling the C++ wrapper directly, synchronous, no GPU
awareness anywhere in this file.**

## The architectural mismatch

Our GPU-native path (`fit_ccd_native()`, `py/specex/specex.py`) is
per-**camera**, not per-bundle: it's a single call that internally spins up
its own `multiprocessing.Pool` of workers spread across however many GPUs
are visible, fits all 20 bundles, and writes one already-merged output file
directly (`write_python_psf()`) -- there's no per-bundle file, no external
merge step. Mapping this onto desispec's "20 ranks, 1 bundle each" model
would mean either (a) fighting the grain of the code to make 20 MPI ranks
somehow share 4 GPUs one bundle at a time (reinventing what
`fit_ccd_native`'s own worker pool already does, badly, at the wrong
layer), or (b) not using those 20 ranks that way at all for the GPU path.
(b) is the right answer, and there's already a precedent for exactly this
in the same codebase.

## The precedent: `gpu_specter` extraction is already wired in exactly this way

`desispec/py/desispec/scripts/proc.py` already branches on
`use_gpu = is_gpu_available() and not args.no_gpu` (`desispec.gpu`, checks
`cupy`/`numba.cuda`) for the spectral **extraction** step, and the GPU
branch looks nothing like the CPU branch's many-small-MPI-groups model:

```python
if args.use_specter:
    extract_subcomm_size = 20                  # CPU, specter: 20 ranks/camera, many groups
elif use_gpu:
    extract_subcomm_size = 2 + 5 * ngpus        # GPU: one small group, sized to the node's GPUs
else:
    extract_subcomm_size = 16                   # CPU, gpu_specter: 16 ranks/camera, many groups

if use_gpu:
    extract_group, num_extract_groups = 0, 1    # ONE GPU group processes every camera, one at a time
else:
    extract_group = rank // extract_subcomm_size
    num_extract_groups = size // extract_subcomm_size   # many CPU groups run in parallel
```

That is: for GPU work, `desi_proc` carves out **one** small subcommunicator
sized to the node's GPU count (not the camera count), and that one group
loops over every camera in the exposure sequentially, calling
`desispec.scripts.extract.main_gpu_specter()`. The many-parallel-MPI-groups
model is CPU-only; GPU work gets a single dedicated group because GPUs are
the scarce, shared resource, not ranks.

**specex's GPU integration should follow this exact pattern**, not
reinvent one.

## Proposed plan

### Phase 1: a `use_gpu` branch in `desispec.scripts.specex`, off by default

Add the same `is_gpu_available()`/`--no-gpu` convention already used for
extraction. When GPU is requested and available:

- Skip the whole bundle-split/`Schedule`/`run_specex()`-per-bundle path
  entirely for that camera.
- One rank (rank 0 of a small GPU-sized subcomm, exactly like the
  extraction precedent -- no need for the other ranks in that subcomm to do
  anything specex-specific, since `fit_ccd_native()` manages its own
  multi-GPU parallelism internally via `multiprocessing`, not MPI) calls
  `fit_ccd_native()` once per camera, looping over every camera in the
  exposure that group owns -- one call per camera, not per bundle.
- No `merge_psf()` step needed for this path: `fit_ccd_native`/
  `write_python_psf` already write one complete, already-merged output
  file, matching the final `fit-psf-<cam>-<expid>.fits` `desi_proc` expects
  directly.
- Argument mapping is close to free: `desispec.scripts.specex.main()`
  already builds a `com` list with `-a`/`--in-psf`/`--out-psf`/
  `--broken-fibers`/`--legendre-deg-wave`/`--fit-continuum` -- the exact
  flag names our own `main()`/`fit_ccd_native()` already accept (by
  design, since our CLI mirrors the C++ one). The GPU branch can reuse
  that same argument-construction logic almost unchanged, just calling
  `fit_ccd_native()` (or a thin wrapper around it) instead of
  shelling out per bundle.

### Phase 2: validate at the desi_proc/MPI layer, not just specex's own harness

specex's own 20-night campaign (`current-status.txt`) validates
`run_night.py`'s orchestration, not `desi_proc`'s. Once Phase 1 lands,
re-run correctness/timing through real `desi_proc --mpi` invocations on a
GPU-allocated node, to catch anything specific to the MPI/`Schedule`
layer (rank-to-GPU binding under SLURM's GPU allocation, interaction with
`desi_proc`'s own logging/error-handling conventions, etc.) that a
standalone `run_night.py` run wouldn't exercise.

### Phase 3: wire the flag through `desi_proc`'s own CLI

`desi_proc` already requests GPU nodes and threads `use_gpu` through for
extraction. The natural endpoint is one GPU-node allocation covering both
extraction (`gpu_specter`) and PSF fitting (this work) in the same nightly
job, rather than requesting GPU nodes twice. A single `--no-gpu`-style
flag (or reusing the existing one) can gate both.

### Phase 4 (later, optional): multi-exposure GPU worker persistence

`run_night.py --worker-mode persistent` (our own fastest configuration,
~18% faster than C++) gets its win partly from keeping a GPU worker pool
alive **across cameras of one exposure** (already compatible with Phase 1
above -- `fit_ccd_native` itself doesn't care whether it's called fresh
each time or from a long-lived process) and partly from JAX's **on-disk**
persistent compilation cache, which already survives independently of
process lifetime -- so even a fresh `desi_proc` invocation per exposure
still benefits from a warm cache after the first exposure of the night,
without needing a long-lived process at all. The *process*-persistence
half of that win (avoiding per-exposure interpreter/JAX-import startup,
not just JIT compilation) would need `desi_proc` itself to run as a
longer-lived service across a whole night's exposures rather than being
invoked fresh per exposure -- a real desispec/workflow architecture
question, well beyond specex's own scope, and not needed for a first
rollout given the disk-cache already captures most of the benefit.

## The forcing-function rename: tried, reverted, then redone (2026-09-13 / 2026-09-11)

Briefly renamed `run_specex` -> `run_specex_cpp` (`py/specex/specex.py`)
the same day this plan was written, so that `../desispec`'s current,
unpatched `from specex.specex import run_specex` would fail loudly
instead of silently keeping the old C++-only path alive -- a forcing
function for whoever eventually does Phase 1 below. In practice this
immediately broke `desi_compute_psf`/`desi_psf_fit` (and therefore
`--backend cpp`/`cpp-direct` in `run_night.py`, and
`testing/full_ccd_campaign.py`'s C++ side) on this branch, which was
still actively needed for ongoing testing/validation work -- so it was
reverted the same day, with the rename kept as the intended plan for
once that testing phase was actually done.

**Redone 2026-09-11**, once this branch's own validation work (item 3's
`26.9` environment testing, the cold-cache OOM fix, all the PSF-shape
investigation) was finished and the branch pushed for review. `run_specex`
is `run_specex_cpp` again, for real this time, across every call site and
doc reference in this repo. `../desispec`'s own copy is deliberately left
untouched -- its `from specex.specex import run_specex` now raises
`ImportError` the moment `desi_compute_psf`/`desi_psf_fit` tries the
C++-only path through this branch, which is the whole point: it forces
whoever picks up Phase 1 below to make an explicit choice (patch that
import to `run_specex_cpp`, or switch to the GPU-native `fit_ccd_native()`
path instead) rather than silently keep working unmodified.

## Open questions this plan doesn't answer yet

- **GPU node allocation for the nightly pipeline** is a NERSC/production-ops
  resourcing decision, not something specex or desispec code alone
  determines -- Phase 3 assumes GPU nodes are actually available to
  `desi_proc`'s SLURM job when this path is used.
- **Which MPI ranks get real GPU affinity** under a GPU-node SLURM
  allocation, and how that interacts with `fit_ccd_native`'s own
  `nvidia-smi`-based GPU auto-detection (which currently assumes it can see
  and freely use every GPU visible to its process, unaware of any prior
  rank-to-GPU binding SLURM/desispec might have already set up) needs a
  real test on a GPU-allocated desi_proc job, not just reasoning about it.
- **The C++ wrapper's own environment portability** (`docs/python-port/how-to-run.md`
  Section 0) is a separate, parallel readiness question -- not blocking
  this plan, since Phase 1 bypasses `run_specex()` entirely for the GPU
  path, but worth tracking if the C++ path needs to keep working
  side-by-side in the same environment.
- This plan does not attempt to estimate desi_proc-layer wall-clock
  savings; that's Phase 2's job, on real data, not a paper estimate here.
