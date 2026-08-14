# Known edge cases

Enumeration of non-standard inputs and failure/divergence modes discovered while
validating the Python/JAX port against real production data. Each entry says what
triggers it, how it presents, and whether it's a Python-port concern or something
upstream/architectural. Full narrative detail and dates are in `porting-notes.md`;
this file is the quick-reference index.

## Input-data edge cases

### Bad/dead amp, masked within an otherwise-complete exposure
- **Trigger**: a CCD amplifier is known-bad; `desispec.scripts.preproc` is run with `--badamps <ampname>`, which masks that quadrant's pixels but still writes a full-size `preproc-<cam>-<expid>.fits.gz`.
- **Presents as**: fit divergence from C++ localized to exactly the masked-amp half of the CCD (e.g. fibers 0-249 for an amp-A mask, 250-499 for amp-B) -- confirmed via full 20-bundle spatial breakdown on `r8@20211028/00106399`/`00106400` and `b8@20221121/00154099`/`00154101`. Whole-camera xrms/yrms can look alarming (1-12px) until broken down spatially; the healthy half of the same exposure matches the normal ~0.02-0.03px baseline.
- **Verdict**: expected, not a bug. Over a masked amp there's no real arc-line signal, so the trace fit there is fundamentally underdetermined; C++ and Python's different priors/regularization/extrapolation naturally diverge from each other in that region. See porting-notes.md 2026-08-13.
- **Distinguish from** "missing preproc file entirely" below -- a masked amp is present data with reduced signal, not absent data.

### Preproc file(s) completely missing for an expid/camera
- **Trigger**: raw-data preprocessing itself failed or was skipped upstream (before specex is ever invoked), OR production deliberately excluded a camera from an expid rather than let it fail downstream.
- **Examples**: `20211028/00106396` -- zero preproc files for the entire expid. `20250822/00307722`, camera z7 -- preproc file for z7 alone doesn't exist (all other 29 cameras present); production appears to have excluded it outright.
- **Presents as**: nothing to fit -- neither pipeline can produce output without input.
- **Verdict**: out of scope for specex entirely, C++ or Python. `run_night.py --backend python` reports this as a clean `SKIPPED` in ~2s; real C++ (`desi_proc`) instead hits the MPI-hang pattern below when the missing camera is expected. Not something a fitting-engine fix can address.

### Private-`SPECPROD` `desi_proc` rerun missing per-night CTE-correction calib files
- **Trigger**: real production pre-generates `calibnight/<night>/ctecorr-<night>.yaml` via full calibration-night processing; a from-scratch single-exposure `desi_proc` rerun (private SPECPROD) never creates it. Any camera the night's characterization says needs CTE correction then fails preprocessing (`RuntimeError: Missing .../ctecorr-<night>.yaml`).
- **Cheap pre-screen**: read production's own `calibnight` yaml directly -- `[]` = clean, a populated list names the exact affected cameras and predicts failure exactly (confirmed on 20230207: z1+z3, 20230805: r6+z1+z3). z1/z3 needed it on both 2023 nights tested -- looks like a persistent hardware property of those spectrographs. All of matterhorn's 2026 calibnight dirs scanned (77 nights) show an empty CTE list -- 2026 data is a safer era to pick random test nights from.
- **Verdict**: an artifact of doing private-SPECPROD reruns from scratch, not a specex or production bug. `--backend python` is structurally immune (reads production's already-generated preproc/PSF files directly, never invokes `desi_proc`) -- confirmed both CTE-gap nights run 30/30 clean in Python.

## Fitting-behavior edge cases

### High-ndead fiber near a fiber-crossing region
- **Trigger**: a fiber has many CCD rows the trace can't directly measure (`ndead`) in a spatial region where neighboring fiber traces cross/overlap, making the local trace geometrically ambiguous.
- **Examples**: `z7@20250822/00307725`, fiber 251 (ndead=2907) -- isolated single-fiber degradation (xrms=1.68px, yrms=3.35px; next-worst fiber in camera ~0.12-0.14px); excluding it, camera yrms drops back to baseline (0.0347px). `z7@20250822/00307722`, bundle 10 (fibers 250-253, ndead 848/5515/28896/3174) -- **real C++ hard-crashes** here (`desi_psf_fit on process 10 failed with return value 1`, no output file written at all); Python's `--trace-prior-ndead-threshold`-gated cross-fiber trace-consensus prior (default threshold 500) automatically activates for these fibers and converges cleanly.
- **Verdict**: Python's ndead-gated prior is a genuine improvement here, not just parity -- C++'s own equivalent trace-consensus mechanism (`trace_prior_deg`) defaults to *off* in real production and is never enabled, so C++ has no equivalent protection in practice. Confirmed by reading `src/specex_psf_fitter.cc` directly.

### Isolated single-fiber trace degradation not near ndead/fiber-crossing
- **Trigger**: not fully pinned down. `z6@20241021/00259030`, fiber 343 (bundle 13) diverges badly (yrms=2.60px vs. camera norm ~0.02-0.03px) despite not having an unusually high ndead or an obviously-broken input.
- **Investigated and falsified**: max-contiguous-dead-row-run as a smarter per-fiber gate (fiber 397 has a longer dead run than 343 but doesn't need the prior); local arc-line-coverage-gap heuristic (fiber 303 has an even bigger gap than 343 but doesn't need the prior either).
- **Verdict**: no simple pre-fit static metric (from weight+linelist data alone) cleanly separates this genuine case from ~42 false positives a naive threshold drop would create. Left open, low priority, current default (`ndead>500`) stays as the aggregate-optimal setting. If revisited: a post-hoc detect-and-refit scheme (flag fibers whose fitted trace looks discontinuous vs. neighbors after a normal fit, selectively re-fit those with the prior on) is the more promising next avenue, not more pre-fit heuristics.

### Bundle-boundary trace weakness, amplified camera-wide on one night
- **Trigger**: normally-bounded (0.05-0.2px) known mechanism at bundle-boundary fiber positions (first/last-of-25 in each bundle) -- but on `r9@20220120/00119496` this was amplified far more severely (0.1-0.7px) and broadly across the entire camera (xrms/yrms 0.1043/0.0629px, worst single camera in the 10-night campaign), correlated with elevated fit iteration counts (876 vs 640 on a clean-night r9) and higher raw chi2.
- **Ruled out as causes**: short ARC exposure time (checked directly -- every arc across all 10 campaign nights is ~5.01s, not a variable), CTE gap (clean), spot yield (normal), trace-correction magnitude (essentially identical to a clean night).
- **Verdict**: known mechanism, amplified; trigger not isolated (interleaved multi-bundle logs made clean per-bundle chi2 comparison impractical without more instrumentation). Isolated to this one night -- 3 additional nights run afterward showed no recurrence. Not blocking.

### Broken/dead fibers (hardware)
- **Trigger**: known persistent hardware fault on specific fibers, tracked per-camera by production and passed via `--broken-fibers`.
- **Example**: `z5@20220120/00119496` carries a genuine 28-fiber broken block (164-191).
- **Verdict**: expected and handled -- always exclude via `--broken-fibers` in any comparison, or a few fibers show fake catastrophic divergence in xrms/yrms (feedback memory: `feedback_broken_fibers.md`). Fiber breakage is a persistent hardware property, so broken-fiber lists can safely be borrowed from a neighboring same-night arc exposure when an expid's own arc log doesn't exist (see `find_cases()` blind spot below).

### GPU OOM silently continuing past a failed bundle
- **Trigger**: per-bundle GPU memory exhaustion during `fit_ccd_native`.
- **Presents as**: `rc==0` at the process level, but the bundle's output is missing/wrong -- looks like a z-band-specific correctness regression unless logs are checked.
- **Verdict**: an infrastructure/resource-sizing issue, not an algorithmic weakness -- always grep logs for `WARNING: Bundle` / `RESOURCE_EXHAUSTED`, never trust `rc==0` alone.

## Infrastructure / orchestration edge cases

### C++ MPI job hangs indefinitely on any single-rank failure
- **Trigger**: any one rank in a `desi_proc -n101` or `desi_compute_psf --mpi -n20` job fails (confirmed 3x across different failure causes).
- **Presents as**: the failed rank logs its error and exits in <1s, but the whole job produces zero further output and sits burning CPU (state `R`, not zombie/D-state) for 10-35+ minutes until manually killed. Must `kill <srun-frontend-pid>`, not `scancel` the whole allocation.
- **Verdict**: a real C++/MPI architectural weakness. `run_night.py --backend python`'s per-camera-subprocess design is structurally immune (no MPI collective between cameras) -- this was a direct motivation for that design.

### `find_cases()` can't locate expids never run through production's normal night pipeline
- **Trigger**: an expid whose original run bypassed the standard automated scripts (`arc*.log` never generated), so `testing/run_night.py`'s log-scraping `find_cases()` (and `testing/select_test_case.py`'s `parse_log_line()`) finds nothing.
- **Example**: `20211028/00106396` has zero arc log entries anywhere in the matterhorn production tree, consistent with its exposure-table "fails psf fitting" comment.
- **Workaround**: manual case dict with explicit input paths + `--broken-fibers` borrowed from a neighboring same-night arc exposure (fiber breakage is a persistent hardware property, safe to borrow). Documented in `how-to-run.md` Section 7.3.

### desispec PR #2732: aggregate rchi2==0 false-positive across multi-arc nights
- **Trigger**: a bad amp exists in only *some* of a night's multiple arc exposures. `desi_proc`'s aggregate failure-detection heuristic (counting fibers with `rchi2==0` as a fit-failure proxy) normally requires ALL arcs to show the same bad-fiber pattern before attributing it to real hardware; a partial pattern across arcs produces a false positive.
- **Verdict**: lives entirely in `desispec`'s orchestration layer (`py/desispec/scripts/specex.py`, a different repo), not in specex's fitting engine (C++ or Python). Confirmed by running the underlying per-camera fits directly against real masked-amp input -- see "Bad/dead amp" above. No specex-side fix needed; useful confirmation that fits behave sanely (and in one case, `z7@20250822/00307722` bundle 10, better than C++, which hard-crashes there) under exactly these input conditions.
