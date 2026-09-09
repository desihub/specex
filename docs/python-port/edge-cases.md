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
- **How it's handled (as of 2026-08-14, per Stephen's directive)**: `fitter.py`'s `find_masked_amp_fibers()` detects fibers via a contiguous-run ndead heuristic (default `ndead>8000` across a run of `>=3` adjacent fibers -- distinct from the existing single-fiber `--trace-prior-ndead-threshold` gate, which handles isolated bad-column fibers that ARE still fittable). Detected fibers are excluded from the fit entirely: the input starting-guess PSF (trace AND PSF-shape, not just trace) is propagated unchanged, and `STATUS` is set to `-1`. A bundle that's fully inside a masked amp short-circuits cleanly (`skip_bundle`) rather than being misreported as a failure. Tunable via `--masked-amp-ndead-threshold` (default 8000 -- a first-pass value calibrated on one real case, explicitly not load-bearing precision).
- **Validated against real C++ ground truth**: STATUS matches real production output on 498/500 fibers for `r8@20211028/00106399`, including the exact mid-bundle boundary (fiber 254 flagged, fiber 255 not) -- C++'s own STATUS=-1/pass-through convention for this case, confirmed by direct inspection, not inferred.
- **Verdict**: the underlying divergence is expected, not a bug -- over a masked amp there's no real arc-line signal, so the trace fit there is fundamentally underdetermined. See porting-notes.md 2026-08-13 (root-cause) and 2026-08-14 (handling + validation).
- **Distinguish from** "missing preproc file entirely" below -- a masked amp is present data with reduced signal, not absent data.

### Preproc file(s) completely missing for an expid/camera
- **Trigger**: raw-data preprocessing itself failed or was skipped upstream (before specex is ever invoked), OR production deliberately excluded a camera from an expid rather than let it fail downstream.
- **Examples**: `20211028/00106396` -- zero preproc files for the entire expid. `20250822/00307722`, camera z7 -- preproc file for z7 alone doesn't exist (all other 29 cameras present); production appears to have excluded it outright.
- **Presents as**: nothing to fit -- neither pipeline can produce output without input.
- **How it's handled (per Stephen's directive, 2026-08-14): ok to crash -- it's the caller's responsibility to ensure input files exist.** `fit_ccd_native`/`main()` reads input files with no surrounding try/except, so a missing preproc file already raises an uncaught exception and a non-zero exit; no change was needed. `run_night.py --backend python` reports this as a clean `SKIPPED` in ~2s (its own pre-flight check, one layer above the crash-is-fine core CLI); real C++ (`desi_proc`) instead hits the MPI-hang pattern below when the missing camera is expected. Not something a fitting-engine fix can address.

### Private-`SPECPROD` `desi_proc` rerun missing per-night CTE-correction calib files
- **Trigger**: real production pre-generates `calibnight/<night>/ctecorr-<night>.yaml` via full calibration-night processing; a from-scratch single-exposure `desi_proc` rerun (private SPECPROD) never creates it. Any camera the night's characterization says needs CTE correction then fails preprocessing (`RuntimeError: Missing .../ctecorr-<night>.yaml`).
- **Cheap pre-screen**: read production's own `calibnight` yaml directly -- `[]` = clean, a populated list names the exact affected cameras and predicts failure exactly (confirmed on 20230207: z1+z3, 20230805: r6+z1+z3). z1/z3 needed it on both 2023 nights tested -- looks like a persistent hardware property of those spectrographs. All of matterhorn's 2026 calibnight dirs scanned (77 nights) show an empty CTE list -- 2026 data is a safer era to pick random test nights from.
- **Verdict**: an artifact of doing private-SPECPROD reruns from scratch, not a specex or production bug. `--backend python` is structurally immune (reads production's already-generated preproc/PSF files directly, never invokes `desi_proc`) -- confirmed both CTE-gap nights run 30/30 clean in Python.

## Fitting-behavior edge cases

### High-ndead fiber near a fiber-crossing region
- **Trigger**: a fiber has many CCD rows the trace can't directly measure (`ndead`) in a spatial region where neighboring fiber traces cross/overlap, making the local trace geometrically ambiguous.
- **Examples**: `z7@20250822/00307725`, fiber 251 (ndead=2907) -- isolated single-fiber degradation (xrms=1.68px, yrms=3.35px; next-worst fiber in camera ~0.12-0.14px); excluding it, camera yrms drops back to baseline (0.0347px). `z7@20250822/00307722`, bundle 10 (fibers 250-253, ndead 848/5515/28896/3174) -- **real C++ hard-crashes** here (`desi_psf_fit on process 10 failed with return value 1`, no output file written at all); Python's `--trace-prior-ndead-threshold`-gated cross-fiber trace-consensus prior (default threshold 500) automatically activates for these fibers and converges cleanly.
- **How it's handled (as of 2026-08-14, per Stephen's directive)**: `io.py`'s `write_python_psf` runs a post-fit trace-crossing QA pass on the merged CCD output (adapted from specex#91's dead-code `trace_psf_qa`), evaluating every adjacent fiber pair's fitted X-position across the wavelength range. On an actual geometric crossing, it flags the pair PLUS their immediate neighbors (`f-1, f, f+1, f+2`) as `STATUS=4` -- deliberately more conservative than specex#91's own pair-only version, per explicit instruction ("just in case"). Never-fit fibers (already `STATUS=-1`) are excluded from the comparison on both sides and never downgraded to 4.
- **Important nuance, confirmed against real C++ output for `z7@20250822/00307725`**: real C++'s STATUS=4 fires here (fibers 250/251 flagged in the actual production file) because **C++'s own trace genuinely crosses** for that pair. Python's independently-fit trace for the same two fibers does **not** cross (checked directly: ~7px clean separation throughout) -- the ndead-gated trace-consensus prior evidently avoids the crossing C++'s unregularized fit falls into. So Python's QA pass correctly does NOT fire here: the residual accuracy gap on fiber 251 is a Python-vs-C++ disagreement, not a self-detectable geometric problem, and a geometry-only QA pass structurally can't (and per explicit decision, shouldn't) flag it.
- **Confirmed no C++ precedent for an ndead-based STATUS flag** (`grep fit_status src/*.cc`): C++'s only STATUS-setting paths are the generic fit-failure enum (bundle-uniform, triggered by solver failure) and specex#91's crossing-specific STATUS=4. Deliberately did not add an ndead trigger to Python's QA pass -- doing so would flag fibers as bad that C++ itself has no mechanism to flag, for a case where Python's fit is objectively better-behaved.
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
- **Old behavior**: `rc==0` at the process level, but the bundle's output is missing/wrong -- looked like a z-band-specific correctness regression unless logs were checked by hand.
- **How it's handled (as of 2026-08-14)**: `fit_bundle_task`'s failure path now tags `is_oom` (`RESOURCE_EXHAUSTED` in the traceback), distinguishing real OOM from any other bug. `fit_ccd_native` automatically resubmits only the OOM-tagged bundles in a follow-up batch at halved worker packing (`--workers-per-gpu`/`--cpu-workers`, floor 1), capped at 2 retry rounds -- bundles are independent tasks, so no need to restart the whole CCD. Non-OOM failures are deliberately NOT retried (a real bug just reproduces itself and wastes GPU time). The CLI now always prints one unambiguous final line (`SPECEX_RESULT: OK N/N bundles` or `FAILED k/N: [ids]`) and **exits non-zero on any unrecovered failure** -- previously it always exited 0 regardless.
- **Validated on a real GPU node** with a deliberately forced extreme-oversubscription stress test: correctly recovered all 12 genuinely-OOM bundles via retry (6 at round 1, 6 more at round 2) while leaving 8 unrelated CUDA-context-corruption failures (an artifact of the artificial stress test, not a realistic scenario) alone and clearly reported rather than futilely retried.
- **Verdict**: an infrastructure/resource-sizing issue, not an algorithmic weakness. Grepping logs for `WARNING: Bundle`/`RESOURCE_EXHAUSTED` is no longer required to catch this -- the exit code and `SPECEX_RESULT` line now do it automatically -- but the gotcha is kept here since older logs/tooling may still rely on it.

## Infrastructure / orchestration edge cases

### Non-OOM crash in a single bundle or camera
- **Trigger**: a genuine bug or unexpected data condition (not GPU OOM) crashes one bundle's worker process.
- **How it's handled**: isolation is structural, not a special case that needed building. `multiprocessing.Pool.starmap` catches per-bundle exceptions internally, so one bad bundle can't take down the rest of its CCD; `run_night.py`'s one-subprocess-per-camera model means one camera can't take down others. A crashed bundle is reported (not silently dropped) via the existing `WARNING: Bundle N failed` log line, tracked in `failed_bundles`, and now (2026-08-14) surfaced via a non-zero exit code and the `SPECEX_RESULT: FAILED k/N: [ids]` summary line -- deliberately NOT retried (retrying a real bug just reproduces the same failure and wastes GPU time; only OOM gets the retry treatment, see above).
- **Verdict**: matches the ask directly -- other cameras/bundles finish, and the one that failed is now loud and unambiguous (exit code + one-line summary) rather than requiring a log grep.

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
