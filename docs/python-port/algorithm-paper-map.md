# Paper-to-code map: Guy et al. (2023) &rarr; specex

This maps the algorithms described in [Guy et al. 2023](../../2209.14482v2.pdf)
("The Spectroscopic Data Processing Pipeline for the Dark Energy
Spectroscopic Instrument", arXiv:2209.14482) to where they're implemented in
this repo, for both the original C++ (`src/`) and the Python/JAX port
(`py/specex/`). The paper describes the *entire* DESI spectroscopic
pipeline; only part of one subsection (&sect;4.3, "Spectrograph point spread
function") describes what `specex` itself does. See the "Out of scope"
section at the end for the paper's other pipeline stages and where (if
anywhere) they live.

For *why* the Python port's numbers occasionally differ from a literal
reading of the paper or the C++ code, `porting-notes.md` (same folder,
chronological session log) has the full investigation and evidence trail for
every deviation noted below -- this document gives the pointer, not the
argument.

## 1. Gauss-Hermite PSF model (paper &sect;4.3.1, Eq. 3)

**What it does:** the PSF at CCD position `(x,y)` for a fiber/wavelength with
center `(xc,yc)` is a 2D Gaussian times a linear combination of Hermite
polynomials (probabilist's convention) in `(x-xc)/sigma_x` and
`(y-yc)/sigma_y`, up to a configurable degree (6 in production). The
coefficients `a_i,j`, `xc`, `yc`, `sigma_x`, `sigma_y` all vary continuously
with wavelength and per fiber, modeled as Legendre polynomials of wavelength
(see section 4 below). The paper also describes an optional power-law
extended-tail term, explicitly **not used in the production pipeline**
("It is however not used in the current version of the pipeline as the
Gauss-Hermite terms appear to be sufficient").

**Python:** `psf.py:16` `GaussHermitePSF.single_pix_value_jnp()` (JAX,
differentiable, the one actually used in the fit) and its NumPy twin
`psf.py:64` `single_pix_value_np()` (dead code, see below). Both integrate
the Gaussian-times-Hermite terms analytically over a pixel (erf-based),
matching the paper's "integral of each component is analytic" property.
Both also carry the optional tail term (`params[param_index:param_index+5]`,
a `t_amp`/`t_core`/`t_xsca`/`t_ysca`/`t_inde` power-law block) -- present as
code, consistent with the paper's own description of it as implemented-but-
unused; nothing in this port's CLI or defaults populates/enables it.

**C++:** `specex_gauss_hermite_psf.cc:31` `GaussHermitePSF::Profile()` and
`specex_gauss_hermite_psf.cc:145` `GaussHermitePSF::PixValue()`.

## 2. Legendre/Hermite polynomial bases (supporting infrastructure for &sect;4.3.1)

**What it does:** the paper doesn't number this as its own algorithm, but
Eq. 3's `a_i,j`, `xc`, `yc`, `sigma_x`, `sigma_y` are each "Legendre
polynomials of the wavelength, per fiber" -- this is the shared basis-
function machinery used throughout the PSF and trace fits.

**Python:** `math.py:4` `legendre_pol()` (NumPy) / `math.py:25`
`legendre_pol_jnp()` (JAX) for the polynomial values themselves;
`math.py:45` `hermite_pol_jnp()` / `math.py:56` `hermite_pol_np()` for the
Hermite polynomials (probabilist's convention, matching the paper's
footnote 6); `math.py:67` `Legendre1DPol` (a 1D fittable/invertible Legendre
polynomial, used for per-fiber wavelength<->Y-coordinate solutions) and
`math.py:127` `SparseLegendre2DPol` (2D Legendre-in-`(x,y)` polynomials, used
for the PSF-shape parameters' spatial variation across the CCD, per the
paper's note that the fit "determines...two dimensional polynomials of the
CCD coordinates").

**C++:** `specex_legendre.cc`/`specex_legendre.h`, `specex_hermite.cc`/
`specex_hermite.h`.

## 3. PSF fit procedure (paper &sect;4.3.2)

**What it does:** a per-bundle (25-fiber block), iterative, non-linear
least-squares fit (paper says "Gauss-Newton"; in practice a Levenberg-
Marquardt-style damped normal-equations solve, see below) with analytic
derivatives, performed in stages: (i) trace coordinates `(xc,yc)` adjusted
from the input PSF, (ii) the Gaussian widths `(sigma_x,sigma_y)` fit, (iii)
the Gauss-Hermite coefficients `a_i,j` fit -- with line intensities fit
simultaneously at every stage. The fit is independent per 25-fiber bundle,
matching DESI's physical pseudo-slit layout.

**Python:** `fitter.py:1051` `class PSF_Fitter`, `fitter.py:1060`
`PSF_Fitter.fit()` is the whole staged procedure (`fit_type` controls which
of trace/psf-shape/full-joint stage runs). The linearized normal-equations
build (the Jacobian/Hessian accumulation the "Gauss-Newton" step needs) is
`fitter.py:336` `_accumulate_bundle_jax()`, with the corresponding forward
model in `fitter.py:298` `_predict_bundle_jax()`. The per-iteration step-size
search inside the loop defaults to a simple coarse grid search
(`line_search='grid'`, the production default) rather than a literal
Numerical-Recipes Brent line search -- see the "Known deviations" note
below.

**C++:** `specex_psf_fitter.cc:2293` `PSF_Fitter::FitEverything()` (top-level
staged driver) and `specex_psf_fitter.cc:1144`
`PSF_Fitter::FitSeveralSpots()` (the actual linear-system build/solve per
stage).

### Known deviations from the paper/C++ here

- **Line search:** C++'s actual step-size search inside the Gauss-Newton
  loop is a Numerical-Recipes-style Brent line search with mode-dependent
  skip logic (not literally described in the paper text, which just says
  "Gauss-Newton"). The Python port's production default (`--line-search
  grid`) replaces this with a coarse 3-point `[0.2, 0.5, 1.0]` search.
  A faithful Brent replica (`fitter.py:10` `_cpp_brent()`, used only via
  `--line-search cpp`) and a continuous-but-not-C++-faithful variant
  (`--line-search brent`) both exist in the code but are **not the default**
  and are kept for reference only -- both were tested and found
  correctness-neutral (no xrms/yrms change) vs. the grid search, so the
  simpler default was kept. See `porting-notes.md`'s line-search
  investigation entries.
- **Mixed float32/float64 precision** for the joint-fit Jacobian
  (`SPECEX_MIXED_PRECISION`, on by default) is a Python-port-only addition
  for GPU memory reasons, with no analog in the paper or C++ (which is
  entirely float64). Validated numerically equivalent (chi2 relative error
  ~2e-6) -- see `porting-notes.md` "Mixed precision, tested exactly as
  directed".

## 4. Individual-spot ("housekeeping") flux/position fits (paper &sect;4.3.2, "line intensities...fit at the same time")

**What it does:** before/alongside the joint bundle fit, each candidate
arc-line spot needs an initial flux (and sometimes position) estimate, both
to seed the joint fit and to gate spot selection by S/N.

**Python:** `fitter.py:989` `_fit_all_spots_batch()` and `fitter.py:1029`
`_get_spot_stats_jax()` are the batched (vmapped) versions actually used.
`fitter.py:963` `_fit_one_spot_jax()` is dead code -- an earlier
per-spot (non-batched) implementation, superseded by the batch versions
above but never removed; it has zero callers in the active package (only
referenced in the already-abandoned `fitter_old.py`/`fitter_old_gpt.py`
drafts).

**C++:** `specex_psf_fitter.cc:1871` `PSF_Fitter::FitIndividualSpotFluxes()`,
`specex_psf_fitter.cc:1935` `PSF_Fitter::FitIndividualSpotPositions()`,
`specex_psf_fitter.cc:1839` `PSF_Fitter::FitOneSpot()` (the direct analog of
Python's now-dead `_fit_one_spot_jax`).

## 5. Arc-line candidate generation and iterative selection (paper &sect;4.3.2, implicit)

**What it does:** for each bundle, candidate arc-line spots are generated
from the known line list (wavelength) and the input trace solution
(position), then iteratively refined/filtered by S/N and other criteria
before being handed to the joint fit. The paper doesn't describe this step
in much algorithmic detail (it's folded into "a list of arc lamp lines...
in Table 4").

**Python:** `fitter.py:719` `generate_bundle_candidates()`,
`fitter.py:796` `select_bundle_spots_iterative()`, `fitter.py:907`
`get_bundle_spots()`, `fitter.py:601` `select_spots_cpp()` (a faithful port
of C++'s own selection-threshold logic), `fitter.py:763`
`fit_candidate_fluxes()`. The line list itself is
`py/specex/data/specex_linelist_desi.txt` (and `..._vacuum.txt`), the
Python-side equivalent of the paper's Table 4 (Appendix B).

**C++:** selection/candidate logic lives inside `specex_psf_fitter.cc`
(no single dedicated file); `read_lamp_lines()` (`io.py:357` on the Python
side) has a direct C++ analog in the line-list-reading code referenced from
`specex_psf_fitter.cc`.

## 6. Production robustness additions (not described in the paper)

These handle real DESI operational data-quality issues the paper doesn't
discuss (it's a methods/performance paper, not an ops-hardening one):

- **Dead-column masking**: `fitter.py:526` `apply_dead_column_mask()` and
  `fitter.py:541` `filter_dead_column_spots()` (the latter's docstring
  already cites its C++ source directly: `specex_psf_fitter.cc:244-258`,
  `InitTmpData`'s `can_measure_flux`/ignore exclusion) -- this one *is* a
  faithful port of existing C++ behavior, just not paper-documented.
- **Masked/dead CCD amplifier detection**: `fitter.py:187`
  `find_masked_amp_fibers()` -- a Python-port-only addition (2026-08-14,
  per `porting-notes.md`), no C++ equivalent. Detects fibers with no real
  data (e.g. a masked bad amp) via a contiguous-run dead-pixel heuristic and
  excludes them from the fit entirely rather than trying to fit noise.
- **Post-fit trace-crossing QA**: the inline pass inside `io.py`'s
  `write_python_psf()` (`io.py:70`, see the large comment around line 192)
  flags neighboring fibers whose fitted traces cross. Adapted from a later
  upstream specex GitHub PR (#91), not from this paper or the original
  C++ baseline this port targets. `py/specex/qa.py`'s `trace_psf_qa()`/
  `specex_psf_qa()` are the *original* PR #91-style functions this was
  adapted from, but they are themselves dead code in this repo -- imported
  by the C++-wrapper path (`run_specex_cpp()`, `specex.py`) but never
  actually called there (the call site is commented out).
- **Trace-consensus prior**: C++ *does* have a trace-prior mechanism
  (`specex_psf_fitter.cc:759-857`, gated by `--trace-prior-deg`, penalizing
  each fiber's high-order trace coefficients for deviating from the
  bundle's cross-fiber mean) -- but it's not mentioned in the paper text,
  and real DESI production runs C++ with it off (`trace_prior_deg=0`).
  The Python port turns an equivalent prior **on by default**
  (`fitter.py:1060`'s `trace_prior_deg`/`trace_prior_weight` params,
  `fitter.py:242` `build_trace_prior_hessian()`), with two real
  corrections found during porting: the penalty weight was dropped from
  C++'s hardcoded `1e8` (found to measurably over-smooth even the fibers
  it's meant to help, per a weight sweep in `porting-notes.md`) to `1e5`,
  and it's gated to only activate for fibers with a high dead-pixel count
  (`ndead > 500` by default) rather than applying blanket-style to every
  fiber, which was found to harm otherwise-healthy bundles.

## Out of scope (paper sections this repo does not implement)

- **&sect;4.4.1, Initial wavelength calibration and trace coordinates fit**
  (the triplet/histogram peak-matching algorithm, inspired by
  Valdes et al. 1995) -- a one-time procedure for calibrating brand-new
  spectrograph hardware. Not present anywhere in this repo (`src/` or
  `py/specex/`); if it exists in the DESI software stack it's a separate
  tool.
- **&sect;4.4.2, Per-exposure adjustment of wavelength calibration and trace
  coordinates** (the `deltaX`/`deltaY` fit that produces DESI's
  `shifted-input-psf-*.fits` files) -- this port's `--in-psf` CLI flag
  *consumes* an already-shifted PSF file as input (see `how-to-run.md`) but
  does not compute the shift itself; that's a separate upstream step.
- **&sect;4.5, Spectral extraction** -- the actual 1D-spectrum forward-model
  extraction that *uses* the PSF this repo fits. A separate codebase
  (`specter`/`desispec`), not `specex`.
- **&sect;4.4.3/4.4.4** (barycentric velocity correction, radial velocity
  measurement), **&sect;4.6-4.14** (fiber flat fielding, sky subtraction,
  stellar model fit, flux calibration, cross-talk correction, co-addition,
  classification/redshift fitting, redshift performance, effective exposure
  time) -- all downstream pipeline stages unrelated to PSF fitting.
- **Eq. 4 / Appendix A, the PSF-stability flux-bias formula** (`delta F/F`)
  -- an analysis formula the paper uses to characterize/validate PSF
  stability over time, not an algorithm either codebase implements as part
  of normal operation.
