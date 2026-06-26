import os
import numpy as np
from specex.io import load_python_psf
from specex.psf import PSF

class DummyOpts:
    def __init__(self):
        self.input_psf_filename = "/dvs_ro/cfs/cdirs/desi/spectro/redux/matterhorn/exposures/20260401/00344649/shifted-input-psf-z8-00344649.fits"

opts = DummyOpts()
psf = load_python_psf(opts.input_psf_filename, opts)

bundle_id = 5
bundle = psf.params_of_bundles[bundle_id]
print(f"Bundle {bundle_id} fibers: {bundle.fiber_min} - {bundle.fiber_max}")
print(f"Param names: {bundle.param_names}")

# Check GHSIGX for middle fiber
fib = 137
wave = 8000.0
rel_idx = fib - psf.fiber_min
sigx = bundle.param_models['GHSIGX'][rel_idx].value(wave)
sigy = bundle.param_models['GHSIGY'][rel_idx].value(wave)
print(f"Fiber {fib} at {wave}A: GHSIGX={sigx:.4f} GHSIGY={sigy:.4f}")
