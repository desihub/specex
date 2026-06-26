import os
from specex.specex import fit_ccd_native

# Setup environment
os.environ["PYTHONPATH"] = "/global/cfs/cdirs/desicollab/users/cdwarner/code/specex/py"

# Params from user
arc_file = "/dvs_ro/cfs/cdirs/desi/spectro/redux/matterhorn/preproc/20260401/00344649/preproc-z8-00344649.fits.gz"
in_psf_file = "/dvs_ro/cfs/cdirs/desi/spectro/redux/matterhorn/exposures/20260401/00344649/shifted-input-psf-z8-00344649.fits"
out_psf_file = "/pscratch/sd/c/cdwarner/specex/pyfit-psf-z8-00344649_05_test.fits"
lamp_lines_file = "/global/cfs/cdirs/desi/users/cdwarner/code/specex/py/specex/data/specex_linelist_desi.txt"

# Bundle 5 only
first_bundle = 5
last_bundle = 5
broken_fibers = [473, 474]

print(f"Running fit_ccd_native for bundle {first_bundle}...")
try:
    fit_ccd_native(
        arc_file=arc_file,
        in_psf_file=in_psf_file,
        out_psf_file=out_psf_file,
        lamp_lines_file=lamp_lines_file,
        first_bundle=first_bundle,
        last_bundle=last_bundle,
        n_gpus=1,
        backend="gpu",
        broken_fibers=broken_fibers
    )
    print("Fit completed successfully.")
except Exception as e:
    import traceback
    traceback.print_exc()
    print(f"Fit failed: {e}")
