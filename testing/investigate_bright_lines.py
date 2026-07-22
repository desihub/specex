"""
Task 21: for the handful of lamp lines where Python's final selection keeps a
"ghost" spot that C++ drops (z8 bundle 5: 8670.325, 9354.8, 7490.9335,
7440.9469, 9356.787, 7516.721 -- see ghost_spots_analysis.txt), compare:
  (a) each pipeline's individual-spot GH-model flux fit (pass-1, pre-selection,
      un-refit candidate positions -- cpp_cp0_pass1.txt vs pyrawspots.txt)
  (b) an independent, model-free raw aperture sum of the actual preproc pixel
      data at that candidate position, to see which pipeline's fitted flux is
      actually closer to "the light that's really there."

Usage:
  python testing/investigate_bright_lines.py
"""
import os
import sys
import numpy as np
import fitsio

CURRENT_DIR = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
OUTDIR = "/pscratch/sd/c/cdwarner/specex/testing/multi"
PREPROC = "/global/cfs/cdirs/desi/spectro/redux/matterhorn/preproc/20260401/00344649/preproc-z8-00344649.fits.gz"

PROBLEM_LINES = [8670.325, 9354.8, 7490.9335, 7440.9469, 9356.787, 7516.721]
WAVE_TOL = 0.02
HSIZE = 3  # matches housekeeping stamp half-size: min(3, h_size_x/y)


def load(path, ncols):
    rows = []
    with open(path) as f:
        for line in f:
            p = line.strip().split(',')
            if len(p) < ncols:
                continue
            rows.append([float(x) for x in p[:ncols]])
    return np.array(rows)


def raw_aperture_flux(data, ivar, xc, yc, hsize=HSIZE):
    """Model-free: sum data within a (2*hsize+1)^2 box, minus a flat local
    background estimated from the box's outer ring, masking bad pixels
    (ivar<=0). Returns (flux, npix_used, bkg_per_pix)."""
    ny, nx = data.shape
    ix, iy = int(np.floor(xc + 0.5)), int(np.floor(yc + 0.5))
    imin, imax = ix - hsize, ix + hsize + 1
    jmin, jmax = iy - hsize, iy + hsize + 1
    if imin < 0 or jmin < 0 or imax > nx or jmax > ny:
        return np.nan, 0, np.nan
    stamp = data[jmin:jmax, imin:imax]
    ivstamp = ivar[jmin:jmax, imin:imax]
    good = ivstamp > 0
    # outer ring = box edge pixels, used for a simple local background estimate
    ring_mask = np.zeros_like(stamp, dtype=bool)
    ring_mask[0, :] = ring_mask[-1, :] = True
    ring_mask[:, 0] = ring_mask[:, -1] = True
    ring_good = ring_mask & good
    bkg = np.median(stamp[ring_good]) if ring_good.sum() > 0 else 0.0
    total = np.sum((stamp - bkg)[good])
    return float(total), int(good.sum()), float(bkg)


def main():
    cpp = load(os.path.join(OUTDIR, "cpp-z8-00344649_05.cpp_cp0_pass1.txt"), 6)
    py = load(os.path.join(OUTDIR, "py-z8-00344649_05.pyrawspots.txt"), 7)

    print(f"Loading preproc image: {PREPROC}")
    fx = fitsio.FITS(PREPROC)
    data = fx['IMAGE'].read().astype(np.float64)
    ivar = fx['IVAR'].read().astype(np.float64)
    print(f"  image shape: {data.shape}")

    hdr = (f"{'fiber':>5} {'wave':>10} {'cpp_flux':>10} {'py_flux':>10} {'ratio':>7} "
           f"{'raw_cpp':>10} {'raw_py':>10} {'cpp_x':>7} {'cpp_y':>7} {'py_x':>7} {'py_y':>7}")
    print(hdr)
    print("-" * len(hdr))

    for line in PROBLEM_LINES:
        cmask = np.abs(cpp[:, 1] - line) < WAVE_TOL
        pmask = np.abs(py[:, 1] - line) < WAVE_TOL
        cpp_rows = {int(r[0]): r for r in cpp[cmask]}
        py_rows = {int(r[0]): r for r in py[pmask]}
        fibers = sorted(set(cpp_rows) & set(py_rows))
        for fib in fibers:
            cr, pr = cpp_rows[fib], py_rows[fib]
            cpp_x, cpp_y, cpp_flux = cr[2], cr[3], cr[4]
            py_x, py_y, py_flux = pr[2], pr[3], pr[4]
            raw_cpp, _, _ = raw_aperture_flux(data, ivar, cpp_x, cpp_y)
            raw_py, _, _ = raw_aperture_flux(data, ivar, py_x, py_y)
            ratio = py_flux / cpp_flux if cpp_flux != 0 else float('nan')
            row = (f"{fib:>5} {line:>10.4f} {cpp_flux:>10.2f} {py_flux:>10.2f} {ratio:>7.3f} "
                   f"{raw_cpp:>10.2f} {raw_py:>10.2f} {cpp_x:>7.2f} {cpp_y:>7.2f} {py_x:>7.2f} {py_y:>7.2f}")
            print(row)


if __name__ == "__main__":
    main()
