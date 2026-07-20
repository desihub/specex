"""
Convert between (fiber, wavelength) and (x, y) CCD pixel position using the
XTRACE/YTRACE Legendre trace model stored in a specex output PSF FITS file.

How the trace is stored (both C++ and Python outputs, same convention):
  - HDU 'XTRACE' and HDU 'YTRACE' are each a (n_fibers, ncoeff) float array.
    Row `fib` of XTRACE holds the Legendre coefficients of X_ccd(wave) for
    fiber `fib`; row `fib` of YTRACE holds the coefficients of Y_ccd(wave).
    `ncoeff = LEGDEG+1` where LEGDEG is the polynomial degree used for the
    *trace* (this is a different, usually higher, degree than the small
    within-bundle joint-fit correction degree discussed elsewhere in
    porting-notes.md -- the trace itself, written per-fiber, is always a
    single 1D Legendre-in-wavelength curve regardless of band).
  - Each HDU's header carries WAVEMIN/WAVEMAX: the wavelength domain the
    Legendre basis is normalized over (rx = 2*(wave-WAVEMIN)/(WAVEMAX-WAVEMIN)-1).
    Evaluate with x(wave) = sum_i coeff[i] * P_i(rx) (standard Legendre polys).
  - Fiber indexing is absolute (0-499 for a full CCD), matching the FITS row
    index directly -- no bundle-relative offset.
  - X is the cross-dispersion axis, Y is the dispersion (wavelength) axis.
    Only Y_vs_W is monotonic in wavelength, so wavelength inversion (given a
    measured y) uses YTRACE only; XTRACE is not invertible for wavelength
    (a given x can repeat at multiple wavelengths depending on trace shape).

Usage:
  python testing/wave_xy_convert.py fiber_wave_to_xy  <psf.fits> <fiber> <wave> [<wave2> ...]
  python testing/wave_xy_convert.py xy_to_wave         <psf.fits> <fiber> <y> [<y2> ...]

Or import and use directly:
  from wave_xy_convert import PSFTrace
  t = PSFTrace('fit-psf-z8-00344649.fits')
  x, y = t.wave_to_xy(fiber=130, wave=8670.33)      # scalar or np.array wave
  wave = t.y_to_wave(fiber=130, y=2015.46)          # scalar or np.array y
"""
import sys
import os
import numpy as np
import fitsio

sys.path.insert(0, os.path.join(os.path.dirname(os.path.dirname(os.path.abspath(__file__))), 'py'))
from specex.math import Legendre1DPol


class PSFTrace:
    def __init__(self, fits_path):
        f = fitsio.FITS(fits_path)
        xt_hdr = f['XTRACE'].read_header()
        yt_hdr = f['YTRACE'].read_header()
        xtrace = f['XTRACE'].read().astype(np.float64)
        ytrace = f['YTRACE'].read().astype(np.float64)
        self.wavemin, self.wavemax = float(xt_hdr['WAVEMIN']), float(xt_hdr['WAVEMAX'])
        assert self.wavemin == float(yt_hdr['WAVEMIN']) and self.wavemax == float(yt_hdr['WAVEMAX']), \
            "XTRACE/YTRACE WAVEMIN/WAVEMAX disagree -- unexpected file"
        self.n_fibers = xtrace.shape[0]
        self.x_pol = {
            fib: Legendre1DPol(deg=xtrace.shape[1] - 1, xmin=self.wavemin, xmax=self.wavemax, coeff=xtrace[fib])
            for fib in range(self.n_fibers)
        }
        self.y_pol = {
            fib: Legendre1DPol(deg=ytrace.shape[1] - 1, xmin=self.wavemin, xmax=self.wavemax, coeff=ytrace[fib])
            for fib in range(self.n_fibers)
        }

    def wave_to_xy(self, fiber, wave):
        """fiber, wave (scalar or np.array) -> (x, y) CCD pixel position."""
        scalar_in = np.isscalar(wave)
        wave = np.atleast_1d(np.asarray(wave, dtype=np.float64))
        x = self.x_pol[fiber].value(wave)
        y = self.y_pol[fiber].value(wave)
        return (float(x[0]), float(y[0])) if scalar_in else (x, y)

    def y_to_wave(self, fiber, y):
        """fiber, measured y pixel (scalar or np.array) -> wavelength, via
        fine-grid inversion of Y_vs_W (Y is the dispersion axis, monotonic)."""
        scalar_in = np.isscalar(y)
        y = np.atleast_1d(np.asarray(y, dtype=np.float64))
        wave = self.y_pol[fiber].invert(y)
        return float(wave[0]) if scalar_in else wave


def main():
    if len(sys.argv) < 5:
        print(__doc__)
        sys.exit(1)
    mode, fits_path, fiber = sys.argv[1], sys.argv[2], int(sys.argv[3])
    vals = [float(v) for v in sys.argv[4:]]
    t = PSFTrace(fits_path)
    if mode == 'fiber_wave_to_xy':
        for w in vals:
            x, y = t.wave_to_xy(fiber, w)
            print(f"fiber={fiber} wave={w:.4f}  ->  x={x:.4f} y={y:.4f}")
    elif mode == 'xy_to_wave':
        for y in vals:
            w = t.y_to_wave(fiber, y)
            print(f"fiber={fiber} y={y:.4f}  ->  wave={w:.4f}")
    else:
        print(f"unknown mode {mode!r}; use fiber_wave_to_xy or xy_to_wave")
        sys.exit(1)


if __name__ == '__main__':
    main()
