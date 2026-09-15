import numpy as np
import sys
import os

# Ensure we can import from py/specex
current_dir = os.getcwd()
sys.path.insert(0, os.path.join(current_dir, 'py'))

from specex.psf import GaussHermitePSF

def test_vectorized_pix_value():
    print("Testing vectorized GaussHermitePSF.pix_value...")
    degree = 6
    psf = GaussHermitePSF(degree)
    
    n_coeff = (degree + 1)**2 - 1
    
    # Test 1: Scalar parity
    params_scalar = np.random.rand(2 + n_coeff)
    params_scalar[0:2] += 0.5 # ensure sx, sy > 0.1
    xc, yc = 10.3, 10.4
    xpix, ypix = 10.0, 10.0
    
    val_scalar = psf.pix_value(xc, yc, xpix, ypix, params_scalar)
    assert isinstance(val_scalar, float)
    print(f"Scalar check passed: {val_scalar}")

    # Test 2: Multiple pixels, single spot -- one value per pixel, shape (n_pix,).
    xpix_arr = np.linspace(5, 15, 5)
    ypix_arr = np.linspace(5, 15, 5)

    val_pixels = psf.pix_value(xc, yc, xpix_arr, ypix_arr, params_scalar)
    assert val_pixels.shape == (5,)

    for k in range(5):
        v = psf.pix_value(xc, yc, xpix_arr[k], ypix_arr[k], params_scalar)
        np.testing.assert_allclose(val_pixels[k], v)
    print("Multiple pixels check passed!")

    # Test 3: Multiple spots, single pixel -- one value per spot, shape (n_spots,).
    xc_arr = np.random.rand(3) + 10
    yc_arr = np.random.rand(3) + 10
    params_arr = np.random.rand(3, 2 + n_coeff)
    params_arr[:, 0:2] += 0.5

    val_spots = psf.pix_value(xc_arr, yc_arr, xpix, ypix, params_arr)
    assert val_spots.shape == (3,)

    for k in range(3):
        v = psf.pix_value(xc_arr[k], yc_arr[k], xpix, ypix, params_arr[k])
        np.testing.assert_allclose(val_spots[k], v)
    print("Multiple spots check passed!")

    # Test 4: Multiple pixels AND multiple spots -- shape (n_spots, n_pix),
    # matching pix_value_jnp's own documented convention.
    val_all = psf.pix_value(xc_arr, yc_arr, xpix_arr, ypix_arr, params_arr)
    assert val_all.shape == (3, 5)

    for j in range(3):
        for i in range(5):
            v = psf.pix_value(xc_arr[j], yc_arr[j], xpix_arr[i], ypix_arr[i], params_arr[j])
            np.testing.assert_allclose(val_all[j, i], v)
    print("Full matrix check passed!")

if __name__ == "__main__":
    test_vectorized_pix_value()
    print("All vectorization tests passed!")
