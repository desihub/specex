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

    # Test 2: Multiple pixels, single spot
    xpix_arr = np.linspace(5, 15, 5)
    ypix_arr = np.linspace(5, 15, 5)
    
    # This should return (5, 1) if we pass xpix_arr, ypix_arr
    # Wait, our current pix_value expects xpix and ypix to be same shape if they are arrays,
    # and it does xpix[:, np.newaxis] which creates (Np, 1) and then broadcasts with xc (Ns,)
    # If xc is scalar, it's atleast_1d -> (1,)
    # So (Np, 1) - (1,) -> (Np, 1). Correct.
    val_pixels = psf.pix_value(xc, yc, xpix_arr, ypix_arr, params_scalar)
    assert val_pixels.shape == (5, 1)
    
    for k in range(5):
        v = psf.pix_value(xc, yc, xpix_arr[k], ypix_arr[k], params_scalar)
        np.testing.assert_allclose(val_pixels[k, 0], v)
    print("Multiple pixels check passed!")

    # Test 3: Multiple spots, single pixel
    xc_arr = np.random.rand(3) + 10
    yc_arr = np.random.rand(3) + 10
    params_arr = np.random.rand(3, 2 + n_coeff)
    params_arr[:, 0:2] += 0.5
    
    val_spots = psf.pix_value(xc_arr, yc_arr, xpix, ypix, params_arr)
    assert val_spots.shape == (1, 3)
    
    for k in range(3):
        v = psf.pix_value(xc_arr[k], yc_arr[k], xpix, ypix, params_arr[k])
        np.testing.assert_allclose(val_spots[0, k], v)
    print("Multiple spots check passed!")

    # Test 4: Multiple pixels AND multiple spots
    val_all = psf.pix_value(xc_arr, yc_arr, xpix_arr, ypix_arr, params_arr)
    assert val_all.shape == (5, 3)
    
    for i in range(5):
        for j in range(3):
            v = psf.pix_value(xc_arr[j], yc_arr[j], xpix_arr[i], ypix_arr[i], params_arr[j])
            np.testing.assert_allclose(val_all[i, j], v)
    print("Full matrix check passed!")

if __name__ == "__main__":
    test_vectorized_pix_value()
    print("All vectorization tests passed!")
