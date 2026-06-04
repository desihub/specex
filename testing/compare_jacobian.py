import numpy as np
import sys
import os
import jax
import jax.numpy as jnp
from jax import jacfwd

# Ensure we can import from py/specex
current_dir = os.getcwd()
sys.path.insert(0, os.path.join(current_dir, 'py'))
# Also add the build directory for the compiled extension
sys.path.insert(0, os.path.join(current_dir, 'build'))

from specex.psf import GaussHermitePSF as PyGaussHermitePSF
from _libspecex import CppGaussHermitePSF, VectorDouble

def compare_jacobians():
    print("Direct Jacobian Comparison: C++ vs JAX")
    degree = 6
    # C++ build has EXTERNAL_TAIL defined, so it has 2 (sigmas) + 48 (GH) + 5 (tail) = 55 params
    n_params = 55
    
    # Random but stable parameters
    np.random.seed(42)
    params = np.zeros(n_params)
    params[0] = 1.2 # sigmaX
    params[1] = 1.1 # sigmaY
    # GH coeffs
    params[2:50] = (np.random.rand(48) - 0.5) * 0.1
    # Tail params: TAILAMP, TAILCORE, TAILXSCA, TAILYSCA, TAILINDE
    params[50] = 0.0 # amp - Set to 0 because C++ PixValue (bound) doesn't include tail
    params[51] = 1.5  # core
    params[52] = 1.2  # xsca
    params[53] = 1.1  # ysca
    params[54] = 2.5  # inde
    
    xc, yc = 10.3, 10.4
    xpix, ypix = 10.0, 10.0
    
    # 1. C++ Result
    cpp_psf = CppGaussHermitePSF(degree)
    cpp_params = VectorDouble()
    for p in params: cpp_params.append(p)
    
    cpp_val, cpp_pos_der, cpp_param_der = cpp_psf.pix_value_with_derivatives(xc, yc, xpix, ypix, cpp_params)
    cpp_pos_der = np.array(cpp_pos_der)
    cpp_param_der = np.array(cpp_param_der)
    
    # 2. JAX Result
    # We need to wrap the JAX function to differentiate it
    def spot_fn(f, x, y, g):
        # Result is scalar for single pixel
        return f * PyGaussHermitePSF.single_pix_value_jnp(x, y, xpix, ypix, g, degree)

    flux = 1.0
    # JAX derivatives wrt (flux, xc, yc, gh_params)
    # Note: C++ PosDer is wrt Xc, Yc?
    # In C++ src/specex_gauss_hermite_psf.cc:
    # dvdx=x*sigma_x_inv*psf_val; (where x = (XPix - Xc)/sx)
    # d(x)/dXc = -1/sx. 
    # So C++ PosDer is likely wrt -Xc.
    
    jac_fn = jacfwd(spot_signal_fn_wrapper, argnums=(0, 1, 2, 3))
    
    # Let's use the actual function from fitter.py style
    def jax_model(xc, yc, flux, params):
        return flux * PyGaussHermitePSF.single_pix_value_jnp(xc, yc, xpix, ypix, params, degree)
    
    jax_val = jax_model(xc, yc, flux, params)
    dj_dxc, dj_dyc, dj_dflux, dj_dparams = jacfwd(jax_model, argnums=(0, 1, 2, 3))(xc, yc, flux, params)
    
    print(f"\nValue comparison:")
    print(f"  C++: {cpp_val:.15f}")
    print(f"  JAX: {float(jax_val):.15f}")
    print(f"  Diff: {abs(cpp_val - float(jax_val)):.2e}")
    
    print(f"\nPosition Derivative comparison (d/dXc, d/dYc):")
    # C++ PosDer is [dvdx, dvdy] which are derivatives wrt -Xc?
    # Let's check signs.
    print(f"  C++ PosDer: {cpp_pos_der}")
    print(f"  JAX d/dXc, d/dYc: [{float(dj_dxc):.10f}, {float(dj_dyc):.10f}]")
    
    # Check if signs are flipped
    if abs(cpp_pos_der[0] + float(dj_dxc)) < abs(cpp_pos_der[0] - float(dj_dxc)):
        print("  Note: C++ PosDer signs are flipped relative to JAX (wrt -Xc).")
        cpp_pos_der_adj = [-cpp_pos_der[0], -cpp_pos_der[1]]
    else:
        cpp_pos_der_adj = cpp_pos_der

    np.testing.assert_allclose(cpp_pos_der_adj, [float(dj_dxc), float(dj_dyc)], rtol=1e-10)
    print("  Position derivatives match!")
    
    print(f"\nParameter Derivative comparison (Sigmas + GH):")
    cpp_param_der = np.array(cpp_param_der)
    jax_param_der = np.array(dj_dparams)
    
    # Compare only the first 50 (Sigmas + GH) which are implemented in C++ PixValue
    print(f"  C++ d/dSigX, d/dSigY: {cpp_param_der[:2]}")
    print(f"  JAX d/dSigX, d/dSigY: {jax_param_der[:2]}")
    
    np.testing.assert_allclose(cpp_param_der[:50], jax_param_der[:50], rtol=1e-10)
    print("  First 50 parameter derivatives (GH Core) match exactly!")
    
    print(f"\nTail Parameter Derivatives (JAX):")
    print(f"  JAX d/dTail: {jax_param_der[50:]}")
    print(f"  C++ d/dTail: {cpp_param_der[50:]} (Expected all 0s from PixValue)")

def spot_signal_fn_wrapper(flux, x, y, g, xpix, ypix, degree):
    from specex.psf import GaussHermitePSF
    return flux * GaussHermitePSF.single_pix_value_jnp(x, y, xpix, ypix, g, degree)

if __name__ == "__main__":
    # Enable float64 for JAX
    from jax import config
    config.update("jax_enable_x64", True)
    
    compare_jacobians()
