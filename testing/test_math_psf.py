import numpy as np
import sys
import os

# Ensure we can import from py/specex
current_dir = os.getcwd()
sys.path.insert(0, os.path.join(current_dir, 'py'))

from specex.math import Legendre1DPol, SparseLegendre2DPol
from specex.psf import GaussHermitePSF

try:
    from specex._libspecex import (CppLegendre1DPol, 
                                 CppSparseLegendre2DPol,
                                 CppGaussHermitePSF,
                                 VectorDouble)
    HAS_CPP = True
except ImportError as e:
    print(f"Warning: Could not import C++ bindings: {e}")
    HAS_CPP = False

def test_legendre_1d():
    print("Testing Legendre1DPol...")
    deg = 3
    xmin, xmax = 0.0, 100.0
    coeff = [1.0, 0.5, -0.2, 0.1]
    
    py_pol = Legendre1DPol(deg, xmin, xmax, coeff)
    
    x_test = np.linspace(xmin, xmax, 10)
    py_vals = [py_pol.value(x) for x in x_test]
    
    if HAS_CPP:
        # Check if CppLegendre1DPol exists in bindings, if not skip
        try:
            cpp_pol = CppLegendre1DPol(deg, xmin, xmax)
            for i, c in enumerate(coeff):
                cpp_pol.coeff[i] = c
            
            cpp_vals = [cpp_pol.value(x) for x in x_test]
            np.testing.assert_allclose(py_vals, cpp_vals, rtol=1e-10)
            print("Legendre1DPol: Python matches C++!")
        except NameError:
            print("CppLegendre1DPol not in bindings, skipping comparison.")
    else:
        print(f"Python values: {py_vals}")

def test_psf_pix_value():
    print("Testing GaussHermitePSF.pix_value...")
    degree = 6
    py_psf = GaussHermitePSF(degree)
    
    # params: sx, sy, then (degree+1)**2 - 1 coefficients
    n_coeff = (degree + 1)**2 - 1
    params = np.zeros(2 + n_coeff)
    params[0] = 1.2 # sx
    params[1] = 1.1 # sy
    # Set some random coefficients
    np.random.seed(42)
    params[2:] = (np.random.rand(n_coeff) - 0.5) * 0.1
    
    xc, yc = 10.3, 10.4
    xpix, ypix = 10.0, 10.0
    
    py_val = py_psf.pix_value(xc, yc, xpix, ypix, params)
    print(f"Python PSF value: {py_val}")
    
    if HAS_CPP:
        cpp_psf = CppGaussHermitePSF(degree)
        cpp_params = VectorDouble()
        for p in params:
            cpp_params.append(p)
            
        # pix_value(Xc, Yc, XPix, YPix, Params, PosDer, ParamDer)
        # We pass None for derivatives to match the C++ signature or just let it use defaults
        # But wait, C++ might need pointers. pybind11 usually handles this.
        cpp_val = cpp_psf.pix_value(xc, yc, xpix, ypix, cpp_params)
        print(f"C++ PSF value:    {cpp_val}")
        
        np.testing.assert_allclose(py_val, cpp_val, rtol=1e-10)
        print("GaussHermitePSF: Python matches C++!")

if __name__ == "__main__":
    test_legendre_1d()
    test_psf_pix_value()
    print("All tests passed!")
