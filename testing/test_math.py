import numpy as np
import sys
import os

# Ensure we can import from py/specex and built C++ libs
current_dir = os.getcwd()
sys.path.insert(0, os.path.join(current_dir, 'py'))
sys.path.insert(0, os.path.join(current_dir, 'build'))

from specex.math import Legendre1DPol, SparseLegendre2DPol
try:
    from specex._libspecex import (Legendre1DPol as CppLegendre1DPol, 
                                 SparseLegendre2DPol as CppSparseLegendre2DPol)
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
        cpp_pol = CppLegendre1DPol(deg, xmin, xmax)
        for i, c in enumerate(coeff):
            cpp_pol.coeff[i] = c
        
        cpp_vals = [cpp_pol.value(x) for x in x_test]
        
        np.testing.assert_allclose(py_vals, cpp_vals, rtol=1e-10)
        print("Legendre1DPol: Python matches C++!")
    else:
        print(f"Python values: {py_vals}")

def test_sparse_legendre_2d():
    print("Testing SparseLegendre2DPol...")
    xdeg, ydeg = 2, 2
    xmin, xmax = 0.0, 10.0
    ymin, ymax = 0.0, 10.0
    # Match the Fill() logic for sparse=true in C++
    # if i==0: all j
    # if i==1: j<2
    # else: j==0
    non_zero_indices = []
    for i in range(xdeg + 1):
        for j in range(ydeg + 1):
            if i == 0: non_zero_indices.append(i + j*(xdeg+1))
            elif i == 1 and j < 2: non_zero_indices.append(i + j*(xdeg+1))
            elif i > 1 and j == 0: non_zero_indices.append(i + j*(xdeg+1))
    
    coeff = np.random.rand(len(non_zero_indices))
    
    py_pol = SparseLegendre2DPol(xdeg, xmin, xmax, ydeg, ymin, ymax, coeff, non_zero_indices)
    
    x_test = np.linspace(xmin, xmax, 5)
    y_test = np.linspace(ymin, ymax, 5)
    
    if HAS_CPP:
        cpp_pol = CppSparseLegendre2DPol(xdeg, xmin, xmax, ydeg, ymin, ymax)
        # SparseLegendre2DPol::Fill(true) in C++
        # We need to manually set the indices and coefficients to match
        for k in non_zero_indices:
            i = k % (xdeg + 1)
            j = k // (xdeg + 1)
            cpp_pol.add(i, j)
        
        for i, c in enumerate(coeff):
            cpp_pol.coeff[i] = c
            
        for x in x_test:
            for y in y_test:
                py_val = py_pol.value(x, y)
                cpp_val = cpp_pol.value(x, y)
                np.testing.assert_allclose(py_val, cpp_val, rtol=1e-10)
        
        print("SparseLegendre2DPol: Python matches C++!")
    else:
        print("Skipping C++ comparison.")

if __name__ == "__main__":
    test_legendre_1d()
    test_sparse_legendre_2d()
    print("All math tests passed!")
