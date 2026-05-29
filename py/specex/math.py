import numpy as np
import jax.numpy as jnp

def legendre_pol(degree, x):
    """
    NumPy version.
    """
    if degree == 0: return np.ones_like(x)
    if degree == 1: return x
    if degree == 2: return 0.5 * (3 * x**2 - 1.)
    if degree == 3: return 0.5 * (5 * x**2 - 3) * x
    if degree == 4: return 0.125 * ((35 * x**2 - 30) * x**2 + 3)
    if degree == 5: return 0.125 * (((63 * x**2 - 70) * x**2 + 15) * x)
    if degree == 6: return (((231 * x**2 - 315) * x**2 + 105) * x**2 - 5) / 16.
    
    p_prev2 = np.ones_like(x)
    p_prev = x
    p_curr = x
    for i in range(2, degree + 1):
        p_curr = ((2 * i - 1) * x * p_prev - (i - 1) * p_prev2) / i
        p_prev2 = p_prev
        p_prev = p_curr
    return p_curr

def legendre_pol_jnp(degree, x):
    """
    JAX version.
    """
    if degree == 0: return jnp.ones_like(x)
    if degree == 1: return x
    if degree == 2: return 0.5 * (3 * x**2 - 1.)
    if degree == 3: return 0.5 * (5 * x**2 - 3) * x
    if degree == 4: return 0.125 * ((35 * x**2 - 30) * x**2 + 3)
    if degree == 5: return 0.125 * (((63 * x**2 - 70) * x**2 + 15) * x)
    
    p_prev2 = jnp.ones_like(x)
    p_prev = x
    p_curr = x
    for i in range(2, degree + 1):
        p_curr = ((2 * i - 1) * x * p_prev - (i - 1) * p_prev2) / i
        p_prev2 = p_prev
        p_prev = p_curr
    return p_curr

class Legendre1DPol:
    def __init__(self, deg=0, xmin=0.0, xmax=0.0, coeff=None):
        self.deg = deg
        self.xmin = xmin
        self.xmax = xmax
        self.coeff = jnp.array(coeff) if coeff is not None else jnp.zeros(deg + 1)

    def monomials(self, x):
        rx = 2 * (x - self.xmin) / (self.xmax - self.xmin) - 1
        return jnp.stack([legendre_pol_jnp(i, rx) for i in range(self.deg + 1)], axis=0)

    def value(self, x):
        m = self.monomials(x)
        return jnp.dot(self.coeff, m)

class SparseLegendre2DPol:
    def __init__(self, xdeg=0, xmin=0.0, xmax=1.0, ydeg=0, ymin=0.0, ymax=1.0, coeff=None, non_zero_indices=None):
        self.xdeg = xdeg
        self.xmin = xmin
        self.xmax = xmax
        self.ydeg = ydeg
        self.ymin = ymin
        self.ymax = ymax
        self.non_zero_indices = non_zero_indices if non_zero_indices is not None else []
        self.coeff = jnp.array(coeff) if coeff is not None else jnp.zeros(len(self.non_zero_indices))

    def monomials(self, x, y):
        rx = 2 * (x - self.xmin) / (self.xmax - self.xmin) - 1
        ry = 2 * (y - self.ymin) / (self.ymax - self.ymin) - 1
        
        mx = [legendre_pol_jnp(i, rx) for i in range(self.xdeg + 1)]
        my = [legendre_pol_jnp(j, ry) for j in range(self.ydeg + 1)]
        
        m = []
        for k in self.non_zero_indices:
            i = k % (self.xdeg + 1)
            j = k // (self.xdeg + 1)
            m.append(mx[i] * my[j])
        return jnp.stack(m, axis=0)

    def value(self, x, y):
        m = self.monomials(x, y)
        return jnp.dot(self.coeff, m)
