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

def hermite_pol_jnp(degree, x):
    if degree == 0: return jnp.ones_like(x)
    if degree == 1: return x
    h_prev2 = jnp.ones_like(x)
    h_prev = x
    h_curr = x
    for i in range(2, degree + 1):
        h_curr = x * h_prev - (i - 1) * h_prev2
        h_prev2 = h_prev; h_prev = h_curr
    return h_curr

def hermite_pol_np(degree, x):
    if degree == 0: return np.ones_like(x)
    if degree == 1: return x
    h_prev2 = np.ones_like(x)
    h_prev = x
    h_curr = x
    for i in range(2, degree + 1):
        h_curr = x * h_prev - (i - 1) * h_prev2
        h_prev2 = h_prev; h_prev = h_curr
    return h_curr

class Legendre1DPol:
    def __init__(self, deg=0, xmin=-1.0, xmax=1.0, coeff=None):
        self.deg = deg
        self.xmin = xmin
        self.xmax = xmax
        # Use numpy for storage to avoid early JAX/GPU allocation
        self.coeff = np.array(coeff) if coeff is not None else np.zeros(deg + 1)


    def monomials(self, x):
        # NumPy, not JAX: this is called from plain-Python hot loops (e.g.
        # gh_params(), once per candidate spot per selection pass) where
        # JAX's eager-mode per-op GPU dispatch overhead (~1ms/op) dominates
        # runtime for what's otherwise a handful of scalar flops. Nothing
        # here needs autodiff -- the JAX-jitted fit machinery in fitter.py
        # operates on raw jnp coefficient arrays directly and never calls
        # this method.
        rx = 2 * (x - self.xmin) / (self.xmax - self.xmin) - 1
        return np.stack([legendre_pol(i, rx) for i in range(self.deg + 1)], axis=0)

    def value(self, x):
        m = self.monomials(x)
        return np.dot(self.coeff, m)

    def derivative(self, x):
        """
        Calculates the first derivative of the Legendre polynomial at x.
        """
        rx = 2 * (x - self.xmin) / (self.xmax - self.xmin) - 1
        drx_dx = 2.0 / (self.xmax - self.xmin)
        
        # dP_n/dx = drx/dx * dP_n/drx
        # For Legendre: (1-x^2) P'_n(x) = n(P_{n-1}(x) - x P_n(x))
        # But it's easier to just use the recurrence for derivatives:
        # P'_n(x) = x P'_{n-1}(x) + n P_{n-1}(x)
        
        p_prev2 = jnp.ones_like(rx); p_prev = rx
        d_prev2 = jnp.zeros_like(rx); d_prev = jnp.ones_like(rx)
        
        derivs = [jnp.zeros_like(rx), jnp.ones_like(rx)]
        for i in range(2, self.deg + 1):
            p_curr = ((2 * i - 1) * rx * p_prev - (i - 1) * p_prev2) / i
            d_curr = d_prev2 + (2 * i - 1) * p_prev
            derivs.append(d_curr)
            p_prev2 = p_prev; p_prev = p_curr
            d_prev2 = d_prev; d_prev = d_curr
            
        d_monomials = jnp.stack(derivs[:self.deg+1], axis=0)
        return jnp.dot(self.coeff, d_monomials) * drx_dx

    def invert(self, y):
        """
        Robust inversion of the polynomial using a fine-grid lookup.
        Finds x such that value(x) == y.
        """
        # Create a grid of 1000 points over the domain
        x_grid = np.linspace(self.xmin, self.xmax, 1000)
        y_grid = self.value(x_grid)
        return np.interp(y, y_grid, x_grid)

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
