import numpy as np
import jax.numpy as jnp

def legendre_pol(degree, x):
    """Evaluate the Legendre polynomial of the given degree at x, NumPy-backed.

    Closed-form expressions are used for degrees 0-6 (avoids recurrence
    overhead in hot loops); higher degrees fall back to the standard
    three-term recurrence.

    Args:
        degree (int): polynomial degree, >= 0.
        x (float or np.ndarray): evaluation point(s), expected in [-1, 1].

    Returns:
        Same type/shape as x: P_degree(x).

    Status: ACTIVE (production default path) -- used by Legendre1DPol.monomials
    and any plain-Python (non-JAX-hot-loop) caller.
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
    """Evaluate the Legendre polynomial of the given degree at x, JAX-backed (jnp).

    Same recurrence as legendre_pol, but using jnp ops so it is traceable/
    differentiable/jittable. Closed-form expressions for degrees 0-5; higher
    degrees use the standard three-term recurrence.

    Args:
        degree (int): polynomial degree, >= 0.
        x (float or jnp.ndarray): evaluation point(s), expected in [-1, 1].

    Returns:
        Same type/shape as x: P_degree(x).

    Status: ACTIVE (production default path) -- used inside the JIT-compiled
    fit machinery (fitter.py) and SparseLegendre2DPol.monomials.
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
    """Evaluate the (physicists') Hermite polynomial of the given degree at x, JAX-backed, via the standard three-term recurrence H_n = x*H_{n-1} - (n-1)*H_{n-2}.

    Args:
        degree (int): polynomial degree, >= 0.
        x (float or jnp.ndarray): evaluation point(s).

    Returns:
        Same type/shape as x: H_degree(x).

    Status: ACTIVE (production default path) -- used by the Gauss-Hermite PSF
    basis (psf.py) inside JIT-compiled code.
    """
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
    """Evaluate the (physicists') Hermite polynomial of the given degree at x, NumPy-backed. Same recurrence as hermite_pol_jnp.

    Args:
        degree (int): polynomial degree, >= 0.
        x (float or np.ndarray): evaluation point(s).

    Returns:
        Same type/shape as x: H_degree(x).

    Status: ACTIVE (production default path) -- used by plain-Python (non-JAX-
    hot-loop) Gauss-Hermite evaluation.
    """
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
    """1D Legendre polynomial with coefficients, defined on [xmin, xmax] (internally rescaled to [-1, 1]). NumPy-backed storage/evaluation (monomials/value/invert); derivative() uses jnp.

    Attributes:
        deg (int): polynomial degree.
        xmin, xmax (float): domain the input x is rescaled from.
        coeff (np.ndarray): shape (deg+1,) coefficient vector.

    Status: ACTIVE (production default path) -- the trace-position wavelength
    basis representation.
    """
    def __init__(self, deg=0, xmin=-1.0, xmax=1.0, coeff=None):
        """Args:
            deg (int): polynomial degree.
            xmin, xmax (float): domain x is rescaled from.
            coeff (array-like or None): initial coefficients, shape (deg+1,);
                zero-initialized if None.

        Status: ACTIVE (production default path).
        """
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
        """Compute the stacked Legendre monomial basis [P_0(rx), ..., P_deg(rx)] at x, where rx is x rescaled to [-1, 1]. Deliberately NumPy (not JAX): called from plain-Python hot loops (e.g. gh_params(), once per candidate spot per selection pass) where JAX's eager-mode per-op dispatch overhead would dominate for what's otherwise a handful of scalar flops. The JIT-compiled fit machinery in fitter.py never calls this -- it operates on raw jnp coefficient arrays directly.

        Args:
            x (float or np.ndarray): evaluation point(s), in the [xmin, xmax]
                domain.

        Returns:
            np.ndarray: shape (deg+1, *x.shape) (or (deg+1,) for scalar x).

        Status: ACTIVE (production default path).
        """
        rx = 2 * (x - self.xmin) / (self.xmax - self.xmin) - 1
        return np.stack([legendre_pol(i, rx) for i in range(self.deg + 1)], axis=0)

    def value(self, x):
        """Evaluate the polynomial (coeff . monomials) at x.

        Args:
            x (float or np.ndarray): evaluation point(s).

        Returns:
            Same shape as x (dot of coeff with monomials(x)).

        Status: ACTIVE (production default path).
        """
        m = self.monomials(x)
        return np.dot(self.coeff, m)

    def derivative(self, x):
        """Compute the first derivative of the polynomial at x, via the Legendre derivative recurrence P'_n = P'_{n-2} + (2n-1)*P_{n-1}, chain-ruled through the [xmin, xmax] -> [-1, 1] rescaling. JAX-backed (jnp), unlike value()/monomials().

        Args:
            x (float or jnp.ndarray): evaluation point(s), in the [xmin, xmax]
                domain.

        Returns:
            Same shape as x: d/dx of the polynomial at x.

        Status: ACTIVE (production default path).
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
        """Find x such that value(x) == y, by linear interpolation over a fine (1000-point) grid of value(x) rather than a closed-form/iterative root solve.

        Args:
            y (float or np.ndarray): target value(s) to invert.

        Returns:
            Same shape as y: x such that value(x) ~= y (interpolated, not exact).

        Status: ACTIVE (production default path).
        """
        # Create a grid of 1000 points over the domain
        x_grid = np.linspace(self.xmin, self.xmax, 1000)
        y_grid = self.value(x_grid)
        return np.interp(y, y_grid, x_grid)

class SparseLegendre2DPol:
    """2D Legendre polynomial in (x, y) with a caller-supplied sparse set of nonzero (i, j) monomial index pairs, rather than the full (xdeg+1)*(ydeg+1) dense product basis. JAX-backed (jnp) throughout.

    Attributes:
        xdeg, ydeg (int): per-axis polynomial degrees.
        xmin, xmax, ymin, ymax (float): domain each axis is rescaled from.
        non_zero_indices (list[int]): flat indices k = i + j*(xdeg+1) of the
            monomial terms actually included (see io.get_sparse_nz).
        coeff (jnp.ndarray): shape (len(non_zero_indices),) coefficient vector.

    Status: ACTIVE (production default path) -- the PSF-shape (Gauss-Hermite
    coefficient) wavelength/position basis representation.
    """
    def __init__(self, xdeg=0, xmin=0.0, xmax=1.0, ydeg=0, ymin=0.0, ymax=1.0, coeff=None, non_zero_indices=None):
        """Args:
            xdeg, ydeg (int): per-axis polynomial degrees.
            xmin, xmax, ymin, ymax (float): domain each axis is rescaled from.
            coeff (array-like or None): initial coefficients, shape
                (len(non_zero_indices),); zero-initialized if None.
            non_zero_indices (list[int] or None): flat sparse monomial indices;
                empty if None.

        Status: ACTIVE (production default path).
        """
        self.xdeg = xdeg
        self.xmin = xmin
        self.xmax = xmax
        self.ydeg = ydeg
        self.ymin = ymin
        self.ymax = ymax
        self.non_zero_indices = non_zero_indices if non_zero_indices is not None else []
        self.coeff = jnp.array(coeff) if coeff is not None else jnp.zeros(len(self.non_zero_indices))

    def monomials(self, x, y):
        """Compute the stacked sparse 2D monomial basis at (x, y): for each k in non_zero_indices, the product of the corresponding 1D Legendre terms P_i(rx)*P_j(ry), where (rx, ry) are (x, y) rescaled to [-1, 1] and (i, j) = (k % (xdeg+1), k // (xdeg+1)).

        Args:
            x, y (float or jnp.ndarray): evaluation point(s), same shape.

        Returns:
            jnp.ndarray: shape (len(non_zero_indices), *x.shape).

        Status: ACTIVE (production default path).
        """
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
        """Evaluate the polynomial (coeff . monomials) at (x, y).

        Args:
            x, y (float or jnp.ndarray): evaluation point(s), same shape.

        Returns:
            Same shape as x/y (dot of coeff with monomials(x, y)).

        Status: ACTIVE (production default path).
        """
        m = self.monomials(x, y)
        return jnp.dot(self.coeff, m)
