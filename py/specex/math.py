import numpy as np

def legendre_pol(degree, x):
    """
    Evaluates the Legendre polynomial of a given degree at points x.
    x should be in the range [-1, 1].
    Matches the recursive definition in specex_legendre.cc.
    """
    if degree == 0:
        return np.ones_like(x)
    if degree == 1:
        return x
    if degree == 2:
        return 0.5 * (3 * x**2 - 1.)
    if degree == 3:
        return 0.5 * (5 * x**2 - 3) * x
    if degree == 4:
        return 0.125 * ((35 * x**2 - 30) * x**2 + 3)
    if degree == 5:
        return 0.125 * (((63 * x**2 - 70) * x**2 + 15) * x)
    if degree == 6:
        return (((231 * x**2 - 315) * x**2 + 105) * x**2 - 5) / 16.
    
    # Recursive case for degree > 6
    p_prev = legendre_pol(degree - 1, x)
    p_prev2 = legendre_pol(degree - 2, x)
    return ((2 * degree - 1) * x * p_prev - (degree - 1) * p_prev2) / degree

class Legendre1DPol:
    def __init__(self, deg=0, xmin=0.0, xmax=0.0, coeff=None):
        self.deg = deg
        self.xmin = xmin
        self.xmax = xmax
        if coeff is not None:
            self.coeff = np.array(coeff)
        else:
            self.coeff = np.zeros(deg + 1)

    def monomials(self, x):
        rx = 2 * (x - self.xmin) / (self.xmax - self.xmin) - 1
        m = np.zeros((self.deg + 1,) + np.shape(x))
        for i in range(self.deg + 1):
            m[i] = legendre_pol(i, rx)
        return m

    def value(self, x):
        m = self.monomials(x)
        return np.dot(self.coeff, m)

class Legendre2DPol:
    def __init__(self, xdeg=0, xmin=0.0, xmax=1.0, ydeg=0, ymin=0.0, ymax=1.0, coeff=None):
        self.xdeg = xdeg
        self.xmin = xmin
        self.xmax = xmax
        self.ydeg = ydeg
        self.ymin = ymin
        self.ymax = ymax
        if coeff is not None:
            self.coeff = np.array(coeff)
        else:
            self.coeff = np.zeros((xdeg + 1) * (ydeg + 1))

    def monomials(self, x, y):
        rx = 2 * (x - self.xmin) / (self.xmax - self.xmin) - 1
        ry = 2 * (y - self.ymin) / (self.ymax - self.ymin) - 1
        
        # This implementation matches the nested loop in specex_legendre.cc
        # index = i + j*(xdeg+1)
        m = np.zeros(((self.xdeg + 1) * (self.ydeg + 1),) + np.shape(x))
        mx = [legendre_pol(i, rx) for i in range(self.xdeg + 1)]
        
        idx = 0
        for j in range(self.ydeg + 1):
            myj = legendre_pol(j, ry)
            for i in range(self.xdeg + 1):
                m[idx] = mx[i] * myj
                idx += 1
        return m

    def value(self, x, y):
        m = self.monomials(x, y)
        return np.dot(self.coeff, m)

class SparseLegendre2DPol:
    def __init__(self, xdeg=0, xmin=0.0, xmax=1.0, ydeg=0, ymin=0.0, ymax=1.0, coeff=None, non_zero_indices=None):
        self.xdeg = xdeg
        self.xmin = xmin
        self.xmax = xmax
        self.ydeg = ydeg
        self.ymin = ymin
        self.ymax = ymax
        self.non_zero_indices = non_zero_indices if non_zero_indices is not None else []
        if coeff is not None:
            self.coeff = np.array(coeff)
        else:
            self.coeff = np.zeros(len(self.non_zero_indices))

    def monomials(self, x, y):
        rx = 2 * (x - self.xmin) / (self.xmax - self.xmin) - 1
        ry = 2 * (y - self.ymin) / (self.ymax - self.ymin) - 1
        
        m = np.zeros((len(self.non_zero_indices),) + np.shape(x))
        for idx, k in enumerate(self.non_zero_indices):
            i = k % (self.xdeg + 1)
            j = k // (self.xdeg + 1)
            m[idx] = legendre_pol(i, rx) * legendre_pol(j, ry)
        return m

    def value(self, x, y):
        m = self.monomials(x, y)
        return np.dot(self.coeff, m)
