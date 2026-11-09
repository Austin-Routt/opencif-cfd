"""Alternative duct expansion (sine in y across the width, cosh in z) to cross-check the band solution."""
import math, numpy as np
from scipy.optimize import brentq
def band_alt(d, a, b, N=20001):
    # u = sum_{m odd} 4a^2/(m^3 pi^3) sin(m pi y/a) [1 - cosh(k(z-b/2))/cosh(k b/2)], k = m pi/a
    m = np.arange(1, N+1, 2, dtype=float); k = m*math.pi/a
    zint = b - 2*np.tanh(k*b/2)/k                      # integral over z in [0,b]
    yint = (1 - np.cos(k*d))/k                         # integral of sin over [0,d]
    return np.sum(4*a**2/(m**3*math.pi**3)*yint*zint)
for a,b,f in [(100,140,4.962120761577921e-4),(100,140,4.8751210307083063e-4),(52,140,4.8518e-4),(100,80,4.8550e-4)]:
    Q = band_alt(a,a,b)
    d = brentq(lambda x: band_alt(x,a,b)/Q - f/2, 1e-9, a/2, xtol=1e-14)
    print(a,b,f, f"delta_alt={d:.10f}")
