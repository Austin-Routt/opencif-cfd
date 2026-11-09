"""The single-gap metrics, applied in the same way to COMSOL's exports and to Python CFD's.

Quarter cell: center half-channel y in [0, wall], wall = w_c/2 at the gap side, symmetry at y = 0; depth z in [0, h/2],
floor at z = 0, symmetry at z = h/2.  The quarter-cell fraction is f_q = Q_side/Q_in (what COMSOL's report calls
f_gap); the fraction of the exact duct law (one wall, full depth, over the full duct) is f_q/2.

  duct_fraction(delta, w_c, h)       fraction of the full-duct flux that a band of width delta at one side wall carries
                                     over the full depth, for developed flow in a w_c x h rectangular duct (Fourier
                                     series in the width coordinate)
  duct_delta(f, w_c, h)              its inverse (the exact duct law)
  skim_native(y, z, U, wall, f_q)    the skim width on an upstream plane: the band at the wall that carries f_q of the
                                     plane's flux, integrating the exported values on their own grid (trapezoid in z,
                                     piecewise-linear in y, inverted exactly)
  skim_resampled(...)                the campaign post-processor's definition: the plane resampled by linear
                                     interpolation onto 400 x 140 points, negative u set to 0, trapezoid integration,
                                     linear interpolation of the cumulative flux (written from its documented method,
                                     for comparison with its reference output; it is not imported)
  trace_dividing(...)                the mid-depth dividing streamline: bisection on the seed height at x_seed between
                                     a path that passes the gap and one that reaches the side channel
"""
import numpy as np
from scipy.optimize import brentq


def duct_fraction(delta, w_c, h, N=400):
    a, b = w_c / 2.0, h / 2.0
    m = 2 * np.arange(N) + 1
    k = m * np.pi / (2 * a)
    depth_int = 2 * b - 2 * np.tanh(k * b) / k               # int_{-b}^{b} [1 - cosh(k z)/cosh(k b)] dz
    c = depth_int / (m ** 3 * k)
    return float(np.sum(c * (1 - np.cos(k * delta))) / np.sum(2 * c))


def duct_delta(f, w_c, h):
    return brentq(lambda d: duct_fraction(d, w_c, h) - f, 1e-9 * w_c, 0.999 * w_c, xtol=1e-14)


def parabolic_delta(f, w_c):
    """Deep-channel limit: f = 3x^2 - 2x^3, x = delta/w_c."""
    return brentq(lambda x: 3 * x ** 2 - 2 * x ** 3 - f, 1e-12, 0.5, xtol=1e-15) * w_c


def _invert_cumulative(s, q, target):
    """s ascending from 0 (distance from the wall), q >= 0 the flux density at s, piecewise linear; returns s* with
    int_0^s* q ds = target exactly for the piecewise-linear q."""
    ds = np.diff(s)
    cum = np.concatenate([[0.0], np.cumsum(0.5 * (q[1:] + q[:-1]) * ds)])
    i = int(np.searchsorted(cum, target) - 1)
    i = min(max(i, 0), len(ds) - 1)
    r = target - cum[i]
    q0, q1, h = q[i], q[i + 1], ds[i]
    g = (q1 - q0) / h
    if abs(g) < 1e-300:
        t = r / q0
    else:
        t = (-q0 + np.sqrt(max(q0 * q0 + 2 * g * r, 0.0))) / g
    return float(s[i] + t)


def skim_native(y, z, U, wall, f_q):
    """y (ny,), z (nz,) ascending; U (ny, nz) streamwise velocity on the upstream plane (quarter cell)."""
    U = np.nan_to_num(U)
    q = np.trapezoid(U, z, axis=1)                          # q(y) = int u dz
    Q = np.trapezoid(q, y)
    s = wall - y[::-1]; qs = q[::-1]
    keep = s >= -1e-12
    s, qs = np.clip(s[keep], 0, None), qs[keep]
    return _invert_cumulative(s, np.clip(qs, 0, None), f_q * Q), float(Q)


def skim_resampled(y, z, u, wall, halfdepth, f_q):
    """Scattered (y, z, u) on the upstream plane; the campaign's resampling definition."""
    from scipy.interpolate import griddata
    yg = np.linspace(0, wall, 400); zg = np.linspace(0, halfdepth, 140)
    U = griddata((y, z), u, (yg[:, None], zg[None, :]), method='linear')
    U = np.nan_to_num(np.clip(U, 0, None))
    q = np.trapezoid(U, zg, axis=1)
    Qtot = np.trapezoid(q, yg)
    s = (wall - yg)[::-1]; qs = q[::-1]; m = s >= 0; s, qs = s[m], qs[m]
    cum = np.concatenate([[0], np.cumsum(0.5 * (qs[1:] + qs[:-1]) * np.diff(s))])
    return float(np.interp(f_q * Qtot, cum, s))


def trace_dividing(xg, yg, Ux, Vy, wall, gap, gap_x0, gap_x1, x_seed, x_report, solid=None):
    """Mid-depth dividing streamline on a tensor grid (xg, yg) with in-plane velocity (Ux, Vy) [nx, ny]; points in the
    solid wall carry 0.  A path is 'captured' when it reaches y > wall + gap + 0.5 (the side channel) and 'passed' when
    it is 5 um downstream of the gap in the center channel.  Returns the wall distance of the dividing streamline at
    x_report, its minimum between x_report and the gap edge, and the path."""
    from scipy.interpolate import RegularGridInterpolator as RGI
    from scipy.integrate import solve_ivp
    fu = RGI((xg, yg), np.nan_to_num(Ux), bounds_error=False, fill_value=None)
    fv = RGI((xg, yg), np.nan_to_num(Vy), bounds_error=False, fill_value=None)

    def trace(y0, smax=900.0):
        def rhs(s, q):
            vel = np.array([fu([q])[0], fv([q])[0]])
            n = np.hypot(*vel)
            return vel / n if n > 1e-14 else np.array([1.0, 0.0])
        def captured(s, q): return q[1] - (wall + gap + 0.5)
        captured.terminal = True
        def passed(s, q): return q[0] - (gap_x1 + 5.0)
        passed.terminal = True
        sol = solve_ivp(rhs, [0, smax], [x_seed, y0], events=[captured, passed], max_step=0.2, rtol=1e-8, atol=1e-10)
        st = 'captured' if len(sol.t_events[0]) else ('passed' if len(sol.t_events[1]) else 'undetermined')
        return st, sol.y

    lo, hi = wall - 6.0, wall - 1e-3
    s_lo, _ = trace(lo); s_hi, _ = trace(hi)
    if not (s_lo == 'passed' and s_hi == 'captured'):
        return dict(error='bracket failed: %s / %s' % (s_lo, s_hi))
    for _ in range(30):
        mid = 0.5 * (lo + hi)
        st, _ = trace(mid)
        if st == 'passed':
            lo = mid
        else:
            hi = mid                                          # captured or trapped counts with the captured side
    st, path = trace(lo)
    xp, yp = path
    up = xp <= gap_x0
    dist = wall - yp
    sel = up & (xp >= x_report)
    return dict(y_seed=float(0.5 * (lo + hi)), delta_at_x_report=float(np.interp(x_report, xp[up], dist[up])),
                delta_min_upstream=float(dist[sel].min()) if sel.any() else None,
                delta_at_gap_edge=float(np.interp(gap_x0, xp[up], dist[up])), path=np.column_stack([xp, yp]))
