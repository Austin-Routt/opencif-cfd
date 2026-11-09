"""Independent single-gap skim-width checker.
Usage: python3 -I singlegap.py <label> <grid_file> <w_c> <h> <f_q> [<label> <file> <w_c> <h> <f_q> ...]
"""
import sys, math
import numpy as np
from scipy.optimize import brentq


def read_grid(fn):
    axes, blocks, names = [], [], []
    cur = None
    with open(fn, encoding="utf-8", errors="replace") as f:
        for line in f:
            s = line.strip()
            if not s:
                continue
            if s.startswith("%"):
                if s.startswith("% Data"):
                    cur = []; blocks.append(cur)
                elif cur is not None and len(cur) == 0 and len(names) < len(blocks):
                    names.append(s[1:].strip())
                continue
            vals = [float(t) for t in s.split()]
            if cur is None:
                axes.append(np.array(vals))
            else:
                cur.extend(vals)
    assert len(axes) == 3, len(axes)
    return axes, [np.array(b) for b in blocks], names


def to_field(flat, ny, nz, order):
    if order == "y_fastest":   # index = iy + ny*iz
        return flat.reshape(nz, ny).T   # [iy, iz]
    else:                      # z fastest: index = iz + nz*iy
        return flat.reshape(ny, nz)     # [iy, iz]


def skim_width(y, q, frac):
    """q(y) piecewise linear on nodes y; return delta with int_{yN-delta}^{yN} q = frac * int q."""
    seg = 0.5 * (q[1:] + q[:-1]) * np.diff(y)
    Q = seg.sum()
    T = frac * Q
    acc = 0.0
    N = len(y) - 1
    for i in range(N, 0, -1):
        a = seg[i - 1]
        if acc + a >= T:
            R = T - acc
            hseg = y[i] - y[i - 1]
            qi, qim = q[i], q[i - 1]
            c = (qim - qi) / (2 * hseg)
            disc = qi * qi + 4 * c * R
            s = 2 * R / (qi + math.sqrt(disc)) if (qi + math.sqrt(disc)) > 0 else math.sqrt(R / c)
            return (y[N] - y[i]) + s, Q
        acc += a
    raise RuntimeError


# ---------- exact rectangular duct (width a in y, depth b in z), no-slip all walls ----------
# u(y,z) = (G/mu) sum_{n odd} 4 b^2/(n^3 pi^3) sin(n pi z/b) [1 - cosh(l_n (y-a/2))/cosh(l_n a/2)], l_n = n pi/b
# band flow 0<=y<=d over full depth:
# Q(d) = (G/mu) sum_{n odd} 8 b^3/(n^4 pi^4) [ d - ( sinh(l(d-a/2)) + sinh(l a/2) ) / (l cosh(l a/2)) ]
NMAX = 4001
def band_flow(d, a, b):
    n = np.arange(1, NMAX + 1, 2, dtype=float)
    l = n * math.pi / b
    t1 = (np.exp(l * (d - a)) - np.exp(-l * d)) / (1 + np.exp(-l * a))   # sinh(l(d-a/2))/cosh(l a/2)
    br = d - (t1 + np.tanh(l * a / 2)) / l
    return np.sum(8 * b ** 3 / (n ** 4 * math.pi ** 4) * br)


def duct_fraction(d, a, b):
    return band_flow(d, a, b) / band_flow(a, a, b)


def duct_delta(target, a, b):
    return brentq(lambda d: duct_fraction(d, a, b) - target, 1e-12, a / 2, xtol=1e-14, rtol=1e-15)


def duct_selftest():
    a, b = 100.0, 140.0
    # closed-form total: a b^3/12 - 16 b^4/pi^5 sum tanh(n pi a/2b)/n^5
    n = np.arange(1, NMAX + 1, 2, dtype=float)
    Qcf = a * b ** 3 / 12 - 16 * b ** 4 / math.pi ** 5 * np.sum(np.tanh(n * math.pi * a / (2 * b)) / n ** 5)
    # brute-force numeric integration of the series velocity
    ny, nz = 2001, 2001
    yy = np.linspace(0, a, ny); zz = np.linspace(0, b, nz)
    U = np.zeros((ny, nz))
    for nn in range(1, 400, 2):
        l = nn * math.pi / b
        fy = 1 - np.cosh(l * (yy - a / 2)) / np.cosh(l * a / 2)
        U += 4 * b ** 2 / (nn ** 3 * math.pi ** 3) * np.outer(fy, np.sin(nn * math.pi * zz / b))
    # Laplacian check at interior point: should be -1
    iy, iz = 700, 900
    dy, dz = yy[1] - yy[0], zz[1] - zz[0]
    lap = (U[iy + 1, iz] - 2 * U[iy, iz] + U[iy - 1, iz]) / dy ** 2 + (U[iy, iz + 1] - 2 * U[iy, iz] + U[iy, iz - 1]) / dz ** 2
    Qnum = np.trapezoid(np.trapezoid(U, zz, axis=1), yy)
    d = 3.0
    m = yy <= d + 1e-12
    Qb = np.trapezoid(np.trapezoid(U[m], zz, axis=1), yy[m])
    print(f"duct self-test: Q_series={band_flow(a,a,b):.10e} Q_closedform={Qcf:.10e} Q_numeric={Qnum:.10e}; "
          f"band(3um) series={band_flow(d,a,b):.8e} numeric={Qb:.8e}; Laplacian at interior={lap:.6f} (expect -1)")


duct_selftest()
args = sys.argv[1:]
res = {}
for k in range(0, len(args), 5):
    label, fn, wc, h, fq = args[k], args[k + 1], float(args[k + 2]), float(args[k + 3]), float(args[k + 4])
    (xa, ya, za), blocks, names = read_grid(fn)
    ny, nz = len(ya), len(za)
    print(f"\n### {label}: x={xa}, ny={ny} y[{ya[0]},{ya[-1]}], nz={nz} z[{za[0]},{za[-1]}], blocks={len(blocks)} sizes={[len(b) for b in blocks]} names={names}")
    assert all(len(b) == ny * nz for b in blocks)
    u_flat = blocks[0]
    print(f"  NaN count per block: {[int(np.isnan(b).sum()) for b in blocks]}")
    for order in ("y_fastest", "z_fastest"):
        try:
            U = to_field(u_flat, ny, nz, order)
        except ValueError:
            continue
        floor = np.nanmax(np.abs(U[:, 0])); sidewall = np.nanmax(np.abs(U[-1, :]))
        print(f"  ordering {order}: max|u| on floor z=0: {floor:.3e}; on wall y=w_c/2: {sidewall:.3e}; max|u| overall {np.nanmax(np.abs(U)):.6e} at (y,z)={np.unravel_index(np.nanargmax(U), U.shape)}")
    U = to_field(u_flat, ny, nz, "y_fastest")
    V = to_field(blocks[1], ny, nz, "y_fastest"); W = to_field(blocks[2], ny, nz, "y_fastest")
    print(f"  (y-fastest) max|v| on y=w_c/2 wall {np.nanmax(np.abs(V[-1,:])):.3e}, max|w| on z=0 floor {np.nanmax(np.abs(W[:,0])):.3e}, max|v| at y=0 {np.nanmax(np.abs(V[0,:])):.3e}, max|w| at z=h/2 {np.nanmax(np.abs(W[:,-1])):.3e}")
    assert abs(ya[-1] - wc / 2) < 1e-9 and abs(za[-1] - h / 2) < 1e-9
    U0 = np.nan_to_num(U)
    q = np.trapezoid(U0, za, axis=1)  # [iy]
    delta, Q = skim_width(ya, q, fq)
    Qplane = Q * 1e-12
    dd = duct_delta(fq / 2, wc, h)
    res[label] = (delta, dd, fq, wc, h)
    print(f"  plane flow (quarter) = {Qplane:.10e} m^3/s; u max {np.nanmax(U):.8e}")
    print(f"  delta = {delta:.10f} um; duct delta = {dd:.10f} um; ratio = {delta/dd:.10f}")
    print(f"  check: duct fraction at dd = {duct_fraction(dd, wc, h):.12e} (target {fq/2:.12e})")
    # sensitivity: delta with linear interp vs simple trapezoid cumulative (same thing) and alternative z-integration (Simpson)
    from scipy.integrate import simpson
    q2 = simpson(U0, x=za, axis=1)
    d2, _ = skim_width(ya, q2, fq)
    print(f"  (info) delta with Simpson in z = {d2:.10f}")

if "RunB_uniform" in res and "RunB_developed" in res:
    du, ddu = res["RunB_uniform"][:2]; dv, ddv = res["RunB_developed"][:2]
    print(f"\nRun B: (dev/duct_dev)/(uni/duct_uni) - 1 = {(dv/ddv)/(du/ddu) - 1:.10e}")
