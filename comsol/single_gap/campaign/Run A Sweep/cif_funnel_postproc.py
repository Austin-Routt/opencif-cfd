#!/usr/bin/env python3
"""
CIF single-gap funnelling post-processor  (Austin Routt, CIF size-cutoff project).

Extracts, from two COMSOL grid exports of one clean single-gap Stokes solution:
    delta_skim_depthavg : depth-averaged near-wall skim thickness (upstream YZ plane)
    delta_skim_mid      : mid-depth skim thickness far upstream  (mid-depth XY plane)
    delta_bulge         : mid-depth dividing-surface closest approach to the centre
    K_inplane           : delta_bulge / delta_skim_mid           (in-plane funnelling)
and compares the depth-averaged skim to the analytic developed-duct skim law.

------------------------------------------------------------------------------
TWO RESULTS THIS SCRIPT SUPPORTS
------------------------------------------------------------------------------
1. K_inplane ~ 1.  The mid-depth dividing streamsurface does NOT bulge toward the
   channel centre at the gap; it THINS toward it (line-sink behaviour), matching
   the analytic Hele-Shaw limit.  delta_bulge is therefore measured with a robust
   FLUX method, not streamline tracing: near the no-slip wall streamlines stall,
   and some poke into the gap then return to the centre channel, which mis-classifies
   the dividing line and spuriously inflates delta_bulge.  K_inplane is a within-plane
   ratio, so it is independent of the fraction convention discussed below.

2. delta_skim matches the exact rectangular-duct skim law to sub-percent -- once the
   PHYSICAL fraction is used (see the factor-of-2 note).

------------------------------------------------------------------------------
THE FACTOR OF 2 (read this before trusting any absolute fraction)
------------------------------------------------------------------------------
The model is a QUARTER cell: y in [0, w_c/2] (wall at w_c/2, symmetry at y=0) and
z in [0, depth/2] (symmetry at depth/2).  COMSOL's report defines

      f_gap = intop_sout / intop_in = Q_sout(quarter) / Q_in(quarter).

The analytic skim law (delta = w_c*sqrt(f/3) and the exact `f_general` below) defines
the fraction as [flux in the near-wall lamina on ONE wall, over the FULL depth] /
[flux through the FULL duct cross-section].  Mapping quarter -> full channel:
the one-wall full-depth lamina flux = 2 x (quarter lamina flux) [z-symmetry], while
the full-duct flux = 4 x Q_in(quarter) [y- and z-symmetry of the developed profile].
Hence

      f_physical  =  f_gap / 2.

So the reported f_gap is TWICE the one-wall/full-duct fraction the law uses.  This
DOES NOT affect delta_skim itself: the extracted layer physically carries
Q_sout(quarter) = f_gap * Q_in(quarter), so `delta_skim_depthavg` is called with the
reported f_gap and returns the correct thickness.  The factor of 2 enters ONLY when
comparing that thickness to the analytic law, which must be evaluated at f_gap/2.
With that correction the depth factor is R ~ 0.95 (the rectangular-duct value at
alpha = depth/w_c ~ 1.4), NOT the ~0.67 you get from comparing against the law at
the un-halved f_gap.

CONSEQUENCE FOR DOWNSTREAM USE: any absolute cutoff/fraction prediction made from the
CFD must use f_physical = f_gap/2.  K_inplane and the delta_skim ~ sqrt(f) scaling are
ratios and are unaffected.

numpy >= 1.24 (uses np.trapezoid), scipy.
"""
import argparse
import numpy as np
from scipy.interpolate import griddata
from scipy.optimize import brentq


# ----------------------------------------------------------------------------
# Analytic developed-duct skim law (exact, and its deep parabolic limit).
# Ported verbatim from cif_cutoff_depth_blocks.py (Blocks 9-10, validated there).
# f_general(delta) = [near-wall one-wall lamina flux] / [full-duct flux], on the exact
# rectangular-duct profile (Fourier y-modes).  Depends on geometry only through
# alpha = h/w_c.  Deep limit: f = 3x^2 - 2x^3, x = delta/w_c.
# ----------------------------------------------------------------------------
def f_general(delta, w_c, h, N=400):
    a, b = w_c / 2.0, h / 2.0
    m = 2 * np.arange(N) + 1                                   # odd modes
    z_int = 2 * b - (4 * a / (m * np.pi)) * np.tanh(m * np.pi * b / (2 * a))
    c = (1.0 / m**3) * z_int * (2 * a / (m * np.pi))
    return np.sum(c * (1 - np.cos(m * np.pi * delta / (2 * a)))) / np.sum(c * 2.0)


def delta_general(f, w_c, h):
    """Invert the exact duct law for the skim thickness at fraction f."""
    return brentq(lambda d: f_general(d, w_c, h) - f, 1e-9 * w_c, 0.999 * w_c)


def delta_parab(f, w_c):
    """Deep-channel parabolic skim: f = 3x^2 - 2x^3, x = delta/w_c.  Small f -> w_c*sqrt(f/3)."""
    return brentq(lambda x: 3 * x**2 - 2 * x**3 - f, 1e-12, 0.5) * w_c


# ----------------------------------------------------------------------------
# COMSOL export parsing.
# ----------------------------------------------------------------------------
def _rows(fn):
    L = [l.rstrip("\n") for l in open(fn)]
    return [l for l in L if l.strip() and not l.lstrip().startswith('%')]


def _vec(l):
    return np.array([float(t) for t in l.split()])


def parse_scatter(fn):
    """Return scattered (x,y,z,u,v,w).  Handles the two COMSOL export layouts seen in
    practice: (a) column 'x y z u v w' rows; (b) Grid block: 3 coordinate lines then
    u,v,w as flattened (n_axisC x n_axisB) blocks with one singleton axis."""
    nc = _rows(fn)
    widths = [len(l.split()) for l in nc[:30]]
    if len(set(widths)) == 1 and widths[0] == 6:                  # (a) columns
        M = np.array([_vec(l) for l in nc])
        return M[:, 0], M[:, 1], M[:, 2], M[:, 3], M[:, 4], M[:, 5]
    a, b, c = _vec(nc[0]), _vec(nc[1]), _vec(nc[2])               # (b) grid block
    data = np.concatenate([_vec(l) for l in nc[3:]])
    axes = [a, b, c]; ns = [len(a), len(b), len(c)]
    sing = int(np.argmin(ns))                                     # singleton = fixed plane
    others = [i for i in range(3) if i != sing]
    A1, A2 = axes[others[0]], axes[others[1]]; n1, n2 = len(A1), len(A2)
    blk = data.reshape(3, n2, n1)                                 # COMSOL: (expr, axis2, axis1)
    G2, G1 = np.meshgrid(A2, A1, indexing='ij')
    flat = np.full((3, n2 * n1), np.nan)
    for k in range(3):
        flat[k] = blk[k].ravel()
    coords = {others[0]: G1.ravel(), others[1]: G2.ravel(),
              sing: np.full(n1 * n2, axes[sing][0])}
    return coords[0], coords[1], coords[2], flat[0], flat[1], flat[2]


def classify_plane(fn):
    """Return 'upstream' (x fixed -> YZ plane) or 'middepth' (z fixed -> XY plane) by
    detecting which coordinate is (near-)constant.  Lets an archive be processed without
    relying on file names."""
    x, y, z, *_ = parse_scatter(fn)
    spans = [np.ptp(x), np.ptp(y), np.ptp(z)]
    fixed = int(np.argmin(spans))
    return 'upstream' if fixed == 0 else ('middepth' if fixed == 2 else 'unknown')


# ----------------------------------------------------------------------------
# Skim extraction.
# ----------------------------------------------------------------------------
def delta_skim_depthavg(upstream_file, w_c, f_gap, depth):
    """Depth-averaged skim thickness from the upstream YZ plane (x fixed).

    Finds delta such that the near-wall band {w_c/2-delta <= y <= w_c/2, all z in the
    quarter} carries flux f_gap * Q_in(quarter) = Q_sout(quarter) -- i.e. the layer the
    gap actually extracts.  Pass the REPORTED f_gap here (not f_gap/2); this returns the
    physical layer thickness.  Compare it to the analytic law at f_gap/2 (module note).
    """
    x, y, z, u, v, w = parse_scatter(upstream_file)
    wall = w_c / 2.0; halfd = depth / 2.0
    yg = np.linspace(0, wall, 400); zg = np.linspace(0, halfd, 140)
    U = griddata((y, z), u, (yg[:, None], zg[None, :]), method='linear')
    U = np.nan_to_num(np.clip(U, 0, None))
    q = np.trapezoid(U, zg, axis=1)                              # q(y) = int u dz
    Qtot = np.trapezoid(q, yg); Qgap = f_gap * Qtot
    s = (wall - yg)[::-1]; qs = q[::-1]; m = s >= 0; s, qs = s[m], qs[m]
    cum = np.concatenate([[0], np.cumsum(0.5 * (qs[1:] + qs[:-1]) * np.diff(s))])
    return float(np.interp(Qgap, cum, s))


def middepth_metrics(mid_file, w_c, gap_x0, gap_x1, x_up=None):
    """Flux-method dividing line y_div(x) in the mid-depth plane.
    Returns (delta_skim_mid, delta_bulge, xs, ydivs, trend, phi)."""
    x, y, z, u, v, w = parse_scatter(mid_file)
    wall = w_c / 2.0
    x0, x1 = x.min(), x.max()
    if x_up is None:
        x_up = x0 + 0.6 * (gap_x0 - x0)                          # upstream sampling station
    xg = np.arange(x0, x1 + 1e-9, 0.5); yg = np.arange(0, y.max() + 1e-9, 0.1)
    XX, YY = np.meshgrid(xg, yg, indexing='ij')
    U = griddata((x, y), u, (XX, YY), method='linear')
    V = griddata((x, y), v, (XX, YY), method='linear')
    from scipy.interpolate import RegularGridInterpolator as RGI
    fu = RGI((xg, yg), np.nan_to_num(U), bounds_error=False, fill_value=0.0)
    fv = RGI((xg, yg), np.nan_to_num(V), bounds_error=False, fill_value=0.0)
    # mid-depth gap inflow per unit depth: int v(x, wall-) dx over the aperture
    xa = np.arange(gap_x0, gap_x1 + 1e-9, 0.1)
    phi = np.trapezoid([float(fv([[xx, wall - 0.05]]).item()) for xx in xa], xa)

    def ydiv(xs):                                               # y where cum flux from wall == phi
        yy = np.arange(wall, wall - 8, -0.02)
        uu = np.array([float(fu([[xs, t]]).item()) for t in yy])
        cum = np.concatenate([[0], np.cumsum(0.5 * (uu[1:] + uu[:-1]) * np.abs(np.diff(yy)))])
        if cum[-1] < phi:
            return np.nan
        return float(np.interp(phi, cum, yy))

    xs = np.arange(x_up, gap_x0 + 1e-9, 2.0)
    # Guard: at low fraction the net gap through-flow is a tiny residual of large
    # antisymmetric recirculation lobes at the gap edges; its discrete value can be
    # dominated by that recirculation (phi <= 0). When that happens the dividing-line
    # flux match is not resolvable and K_inplane is reported as indeterminate (NaN).
    # This is a post-processing limit at small f, NOT a model issue -- the depth-averaged
    # delta_skim (which drives the validation) is unaffected, and K_inplane is established
    # at the higher fractions and in the shallower runs.
    if phi <= 0:
        nan = float('nan')
        return nan, nan, xs, np.full(len(xs), np.nan), \
            "indeterminate (gap through-flow below recirculation floor at this low f)", phi
    yd = np.array([ydiv(t) for t in xs])
    d_mid = wall - yd[0]                                         # mid-depth skim (upstream)
    d_bulge = wall - np.nanmin(yd)                              # closest approach to centre
    if np.isnan(d_mid) or d_mid < 1e-6:
        nan = float('nan')
        return nan, nan, xs, yd, \
            "indeterminate (dividing line unresolved at this low f)", phi
    trend = ("THINS toward gap (line-sink, K_inplane<=1)" if yd[-1] > yd[0]
             else "BULGES toward centre (K_inplane>1)")
    return d_mid, d_bulge, xs, yd, trend, phi


def analyze(upstream_file, middepth_file, w_c, f_gap, depth, gap_x0, gap_x1):
    """Return all per-run metrics as a dict.  f_gap is the reported COMSOL fraction."""
    f_phys = f_gap / 2.0                                         # one-wall/full-duct fraction
    dsk_da = delta_skim_depthavg(upstream_file, w_c, f_gap, depth)
    dsk_md, dbulge, xs, yd, trend, phi = middepth_metrics(middepth_file, w_c, gap_x0, gap_x1)
    dp = delta_parab(f_phys, w_c)                               # deep parabolic at f_phys
    dg = delta_general(f_phys, w_c, depth)                      # exact duct law at f_phys
    return dict(f_gap=f_gap, f_phys=f_phys, alpha=depth / w_c,
                delta_skim_da=dsk_da, delta_skim_mid=dsk_md, delta_bulge=dbulge,
                K_inplane=dbulge / dsk_md, trend=trend,
                delta_parab=dp, delta_general=dg,
                R_depth=dsk_da / dp,            # depth correction vs deep parabolic (~0.95)
                validation=dsk_da / dg)         # CFD vs exact duct law (~1.00)


def main():
    p = argparse.ArgumentParser(description="Per-run CIF single-gap funnelling metrics.")
    p.add_argument("--upstream", required=True, help="upstream YZ-plane export (x fixed)")
    p.add_argument("--middepth", required=True, help="mid-depth XY-plane export (z fixed)")
    p.add_argument("--wc", type=float, required=True, help="channel width w_c [um]")
    p.add_argument("--f", type=float, required=True,
                   help="REPORTED f_gap = intop_sout/intop_in from the COMSOL report")
    p.add_argument("--depth", type=float, required=True, help="channel depth h [um]")
    p.add_argument("--gap_x0", type=float, required=True)
    p.add_argument("--gap_x1", type=float, required=True)
    p.add_argument("--label", default="run")
    a = p.parse_args()
    r = analyze(a.upstream, a.middepth, a.wc, a.f, a.depth, a.gap_x0, a.gap_x1)
    print("=" * 66)
    print(f"  CIF funnelling post-process  [{a.label}]")
    print(f"  w_c={a.wc} um  depth={a.depth} um  alpha={r['alpha']:.2f}")
    print(f"  f_gap (reported) = {r['f_gap']:.4e}   f_physical = f_gap/2 = {r['f_phys']:.4e}")
    print("=" * 66)
    print(f"  analytic exact duct  delta_general(f_phys) = {r['delta_general']:.4f} um")
    print(f"  analytic deep parab  w_c*sqrt(f_phys/3)    = {r['delta_parab']:.4f} um")
    print(f"  delta_skim (depth-averaged)                = {r['delta_skim_da']:.4f} um")
    print(f"  delta_skim (mid-depth, upstream)           = {r['delta_skim_mid']:.4f} um")
    print(f"  delta_bulge(mid-depth, closest)            = {r['delta_bulge']:.4f} um")
    print(f"  mid-depth dividing line: {r['trend']}")
    print("  " + "-" * 62)
    print(f"  K_inplane = delta_bulge/delta_skim_mid     = {r['K_inplane']:.3f}"
          "   (expect ~1: no funnelling)")
    print(f"  depth factor R = delta_skim_da/delta_parab = {r['R_depth']:.3f}"
          "   (rect-duct R at this alpha)")
    print(f"  VALIDATION  delta_skim_da/delta_general    = {r['validation']:.3f}"
          "   (~1.00: CFD == exact duct law)")
    print("=" * 66)


if __name__ == "__main__":
    main()
