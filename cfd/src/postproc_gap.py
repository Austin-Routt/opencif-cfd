#!/usr/bin/env python3
"""Post-process a single-gap solve (stokes_gap.py) and, with the same routines, COMSOL's exports of the same run.

usage: python3 postproc_gap.py RUN_DIR --comsol_upstream FILE --comsol_middepth FILE --comsol_fq F [--x_up 75]

Python CFD's field is probed on the grids of COMSOL's own exports (read from the files), so both codes are
measured on identical points with identical routines (skim.py):
  skim width   the depth-averaged band at the gap-side wall that carries the fraction f_q of the upstream plane's flux
               (f_q = Q_side/Q_in of the quarter cell): on the exported grid (native) and with the campaign's resampling
  duct law     the exact rectangular-duct skim width at f_q/2
  trace        the mid-depth dividing streamline (seeded at x = 60 um, reported at x = x_up)
  field        the streamwise velocity on the upstream plane, each code normalized by its own plane-mean velocity
  gap mouth    v at 0.05 um inside the wall line across the opening (the inflow and outflow lobes)
Writes RUN_DIR/metrics.json, RUN_DIR/upstream_x75.txt.gz and RUN_DIR/middepth.txt.gz (x y z u v w, U_in = 1).
"""
import sys, os, json, time, argparse, gzip
import numpy as np
from skfem import MeshTet, Basis, ElementTetP2, ElementVector

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import comsol_io as cio                                                   # noqa: E402
import skim                                                               # noqa: E402


def tensor_finder(m, xs, ys, zs):
    """Point location for a tet mesh made from a tensor hex grid: the hex cell by searchsorted, then the barycentric
    test over that cell's tets and its neighbours' (skfem's own finder can exhaust memory on graded meshes)."""
    p, t = m.p, m.t
    cen = p[:, t].mean(axis=1)
    nx, ny, nz = len(xs) - 1, len(ys) - 1, len(zs) - 1
    ix = np.clip(np.searchsorted(xs, cen[0]) - 1, 0, nx - 1)
    iy = np.clip(np.searchsorted(ys, cen[1]) - 1, 0, ny - 1)
    iz = np.clip(np.searchsorted(zs, cen[2]) - 1, 0, nz - 1)
    hexid = (ix * ny + iy) * nz + iz
    order = np.argsort(hexid, kind='stable')
    hs = hexid[order]
    starts = np.searchsorted(hs, np.arange(nx * ny * nz))
    ends = np.searchsorted(hs, np.arange(nx * ny * nz), side='right')
    A = np.stack([p[:, t[1]] - p[:, t[0]], p[:, t[2]] - p[:, t[0]], p[:, t[3]] - p[:, t[0]]], axis=1)
    Ainv = np.linalg.inv(np.transpose(A, (2, 0, 1)))
    p0 = p[:, t[0]].T

    def finder(x, y, z, _search_all=False):
        x = np.asarray(x); y = np.asarray(y); z = np.asarray(z)
        jx = np.clip(np.searchsorted(xs, x) - 1, 0, nx - 1)
        jy = np.clip(np.searchsorted(ys, y) - 1, 0, ny - 1)
        jz = np.clip(np.searchsorted(zs, z) - 1, 0, nz - 1)
        out = np.full(x.shape, -1, dtype=np.int64)
        best = np.full(x.shape, -np.inf)
        P = np.stack([x, y, z], axis=1)
        for (ox, oy, oz) in [(0, 0, 0), (0, -1, 0), (0, 1, 0), (-1, 0, 0), (1, 0, 0), (0, 0, -1), (0, 0, 1),
                             (-1, -1, 0), (1, -1, 0), (-1, 1, 0), (1, 1, 0)]:
            hid = (np.clip(jx + ox, 0, nx - 1) * ny + np.clip(jy + oy, 0, ny - 1)) * nz + np.clip(jz + oz, 0, nz - 1)
            for k in range(6):
                s0, e0 = starts[hid], ends[hid]
                has = s0 + k < e0
                cand = order[np.minimum(s0 + k, len(order) - 1)]
                d = P - p0[cand]
                lam = np.einsum('nij,nj->ni', Ainv[cand], d)
                l0 = 1 - lam.sum(axis=1)
                score = np.minimum(np.min(lam, axis=1), l0)
                better = has & (score > best)
                out[better] = cand[better]; best[better] = score[better]
            if (best > -1e-9).all():
                break
        if (best < -1e-6).any():
            raise ValueError('%d points outside the mesh (worst barycentric %.2e)' % ((best < -1e-6).sum(), best.min()))
        return out
    return finder


def load(run_dir):
    d = np.load(os.path.join(run_dir, 'solution.npz'))
    meta = json.load(open(os.path.join(run_dir, 'meta.json')))
    m = MeshTet(d['p'], d['t'])
    fnd = tensor_finder(m, d['xs'], d['ys'], d['zs'])
    m.element_finder = lambda mapping=None: fnd
    bu = Basis(m, ElementVector(ElementTetP2()), intorder=4)
    return m, bu, d['u'], meta


def probe(bu, u, pts, chunk=40000):
    out = np.empty((3, pts.shape[1]))
    for i in range(0, pts.shape[1], chunk):
        out[:, i:i + chunk] = (bu.probes(pts[:, i:i + chunk]) @ u).reshape(3, -1)
    return out


def write_plane(path, X, Y, Z, V):
    M = np.column_stack([X.ravel(), Y.ravel(), Z.ravel(), V[0].ravel(), V[1].ravel(), V[2].ravel()])
    with gzip.open(path, 'wt') as fh:
        fh.write('% x y z (um)  u v w (units of U_in, the mean inlet velocity of the quarter cell); NaN in the solid\n')
        np.savetxt(fh, M, fmt='%.9g')


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('run_dir')
    ap.add_argument('--comsol_upstream', required=True)
    ap.add_argument('--comsol_middepth', required=True)
    ap.add_argument('--comsol_fq', type=float, required=True, help="COMSOL's achieved Q_sout/Q_in")
    ap.add_argument('--x_up', type=float, default=75.0)
    ap.add_argument('--x_seed', type=float, default=60.0)
    a = ap.parse_args()
    t0 = time.time()
    m, bu, u, meta = load(a.run_dir)
    G = meta['G']
    wc2, gap, ws, h2 = G['w_c'] / 2, G['gap'], G['w_s'], G['depth'] / 2
    x0, x1, Lx = G['L_entry'], G['L_entry'] + G['L_gap_x'], G['L_x']
    eps = 1e-7
    out = dict(meta=meta)

    # --- mass balance from the solution (fine trapezoid sums of probed values)
    def flux_x(xpos, ylo, yhi):
        yy = np.linspace(ylo + eps, yhi - eps, 401); zz = np.linspace(eps, h2 - eps, 281)
        Y, Z = np.meshgrid(yy, zz, indexing='ij')
        V = probe(bu, u, np.vstack([np.full(Y.size, xpos), Y.ravel(), Z.ravel()]))[0].reshape(Y.shape)
        return float(np.trapezoid(np.trapezoid(V, zz, axis=1), yy))
    Q_in = flux_x(eps, 0, wc2); Q_c = flux_x(Lx - eps, 0, wc2); Q_s = flux_x(Lx - eps, wc2 + gap, wc2 + gap + ws)
    out['mass'] = dict(Q_in=Q_in, Q_center_out=Q_c, Q_side_out=Q_s, imbalance=(Q_in - Q_c - Q_s) / Q_in,
                       f_q_achieved=Q_s / Q_in, Q_in_nominal=wc2 * h2)
    f_q = Q_s / Q_in

    # --- upstream plane on COMSOL's grid
    C_up = cio.read_plane(a.comsol_upstream)
    ys, zs, Cg = cio.as_grid(C_up, 1, 2)
    if abs(np.unique(C_up['x'])[0] - a.x_up) > 1e-6:
        raise ValueError('COMSOL upstream plane is not at x = %g' % a.x_up)
    YY, ZZ = np.meshgrid(ys, zs, indexing='ij')
    pts = np.vstack([np.full(YY.size, a.x_up), np.minimum(YY.ravel(), wc2 - eps), np.clip(ZZ.ravel(), eps, h2 - eps)])
    Vup = probe(bu, u, pts).reshape(3, *YY.shape)
    Vup[:, np.isclose(YY, wc2)] = 0.0; Vup[:, np.isclose(ZZ, 0.0)] = 0.0          # walls (probes sat eps inside)
    write_plane(os.path.join(a.run_dir, 'upstream_x75.txt.gz'), np.full(YY.shape, a.x_up), YY, ZZ, Vup)
    d_py, Qp_py = skim.skim_native(ys, zs, Vup[0], wc2, f_q)
    d_py_rs = skim.skim_resampled(YY.ravel(), ZZ.ravel(), Vup[0].ravel(), wc2, h2, f_q)
    d_c, Qp_c = skim.skim_native(ys, zs, Cg['u'], wc2, a.comsol_fq)
    d_c_rs = skim.skim_resampled(C_up['y'], C_up['z'], C_up['u'], wc2, h2, a.comsol_fq)
    law_py = skim.duct_delta(f_q / 2, G['w_c'], G['depth'])
    law_c = skim.duct_delta(a.comsol_fq / 2, G['w_c'], G['depth'])
    A_plane = wc2 * h2
    up_py = Vup[0] / (Qp_py / A_plane); up_c = np.nan_to_num(Cg['u']) / (Qp_c / A_plane)
    diff = up_py - up_c
    out['skim'] = dict(f_q_python=f_q, f_q_comsol=a.comsol_fq,
                       python=dict(native=d_py, resampled=d_py_rs, duct_law=law_py, native_over_law=d_py / law_py),
                       comsol=dict(native=d_c, resampled=d_c_rs, duct_law=law_c, native_over_law=d_c / law_c),
                       python_over_comsol_native=d_py / d_c, python_over_comsol_resampled=d_py_rs / d_c_rs,
                       python_over_comsol_at_equal_law=(d_py / law_py) / (d_c / law_c),
                       upstream_field=dict(max_abs_diff=float(np.abs(diff).max()),
                                           rel_l2=float(np.sqrt((diff ** 2).mean() / (up_c ** 2).mean())),
                                           peak_python=float(up_py.max()), peak_comsol=float(up_c.max())))
    print('skim: python %.5f (resampled %.5f), comsol %.5f (%.5f); law %.5f / %.5f (%.0f s)'
          % (d_py, d_py_rs, d_c, d_c_rs, law_py, law_c, time.time() - t0), flush=True)

    # --- mid-depth plane on COMSOL's grid (z = h/2), the solid masked
    C_mid = cio.read_plane(a.comsol_middepth)
    xg, yg, Cm = cio.as_grid(C_mid, 0, 1)
    XX, YY = np.meshgrid(xg, yg, indexing='ij')
    solid = (YY > wc2 + 1e-9) & (YY < wc2 + gap - 1e-9) & ~((XX > x0 + 1e-9) & (XX < x1 - 1e-9))
    ingap = (XX >= x0 - 1e-9) & (XX <= x1 + 1e-9) & (YY > wc2) & (YY < wc2 + gap)
    xq = np.where(ingap, np.clip(XX, x0 + eps, x1 - eps), XX)
    yq = np.where(np.isclose(YY, wc2) & ~ingap, wc2 - eps, YY)
    Vm = np.full((3, *XX.shape), np.nan)
    ok = ~solid
    Vm[:, ok] = probe(bu, u, np.vstack([xq[ok], yq[ok], np.full(ok.sum(), h2 - eps)]))
    write_plane(os.path.join(a.run_dir, 'middepth.txt.gz'), XX, YY, np.full(XX.shape, h2), Vm)
    print('mid-depth plane probed (%.0f s)' % (time.time() - t0), flush=True)
    Uc = np.where(solid, 0.0, np.nan_to_num(Cm['u'])); Vc = np.where(solid, 0.0, np.nan_to_num(Cm['v']))
    Up = np.where(solid, 0.0, np.nan_to_num(Vm[0])); Vp = np.where(solid, 0.0, np.nan_to_num(Vm[1]))
    tr_py = skim.trace_dividing(xg, yg, Up, Vp, wc2, gap, x0, x1, a.x_seed, a.x_up)
    tr_c = skim.trace_dividing(xg, yg, Uc, Vc, wc2, gap, x0, x1, a.x_seed, a.x_up)
    for nm, tr in (('python', tr_py), ('comsol', tr_c)):
        if 'path' in tr:
            np.savetxt(os.path.join(a.run_dir, 'dividing_streamline_%s.txt' % nm), tr.pop('path'),
                       header='x y (um) of the mid-depth dividing streamline traced on the %s field' % nm)
    # gap-mouth lobes at 0.05 um inside the wall line, in units of each code's plane-mean velocity at x_up
    jy = int(np.argmin(np.abs(yg - (wc2 - 0.05)))); sel = (xg >= x0) & (xg <= x1)
    lob = {}
    for nm, V, ub in (('python', Vp, Qp_py / A_plane), ('comsol', Vc, Qp_c / A_plane)):
        vv = V[sel, jy] / ub
        lob[nm] = dict(inflow_max=float(vv.max()), outflow_min=float(vv.min()),
                       net_per_um=float(np.trapezoid(vv, xg[sel]) / (x1 - x0)))
    out['middepth'] = dict(trace_python=tr_py, trace_comsol=tr_c, gap_mouth=lob, y_lobes=float(yg[jy]))
    print('traces: python %s, comsol %s (%.0f s)' % (tr_py.get('delta_at_x_report'), tr_c.get('delta_at_x_report'),
                                                      time.time() - t0), flush=True)
    json.dump(out, open(os.path.join(a.run_dir, 'metrics.json'), 'w'), indent=1)
    print('done (%.0f s)' % (time.time() - t0))


if __name__ == '__main__':
    main()
