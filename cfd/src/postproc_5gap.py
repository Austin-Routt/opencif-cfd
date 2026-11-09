#!/usr/bin/env python3
"""Fluxes, F, the probe pressures, the flow through each gap and the two field planes of a five-gap solve
(stokes_5gap.py), on the points of COMSOL's exports.

usage: python3 postproc_5gap.py RUN_DIR

  Q_in (x = 0), Q_c and Q_s (x = L_total), F = 1 - Q_c/Q_in, in m^3/s for the half-width model at full depth
  pressures at (38 + 49k, 0, 70) and (38 + 49k, w_c/2 + gap + w_s/2, 70), k = 0..4, and their difference per gap
  the flow through each gap: v integrated over the gap mid-plane y = w_c/2 + gap/2 on COMSOL's export grid
      (x every 0.25 um, z every 1 um), doubled for the full depth
  the mid-depth plane z = 70 um on COMSOL's export grid (x and y every 0.5 um), written to middepth_z70.txt.gz,
      and the gap mid-plane, written to gapplane.txt.gz (x y z u v w p in um, m/s and Pa; NaN in the solid)
Writes RUN_DIR/metrics.json.  The comparison with COMSOL is compare.py's.
"""
import sys, os, json, argparse, gzip, time
import numpy as np
from skfem import MeshTet, Basis, ElementTetP2, ElementTetP1, ElementVector
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from postproc_gap import tensor_finder                                        # noqa: E402
import fivegap                                                                # noqa: E402


def write_plane(path, cols, header):
    with gzip.open(path, 'wt') as fh:
        fh.write('% ' + header + '\n')
        np.savetxt(fh, np.column_stack(cols), fmt='%.10g')


def main():
    ap = argparse.ArgumentParser(); ap.add_argument('run_dir'); a = ap.parse_args()
    t0 = time.time()
    d = np.load(os.path.join(a.run_dir, 'solution.npz')); meta = json.load(open(os.path.join(a.run_dir, 'meta.json')))
    G = meta['G']; wc2, gap, ws, h2 = G['w_c'] / 2, G['gap'], G['w_s'], G['depth'] / 2; Lx = meta['Lx']
    m = MeshTet(d['p'], d['t']); fnd = tensor_finder(m, d['xs'], d['ys'], d['zs']); m.element_finder = lambda mapping=None: fnd
    bu = Basis(m, ElementVector(ElementTetP2()), intorder=4); bp = bu.with_element(ElementTetP1())
    u, pr = d['u'], d['pr']
    scale_u = meta['p_in_Pa'] / meta['mu_Pa_s'] * 1e-6     # m/s per unit of the solved velocity (lengths in um)
    eps = 1e-7

    def probe_u(pts, chunk=40000):
        out = np.empty((3, pts.shape[1]))
        for i in range(0, pts.shape[1], chunk):
            out[:, i:i + chunk] = (bu.probes(pts[:, i:i + chunk]) @ u).reshape(3, -1)
        return out

    def probe_p(pts, chunk=40000):
        out = np.empty(pts.shape[1])
        for i in range(0, pts.shape[1], chunk):
            out[i:i + chunk] = bp.probes(pts[:, i:i + chunk]) @ pr
        return out

    def flux_x(xpos, ylo, yhi):
        yy = np.linspace(ylo + eps, yhi - eps, 321); zz = np.linspace(eps, h2 - eps, 281)
        Y, Z = np.meshgrid(yy, zz, indexing='ij')
        V = probe_u(np.vstack([np.full(Y.size, xpos), Y.ravel(), Z.ravel()]))[0].reshape(Y.shape)
        return 2.0 * float(np.trapezoid(np.trapezoid(V, zz, axis=1), yy)) * scale_u * 1e-12   # m^3/s, full depth
    Q_in = flux_x(eps, 0, wc2); Q_c = flux_x(Lx - eps, 0, wc2); Q_s = flux_x(Lx - eps, wc2 + gap, wc2 + gap + ws)
    xs_probe = np.array([G['L_post'] + G['L_gap_x'] / 2 + k * (G['L_post'] + G['L_gap_x']) for k in range(G['n_gaps'])])
    y_s = wc2 + gap + ws / 2
    Pp = probe_p(np.vstack([xs_probe, np.full(5, eps), np.full(5, h2 - eps)])) * meta['p_in_Pa']
    Ps = probe_p(np.vstack([xs_probe, np.full(5, y_s), np.full(5, h2 - eps)])) * meta['p_in_Pa']
    out = dict(meta=meta, Q_in=Q_in, Q_c=Q_c, Q_s=Q_s, F=1 - Q_c / Q_in, F_split=Q_s / (Q_c + Q_s),
               imbalance=(Q_in - Q_c - Q_s) / Q_in, x_probe=xs_probe.tolist(), y_side_probe=y_s,
               p_c=Pp.tolist(), p_s=Ps.tolist(), dp_gap=(Pp - Ps).tolist())
    print('fluxes and probes (%.0f s): F %.6f' % (time.time() - t0, out['F']), flush=True)

    # gap mid-plane on COMSOL's grid: x every 0.25 um, z every 1 um (here 0..h2, mirrored for the full depth)
    yg = wc2 + gap / 2
    xg = np.round(np.arange(0, Lx + 1e-9, 0.25), 9); zg = np.arange(0, h2 + 1e-9, 1.0)
    XX, ZZ = np.meshgrid(xg, zg, indexing='ij')
    ingap = fivegap.in_gap(XX, G)
    Vg = np.full((3, *XX.shape), np.nan); Pg = np.full(XX.shape, np.nan)
    xq = XX.copy()
    for k in range(G['n_gaps']):                          # the gap edges lie on the post faces: probe a hair inside
        x0 = G['L_post'] + k * (G['L_post'] + G['L_gap_x']); x1 = x0 + G['L_gap_x']
        xq = np.where(np.isclose(XX, x0), x0 + eps, np.where(np.isclose(XX, x1), x1 - eps, xq))
    pts = np.vstack([xq[ingap], np.full(ingap.sum(), yg), np.clip(ZZ[ingap], eps, h2 - eps)])
    Vg[:, ingap] = probe_u(pts) * scale_u
    Pg[ingap] = probe_p(pts) * meta['p_in_Pa']
    Vg[:, ingap & (np.isclose(ZZ, 0.0))] = 0.0
    q_gap = fivegap.gap_flows(xg, zg, Vg[1], G, depth_factor=2.0)
    out['Q_gap'] = q_gap; out['Q_gap_sum_over_Qs'] = sum(q_gap) / Q_s
    Qarr = Q_in - np.concatenate([[0], np.cumsum(q_gap)[:-1]])
    out['f_gap'] = [float(q / qa) for q, qa in zip(q_gap, Qarr)]
    write_plane(os.path.join(a.run_dir, 'gapplane.txt.gz'),
                [XX.ravel(), np.full(XX.size, yg), ZZ.ravel(), Vg[0].ravel(), Vg[1].ravel(), Vg[2].ravel(), Pg.ravel()],
                'x y z (um), u v w (m/s), p (Pa); half depth (z = 0 to %g, symmetric about z = %g); NaN in the solid' % (h2, h2))
    print('gap plane (%.0f s): sum of gap flows / Q_s %.5f' % (time.time() - t0, out['Q_gap_sum_over_Qs']), flush=True)

    # mid-depth plane on COMSOL's grid
    Ly = wc2 + gap + ws
    xm = np.arange(0, Lx + 1e-9, 0.5); ym = np.arange(0, Ly + 1e-9, 0.5)
    XX, YY = np.meshgrid(xm, ym, indexing='ij')
    solid = fivegap.in_solid(XX, YY, G)
    Vm = np.full((3, *XX.shape), np.nan); Pm = np.full(XX.shape, np.nan)
    ok = ~solid
    xq = np.clip(XX, eps, Lx - eps); yq = np.clip(YY, eps, Ly - eps)
    # points on the post faces: nudge into the fluid
    band = (YY > wc2 - 1e-9) & (YY < wc2 + gap + 1e-9)
    for k in range(G['n_gaps']):
        x0 = G['L_post'] + k * (G['L_post'] + G['L_gap_x']); x1 = x0 + G['L_gap_x']
        xq = np.where(band & np.isclose(XX, x0), x0 + eps, np.where(band & np.isclose(XX, x1), x1 - eps, xq))
    yq = np.where(np.isclose(YY, wc2) & ~fivegap.in_gap(XX, G), wc2 - eps, yq)
    yq = np.where(np.isclose(YY, wc2 + gap) & ~fivegap.in_gap(XX, G), wc2 + gap + eps, yq)
    pts = np.vstack([xq[ok], yq[ok], np.full(ok.sum(), h2 - eps)])
    Vm[:, ok] = probe_u(pts) * scale_u
    Pm[ok] = probe_p(pts) * meta['p_in_Pa']
    write_plane(os.path.join(a.run_dir, 'middepth_z70.txt.gz'),
                [XX.ravel(), YY.ravel(), np.full(XX.size, h2), Vm[0].ravel(), Vm[1].ravel(), Vm[2].ravel(), Pm.ravel()],
                'x y z (um), u v w (m/s), p (Pa) on the mid-depth symmetry plane; NaN in the solid')
    print('mid-depth plane (%.0f s)' % (time.time() - t0), flush=True)
    json.dump(out, open(os.path.join(a.run_dir, 'metrics.json'), 'w'), indent=1)
    print(json.dumps({k: v for k, v in out.items() if k != 'meta'}))


if __name__ == '__main__':
    main()
