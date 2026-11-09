#!/usr/bin/env python3
"""The two-code comparison, from the archived files only: COMSOL's exported tables and planes (comsol/) and the second
code's run records and planes (cfd/runs/).

usage: python3 compare.py ARCHIVE_ROOT        writes cfd/results/comparison.json and cfd/results/COMPARISON.md

Numerical uncertainty follows ASME V&V 20-2009, Section 2-4:
  COMSOL, two meshes (unstructured, not geometrically similar): GCI = Fs |F1 - F2| / (r^p - 1) with Fs = 3, p = 2 and
      r = (N1/N2)^(1/3) from the element counts; u_num = GCI/1.15 (eq. 2-4-14)
  Python CFD, three meshes (structured tensor grids, refined by a scale factor but not geometrically similar):
      observed order p from eqs. 2-4-5 to 2-4-7 with h = N^(-1/3) (the representative cell size of a fixed volume);
      F_ext from eq. 2-4-8; GCI_fine = Fs |F2 - F1| / (r21^p - 1) with Fs = 1.25 (three grids) and, for comparison,
      Fs = 3; u_num = GCI/1.15
The difference between the codes is reported beside these, as supporting evidence; V&V 20 does not count it as u_num.
"""
import sys, os, json, gzip
import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
import comsol_io as cio                                                       # noqa: E402
import fivegap                                                                # noqa: E402

FIVE = {   # COMSOL five-gap models: folder, the five table files, mid-depth plane, gap plane, element count
    'V1_coarse': ('comsol/five_gap/V1_coarse', ['table%d.txt' % i for i in range(1, 6)],
                  'V1_coarse_middepth_z70.txt', 'V1_coarse_gapplane_y35p5.txt', 115904),
    'V1_extra_coarse': ('comsol/five_gap/V1_extra_coarse', ['Table%d.txt' % i for i in range(1, 6)],
                        'V1_xcoarse_middepth_z70.txt', 'V1_xcoarse_gapplane_y35p5.txt', 24783),
    'V3_coarse': ('comsol/five_gap/V3_coarse_resolved', ['Table%d_after.txt' % i for i in range(1, 6)],
                  'V3_coarse_middepth_z70_after.txt', 'V3_coarse_gapplane_y35p5_after.txt', 104705),
    'V3_coarse_stored': ('comsol/five_gap/V3_coarse_stored_May7', ['Table%d_before.txt' % i for i in range(1, 6)],
                         'V3_coarse_middepth_z70_before.txt', 'V3_coarse_gapplane_y35p5_before.txt', 104206),
    'V3_extra_coarse': ('comsol/five_gap/V3_extra_coarse', ['Table%d.txt' % i for i in range(1, 6)],
                        'V3_xcoarse_middepth_z70.txt', 'V3_xcoarse_gapplane_y35p5.txt', 21448),
}
PY5 = {'V1': ['V1_s3', 'V1_s2', 'V1_s1.7'], 'V3': ['V3_s3', 'V3_s2', 'V3_s1.7']}       # coarse to fine
DINH = {'V1': (0.505, 0.033), 'V3': (0.649, 0.036)}        # Dinh et al. 2024, Fig. 4: mean and SD of F
U_INPUT_OVER_UD = {'V1': 0.031, 'V3': 0.015}               # the dissertation's input uncertainty, as a fraction of u_D


def comsol_five(root, key):
    d, tabs, mid, gp, nel = FIVE[key]
    t = [cio.read_table(os.path.join(root, d, f)) for f in tabs]
    Q_in, Q_c, Q_s = t[0][0], t[1][0], t[2][0]
    G = fivegap.GEOM[key[:2]]
    r = dict(elements=nel, Q_in=Q_in, Q_c=Q_c, Q_s=Q_s, F=1 - Q_c / Q_in, imbalance=(Q_in - Q_c - Q_s) / Q_in,
             p_c=t[3], p_s=t[4], dp_gap=[a - b for a, b in zip(t[3], t[4])])
    P = cio.read_plane(os.path.join(root, d, gp))
    xs, zs, g = cio.as_grid(P, 0, 2)
    q = fivegap.gap_flows(xs, zs, g['v'], G)
    Qarr = Q_in - np.concatenate([[0], np.cumsum(q)[:-1]])
    r.update(Q_gap=q, Q_gap_sum_over_Qs=sum(q) / Q_s, f_gap=[float(a / b) for a, b in zip(q, Qarr)])
    M = cio.read_plane(os.path.join(root, d, mid))
    U = np.sqrt(M['u'] ** 2 + M['v'] ** 2 + M['w'] ** 2)
    r['middepth_w'] = dict(max_abs_w_over_max_speed=float(np.nanmax(np.abs(M['w'])) / np.nanmax(U)),
                           rms_w_over_rms_speed=float(np.sqrt(np.nanmean(M['w'] ** 2) / np.nanmean(U ** 2))))
    return r, M


def plane_diff(A, B, G):
    """In-plane velocity difference of two mid-depth planes on the same grid, over the points that are fluid in both
    and outside the wall."""
    xa, ya, ga = cio.as_grid(A, 0, 1); xb, yb, gb = cio.as_grid(B, 0, 1)
    assert np.allclose(xa, xb) and np.allclose(ya, yb)
    X, Y = np.meshgrid(xa, ya, indexing='ij')
    ok = ~np.isnan(ga['u']) & ~np.isnan(gb['u']) & ~fivegap.in_solid(X, Y, G)
    du = ga['u'][ok] - gb['u'][ok]; dv = ga['v'][ok] - gb['v'][ok]
    sb = np.sqrt(gb['u'][ok] ** 2 + gb['v'][ok] ** 2)
    return dict(points=int(ok.sum()), rel_l2=float(np.sqrt(np.mean(du ** 2 + dv ** 2) / np.mean(sb ** 2))),
                max_abs_over_max_speed=float(np.sqrt(du ** 2 + dv ** 2).max() / sb.max()))


def read_gz_plane(path, cols=7):
    with gzip.open(path, 'rt') as fh:
        M = np.loadtxt(fh, comments='%')
    out = dict(axes=None, x=M[:, 0], y=M[:, 1], z=M[:, 2], u=M[:, 3], v=M[:, 4], w=M[:, 5])
    if cols == 7:
        out['p'] = M[:, 6]
    return out


def two_grid(fine, coarse, Fs=3.0, p=2.0, k=1.15, key='F'):
    r = (fine['elements'] / coarse['elements']) ** (1 / 3)
    d = fine[key] - coarse[key]
    gci = Fs * abs(d) / (r ** p - 1)
    return dict(r=r, delta=d, GCI=gci, u_num=gci / k, factor_u_over_delta=Fs / (r ** p - 1) / k,
                extrapolated=fine[key] + d / (r ** p - 1))


def three_grid(N, phi, Fs=1.25, k=1.15):
    """N, phi ordered coarse -> fine; h = N^(-1/3).  V&V 20 eqs. 2-4-5 to 2-4-8 and 2-4-12; the observed order is the
    root of eq. 2-4-5 found by bisection on [0.05, 10] (the fixed-point iteration of the standard can land on a
    spurious root when r21 != r32)."""
    h = np.asarray(N, float) ** (-1 / 3)
    h3, h2, h1 = h; p3, p2, p1 = phi
    r21, r32 = h2 / h1, h3 / h2
    e21, e32 = p2 - p1, p3 - p2
    s = np.sign(e32 / e21)
    if s < 0:     # oscillatory convergence: no extrapolation; report the spread of the three values
        return dict(r21=r21, r32=r32, monotone=False, p=None, spread=float(max(phi) - min(phi)),
                    spread_rel=float((max(phi) - min(phi)) / abs(p1)))
    def g(pp):
        return pp - abs(np.log(abs(e32 / e21)) + np.log((r21 ** pp - s) / (r32 ** pp - s))) / np.log(r21)
    grid = np.linspace(0.05, 10, 4000); gv = np.array([g(x) for x in grid])
    idx = [i for i in range(len(grid) - 1) if np.sign(gv[i]) != np.sign(gv[i + 1])]
    roots = [grid[i] for i in idx]
    if not roots:
        return dict(r21=r21, r32=r32, monotone=bool(s > 0), p=None)
    from scipy.optimize import brentq
    pp = brentq(g, grid[idx[0]], grid[idx[0] + 1])
    ext = (r21 ** pp * p1 - p2) / (r21 ** pp - 1)
    gci = Fs * abs(e21) / (r21 ** pp - 1)
    gci3 = 3.0 * abs(e21) / (r21 ** pp - 1)
    return dict(r21=r21, r32=r32, monotone=bool(s > 0), p=pp, other_roots=[float(x) for x in roots[1:]],
                extrapolated=ext, GCI_fine=gci, u_num=gci / k, GCI_fine_Fs3=gci3, u_num_Fs3=gci3 / k)


def main():
    root = sys.argv[1]
    runs = os.path.join(root, 'cfd', 'runs')
    out = dict(five_gap=dict(comsol={}, python={}, grid={}, compare={}), single_gap={})
    planes = {}
    for key in FIVE:
        out['five_gap']['comsol'][key], planes[key] = comsol_five(root, key)
    C = out['five_gap']['comsol']
    out['five_gap']['grid']['comsol'] = {'V1': two_grid(C['V1_coarse'], C['V1_extra_coarse']),
                                         'V3': two_grid(C['V3_coarse'], C['V3_extra_coarse']),
                                         'V3_stored': two_grid(C['V3_coarse_stored'], C['V3_extra_coarse'])}
    out['five_gap']['grid']['comsol_Q_in'] = {
        'V1': two_grid(C['V1_coarse'], C['V1_extra_coarse'], key='Q_in'),
        'V3': two_grid(C['V3_coarse'], C['V3_extra_coarse'], key='Q_in')}
    out['five_gap']['mesh_realization_V3'] = dict(
        F_resolved=C['V3_coarse']['F'], F_stored_solution=C['V3_coarse_stored']['F'],
        F_preprint_tables=0.6280649, difference=C['V3_coarse']['F'] - C['V3_coarse_stored']['F'],
        field=plane_diff(planes['V3_coarse'], planes['V3_coarse_stored'], fivegap.GEOM['V3']))
    for v in ('V1', 'V3'):
        out['five_gap']['compare'].setdefault(v, {})['comsol_coarse_vs_extra_coarse_field'] = plane_diff(
            planes[v + '_coarse'], planes[v + '_extra_coarse'], fivegap.GEOM[v])
    for v, names in PY5.items():
        recs = []
        for nm in names:
            fn = os.path.join(runs, nm, 'metrics.json')
            if not os.path.exists(fn):
                recs = None; break
            M = json.load(open(fn)); recs.append(M)
        if not recs:
            continue
        out['five_gap']['python'][v] = [dict(run=nm, scale=M['meta']['scale'], n_tets=M['meta']['n_tets'],
                                             n_dofs=M['meta']['n_dofs'], minres_residual=M['meta']['residual'],
                                             seconds=M['meta']['seconds'], Q_in=M['Q_in'], Q_c=M['Q_c'], Q_s=M['Q_s'],
                                             F=M['F'], imbalance=M['imbalance'], dp_gap=M['dp_gap'], p_c=M['p_c'],
                                             p_s=M['p_s'], Q_gap=M['Q_gap'], f_gap=M['f_gap'],
                                             Q_gap_sum_over_Qs=M['Q_gap_sum_over_Qs'])
                                        for nm, M in zip(names, recs)]
        P = out['five_gap']['python'][v]
        N = [r['n_tets'] for r in P]
        g = dict(F=three_grid(N, [r['F'] for r in P]), Q_in=three_grid(N, [r['Q_in'] for r in P]),
                 F_by_scale_parameter=None)
        # the package's earlier definition (h = mesh scale), for the record
        hs = [r['scale'] for r in P]
        g['F_by_scale_parameter'] = three_grid([s ** -3 for s in hs], [r['F'] for r in P])
        out['five_gap']['grid'][v + '_python'] = g
        cc = C[v + '_coarse']; fine = P[-1]; tg = out['five_gap']['grid']['comsol'][v]
        cmp_ = out['five_gap']['compare'][v]
        cmp_.update(
            F_comsol=cc['F'], F_comsol_extrapolated=tg['extrapolated'], u_num_comsol=tg['u_num'],
            F_python_finest=fine['F'], F_python_extrapolated=g['F']['extrapolated'], u_num_python=g['F']['u_num'],
            dF_python_finest_minus_comsol=fine['F'] - cc['F'],
            dF_python_ext_minus_comsol=g['F']['extrapolated'] - cc['F'],
            dF_extrapolations=g['F']['extrapolated'] - tg['extrapolated'],
            Q_in_ratio_python_finest=fine['Q_in'] / cc['Q_in'],
            Q_in_comsol_extrapolated=out['five_gap']['grid']['comsol_Q_in'][v]['extrapolated'],
            Q_in_python_extrapolated=g['Q_in'].get('extrapolated'),
            Q_in_python_finest=fine['Q_in'],
            dp_rel_diff=[a / b - 1 for a, b in zip(fine['dp_gap'], cc['dp_gap'])],
            Q_gap_rel_diff=[a / b - 1 for a, b in zip(fine['Q_gap'], cc['Q_gap'])],
            f_gap_python=fine['f_gap'], f_gap_comsol=cc['f_gap'])
        mid = os.path.join(runs, PY5[v][-1], 'middepth_z70.txt.gz')
        if os.path.exists(mid):
            cmp_['field_python_finest_vs_comsol_coarse'] = plane_diff(read_gz_plane(mid), planes[v + '_coarse'],
                                                                      fivegap.GEOM[v])
    # V&V 20 validation comparison with Dinh et al. (COMSOL's coarse F and its two-grid u_num)
    for v in ('V1', 'V3'):
        cc = C[v + '_coarse']; tg = out['five_gap']['grid']['comsol'][v]
        mu, sd = DINH[v]
        u_in = U_INPUT_OVER_UD[v] * sd
        un = tg['u_num']
        u_val = np.sqrt(un ** 2 + u_in ** 2 + sd ** 2)
        out['five_gap']['compare'].setdefault(v, {})['validation_comsol'] = dict(
            F=cc['F'], E=cc['F'] - mu, u_D=sd, u_input=u_in, u_num=un, u_val=u_val, E_over_uD=(cc['F'] - mu) / sd,
            E_over_uval=(cc['F'] - mu) / u_val, uval_excess_pct=(u_val / sd - 1) * 100, uD_over_unum=sd / un,
            E_over_sem={n: (cc['F'] - mu) / (sd / np.sqrt(n)) for n in (3, 4)})
    # pressure-difference half-lives (gaps per halving), fitted on all five probes and on the four upstream ones
    def half_lives(dp):
        dp = np.asarray(dp); n = np.arange(5)
        return dict(all5=float(np.log(2) / -np.polyfit(n, np.log(dp), 1)[0]),
                    first4=float(np.log(2) / -np.polyfit(n[:4], np.log(dp[:4]), 1)[0]),
                    step_ratios=(dp[1:] / dp[:-1]).tolist())
    out['five_gap']['half_lives'] = {k: half_lives(v['dp_gap']) for k, v in C.items()}
    for v in PY5:
        if v in out['five_gap']['python']:
            out['five_gap']['half_lives'][v + '_python_finest'] = half_lives(out['five_gap']['python'][v][-1]['dp_gap'])
    # COMSOL single gap at fraction 5e-4: the campaign's exports and the November developed-inlet solve, measured with
    # the same routines as Python CFD (skim.py)
    import skim
    camp = os.path.join(root, 'comsol', 'single_gap', 'campaign')
    SGC = {'B_uniform': (os.path.join(camp, 'Run B Sweep', '5e-4', 'upstream_grid_frac5e-4 v3.txt'), 100., 140.,
                         4.962120761577921e-4),
           'B_developed': (os.path.join(root, 'comsol', 'single_gap', 'runB_5e-4_developed',
                                        'upstream_grid_frac5e-4 v3.txt'), 100., 140., 4.8751210307083063e-4),
           'A_developed': (os.path.join(camp, 'Run A Sweep', '5e-4', 'upstream_official_5e-4.txt'), 52., 140., 4.8518e-4),
           'C_developed': (os.path.join(camp, 'Run C Sweep', '5e-4', 'upstream_official_5e-4.txt'), 100., 80., 4.8550e-4)}
    sgc = {}
    for k, (fn, wc, h, fq) in SGC.items():
        P = cio.read_plane(fn); ys, zs, g = cio.as_grid(P, 1, 2)
        dn, _ = skim.skim_native(ys, zs, g['u'], wc / 2, fq)
        dr = skim.skim_resampled(P['y'], P['z'], P['u'], wc / 2, h / 2, fq)
        law = skim.duct_delta(fq / 2, wc, h)
        sgc[k] = dict(f_q=fq, native=dn, resampled=dr, duct_law=law, native_over_law=dn / law,
                      resampled_over_law=dr / law)
    sgc['inlet_effect_B_native'] = sgc['B_developed']['native_over_law'] / sgc['B_uniform']['native_over_law'] - 1
    sgc['inlet_effect_B_resampled'] = (sgc['B_developed']['resampled_over_law'] /
                                       sgc['B_uniform']['resampled_over_law'] - 1)
    # the upstream plane against the exact developed duct profile (developed inlets: an exact Stokes solution of the
    # straight duct, perturbed only by the gap 25 um downstream), each field over its own plane mean
    from stokes_gap import duct_profile
    def profile_err(ys, zs, U, wc, h):
        Y, Z = np.meshgrid(ys, zs, indexing='ij')
        ua = duct_profile(Y, Z, wc / 2, h / 2)
        Q = np.trapezoid(np.trapezoid(np.nan_to_num(U), zs, axis=1), ys) / ((wc / 2) * (h / 2))
        d = np.nan_to_num(U) / Q - ua
        near = (wc / 2 - Y) <= 0.1 * wc                  # within 0.1 w_c of the gap-side wall, where the skim lies
        return dict(rel_l2=float(np.sqrt((d ** 2).mean() / (ua ** 2).mean())), max_abs_over_peak=float(np.abs(d).max() / ua.max()),
                    rel_l2_near_wall=float(np.sqrt((d[near] ** 2).mean() / (ua[near] ** 2).mean())))
    for k, (fn, wc, h, fq) in SGC.items():
        P = cio.read_plane(fn); ys, zs, g = cio.as_grid(P, 1, 2)
        sgc[k]['profile_vs_duct'] = profile_err(ys, zs, g['u'], wc, h)
    out['single_gap']['comsol'] = sgc
    # single gap
    sg = {}
    for d in sorted(os.listdir(runs)):
        fn = os.path.join(runs, d, 'metrics.json')
        if d[0] in 'ABC' and os.path.exists(fn):
            M = json.load(open(fn))
            sg[d] = dict(run=M['meta']['run'], inlet=M['meta']['inlet'], scale=M['meta']['scale'],
                         n_tets=M['meta']['n_tets'], n_dofs=M['meta']['n_dofs'], residual=M['meta']['residual'],
                         f_q=M['mass']['f_q_achieved'], imbalance=M['mass']['imbalance'], skim=M['skim'],
                         middepth=M['middepth'])
            up = os.path.join(runs, d, 'upstream_x75.txt.gz')
            if os.path.exists(up):
                Pu = read_gz_plane(up, cols=6); ys, zs, g = cio.as_grid(Pu, 1, 2)
                Gm = M['meta']['G']
                sg[d]['profile_vs_duct'] = profile_err(ys, zs, g['u'], Gm['w_c'], Gm['depth'])
    out['single_gap']['runs'] = sg
    ie = {}
    for s in ('2', '1.414', '1'):
        u_, d_ = sg.get('B_uni_s' + s), sg.get('B_dev_s' + s)
        if u_ and d_:
            ie[s] = d_['skim']['python']['native_over_law'] / u_['skim']['python']['native_over_law'] - 1
    out['single_gap']['inlet_effect_B_python'] = ie
    conv = {}
    for fam in ('B_uni', 'B_dev', 'A_dev', 'C_dev'):
        rs = [sg.get('%s_s%s' % (fam, s)) for s in ('2', '1.414', '1')]
        if all(rs):
            conv[fam] = dict(delta=three_grid([r['n_tets'] for r in rs], [r['skim']['python']['native'] for r in rs]),
                             trace=three_grid([r['n_tets'] for r in rs],
                                              [r['middepth']['trace_python']['delta_at_x_report'] for r in rs]))
    out['single_gap']['grid_python'] = conv
    os.makedirs(os.path.join(root, 'cfd', 'results'), exist_ok=True)
    json.dump(out, open(os.path.join(root, 'cfd', 'results', 'comparison.json'), 'w'), indent=1, default=float)
    print('written', os.path.join(root, 'cfd', 'results', 'comparison.json'))


if __name__ == '__main__':
    main()
