#!/usr/bin/env python3
"""The figures of the two-code comparison (COMSOL and Python CFD), in the dissertation's style (pubstyle.py: no titles
or text inside the axes; axis labels, ticks, legends and panel letters only; the colors of Chapters 1 and 2, one color
per entity: V1 blue, V3 vermillion; runs B blue, A vermillion, C green). Everything is drawn from
cfd/results/comparison.json and the archived planes.

usage: python3 figures.py ARCHIVE_ROOT       writes cfd/figures/*.pdf and *.png
  F5-1c_two_code_convergence   F against the representative element size h = (V/N)^(1/3) for V1 (a) and V3 (b): COMSOL's
                               two meshes, Python CFD's three, and each code's extrapolation at h = 0 with its u_num
  F4-5d_skim_two_code          single gap: skim width at x = 75 um over the exact duct law, against the first cell at
                               the gap-side wall, for both codes and both inlets
  F4-1c_gap_passage            COMSOL's mid-depth speed in the gap passage of Run B (fine export), with the dividing
                               streamline traced on each code's mid-depth field (COMSOL's in vermillion, as in F4-1a and b)
  two_code_five_gap            (a) the center-to-side pressure difference per gap, (b) the relative differences of the
                               pressure difference and of the flow through each gap, V1 and V3
"""
import sys, os, json, gzip
import numpy as np
HERE = os.path.dirname(os.path.abspath(__file__)); sys.path.insert(0, HERE)
import pubstyle as ps                                                          # noqa: E402
import comsol_io as cio                                                        # noqa: E402
import matplotlib.pyplot as plt                                                # noqa: E402
import matplotlib.patheffects as pe                                            # noqa: E402
from matplotlib.lines import Line2D                                            # noqa: E402
from matplotlib.patches import Rectangle                                       # noqa: E402

ps.apply()
VOL_HALF = {'V1': (26 * 272 + 5 * 22 * 19 + 44 * 272) * 70.0, 'V3': (26 * 272 + 5 * 22 * 19 + 57 * 272) * 70.0}  # um^3
COL = ps.VCOL
SOLID = '#d9d9d9'                                                              # the posts, as in F4-1a and b


def fig_f51c(R, out):
    fig, axs = plt.subplots(2, 1, figsize=(ps.HALF_W, 4.2), sharex=True)
    for ax, v, L in zip(axs, ('V1', 'V3'), 'ab'):
        C = R['five_gap']['comsol']; g = R['five_gap']['grid']
        hc = [(2 * VOL_HALF[v] / C[k]['elements']) ** (1 / 3) for k in (v + '_extra_coarse', v + '_coarse')]
        Fc = [C[k]['F'] for k in (v + '_extra_coarse', v + '_coarse')]
        P = R['five_gap']['python'][v]
        hp = [(VOL_HALF[v] / r['n_tets']) ** (1 / 3) for r in P]; Fp = [r['F'] for r in P]
        tg = g['comsol'][v]; tp = g[v + '_python']['F']
        ax.plot(hp, Fp, 'o-', color=COL[v], mfc=COL[v])
        ax.plot(hc, Fc, 's--', color=COL[v], mfc='white')
        # V&V 20's u_num is the uncertainty of the finest solution, so the bars sit on the finest points
        ax.errorbar([hp[-1]], [Fp[-1]], yerr=[tp['u_num']], fmt='none', ecolor=COL[v], capsize=2, lw=0.8)
        ax.errorbar([hc[-1]], [Fc[-1]], yerr=[tg['u_num']], fmt='none', ecolor=COL[v], capsize=2, lw=0.8)
        ax.plot([0.0], [tp['extrapolated']], 'o', color=COL[v], mfc=COL[v], ms=3)
        ax.plot([0.0], [tg['extrapolated']], 's', color=COL[v], mfc='white', ms=5)
        ax.set_ylabel('filtrate fraction $F$ (%s)' % v)
        ax.set_xlim(-0.3, 5.9)
        ps.letter(ax, L)
    axs[-1].set_xlabel('representative element size $h = (V/N)^{1/3}$ (µm)')
    # the legend explains the marks of both panels (V1 blue, V3 vermillion), so its keys are drawn in black
    K = ps.BLACK
    eb = axs[0].errorbar([np.nan], [np.nan], yerr=[np.nan], fmt='none', ecolor=K, capsize=2, lw=0.8)
    h = [Line2D([], [], color=K, marker='o', ls='-', mfc=K), Line2D([], [], color=K, marker='s', ls='--', mfc='white'),
         Line2D([], [], color=K, marker='o', ls='none', mfc=K, ms=3),
         Line2D([], [], color=K, marker='s', ls='none', mfc='white', ms=5), eb]
    l = ['Python CFD, three meshes', 'COMSOL, two meshes', 'Python CFD, extrapolated (observed order)',
         'COMSOL, extrapolated (order 2)', '$\\pm u_\\mathrm{num}$ of the finest mesh']
    fig.legend(h, l, loc='upper center', bbox_to_anchor=(0.55, 0.0), ncol=1, fontsize=7)
    fig.subplots_adjust(hspace=0.12)
    ps.save(fig, os.path.join(out, 'F5-1c_two_code_convergence'))


def fig_f45d(R, out):
    fig, ax = plt.subplots(figsize=(ps.HALF_W + 0.6, 2.6))
    runs = R['single_gap']['runs']; cs = R['single_gap']['comsol']
    spec = [('B_uni', 'B_uniform', 1.40, 'B, uniform inlet', 0.1), ('B_dev', 'B_developed', 1.40, 'B, developed inlet', 0.1),
            ('A_dev', 'A_developed', 2.69, 'A, developed inlet', 0.05), ('C_dev', 'C_developed', 0.80, 'C, developed inlet', 0.1)]
    for fam, ck, alpha, lab, hw_c in spec:
        pts = sorted([(0.1 * r['scale'], r['skim']['python']['native_over_law']) for k, r in runs.items()
                      if k.startswith(fam + '_s')])
        if not pts:
            continue
        col = ps.ACOL[alpha]; ls, mk = ps.ASTY[alpha]
        fill = 'white' if fam == 'B_uni' else col
        ax.plot([p[0] for p in pts], [p[1] for p in pts], ls=ls, marker=mk, color=col, mfc=fill,
                label='%s, Python CFD' % lab)
        ax.plot([hw_c], [cs[ck]['native_over_law']], marker='D', ls='none', color=col, mfc=fill, ms=4.5,
                label='%s, COMSOL' % lab)
    ax.axhline(1.0, color=ps.OI['grey'], lw=0.8, ls=':', label='exact rectangular-duct law')
    ax.set_xlabel('first cell at the gap-side wall (µm)')
    ax.set_ylabel('skim width / duct law')
    ax.set_xlim(0.03, 0.22)
    ax.legend(loc='upper center', bbox_to_anchor=(0.45, -0.22), ncol=2, fontsize=7, columnspacing=0.8,
              handlelength=1.8)
    ps.save(fig, os.path.join(out, 'F4-5d_skim_two_code'))


def fig_f41c(root, R, out, py_run):
    P = cio.read_plane(os.path.join(root, 'comsol', 'single_gap', 'runB_5e-4_passage', 'RunB_5e-4_gap_passage_z70.txt'))
    xs, ys, g = cio.as_grid(P, 0, 1)
    U_in = 0.095238                                                   # m/s, Run B's U_in (its report)
    S = np.sqrt(g['u'] ** 2 + g['v'] ** 2) / U_in
    fig, ax = plt.subplots(figsize=(ps.HALF_W, 2.3))
    im = ax.pcolormesh(xs, ys, S.T, shading='auto', cmap=ps.CMAP, rasterized=True)
    U = np.nan_to_num(g['u']).T; V = np.nan_to_num(g['v']).T
    ax.streamplot(xs, ys, U, V, color='white', linewidth=0.4, density=1.4, arrowsize=0.5)
    # the posts (y > 50 um outside the gap, 100 < x < 120 um), filled as in F4-1a and b
    for xa, xb in ((xs[0] - 1.0, 100.0), (120.0, xs[-1] + 1.0)):
        ax.add_patch(Rectangle((xa, 50.0), xb - xa, ys[-1] + 1.0 - 50.0, fc=SOLID, ec=ps.BLACK, lw=0.6, zorder=3))
    # COMSOL's dividing streamline in vermillion, as in F4-1a and b; Python CFD's in blue, dashed, with a thin white
    # edge so that it shows on the dark end of the color map
    halo = [pe.Stroke(linewidth=2.2, foreground='white'), pe.Normal()]
    for nm, ls, col, lw, fx, lab in (('comsol', '-', ps.OI['vermillion'], 1.3, None, 'dividing streamline, COMSOL'),
                                     ('python', '--', ps.OI['blue'], 1.0, halo, 'dividing streamline, Python CFD')):
        fn = os.path.join(root, 'cfd', 'runs', py_run, 'dividing_streamline_%s.txt' % nm)
        if os.path.exists(fn):
            p = np.loadtxt(fn)
            # the dividing streamline ends at the stagnation point on the downstream post (x = 120 um); beyond it the
            # tracer's path is not a streamline, so it is not drawn
            end = np.nonzero(p[:, 0] >= 120.0 - 1e-3)[0]
            p = p[:end[0] + 1] if len(end) else p
            ax.plot(p[:, 0], p[:, 1], ls=ls, color=col, lw=lw, label=lab, zorder=4, path_effects=fx)
    ax.set_xlim(xs[0], xs[-1]); ax.set_ylim(ys[0], ys[-1]); ax.set_aspect('equal')
    ax.set_xlabel('$x$ (µm)'); ax.set_ylabel('$y$ (µm)')
    cb = fig.colorbar(im, ax=ax, fraction=0.035, pad=0.03); cb.set_label('speed / $U_\\mathrm{in}$')
    ax.legend(loc='upper center', bbox_to_anchor=(0.45, -0.3), ncol=1, fontsize=7)
    ps.save(fig, os.path.join(out, 'F4-1c_gap_passage'))


def fig_five(R, out):
    from matplotlib.ticker import FixedLocator, NullFormatter, FormatStrFormatter
    fig, axs = plt.subplots(1, 2, figsize=(ps.TEXT_W, 2.5))
    k = np.arange(1, 6)
    for v in ('V1', 'V3'):
        C = R['five_gap']['comsol'][v + '_coarse']; X = R['five_gap']['comsol'][v + '_extra_coarse']
        P = R['five_gap']['python'][v][-1]
        axs[0].semilogy(k, C['dp_gap'], 's--', color=COL[v], mfc='white', ms=5, label='%s, COMSOL' % v)
        axs[0].semilogy(k, P['dp_gap'], 'o-', color=COL[v], mfc=COL[v], ms=3, label='%s, Python CFD' % v)
        dp = 100 * (np.array(P['dp_gap']) / np.array(C['dp_gap']) - 1)
        dq = 100 * (np.array(P['Q_gap']) / np.array(C['Q_gap']) - 1)
        dm = 100 * (np.array(C['dp_gap']) / np.array(X['dp_gap']) - 1)
        axs[1].plot(k, dp, 'o-', color=COL[v], mfc=COL[v], ms=3.5, label='%s, pressure difference, Python CFD vs COMSOL' % v)
        axs[1].plot(k, dq, '^-', color=COL[v], mfc='white', ms=4, label='%s, gap flow, Python CFD vs COMSOL' % v)
        axs[1].plot(k, dm, ':', marker='s', color=COL[v], mfc='white', ms=3.5, lw=0.9,
                    label='%s, pressure difference, COMSOL coarse vs extra coarse' % v)
    axs[0].set_ylabel('center-to-side pressure difference (Pa)')
    axs[0].yaxis.set_major_locator(FixedLocator([1, 2, 5, 10, 20, 50]))
    axs[0].yaxis.set_major_formatter(FormatStrFormatter('%g')); axs[0].yaxis.set_minor_formatter(NullFormatter())
    axs[1].axhline(0, color=ps.OI['grey'], lw=0.6)
    axs[1].set_ylabel('relative difference (%)')
    for ax, L in zip(axs, 'ab'):
        ax.set_xlabel('gap'); ax.set_xticks(k); ax.xaxis.set_minor_locator(FixedLocator([])); ps.letter(ax, L)
    axs[0].legend(loc='upper right', fontsize=7, ncol=1)
    axs[1].legend(loc='upper center', bbox_to_anchor=(0.5, -0.24), ncol=1, fontsize=7)
    fig.subplots_adjust(wspace=0.32)
    ps.save(fig, os.path.join(out, 'two_code_five_gap'))


def main():
    root = sys.argv[1]
    R = json.load(open(os.path.join(root, 'cfd', 'results', 'comparison.json')))
    out = os.path.join(root, 'cfd', 'figures'); os.makedirs(out, exist_ok=True)
    fig_f51c(R, out); fig_f45d(R, out); fig_five(R, out)
    finest = [k for k in ('B_uni_s1', 'B_uni_s1.414') if k in R['single_gap']['runs']]
    if finest:
        fig_f41c(root, R, out, finest[0])
    print('figures written to', out)


if __name__ == '__main__':
    main()
