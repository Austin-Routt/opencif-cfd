"""figures_kb.py: dissertation figures computed from the frozen appendices and the knowledge-base data alone, in the
journal style of tools/pubstyle.py.

Usage: python figures_kb.py OUTLINE_DIR

Writes OUTLINE_DIR/figures/<id>_<name>.{pdf,png} and OUTLINE_DIR/figures/FIGURES_KB_CHECK.txt, the numbers each figure
rests on, recomputed here and compared with the values that Appendix A, Appendix B, the glossary and the preprint print.

  F3-3  skim_cubic      the skim cubic, its square law and the DLD first-streamline cubic (App. A R1, R3, R4), and the
                        error of the square law over the range of the calibration data
  F3-4  kappa_app       the apparent prefactor of each calibration point against u = D/w_c (Gifford et al.'s Fig. 3 as
                        digitized by the author, kb/kb_data.json), with the plateau, its interval, the plateau edge, the foot
                        and the steric value
  F3-5  depth_factor    K_depth(alpha) = sqrt(3/c2(alpha)) (App. A R9, eq:c2), with the single-gap Stokes values of
                        Runs C, B and A and the three calibration widths at the assumed depth
  F3-6  row_count       N against D at fixed F for the canonical taper (App. A R13(b)), with the region below the plateau
                        and, for each F, the cutoffs whose outlet rows pass the calibration cap (lighter line)
  F3-7  duct_field      the rectangular-duct velocity (App. A eq:duct) over two cross-sections, alpha = 0.80 and 2.69,
                        and the depth-averaged velocity across the width against the parabola of the same mean
  F5-1  outlet_slab     the 5-gap outlet model of the COMSOL reports (kb/comsol_5gap_reports.json) and the decay of the
                        center-to-side pressure difference along it
"""
import json
import math
import os
import sys

import numpy as np

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import pubstyle as P                                   # noqa: E402
import matplotlib.pyplot as plt                        # noqa: E402
from matplotlib.patches import Rectangle, FancyArrowPatch   # noqa: E402

LOG = []
C1, C2, C3 = P.OI['blue'], P.OI['vermillion'], P.OI['green']
SHADE = '#e6e6e6'


def note(line):
    LOG.append(line)
    print(line)


def check(label, got, want, tol):
    ok = abs(got - want) <= tol
    note('%-58s computed %-12.6g printed %-10s %s' % (label, got, want, 'ok' if ok else 'MISMATCH'))
    return ok


def save(fig, out, stem):
    P.save(fig, os.path.join(out, stem))
    note('wrote %s.{pdf,png}' % stem)


# ---------------------------------------------------------------- Appendix A, Result 9 (the duct field)
ODD = np.arange(1, 4_000_002, 2, dtype=float)       # two million odd terms; the tail beyond is below 1e-7 relative


def c_m(alpha, m):
    return alpha - 2.0 / (m * math.pi) * np.tanh(m * math.pi * alpha / 2.0)


def c2(alpha):
    c = c_m(alpha, ODD)
    return math.pi ** 2 / 4.0 * np.sum(c / ODD ** 2) / np.sum(c / ODD ** 4)


def k_depth(alpha):
    return math.sqrt(3.0 / c2(alpha))


def duct_u(w, h, ny=401, nz=401, terms=399):
    """App. A eq:duct on -w/2 <= y <= w/2, 0 <= z <= h, divided by its mean over the section."""
    y = np.linspace(-w / 2.0, w / 2.0, ny)
    z = np.linspace(0.0, h, nz)
    Y, Z = np.meshgrid(y + w / 2.0, z - h / 2.0)
    u = np.zeros_like(Y)
    for m in range(1, terms + 1, 2):
        k = m * math.pi / w
        ratio = np.exp(k * (np.abs(Z) - h / 2.0)) * (1.0 + np.exp(-2.0 * k * np.abs(Z))) / (1.0 + math.exp(-k * h))
        u += (1.0 / m ** 3) * (1.0 - ratio) * np.sin(k * Y)
    return y, z, u / u.mean()


def depth_avg(alpha, x, terms=20001):
    """Depth-averaged velocity across the width, divided by its mean (App. A eq:cm with eq:duct)."""
    m = np.arange(1, terms + 1, 2, dtype=float)
    c = c_m(alpha, m)
    ub = np.sum((c / m ** 3)[None, :] * np.sin(np.outer(x, m) * math.pi), axis=1)
    mean = np.sum(c / m ** 3 * 2.0 / (m * math.pi))
    return ub / mean


# ---------------------------------------------------------------- F3-3
def fig_cubic(out, kb):
    fstar = np.array(kb['FIG3'])[:, 1] * 1e-4
    lo, hi = fstar.min(), fstar.max()
    xlo, xhi = math.sqrt(lo / 3.0), math.sqrt(hi / 3.0)
    note('F3-3 calibration range of f* %.3g to %.3g, so x = sqrt(f/3) from %.4f to %.4f' % (lo, hi, xlo, xhi))
    note('F3-3 relative error of 3x^2 at the top of that range: %.3f%% (2x/(3 - 2x))' % (100 * 2 * xhi / (3 - 2 * xhi)))
    check('F3-3 cubic at x = 1 (all the flow)', 3 - 2, 1.0, 0)
    check('F3-3 cubic at x = 1/2 (half the flow, by symmetry)', 3 * 0.25 - 2 * 0.125, 0.5, 0)
    P.apply()
    fig = plt.figure(figsize=(P.TEXT_W, 2.55))
    ax = fig.add_axes([0.075, 0.17, 0.38, 0.77])
    axe = fig.add_axes([0.6, 0.17, 0.38, 0.77])
    x = np.linspace(0, 1, 401)
    ax.plot(x, 3 * x ** 2 - 2 * x ** 3, color='k', lw=1.1, label=r'CIF skim cubic, $3x^2 - 2x^3$')
    ax.plot(x, 3 * x ** 2, color=C1, lw=1.0, ls='--', label=r'square law, $3x^2$')
    yy = np.linspace(0.05, 0.95, 10)
    ax.plot(yy, 3 * yy ** 2 - 2 * yy ** 3, 'o', ms=3.8, mfc='white', mec=C2, mew=0.9,
            label=r'DLD cubic, $3y^2 - 2y^3$')
    ax.set_xlim(0, 1)
    ax.set_ylim(0, 1.05)
    ax.set_xlabel(r'lamina width, $x = \delta/w_\mathrm{c}$ (DLD: $y = b/g$)')
    ax.set_ylabel('fraction of the channel flow')
    ax.legend(loc='lower right', handlelength=1.8)
    P.letter(ax, 'a')
    xs = np.geomspace(1e-3, 0.5, 200)
    axe.set_xscale('log')
    axe.set_yscale('log')
    axe.axvspan(xlo, xhi, color=SHADE, lw=0)
    axe.plot(xs, 100 * 2 * xs / (3 - 2 * xs), color='k', lw=1.1)
    axe.set_xlim(1e-3, 0.5)
    axe.set_ylim(0.05, 100)
    axe.set_xlabel(r'lamina width, $x$')
    axe.set_ylabel(r'error of $3x^2$ (% of $f$)')
    P.letter(axe, 'b')
    save(fig, out, 'F3-3_skim_cubic')


# ---------------------------------------------------------------- F3-4
def fig_kappa(out, kb):
    a = np.array(kb['FIG3'])
    D, f, s, w = a.T
    f, s = f * 1e-4, s * 1e-4
    u = D / w
    k = D / (w * np.sqrt(f / 3.0))
    sk = 0.5 * k * s / f
    pl, ft = u >= 0.035, u < 0.012
    kp = k[pl].mean()
    check('F3-4 plateau mean of kappa_app (u >= 0.035, n = %d)' % pl.sum(), kp, 4.9137, 5e-5)
    check('F3-4 plateau CV (%)', 100 * k[pl].std(ddof=1) / kp, 4.0, 0.05)
    check('F3-4 foot mean (u < 0.012, n = %d)' % ft.sum(), k[ft].mean(), 1.34, 0.005)
    P.apply()
    fig = plt.figure(figsize=(P.TEXT_W, 2.9))
    ax = fig.add_axes([0.08, 0.15, 0.9, 0.81])
    ax.fill_between([0.035, 0.0705], 4.75, 5.08, color=SHADE, lw=0, zorder=0)
    ax.plot([0.035, 0.0705], [4.9138] * 2, color='k', lw=0.9, zorder=1)
    ax.axvline(0.035, color='k', lw=0.6, ls='--', zorder=1)
    ax.axhline(2.0, color='k', lw=0.6, ls=':', zorder=1)
    ax.plot([0.004, 0.012], [1.344] * 2, color='k', lw=0.9, zorder=1)
    for wc, col, mk, fill in ((100, C1, 'o', True), (125, C2, 's', False), (150, C3, '^', True)):
        sel = w == wc
        o = np.argsort(u[sel])
        ax.errorbar(u[sel][o], k[sel][o], yerr=sk[sel][o], fmt=mk, ms=4.2, color=col,
                    mfc=col if fill else 'white', mec=col, mew=0.9, ecolor=col, elinewidth=0.7, capsize=1.6,
                    capthick=0.7, lw=0, zorder=3, label=r'$w_\mathrm{c}$ = %d μm' % wc)
    ax.set_xlim(0.004, 0.0705)
    ax.set_ylim(0.0, 6.0)
    ax.set_xlabel(r'cutoff over center width, $u = D/w_\mathrm{c}$')
    ax.set_ylabel(r'apparent prefactor, $\kappa_\mathrm{app}$')
    ax.legend(loc='upper left')
    save(fig, out, 'F3-4_kappa_app')


# ---------------------------------------------------------------- F3-5
def fig_depth(out, kb):
    for al, c2p, kp in ((0.80, 3.8054, 0.8879), (1.40, 3.3240, 0.9500), (2.69, 3.1273, 0.9794)):
        check('F3-5 c2(%.2f)' % al, c2(al), c2p, 6e-5)
        check('F3-5 K_depth(%.2f)' % al, k_depth(al), kp, 6e-5)
    for al, kp in ((0.64, 0.849), (0.8 / 1.5, 0.811)):
        check('F3-5 K_depth(%.4f) (calibration width)' % al, k_depth(al), kp, 6e-4)
    cfd = {r: (v[2], v[3]) for r, v in kb['CFD_RUNS'].items()}
    printed = {'C': 0.891, 'B': 0.945, 'A': 0.981}
    for r in 'CBA':
        note('F3-5 run %s alpha %.2f R_depth per row %s (glossary %.3f)' % (r, cfd[r][0], [x[4] for x in cfd[r][1]],
                                                                           printed[r]))
    worst = max(abs(k_depth(cfd[r][0]) - printed[r]) / printed[r] for r in 'CBA')
    check('F3-5 largest |K_th - K_cfd| / K_cfd (%)', 100 * worst, 0.5301, 0.002)
    al = np.geomspace(0.3, 10, 160)
    K = np.array([k_depth(x) for x in al])
    P.apply()
    fig = plt.figure(figsize=(P.HALF_W * 1.45, 2.8))
    ax = fig.add_axes([0.15, 0.16, 0.82, 0.8])
    ax.set_xscale('log')
    ax.axhline(1.0, color='k', lw=0.6, ls=':')
    ax.plot(al, K, color='k', lw=1.1, label=r'analytic, $\sqrt{3/c_2(\alpha)}$')
    cal = ((150, 80 / 150), (125, 0.64), (100, 0.80))
    ax.plot([c[1] for c in cal], [k_depth(c[1]) for c in cal], 'o', ms=4.5, mfc='white', mec='0.35', mew=1.0,
            label='calibration widths')
    for r in 'CBA':
        a_ = round(cfd[r][0], 2)
        ax.plot([cfd[r][0]], [printed[r]], P.ASTY[a_][1], ms=4.5, color=P.ACOL[a_], mec='white', mew=0.5,
                label=r'Stokes solve, run %s ($\alpha$ = %.2f)' % (r, a_))
    ax.set_xlim(0.3, 10)
    ax.set_ylim(0.6, 1.03)
    ax.set_xticks([0.3, 0.5, 1, 2, 3, 5, 10])
    ax.set_xticklabels(['0.3', '0.5', '1', '2', '3', '5', '10'])
    ax.set_xlabel(r'aspect ratio, $\alpha = h/w_\mathrm{c}$')
    ax.set_ylabel(r'depth factor, $K_\mathrm{depth}$')
    ax.legend(loc='lower right')
    save(fig, out, 'F3-5_depth_factor')


# ---------------------------------------------------------------- F3-6
def fig_rows(out):
    kap, w0 = 4.9138, 100.0
    D = np.geomspace(1.0, 20.0, 120)
    P.apply()
    fig = plt.figure(figsize=(P.HALF_W * 1.45, 2.8))
    ax = fig.add_axes([0.15, 0.16, 0.82, 0.8])
    ax.set_xscale('log')
    ax.set_yscale('log')
    ax.axvspan(1.0, 0.035 * w0, color=SHADE, lw=0)
    ucap = kap * math.sqrt(5.76e-4 / 3.0)            # the cap of f_gap, reached first at the narrowest (outlet) row
    for F in (0.3, 0.5, 0.7):
        Dcap = ucap * w0 * math.sqrt(1.0 - F)        # the canonical taper's outlet width is w0 sqrt(1 - F)
        Nf = kap ** 2 * F * w0 ** 2 / (6 * D ** 2)
        inside = D <= Dcap
        ax.plot(D[inside], Nf[inside], color=P.FCOL[F], lw=1.1, ls=P.FSTY[F][0], label='$F$ = %.1f' % F)
        ax.plot(D[D >= D[inside][-1]], Nf[D >= D[inside][-1]], color=P.FCOL[F], lw=0.9, ls=P.FSTY[F][0], alpha=0.35)
        ax.plot([Dcap], [kap ** 2 * F * w0 ** 2 / (6 * Dcap ** 2)], 'o', ms=4.5, mfc='white', mec=P.FCOL[F], mew=1.0,
                zorder=5)
        note('F3-6 F = %.1f: the outlet rows pass the cap f_gap 5.76e-4 for D above %.3f um (lighter line beyond)'
             % (F, Dcap))
    Dm = np.array([4.0, 5.0, 7.0, 10.0, 14.0, 20.0])
    xm = Dm / (kap * w0)
    Nm = np.ceil(np.log(1 - 0.5) / np.log(1 - 2 * (3 * xm ** 2 - 2 * xm ** 3)))
    r = (kap ** 2 * 0.5 * w0 ** 2 / 6) / Nm / Dm ** 2
    note('F3-6 constant width over canonical taper at F = 0.5, N D^2 ratio: %s (closed forms: -ln(0.5)/0.5 = %.4f)'
         % (np.round(1 / r, 3).tolist(), -math.log(0.5) / 0.5))
    check('F3-6 halving D (10 to 5 um) multiplies N by', (kap ** 2 * 0.5 * w0 ** 2 / 150) /
          (kap ** 2 * 0.5 * w0 ** 2 / 600), 4.0, 1e-12)
    ax.set_xlim(1.0, 20.0)
    ax.set_ylim(20, 4e4)
    ax.set_xticks([1, 2, 5, 10, 20])
    ax.set_xticklabels(['1', '2', '5', '10', '20'])
    ax.set_xlabel(r'cutoff, $D$ (μm)')
    ax.set_ylabel(r'number of gap rows, $N$')
    ax.legend(loc='upper right')
    save(fig, out, 'F3-6_row_count')


# ---------------------------------------------------------------- F3-7
def fig_duct(out):
    P.apply()
    fig = plt.figure(figsize=(P.TEXT_W, 4.9))
    axa = fig.add_axes([0.075, 0.555, 0.40, 0.41])
    axb = fig.add_axes([0.55, 0.555, 0.164, 0.41])
    cax = fig.add_axes([0.73, 0.555, 0.012, 0.41])
    axc = fig.add_axes([0.075, 0.085, 0.40, 0.37])
    axd = fig.add_axes([0.585, 0.085, 0.40, 0.37])
    vmax = 0.0
    fields = []
    for w, h in ((100.0, 80.0), (52.0, 140.0)):
        y, z, u = duct_u(w, h)
        fields.append((y, z, u))
        vmax = max(vmax, float(u.max()))
    for ax, (y, z, u), lab in ((axa, fields[0], 'a'), (axb, fields[1], 'b')):
        im = ax.pcolormesh(y, z, u, cmap=P.CMAP, shading='gouraud', vmin=0.0, vmax=vmax, rasterized=True)
        ax.contour(y, z, u, levels=[0.5, 1.0, 1.5, 2.0], colors='white', linewidths=0.6)
        ax.set_aspect('equal')
        ax.set_xlabel(r'$y$ (μm)')
        P.letter(ax, lab)
    axa.set_ylabel(r'$z$ (μm)')
    note('F3-7 peak u/u_mean: alpha 0.80 %.4f, alpha 2.69 %.4f' % (fields[0][2].max(), fields[1][2].max()))
    cb = fig.colorbar(im, cax=cax)
    cb.set_label(r'$u/\bar{u}$')
    cb.ax.minorticks_off()
    cb.outline.set_linewidth(0.6)
    x = np.linspace(0, 1, 401)
    par = 6 * x * (1 - x)
    axc.plot(x, par, color=P.OI['grey'], lw=0.9, ls=':', label='parabola')
    axd.axhline(0.0, color=P.OI['grey'], lw=0.9, ls=':')
    for al in (0.80, 1.40, 2.69):
        col, ls = P.ACOL[al], P.ASTY[al][0]
        ub = depth_avg(al, x)
        axc.plot(x, ub, color=col, lw=1.1, ls=ls, label=r'$\alpha$ = %.2f' % al)
        axd.plot(x, ub - par, color=col, lw=1.1, ls=ls)
        m = np.arange(1, 20002, 2, dtype=float)
        cm_ = c_m(al, m)
        slope0 = np.sum(cm_ / m ** 3 * m * math.pi) / np.sum(cm_ / m ** 3 * 2.0 / (m * math.pi))
        check('F3-7 wall slope of the depth average, c2(%.2f) = ub\'(0)/2' % al, slope0 / 2.0, c2(al), 2e-4)
    axc.set_xlim(0, 1)
    axc.set_ylim(0, 1.6)
    axc.set_xlabel(r'$y/w_\mathrm{c} + 1/2$')
    axc.set_ylabel(r'depth-averaged $u/\bar{u}$')
    axc.legend(loc='lower center', ncol=2)
    P.letter(axc, 'c')
    axd.set_xlim(0, 1)
    axd.set_ylim(-0.1, 0.1)
    axd.set_xlabel(r'$y/w_\mathrm{c} + 1/2$')
    axd.set_ylabel('difference from the parabola')
    P.letter(axd, 'd')
    save(fig, out, 'F3-7_duct_field')


# ---------------------------------------------------------------- F5-1
def fig_slab(out, cs):
    p = cs['common_parameters_um']
    fits = {}
    for v in ('V1', 'V3'):
        r = cs[v]
        F = r['Q_s_out_top_m3s'] / r['Q_in_m3s']
        dp = np.array(r['p_center_Pa']) - np.array(r['p_side_Pa'])
        i = np.arange(1, 6)
        sl, b = np.polyfit(i, np.log(dp), 1)
        yh = sl * i + b
        R2 = 1 - ((np.log(dp) - yh) ** 2).sum() / ((np.log(dp) - np.log(dp).mean()) ** 2).sum()
        fits[v] = (dp, sl, b, R2, F)
        check('F5-1 %s F = Q_side / Q_in' % v, F, {'V1': 0.5070, 'V3': 0.6281}[v], 5e-5)
        check('F5-1 %s pressure half-life (gaps)' % v, math.log(2) / -sl, {'V1': 0.944, 'V3': 1.041}[v], 5e-4)
        check('F5-1 %s R^2 of the semi-log fit' % v, R2, {'V1': 0.987, 'V3': 0.986}[v], 5e-4)
        note('F5-1 %s outlet sum Q_c + Q_s - Q_in = %.3g m3/s' % (v, r['Q_c_out_m3s'] + r['Q_s_out_top_m3s'] -
                                                                r['Q_in_m3s']))
    P.apply()
    fig = plt.figure(figsize=(P.TEXT_W, 3.9))
    ax = fig.add_axes([0.06, 0.5, 0.9, 0.48])
    ax2 = fig.add_axes([0.1, 0.1, 0.86, 0.33])
    wc2, g, ws = p['w_c'] / 2, p['gap'], cs['V1']['w_s_um']
    L, Lp, Lg, n = p['L_total'], p['L_post'], p['L_gap_x'], p['n_gaps']
    top = wc2 + g + ws
    ax.plot([0, L], [0, 0], color='k', lw=0.6, ls=(0, (4, 3)))
    ax.plot([0, L], [top, top], color='k', lw=1.0)
    ax.plot([0, 0], [wc2 + g, top], color='k', lw=1.0)
    for kk in range(n + 1):
        ax.add_patch(Rectangle((kk * (Lp + Lg), wc2), Lp, g, fc='#bdbdbd', ec='k', lw=0.8))
    for yy, lab in ((wc2 / 2, 'center outlet, $p$ = 0'), (wc2 + g + ws / 2, 'side outlet, $p$ = 0')):
        ax.add_patch(FancyArrowPatch((L + 4, yy), (L + 34, yy), arrowstyle='-|>', mutation_scale=7, lw=0.8, color='k'))
        ax.text(L + 38, yy, lab, fontsize=7, va='center')
    ax.add_patch(FancyArrowPatch((-34, wc2 / 2), (-4, wc2 / 2), arrowstyle='-|>', mutation_scale=7, lw=0.8, color='k'))
    ax.text(-38, wc2 / 2, r'inlet, $p_\mathrm{in}$', fontsize=7, va='center', ha='right')
    ax.text(L / 2, -3.5, 'symmetry plane', fontsize=7, ha='center', va='top')
    for xs in cs['V1']['x_um']:
        ax.plot([xs], [wc2 * 0.35], marker='x', ms=3.5, mew=0.8, color='k')
        ax.plot([xs], [cs['V1']['side_y_um']], marker='x', ms=3.5, mew=0.8, color='k')

    def dim(x0, x1, yv, txt):
        ax.add_patch(FancyArrowPatch((x0, yv), (x1, yv), arrowstyle='<|-|>', mutation_scale=4, lw=0.5, color='k'))
        ax.text(0.5 * (x0 + x1), yv + 1.5, txt, fontsize=7, ha='center', va='bottom')
    dim(Lp, Lp + Lg, wc2 + g + 4.5, '%d' % Lg)
    dim(0, Lp, wc2 + g + 4.5, '%d' % Lp)
    xd = 112.0
    ax.add_patch(FancyArrowPatch((xd, 0), (xd, wc2), arrowstyle='<|-|>', mutation_scale=4, lw=0.5, color='k'))
    ax.text(xd + 2, wc2 / 2, '%d' % wc2, fontsize=7, ha='left', va='center')
    ax.add_patch(FancyArrowPatch((xd, wc2 + g), (xd, top), arrowstyle='<|-|>', mutation_scale=4, lw=0.5,
                                 color='k'))
    ax.text(xd + 2, wc2 + g + 0.8 * ws, '%d' % ws, fontsize=7, ha='left', va='center')
    ax.set_xlim(-80, L + 120)
    ax.set_ylim(-10, top + 6)
    ax.set_aspect('equal')
    ax.axis('off')
    P.letter(ax, 'a', dx=0.02, dy=0.98)
    i = np.arange(1, 6)
    for v, mk, ls in (('V1', 'o', '-'), ('V3', 's', '--')):
        col = P.VCOL[v]
        dp, sl, b, R2, F = fits[v]
        xx = np.linspace(1, 5, 50)
        ax2.plot(xx, np.exp(sl * xx + b) / dp[0], color=col, lw=0.9, ls=ls)
        ax2.plot(i, dp / dp[0], mk, ms=4.2, mfc='white' if v == 'V3' else col, mec=col, mew=0.9,
                 label=v)
        note('F5-1 %s: half-life %.3f gaps, R^2 %.3f (caption)' % (v, math.log(2) / -sl, R2))
    ax2.set_yscale('log')
    ax2.set_xlim(0.7, 5.3)
    ax2.set_xticks(i)
    ax2.minorticks_off()
    ax2.set_yticks([1.0, 0.5, 0.2, 0.1, 0.05])
    ax2.set_yticklabels(['1', '0.5', '0.2', '0.1', '0.05'])
    ax2.set_ylim(0.04, 1.3)
    ax2.set_xlabel('gap number')
    ax2.set_ylabel(r'$\Delta p_i / \Delta p_1$')
    ax2.legend(loc='lower left')
    P.letter(ax2, 'b')
    save(fig, out, 'F5-1_outlet_slab')


def main():
    root = sys.argv[1]
    out = os.path.join(root, 'figures')
    os.makedirs(out, exist_ok=True)
    kb = json.load(open(os.path.join(root, 'kb', 'kb_data.json'), encoding='utf-8'))
    cs = json.load(open(os.path.join(root, 'kb', 'comsol_5gap_reports.json'), encoding='utf-8'))
    fig_cubic(out, kb)
    fig_kappa(out, kb)
    fig_depth(out, kb)
    fig_rows(out)
    fig_duct(out)
    fig_slab(out, cs)
    bad = [line for line in LOG if 'MISMATCH' in line]
    note('%d checks, %d mismatches' % (sum('computed' in line for line in LOG), len(bad)))
    with open(os.path.join(out, 'FIGURES_KB_CHECK.txt'), 'w', encoding='utf-8') as fh:
        fh.write('\n'.join(LOG) + '\n')


if __name__ == '__main__':
    main()
