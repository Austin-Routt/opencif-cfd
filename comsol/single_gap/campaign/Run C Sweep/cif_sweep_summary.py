#!/usr/bin/env python3
"""
CIF f-sweep verification -- self-contained, geometry-agnostic archive reproducer.

Works for ANY of the aspect-ratio runs (Run B alpha=1.4, Run C alpha=0.8, Run A/Dinh
alpha=2.69): point it at a run's archive root and it auto-detects the geometry from the
exports. For each fraction subfolder it:
  * finds the two COMSOL field exports and classifies them (upstream YZ vs mid-depth XY)
    by which coordinate is fixed -- no reliance on file names;
  * auto-detects w_c and depth from the upstream export extents
    (y spans [0, w_c/2]; z spans [0, depth/2]);
  * reads the ACHIEVED f_gap straight from that run's COMSOL report (.docx),
    column intop_sout/intop_in -- no hand-transcribed numbers;
  * computes the per-run metrics via cif_funnel_postproc.analyze.

Reports:
  CHECK 1  delta_skim ~ sqrt(f)         -> extraction-method validation (a ratio; the
                                           factor-of-2 fraction convention cancels)
  CHECK 2  K_inplane ~ 1 at every f     -> no in-plane funnelling at any operating point
  CHECK 3  depth factor R f-independent -> R = delta_skim/delta_parab(f_gap/2),
           and VALIDATION delta_skim/delta_general(f_gap/2) ~ 1.00 at every point
           (the CFD reproduces the exact rectangular-duct law at this aspect ratio)
  ROBUSTNESS  K_inplane vs upstream sampling station + a monotonicity (bulge) test.

The streamwise gap window (GAP_X0, GAP_X1) is fixed because the x-layout is identical
across runs; only the cross-section (w_c, depth) changes. The factor-of-2 (report f_gap =
2 x physical one-wall/full-duct fraction) is handled in cif_funnel_postproc.analyze.

Usage:   python cif_sweep_summary.py [ARCHIVE_ROOT]   (default: this script's directory)
Requires cif_funnel_postproc.py in the same directory; numpy>=1.24, scipy.
"""
import os
import re
import sys
import glob
import zipfile
import numpy as np
from cif_funnel_postproc import analyze, classify_plane, middepth_metrics, parse_scatter

GAP_X0, GAP_X1 = 100.0, 120.0          # streamwise gap window (same x-layout for B/A/C)


def detect_geometry(upstream_file):
    """w_c, depth from the upstream YZ export: y spans [0, w_c/2], z spans [0, depth/2]."""
    x, y, z, *_ = parse_scatter(upstream_file)
    return float(round(2 * max(y))), float(round(2 * max(z)))


def read_fgap_from_docx(path):
    """f_gap = intop_sout/intop_in, read from the report's Global Evaluation table.
    The table row is [Q_in, Q_cout, Q_sout, mass_balance, f_gap]; f_gap is its last entry
    and the only one in the physical fraction range (1e-6, 1e-2). We locate the row
    directly -- a run of >=5 scientific-notation numbers whose cells sit close together --
    rather than trusting a heading, since the heading text can recur in the document.

    If no such row exists (a report exported without the Global Evaluation), fall back to
    the design 'frac' parameter from the Parameters table.  The achieved f_gap differs from
    frac only by the side-outlet area factor (a few %), and the validation ratio is
    sqrt(f)-independent, so this substitution does not change any conclusion -- but it
    prevents grabbing an unrelated number (e.g. a solver tolerance).  A warning is emitted."""
    xml = zipfile.ZipFile(path).read('word/document.xml').decode('utf-8', 'ignore')
    txt = re.sub(r'<[^>]+>', ' ', xml)
    toks = list(re.finditer(r'-?\d\.?\d*[eE][-+\u2212]?\d+', txt))
    groups, cur = [], []
    for m in toks:
        if cur and m.start() - cur[-1].end() > 500:
            groups.append(cur); cur = []
        cur.append(m)
    if cur:
        groups.append(cur)

    def _f(s):                                   # COMSOL uses U+2212 (minus sign) in docx
        return float(s.replace('\u2212', '-'))

    # 1) Global Evaluation row: a cluster of >=5 numbers, last one in the fraction range.
    for g in reversed(groups):
        if len(g) >= 5 and 1e-6 < abs(_f(g[-1].group())) < 1e-2:
            return abs(_f(g[-1].group()))
    # 2) Fallback: the design 'frac' parameter (always present in the Parameters table).
    #    Collapse whitespace first -- COMSOL puts the name and value in separate table cells
    #    with newlines between them, so they are not adjacent in the tag-stripped text.
    flat = re.sub(r'\s+', ' ', txt)
    m = re.search(r'\bfrac\b\s+(\d\.?\d*[eE][-+\u2212]?\d+)', flat)
    if m and 1e-6 < abs(_f(m.group(1))) < 1e-2:
        val = abs(_f(m.group(1)))
        sys.stderr.write(f"  [warn] {os.path.basename(path)}: no Global Evaluation table found; "
                         f"using design frac={val:.3e} (achieved f_gap unavailable). "
                         f"Re-export the report with the f_gap evaluation for full consistency.\n")
        return val
    # 3) last resort: last in-range token (legacy behaviour).
    for m in reversed(toks):
        if 1e-6 < abs(_f(m.group())) < 1e-2:
            return abs(_f(m.group()))
    raise ValueError(f"no f_gap-like value found in {path}")


def find_run(folder):
    """Return (f_gap, upstream_txt, middepth_txt, w_c, depth) for one fraction subfolder."""
    txts = sorted(glob.glob(os.path.join(folder, '*.txt')))
    docs = sorted(glob.glob(os.path.join(folder, '*.docx')))
    if len(txts) < 2 or not docs:
        raise FileNotFoundError(f"{folder}: need >=2 .txt and a .docx report")
    up = md = None
    for t in txts:
        kind = classify_plane(t)
        if kind == 'upstream' and up is None:
            up = t
        elif kind == 'middepth' and md is None:
            md = t
    if up is None or md is None:
        raise RuntimeError(f"{folder}: could not classify exports as upstream + middepth")
    w_c, depth = detect_geometry(up)
    return read_fgap_from_docx(docs[0]), up, md, w_c, depth


def main():
    root = sys.argv[1] if len(sys.argv) > 1 else os.path.dirname(os.path.abspath(__file__))
    subs = sorted(d for d in glob.glob(os.path.join(root, '*'))
                  if os.path.isdir(d) and glob.glob(os.path.join(d, '*.docx')))
    if not subs:
        sys.exit(f"no fraction subfolders (each with a .docx report) found under {root}")

    rows = []
    for folder in subs:
        lab = os.path.basename(folder)
        try:
            f_gap, up, md, w_c, depth = find_run(folder)
            r = analyze(up, md, w_c, f_gap, depth, GAP_X0, GAP_X1)
            r.update(label=lab, _md=md, _wc=w_c, _depth=depth)
            rows.append(r)
        except Exception as e:
            print(f"  [skip {lab}] {e}")
    if not rows:
        sys.exit("no runs processed")

    wcs = {r['_wc'] for r in rows}; deps = {r['_depth'] for r in rows}
    wc = sorted(wcs)[0]
    geo = (f"w_c={wc:.0f} depth={sorted(deps)[0]:.0f} alpha={sorted(deps)[0]/wc:.2f}"
           if len(wcs) == 1 and len(deps) == 1
           else f"MIXED geometry across subfolders: w_c={wcs}, depth={deps}  (check archive!)")
    print(f"\n  auto-detected geometry: {geo}")
    print(f"\n{'frac':<6}{'f_gap':>11}{'f_phys':>11}{'d_skim':>8}{'d_bulge':>9}"
          f"{'K_in':>7}{'R_depth':>9}{'valid':>8}")
    print("  " + "-" * 67)
    for r in rows:
        print(f"{r['label']:<6}{r['f_gap']:>11.4e}{r['f_phys']:>11.4e}"
              f"{r['delta_skim_da']:>8.3f}{r['delta_bulge']:>9.3f}"
              f"{r['K_inplane']:>7.3f}{r['R_depth']:>9.3f}{r['validation']:>8.3f}")

    f = np.array([r['f_phys'] for r in rows])
    d = np.array([r['delta_skim_da'] for r in rows])
    K = np.array([r['K_inplane'] for r in rows])
    R = np.array([r['R_depth'] for r in rows])
    valid = np.array([r['validation'] for r in rows])

    sq = np.sqrt(f); m = np.sum(d * sq) / np.sum(sq * sq); pred = m * sq
    ss_tot = np.sum((d - d.mean()) ** 2)
    R2 = (1 - np.sum((d - pred) ** 2) / ss_tot) if ss_tot > 0 else float('nan')
    R2s = f"{R2:.5f}" if np.isfinite(R2) else "n/a (need >=2 distinct fractions)"
    rat = d / sq; deep_const = wc / np.sqrt(3.0)
    print(f"\n  CHECK 1  delta_skim/sqrt(f_phys) spread "
          f"{100 * (rat.max() - rat.min()) / rat.mean():.2f}%   "
          f"fit d = {m:.3f}*sqrt(f_phys), R^2 = {R2s}")
    print(f"           (deep-limit constant w_c/sqrt(3) = {deep_const:.3f}; "
          f"m/that = {m / deep_const:.3f} = mean depth factor R)")
    print(f"  CHECK 2  K_inplane  {', '.join(f'{x:.3f}' for x in K)}   "
          f"(~1 => no funnelling at any f)")
    print(f"  CHECK 3  R_depth    {', '.join(f'{x:.3f}' for x in R)}   "
          f"spread {100 * (R.max() - R.min()) / R.mean():.2f}%  (f-independent)")
    print(f"           VALIDATION {', '.join(f'{x:.3f}' for x in valid)}   "
          f"(delta_skim / exact duct law at f_gap/2; ~1.00 => CFD reproduces the law)")

    print("\n  ROBUSTNESS  K_inplane vs upstream sampling station (rules out a hidden bulge):")
    for r in rows:
        ks = []
        for xup in [None, 75, 70, 65, 62]:
            dmid, db, *_ = middepth_metrics(r['_md'], r['_wc'], GAP_X0, GAP_X1, x_up=xup)
            ks.append(db / dmid)
        dmid, db, xs, yd, tr, phi = middepth_metrics(r['_md'], r['_wc'], GAP_X0, GAP_X1,
                                                     x_up=GAP_X0 - 37)
        sw = r['_wc'] / 2 - yd
        dmax = max(0.0, np.nanmax(np.diff(sw[~np.isnan(sw)])))
        print(f"    f={r['label']:<5} K at x_up[auto,75,70,65,62]="
              f"{[f'{x:.3f}' for x in ks]}  max bulge {dmax:.4f} um")


if __name__ == "__main__":
    main()
