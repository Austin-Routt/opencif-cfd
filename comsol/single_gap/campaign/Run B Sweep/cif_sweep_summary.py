#!/usr/bin/env python3
"""
CIF f-sweep verification (Run B geometry) -- self-contained archive reproducer.

Point this at the Run B Sweep archive root.  For each fraction subfolder it:
  * finds the two COMSOL field exports and classifies them (upstream YZ vs mid-depth XY)
    by which coordinate is fixed -- no reliance on file names;
  * reads the ACHIEVED f_gap straight from that run's COMSOL report (.docx),
    column intop_sout/intop_in -- no hand-transcribed numbers;
  * computes the per-run metrics via cif_funnel_postproc.analyze.

It then reports:
  CHECK 1  delta_skim ~ sqrt(f)         -> extraction-method validation (a ratio; the
                                           factor-of-2 fraction convention cancels)
  CHECK 2  K_inplane ~ 1 at every f     -> no in-plane funnelling at any operating point
  CHECK 3  depth factor R f-independent -> R = delta_skim/delta_parab(f_gap/2) ~ 0.95,
           and VALIDATION delta_skim/delta_general(f_gap/2) ~ 1.00 at every point
           (the CFD reproduces the exact rectangular-duct law)
  ROBUSTNESS  K_inplane vs upstream sampling station + a monotonicity (bulge) test.

See the factor-of-2 note in cif_funnel_postproc.py: the report's f_gap is twice the
one-wall/full-duct fraction the analytic law uses, so the law is evaluated at f_gap/2.

Usage:   python cif_sweep_summary.py [ARCHIVE_ROOT]   (default: this script's directory)
Requires cif_funnel_postproc.py in the same directory; numpy>=1.24, scipy.
"""
import os
import re
import sys
import glob
import zipfile
import numpy as np
from cif_funnel_postproc import analyze, classify_plane, middepth_metrics

# fixed geometry of the Run B family
W_C, DEPTH, GAP_X0, GAP_X1 = 100.0, 140.0, 100.0, 120.0


def read_fgap_from_docx(path):
    """f_gap = intop_sout/intop_in, read from the report's Global Evaluation table.

    The table row is [Q_in, Q_cout, Q_sout, mass_balance, f_gap]; f_gap is its last
    entry and the only one in the physical fraction range (1e-6, 1e-2).  We locate the
    row directly -- a run of >=5 scientific-notation numbers whose cells sit close
    together (~hundreds of chars apart after tag-stripping) -- rather than trusting a
    section heading, since the heading text can recur elsewhere in the document."""
    xml = zipfile.ZipFile(path).read('word/document.xml').decode('utf-8', 'ignore')
    txt = re.sub(r'<[^>]+>', ' ', xml)
    toks = list(re.finditer(r'-?\d\.?\d*[eE][-+]?\d+', txt))
    groups, cur = [], []
    for m in toks:                                # cluster numbers into table rows
        if cur and m.start() - cur[-1].end() > 500:
            groups.append(cur); cur = []
        cur.append(m)
    if cur:
        groups.append(cur)
    for g in reversed(groups):                    # last >=5-number row ending in a fraction
        if len(g) >= 5 and 1e-6 < abs(float(g[-1].group())) < 1e-2:
            return abs(float(g[-1].group()))
    for m in reversed(toks):                       # fallback: last fraction-range number
        if 1e-6 < abs(float(m.group())) < 1e-2:
            return abs(float(m.group()))
    raise ValueError(f"no f_gap-like value found in {path}")


def find_run(folder):
    """Return (f_gap, upstream_txt, middepth_txt) for one fraction subfolder."""
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
    return read_fgap_from_docx(docs[0]), up, md


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
            f_gap, up, md = find_run(folder)
            r = analyze(up, md, W_C, f_gap, DEPTH, GAP_X0, GAP_X1)
            r['label'] = lab
            r['_md'] = md
            rows.append(r)
        except Exception as e:
            print(f"  [skip {lab}] {e}")

    if not rows:
        sys.exit("no runs processed")

    print(f"\n{'frac':<6}{'f_gap':>11}{'f_phys':>11}{'d_skim':>8}{'d_bulge':>9}"
          f"{'K_in':>7}{'R_depth':>9}{'valid':>8}")
    print("  " + "-" * 67)
    for r in rows:
        print(f"{r['label']:<6}{r['f_gap']:>11.4e}{r['f_phys']:>11.4e}"
              f"{r['delta_skim_da']:>8.3f}{r['delta_bulge']:>9.3f}"
              f"{r['K_inplane']:>7.3f}{r['R_depth']:>9.3f}{r['validation']:>8.3f}")

    f = np.array([r['f_phys'] for r in rows])           # physical fraction
    d = np.array([r['delta_skim_da'] for r in rows])
    K = np.array([r['K_inplane'] for r in rows])
    R = np.array([r['R_depth'] for r in rows])
    valid = np.array([r['validation'] for r in rows])

    # CHECK 1: delta_skim = m*sqrt(f_phys); R^2 and the spread of the ratio
    sq = np.sqrt(f); m = np.sum(d * sq) / np.sum(sq * sq); pred = m * sq
    ss_tot = np.sum((d - d.mean()) ** 2)
    R2 = (1 - np.sum((d - pred) ** 2) / ss_tot) if ss_tot > 0 else float('nan')
    R2s = f"{R2:.5f}" if np.isfinite(R2) else "n/a (need >=2 distinct fractions)"
    rat = d / sq
    deep_const = W_C / np.sqrt(3.0)                     # delta/sqrt(f) in the deep limit
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
            dmid, db, *_ = middepth_metrics(r['_md'], W_C, GAP_X0, GAP_X1, x_up=xup)
            ks.append(db / dmid)
        dmid, db, xs, yd, tr, phi = middepth_metrics(r['_md'], W_C, GAP_X0, GAP_X1,
                                                     x_up=GAP_X0 - 37)
        sw = W_C / 2 - yd
        dmax = max(0.0, np.nanmax(np.diff(sw[~np.isnan(sw)])))
        print(f"    f={r['label']:<5} K at x_up[auto,75,70,65,62]="
              f"{[f'{x:.3f}' for x in ks]}  max bulge {dmax:.4f} um")


if __name__ == "__main__":
    main()
