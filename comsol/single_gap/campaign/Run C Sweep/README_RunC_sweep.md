# Run C Sweep — shallow aspect ratio (α = 0.8) confirmation

This folder is the **shallow-channel** companion to the Run B Sweep. Same single-gap
quarter-cell CFD, same post-processing, one geometry change: the channel depth is reduced
from 140 µm to **80 µm** (w_c = 100, gap = 20, w_s = 60 unchanged), giving aspect ratio
**α = depth/w_c = 0.80**. Run B established the method at α = 1.4; Run C is a *pre-registered,
out-of-sample test* of the analytic skim law at a genuinely different — and more
Hele-Shaw-like — aspect ratio.

## Result (three checks, all pass)

1. **K_inplane = 1.000 at every fraction and every upstream station** (max bulge 0.0000 µm).
   The mid-depth dividing streamsurface thins monotonically toward the gap; it never bulges
   toward the channel centre. α = 0.8 is the shallow regime where in-plane funnelling, if it
   existed, would be *most* likely to appear — and it does not. The funnelling-enhancement
   hypothesis is refuted at a second aspect ratio, the one most favourable to it.

2. **δ_skim ∝ √f_phys** to R² = 0.99997 across f = 2/5/9 ×10⁻⁴ (spread 0.35%).

3. **δ_skim matches the exact rectangular-duct law to ≤0.2%** — validation = 1.000–1.002,
   with R_depth = δ_skim/δ_parab = 0.891 (spread 0.20%) against the **pre-registered
   prediction R(0.8) = 0.890**. The prediction was written down before the run; the CFD hit
   it to the third decimal.

## Pre-registered targets vs CFD

| frac  | f_phys   | predicted δ_skim | CFD δ_skim | predicted R(0.8) | CFD validation |
|-------|----------|------------------|------------|------------------|----------------|
| 2e-4  | 1.0e-4   | 0.515 µm         | 0.508 µm   | 0.890            | 1.002          |
| 5e-4  | 2.5e-4   | 0.816 µm         | 0.804 µm   | 0.891            | 1.000          |
| 9e-4  | 4.5e-4   | 1.097 µm         | 1.081 µm   | 0.892            | 1.001          |

(CFD δ_skim uses the achieved f_gap = 0.971·frac, so each entry sits a touch below the
input-frac prediction; the validation ratio, which is what the test turns on, is insensitive
to this and lands on 1.00.)

## Two conventions to keep straight

- **Factor of 2 (flow):** the CFD's reported `f_gap = Q_sout(quarter)/Q_in(quarter)` is
  **twice** the physical one-wall/full-duct fraction the √(f/3) law uses. The post-processor
  forms `f_phys = f_gap/2` and compares δ_skim against the exact duct law at f_phys. This is
  bookkeeping, not a model error; the √f and K_inplane checks are ratios and are unaffected.

- **Achieved vs design fraction:** the measured gap fraction is `f_gap = Q_sout/Q_in =
  0.971 × frac` for every point in this run — a *geometric* factor, constant across fractions:
  the side-outlet boundary area is 2.9% under the nominal `w_s·depth/2`, so the velocity BC
  delivers 0.971·frac·Q0. (Run B's factor was 0.992; it is larger here because depth 80 makes
  the same edge effect a bigger relative fraction.) The achieved fraction is the physically
  correct band-flux target and is read straight from each report's Global Evaluation row
  `[Q_in, Q_cout, Q_sout, mass_bal, f_gap]`. All three reports carry that row, so the summary
  uses the achieved fraction throughout and runs warning-free. (Because the validation ratio
  is √f-independent, using the design frac instead would change no check — the summary script
  falls back to it automatically, with a warning, if a report is ever exported without the
  evaluation.)

## Folder layout

```
Run C Sweep/
  2e-4/   Report Run C - Sweep 2e-4.docx   upstream_grid_frac2e-4.txt   middepth_grid_frac2e-4.txt
  5e-4/   Report Run C - Sweep 5e-4.docx   upstream_grid_frac5e-4.txt   middepth_grid_frac5e-4.txt
  9e-4/   Report Run C - Sweep 9e-4.docx   upstream_grid_frac9e-4.txt   middepth_grid_frac9e-4.txt
  cif_funnel_postproc.py      per-point analyzer (δ_skim, δ_bulge, K_inplane, R, validation)
  cif_sweep_summary.py        sweep driver — auto-discovers the fraction subfolders
  README_RunC_sweep.md        this file
  REFERENCE_OUTPUT_RunC.txt   known-good console output to diff against
```

The two `.txt` per fraction are the COMSOL grid exports: upstream YZ plane at x = 75
(y ∈ [0, 50], z ∈ [0, 40]) and mid-depth XY plane at z = 40 (x ∈ [60, 140]). File names do
not matter — the scripts classify each export by which coordinate is fixed.

## Reproduce

```
python cif_sweep_summary.py "Run C Sweep"
```

Geometry (w_c = 100, depth = 80, α = 0.80) is auto-detected from the upstream export
(w_c = 2·max y, depth = 2·max z), so the same script serves Runs B, A, and C unchanged.
Requires `cif_funnel_postproc.py` in the same folder; numpy ≥ 1.24, scipy. Compare the
console output to `REFERENCE_OUTPUT_RunC.txt`.

## What this settles — and what it doesn't

**Settles:** the exact rectangular-duct skim law δ_skim = w_c·√(f_phys/3)·R(α) is correct at
α = 0.8 (to ≤0.2%), the √f scaling holds, and there is no in-plane funnelling (K_inplane = 1)
even in the shallow, near-Hele-Shaw regime. Together with Run B (α = 1.4) this brackets the
depth correction R(α) on the shallow side and confirms it is a real, predictable geometric
factor — not a fitted fudge.

**Does not settle:** the absolute particle size cut-off. δ_skim is the *hydrodynamic* skim
depth (here ~0.5–1.1 µm). The empirical cut-off prefactor κ ≈ 5 lives in the steps *downstream*
of δ_skim — the steric capture criterion (κ_steric ∈ {1,2}), the effective fraction set by the
device's resistance network (f_eff several× nominal), and cumulative finite-particle
displacement over many gaps — none of which this single-gap, point-tracer CFD measures.
Run C confirms the skim hydrodynamics; the cut-off closure is separate work.
