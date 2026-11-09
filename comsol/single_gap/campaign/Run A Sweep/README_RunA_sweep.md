# Run A Sweep — deep aspect ratio (α = 2.69), Dinh geometry confirmation

This folder is the **deep-channel** companion to the Run B and Run C sweeps, and it is the
one run that sits on a *real device*: the geometry is the narrowest gap of Dinh et al. 2024
(Lab Chip, d3lc00842h) — middle channel **w_c = 52 µm** at the outlet, bilateral **19 µm**
gaps, depth **140 µm** — giving aspect ratio **α = depth/w_c = 2.69**. Same single-gap
quarter-cell CFD and post-processing as the other runs.

Run B established the method at α = 1.4 and Run C tested the shallow side at α = 0.8. Run A is
the *pre-registered, out-of-sample* test at the **deep** extreme — the opposite end of the
aspect-ratio range from Run C — and it closes the bracket on the depth correction R(α).

## Result (three checks, all pass)

1. **K_inplane = 1.000 at every resolvable fraction and every upstream station** (max bulge
   0.0000 µm). The mid-depth dividing streamsurface thins monotonically toward the gap; it
   never bulges toward the channel centre. Together with Run C (shallow) this brackets the
   full aspect-ratio range α = 0.8 → 1.4 → 2.69 — and there is **no in-plane funnelling
   anywhere in it.** The funnelling-enhancement hypothesis is refuted at both extremes.
   (The lowest fraction, 2e-4, returns K_inplane = n/a for a documented low-f
   post-processing reason — see *The 2e-4 mid-depth point* below. It is not a model issue and
   does not affect the conclusion, which is established at 5e-4, 9e-4, and all of Runs B/C.)

2. **δ_skim ∝ √f_phys** to R² = 0.99993 across f = 2/5/9 ×10⁻⁴ (spread 0.71%).

3. **δ_skim matches the exact rectangular-duct law to ≤0.3%** — validation = 0.998–1.003,
   with R_depth = δ_skim/δ_parab = 0.983 (spread 0.55%) against the **pre-registered
   prediction R(2.69) = 0.980**. The prediction was written down before the run; the CFD hit
   it to the third decimal at the deep extreme.

## Pre-registered targets vs CFD

| frac  | f_phys (input) | predicted δ_skim | CFD δ_skim | predicted R(2.69) | CFD validation |
|-------|----------------|------------------|------------|-------------------|----------------|
| 2e-4  | 1.0e-4         | 0.295 µm         | 0.290 µm   | 0.980             | 0.998          |
| 5e-4  | 2.5e-4         | 0.467 µm         | 0.461 µm   | 0.980             | 1.003          |
| 9e-4  | 4.5e-4         | 0.627 µm         | 0.619 µm   | 0.980             | 1.002          |

(CFD δ_skim uses the achieved f_gap = 0.970·frac, so each entry sits a touch below the
input-frac prediction; the validation ratio, which is what the test turns on, is insensitive
to this and lands on 1.00.)

## How Run A compares to the other two aspect ratios

| run | α    | R(α) predicted | R(α) CFD | √f R²    | K_inplane            | validation  |
|-----|------|----------------|----------|----------|----------------------|-------------|
| C   | 0.80 | 0.890          | 0.891    | 0.99997  | 1.000 (all f)        | 1.000–1.002 |
| B   | 1.40 | 0.951          | ~0.95    | —        | 1.000 (all f)        | ~1.00       |
| A   | 2.69 | 0.980          | 0.983    | 0.99993  | 1.000 (5e-4, 9e-4)   | 0.998–1.003 |

R(α) rises monotonically toward 1 as the channel deepens, exactly as the analytic duct law
predicts. Three independent aspect ratios, three pre-registered hits. The depth correction is
a real, predictable geometric factor — not a fitted fudge.

## Two conventions to keep straight

- **Factor of 2 (flow):** the CFD's reported `f_gap = Q_sout(quarter)/Q_in(quarter)` is
  **twice** the physical one-wall/full-duct fraction the √(f/3) law uses. The post-processor
  forms `f_phys = f_gap/2` and compares δ_skim against the exact duct law at f_phys. This is
  bookkeeping, not a model error; the √f and K_inplane checks are ratios and are unaffected.
  For Dinh's geometry this is doubly confirmed: Dinh defines the gap fraction f*gap *per gap*
  with bilateral gaps, so Dinh's f*gap equals the physical one-wall fraction = f_gap(CFD)/2.

- **Achieved vs design fraction:** the measured gap fraction is `f_gap = Q_sout/Q_in =
  0.970 × frac` for every point in this run — a *geometric* factor, constant across all three
  fractions to four digits (0.9703 / 0.9704 / 0.9704): the side-outlet boundary area is ~3%
  under the nominal `w_s·depth/2`, so the velocity BC delivers 0.970·frac·Q0. (Run B's factor
  was 0.992, Run C's 0.971.) The achieved fraction is the physically correct band-flux target
  and is read straight from each report's Global Evaluation row `[Q_in, Q_cout, Q_sout,
  mass_bal, f_gap]`. All three reports carry that row, so the summary uses the achieved
  fraction throughout and runs warning-free.

## The 2e-4 mid-depth point (why K_inplane = n/a there)

Run A is the most demanding aspect ratio for the *mid-depth* (K_inplane) measurement because
the deep channel makes the mid-plane skim very thin (δ_skim_mid ≈ 0.16 µm at 2e-4). The
bulge check needs the net flux through the gap at mid-depth, but along the gap aperture the
in-plane velocity is dominated by large **antisymmetric recirculation lobes** (±1.3×10⁻² at
the gap edges — flow turning in on the upstream side and out on the downstream side), while
the actual net through-flow at the gap centre is tiny (+2.6×10⁻⁵ at 2e-4). At 5e-4 and 9e-4
the through-flow is large enough to survive that subtraction; at the lowest fraction it is
~5× smaller and the imperfect discrete cancellation of the recirculation flips the residual
negative, collapsing the dividing line onto the wall. The analyzer detects this (φ_gap ≤ 0)
and reports K_inplane as **n/a (mid-depth flux indeterminate at low f)** rather than a
spurious value.

This is a post-processing limit at small f, **not** a model or resolution problem:

- the **depth-averaged** δ_skim — which drives the validation and R_depth, i.e. the actual
  publishable result — is computed from the well-resolved upstream u-field and is unaffected
  (δ_skim = 0.290 µm, validation = 0.998, fully export-converged); and
- K_inplane = 1.000 is established robustly at Run A's 5e-4 and 9e-4 points and at *every*
  Run B and Run C fraction — a dozen-plus points spanning the whole α range. The 2e-4
  mid-depth failure is not evidence of a bulge; it is simply below the recirculation floor.

(A mass-balance flux estimator recovers a positive φ_gap for this point, but it breaks the
shallow-run cases, so it is not a clean universal swap; the validated gap-flux method is kept
and the indeterminate point is flagged rather than forced.)

## Dinh operating point

Because Run A *is* the Dinh outlet geometry, the sweep brackets the paper's real operating
point. Dinh's main condition is f*gap = 1.04×10⁻⁴ per gap, which in CFD units is frac =
2.08×10⁻⁴ — so the **2e-4 point sits essentially on it** (within 4%; use frac = 2.08e-4 to
land on it exactly). The CFD lamina there is δ_skim ≈ 0.30 µm. Dinh's SI Fig S1 puts the
*measured* bead cut-off between 1.0 µm (passes, recovers ~fluid) and 4.7 µm (retained), and
the empirical κ ≈ 5 law gives d_crit = 1.53 µm — inside that bracket. The gap between the
0.30 µm hydrodynamic skim and the ~1.5 µm measured cut-off (a factor of ≥1.7–3.3×, depending
on the steric convention κ_steric ∈ {1,2}) is the part that lives **downstream** of δ_skim:
effective fraction (resistance network) and cumulative finite-particle displacement. Run A
validates the skim hydrodynamics on the real geometry; that downstream factor is the
remaining open piece of the cut-off closure.

## Folder layout

```
Run A Sweep/
  2e-4/   Report Run A - Sweep 2e-4.docx   upstream_grid_frac2e-4.txt   middepth_grid_frac2e-4.txt
  5e-4/   Report Run A - Sweep 5e-4.docx   upstream_grid_frac5e-4.txt   middepth_grid_frac5e-4.txt
  9e-4/   Report Run A - Sweep 9e-4.docx   upstream_grid_frac9e-4.txt   middepth_grid_frac9e-4.txt
  cif_funnel_postproc.py      per-point analyzer (δ_skim, δ_bulge, K_inplane, R, validation)
  cif_sweep_summary.py        sweep driver — auto-discovers the fraction subfolders
  README_RunA_sweep.md        this file
  REFERENCE_OUTPUT_RunA.txt   known-good console output to diff against
```

The two `.txt` per fraction are the COMSOL grid exports: upstream YZ plane at x = 75
(y ∈ [0, 26], z ∈ [0, 70]) and mid-depth XY plane at z = 70 (x ∈ [40, 140]). File names do
not matter — the scripts classify each export by which coordinate is fixed. Add your `.mph`
files to each subfolder if you want the solved models alongside the exports.

## Reproduce

```
python cif_sweep_summary.py "Run A Sweep"
```

Geometry (w_c = 52, depth = 140, α = 2.69) is auto-detected from the upstream export
(w_c = 2·max y, depth = 2·max z), so the same script serves Runs B, C, and A unchanged.
Requires `cif_funnel_postproc.py` in the same folder; numpy ≥ 1.24, scipy. Compare the
console output to `REFERENCE_OUTPUT_RunA.txt`.

## What this settles — and what it doesn't

**Settles:** the exact rectangular-duct skim law δ_skim = w_c·√(f_phys/3)·R(α) is correct at
α = 2.69 (to ≤0.3%), the √f scaling holds, and there is no in-plane funnelling (K_inplane = 1)
at the deep extreme. With Runs C (0.8) and B (1.4) this confirms R(α) across the full
aspect-ratio range on three pre-registered points, on real Dinh geometry at the deep end.

**Does not settle:** the absolute particle size cut-off. δ_skim is the *hydrodynamic* skim
depth (~0.3–0.6 µm here). The empirical cut-off prefactor κ ≈ 5 lives in the steps
*downstream* of δ_skim — the steric capture criterion (κ_steric ∈ {1,2}), the effective
fraction set by the device's resistance network (f_eff several× nominal), and cumulative
finite-particle displacement over many gaps — none of which this single-gap, point-tracer CFD
measures. On the Dinh point specifically, that downstream factor must supply ≥1.7–3.3× to
reach the measured 1.0–4.7 µm cut-off. Run A confirms the skim hydrodynamics; the cut-off
closure is separate work.
