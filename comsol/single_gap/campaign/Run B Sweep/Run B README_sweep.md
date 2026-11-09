# Run B f-Sweep — CIF single-gap funnelling CFD (archive)

Single-gap Stokes (creeping-flow) solves of one CIF unit cell, sweeping the per-gap
fraction. Geometry: deep "Run B" family — channel width `w_c = 100 µm`, depth
`h = 140 µm` (aspect α = 1.4), gap aperture 20 µm. Quarter-cell model with symmetry at
`y = 0` and `z = 70 µm`. COMSOL 6.3.

## Result (three findings)

Across `f_gap = 2–9×10⁻⁴` (five points, a 4.5× range):

1. **K_inplane = 1.000 at every fraction.** The mid-depth dividing streamsurface does
   **not** bulge toward the channel centre at the gap; it *thins* toward it (line-sink
   behaviour). There is no in-plane funnelling, at any operating point. This refutes the
   funnelling-enhancement hypothesis (κ ≈ 5 is **not** an in-plane convergence effect).

2. **δ_skim ∝ √f to R² = 0.99992.** The leading-order skim law `δ/w_c = √(f/3)` holds
   across the fraction range, confirmed by interior points, not just the endpoints.

3. **δ_skim reproduces the exact rectangular-duct skim law to < 1%**, at the physical
   fraction `f_gap/2` (see the factor-of-2 note below), with depth factor **R ≈ 0.95** —
   the analytic rectangular-duct value at α = 1.4. The 3-D CFD validates the
   developed-duct depth correction; there is no "offset."

> **⚠ Factor of 2 in the reported fraction.** The model is a quarter cell, and COMSOL
> reports `f_gap = intop_sout/intop_in = Q_sout(quarter)/Q_in(quarter)`. The analytic
> skim law (and the empirical κ work) define the fraction as
> [one-wall, full-depth lamina flux] / [full-duct flux], which equals **f_gap / 2**:
> the one-wall full-depth lamina is 2× the quarter lamina (z-symmetry), while the full
> duct is 4× the quarter inlet (y- and z-symmetry). So **the reported f_gap is twice the
> physical one-wall/full-duct fraction.** This does *not* affect δ_skim, K_inplane, or
> the √f scaling (δ_skim is the layer that physically carries Q_sout; the other two are
> ratios). It *does* mean any absolute cutoff/fraction prediction taken from this CFD
> must use **f_physical = f_gap/2**. The depth comparison above is correct only at
> f_gap/2; comparing δ_skim to the law at the un-halved f_gap gives a spurious R ≈ 0.67.

## Archive layout

```
Run B Sweep/
├── README_sweep.md            ← this file
├── cif_funnel_postproc.py     ← per-run extractor (δ_skim, δ_bulge, K_inplane, R, validation)
├── cif_sweep_summary.py       ← reproduces the whole table:  python cif_sweep_summary.py
├── REFERENCE_OUTPUT.txt       ← known-good output of the summary, to diff against
├── CIF_RunB_Sweep_StepByStep.md ← how the COMSOL runs were generated (procedure of record)
├── 2e-4/ … 9e-4/              ← one folder per fraction, each containing:
│     • the .mph model         (COMSOL 6.3; holds that fraction's solve)
│     • two field exports .txt (upstream YZ plane at x = 75; mid-depth XY plane at z = 70)
│     • the COMSOL report .docx (Global Evaluation table: Q_in, Q_cout, Q_sout, mass-balance, f_gap)
```

The `5e-4` folder is the **canonical single-point Run B**; the sweep is the robustness
layer around it, not a replacement.

## Reproduce

```
python cif_sweep_summary.py            # run from inside Run B Sweep/
# or:  python cif_sweep_summary.py /path/to/"Run B Sweep"
```

Requires `cif_funnel_postproc.py` beside it; `numpy >= 1.24`, `scipy`. The summary script
discovers each fraction subfolder, classifies the two exports by which coordinate is fixed
(so file names do not matter), reads the achieved `f_gap` directly from each `.docx`
report, computes every metric, and prints the table plus the three checks and the
robustness block. Output should match `REFERENCE_OUTPUT.txt` to the last digit shown.

## Per-point CFD verification (from each report's Global Evaluation table)

| frac | mass balance | achieved f_gap | Q_in (m³/s) | Q_sout (m³/s) | Q_sout / Q_sout(5e-4) |
|---|---|---|---|---|---|
| 2e-4 | −1.13×10⁻¹⁵ | 1.9848×10⁻⁴ | 3.2749×10⁻¹⁰ | 6.5002×10⁻¹⁴ | 0.400 |
| 3e-4 | −2.09×10⁻¹⁵ | 2.9773×10⁻⁴ | 3.2749×10⁻¹⁰ | 9.7502×10⁻¹⁴ | 0.600 |
| 5e-4 | 3.09×10⁻¹⁵ | 4.9621×10⁻⁴ | 3.2749×10⁻¹⁰ | 1.6250×10⁻¹³ | 1.000 |
| 7e-4 | 1.80×10⁻¹⁵ | 6.9470×10⁻⁴ | 3.2749×10⁻¹⁰ | 2.2751×10⁻¹³ | 1.400 |
| 9e-4 | 3.68×10⁻¹⁷ | 8.9318×10⁻⁴ | 3.2749×10⁻¹⁰ | 2.9251×10⁻¹³ | 1.800 |

Each fraction lands ~0.76% under its nominal target (consistent BC discretisation). Q_in
is identical across all five (only the side fraction varies), Q_sout scales linearly with
f_gap (last column matches 0.400/0.600/1.000/1.400/1.800 exactly), and mass is conserved
to ≤ 10⁻¹⁵ everywhere.

## Results table (matches `REFERENCE_OUTPUT.txt`)

`f_phys = f_gap/2`; `R_depth = δ_skim / [w_c√(f_phys/3)]`;
`valid = δ_skim / δ_general(f_phys)` (δ_general = exact duct law).

| frac | f_gap | f_phys | δ_skim | δ_bulge | K_inplane | R_depth | valid |
|---|---|---|---|---|---|---|---|
| 2e-4 | 1.985×10⁻⁴ | 9.924×10⁻⁵ | 0.543 | 0.484 | 1.000 | 0.942 | 0.991 |
| 3e-4 | 2.977×10⁻⁴ | 1.489×10⁻⁴ | 0.667 | 0.599 | 1.000 | 0.944 | 0.993 |
| 5e-4 | 4.962×10⁻⁴ | 2.481×10⁻⁴ | 0.864 | 0.781 | 1.000 | 0.947 | 0.995 |
| 7e-4 | 6.947×10⁻⁴ | 3.474×10⁻⁴ | 1.023 | 0.928 | 1.000 | 0.947 | 0.995 |
| 9e-4 | 8.932×10⁻⁴ | 4.466×10⁻⁴ | 1.160 | 1.055 | 1.000 | 0.947 | 0.995 |

- **CHECK 1 — δ_skim ∝ √f:** ratio constant to 0.75%; fit δ_skim = 54.82·√f_phys,
  **R² = 0.99992**. The fit constant ÷ the deep-limit constant w_c/√3 = 57.735 gives
  0.950 = the mean depth factor R — an independent route to the same R.
- **CHECK 2 — K_inplane = 1.000** at all five fractions.
- **CHECK 3 — R_depth = 0.942 … 0.947** (spread 0.53%): fraction-independent, as a
  geometric factor must be, and equal to the analytic R(α = 1.4). **Validation =
  0.991 … 0.995:** δ_skim equals the exact rectangular-duct law (evaluated at f_phys) to
  < 1% at every point.

## Why K_inplane = 1.000 is a real null, not a degeneracy

δ_bulge is the dividing line's closest approach to centre, taken as a `max` flux-width
over a window that starts at the upstream station; δ_skim_mid is the value *at* that
station. So δ_bulge ≥ δ_skim_mid by construction, and **K = 1.000 only when the line never
bulges outward** — genuine funnelling would read K > 1. Two stress tests confirm the null:

1. **Robust to sampling location.** K = 1.000 at every upstream station from x = 62 to
   x = 84, at all five fractions (see the robustness block).
2. **Monotonic thinning, zero bulge.** The dividing-line skim width is flat far upstream,
   then decreases monotonically into the gap; the largest downstream *increase* in skim
   width is **0.0000 µm** at every fraction.

This is the line-sink behaviour predicted by the analytic Hele-Shaw limit, here confirmed
in the full 3-D deep geometry.

## What the sweep does and does not settle

- **Settles:** no in-plane funnelling at any fraction (K_inplane = 1); √f scaling; and that
  δ_skim matches the exact developed-duct law (depth factor R ≈ 0.95) at f_phys.
- **Does not settle:** the absolute size cutoff. With K_inplane = 1 and K_depth = R ≈ 0.95,
  the geometric prefactor is κ_steric·1·0.95 ≈ 1.9 (κ_steric = 2), still well short of the
  empirical κ ≈ 5. The remaining enhancement is attributable to the *effective* fraction
  (f_eff several × nominal, from the resistance-network argument that an isolated unit cell
  over-extracts) and likely cumulative finite-particle displacement — **not** funnelling.
  Carry `f_physical = f_gap/2` into any of that analysis.

## Correction note (record of a fixed error)

An earlier version of this archive reported `K_depth ≈ 0.67` and described a "~25%
under-development offset," with a prediction that a developed-inlet sampling plane would
raise K_depth toward the analytic R. A confirmatory developed-inlet run was performed and
**falsified that prediction**: δ_skim was unchanged (0.864 → 0.858 µm), proving the flow
at x = 75 was already developed for the near-wall layer (the skim develops over a far
shorter length than the centreline). The recheck that followed found the real cause — the
factor-of-2 fraction convention above. δ_skim was correct all along; the "offset" was an
artifact of comparing it to the law at the un-halved f_gap. The funnelling result
(K_inplane = 1) and the √f scaling were never affected, as both are ratios.
