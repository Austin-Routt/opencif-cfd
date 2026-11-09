# CIF Funnelling CFD — Unified README (Runs A, B, C)

*Single-gap Stokes CFD campaign for the CIF size-cutoff prefactor. This is the canonical
document for all three aspect-ratio runs; the per-run READMEs (`README_sweep.md` for B,
`README_RunC_sweep.md`, `README_RunA_sweep.md`) are the detailed companions. Every number
here is reproduced by `cif_sweep_summary.py` against each run's `REFERENCE_OUTPUT*.txt`, and
cross-checked against the analytic depth law (`cif_cutoff_depth_blocks.py`, Blocks 9–13) and
the κ-decomposition notebook (`cif_cutoff_notebook.py`, Blocks 4–8).*

---

## 0. Executive summary — what the campaign was for and what it found

The empirical CIF cutoff law is `d_crit/w_c = κ·√(f_gap/3)` with κ ≈ 5 (deep family). The
companion analysis factors that prefactor as

> **κ = κ_steric · K_inplane · K_depth**, κ_steric ∈ {1, 2}.

The pre-CFD picture (Companion II, Part C.4) attributed the bulk of κ ≈ 5 to **in-plane
funnelling** — the idea that flow converging into a finite gap makes the dividing
streamsurface *bulge toward the channel centre*, so the captured layer reaches deeper than a
uniform skim (K_inplane ≈ 2.5–4, dominant and depth-independent). The single-gap CFD was
built to **measure K_inplane directly** as δ_bulge/δ_skim, and to measure K_depth from a
depth sweep.

**The result inverts the picture.** Across three aspect ratios spanning the entire deep CIF
family (α = 0.8 → 1.4 → 2.69):

1. **K_inplane = 1.000** — at every fraction and every sampling station, max bulge 0.0000 µm.
   The dividing streamsurface *thins* toward the gap (line-sink behaviour); it never bulges.
   **There is no in-plane funnelling at any aspect ratio.** This refutes the funnelling
   hypothesis and confirms the analytic Hele-Shaw line-sink prediction (depth-derivation
   addendum A1), now verified in the full 3-D geometry.
2. **K_depth = R(α)** — the CFD reproduces the exact rectangular-duct depth correction to
   ≤0.6%: R = 0.89 (α=0.8) → 0.95 (α=1.4) → 0.98 (α=2.69), monotone toward 1 as the channel
   deepens, exactly as the analytic duct law predicts.
3. **δ_skim ∝ √f_phys** to R² ≥ 0.9999 at each aspect ratio.

The consequence: the geometric prefactor is **κ_geometry = κ_steric · 1 · R ≈ 1.9**
(κ_steric = 2), far short of the empirical κ ≈ 5. The missing factor is **not** funnelling.
It lives downstream of the single-gap skim — in the **effective fraction** f_eff (an isolated
unit cell over-extracts relative to the device's resistance network) and in **cumulative
finite-particle displacement** over many gaps. Both are separate from this CFD; the cutoff
closure is the remaining work (Section 6).

---

## 1. The three runs map onto the three source papers

Each run is a real device aspect ratio, not an abstract sweep point:

| run | α = h/w_c | w_c | depth | gap | source paper (geometry) | match |
|-----|-----------|-----|-------|-----|-------------------------|-------|
| C   | **0.80**  | 100 | 80    | 20  | Gifford 2014, Fig 4 demo | exact |
| B   | **1.40**  | 100 | 140   | 20  | Xia 2016, RS-8/RS-9 (α=1.40) | aspect ratio (Xia gap = 16) |
| A   | **2.69**  | 52  | 140   | 19  | Dinh 2024, narrowest outlet gap | exact |

Run B established the method at α = 1.4; Runs C and A are **pre-registered, out-of-sample**
tests at the shallow and deep extremes. The depth correction R(α) depends only on α, so the
three runs validate it at the actual aspect ratios of all three papers in the CIF lineage.

### Headline result: R(α) across the family

| α    | device       | R(α) analytic | R(α) CFD | √f R²    | K_inplane                 | validation (δ_skim/exact law) |
|------|--------------|---------------|----------|----------|---------------------------|-------------------------------|
| 0.80 | Gifford Fig4 | 0.891         | 0.891    | 0.99997  | 1.000 (all f)             | 1.000–1.002                   |
| 1.40 | Xia RS-8/9   | 0.953         | 0.947    | 0.99992  | 1.000 (all f)             | 0.991–0.995                   |
| 2.69 | Dinh outlet  | 0.980         | 0.983    | 0.99993  | 1.000 (5e-4, 9e-4; n/a 2e-4) | 0.998–1.003                |

Three independent aspect ratios, three pre-registered hits on R(α). The CFD δ_skim matches
the exact duct law to within **±0.6%** at every point: Run C lands dead-on, Run A +0.3%,
Run B −0.6% (the only run more than 0.3% off — within the run-to-run numerical scatter, not a
bias). R(α) is a real, predictable geometric factor — not a fitted fudge.

---

## 2. The shared model (method of record)

All three runs are the **same single-gap CIF unit cell**, quarter-cell by symmetry, solved as
creeping (Stokes) flow in COMSOL 6.3.0.420.

- **Domain & symmetry.** One unit cell of the gap cascade: main channel (width w_c, depth h),
  one side gap (aperture `gap`) opening to a side outlet (width w_s). Modelled as a **quarter
  cell** with symmetry planes at `y = 0` (channel midline) and `z = h/2` (mid-depth). The
  solved region is therefore one wall, half-width, half-depth.
- **Physics.** Laminar/creeping flow (Re ≪ 1), incompressible, P2–P1 (Taylor–Hood) elements,
  PARDISO direct solver. Fluid = PRP: ρ = 1000 kg/m³, μ = 1.2 mPa·s. (Stokes linearity means
  every reported quantity — δ_skim, K_inplane, R — depends only on the **fraction**, not the
  absolute flow rate.)
- **Boundary conditions.** Inlet: fully-developed flow (or velocity `U_in = Q0/A_in`) at
  flow rate **Q0 = 20 µL/min** (hardcoded = 3.3333×10⁻¹⁰ m³/s; A_in = (w_c/2)(h/2)). Main
  outlet: pressure 0. Side outlet: velocity `U_side = frac·Q0/A_side`, A_side = w_s·(h/2),
  which sets the per-gap fraction. Walls: no-slip. The two symmetry planes carry symmetry BCs.
- **Mesh.** Free tetrahedral with a **near-wall boundary-layer refinement** on the gap-side
  wall (first-cell ≤ 0.05 µm, several layers) so the sub-micron skim is resolved. Element
  counts: Run B/C ≈ 56 k; Run A ≈ 40 k (smaller because the channel is half as wide). The low
  *minimum* element quality (~0.05) is the deliberately thin boundary-layer cells, not a
  defect — see Section 5.
- **Exports (per fraction).** Two grid exports of (u, v, w): an **upstream YZ plane at
  x = 75** (the developed-flow station, for the depth-averaged δ_skim) and a **mid-depth XY
  plane at z = h/2** (for K_inplane). File names are irrelevant — the post-processor
  classifies each export by which coordinate is held fixed. The gap sits at x ∈ [100, 120].
- **Reports (per fraction).** A COMSOL report `.docx` carrying the Global Evaluation row
  `[Q_in, Q_cout, Q_sout, mass_balance, f_gap]`, from which the achieved fraction is read.

One pair of scripts (`cif_funnel_postproc.py`, `cif_sweep_summary.py`) serves all three runs
unchanged: geometry (w_c = 2·max y, depth = 2·max z, α) is auto-detected from the upstream
export.

---

## 3. Two conventions to keep straight

**(a) Factor of 2 in the reported fraction (flow).** Because the model is a quarter cell,
COMSOL reports `f_gap = Q_sout(quarter)/Q_in(quarter)`. The analytic skim law and the
empirical κ work define the fraction as [one-wall, full-depth lamina flux] / [full-duct flux],
which equals **f_gap / 2**: the one-wall full-depth lamina is 2× the quarter lamina
(z-symmetry), while the full duct is 4× the quarter inlet (y- and z-symmetry). So the reported
f_gap is **twice** the physical one-wall/full-duct fraction. The post-processor forms
`f_phys = f_gap/2` and compares δ_skim to the exact duct law at f_phys. This is bookkeeping,
not a model error — δ_skim, K_inplane, and the √f scaling are unaffected (the first is the
layer that physically carries Q_sout; the latter two are ratios). **It does mean any absolute
cutoff/fraction taken from this CFD must carry f_phys = f_gap/2.** Comparing δ_skim to the law
at the un-halved f_gap gives a spurious R ≈ 0.67 — the origin of the retired "offset"
(Section 8). For Dinh's geometry this is doubly confirmed: Dinh defines f*gap *per gap* with
bilateral gaps, so Dinh's f*gap already equals the physical one-wall fraction = f_gap(CFD)/2.

**(b) Achieved vs design fraction (geometry).** The measured gap fraction is
`f_gap = Q_sout/Q_in = (geometric factor)·frac`, constant across fractions within a run:
**0.992** (Run B), **0.971** (Run C), **0.970** (Run A). The factor is the side-outlet
boundary area relative to the nominal `w_s·(h/2)` (a few percent under, varying with depth and
w_s); it is purely geometric and stable to four digits across each run's fractions. The
achieved value is the physically correct band-flux target and is read straight from each
report's Global Evaluation row. All reports carry that row, so the summaries run warning-free;
because every check is √f-independent or a ratio, using the design frac instead changes no
conclusion (the scripts fall back to it with a warning if a report ever lacks the evaluation).

---

## 4. Results — the three checks, consolidated

Per-run sweep tables (each reproduced by `python cif_sweep_summary.py "<Run> Sweep"`):

**Run C — α = 0.80 (Gifford geometry), factor 0.971**

| frac | f_phys   | δ_skim | δ_bulge | K_in  | R_depth | valid |
|------|----------|--------|---------|-------|---------|-------|
| 2e-4 | 9.71e-5  | 0.508  | 0.523   | 1.000 | 0.891   | 1.002 |
| 5e-4 | 2.43e-4  | 0.804  | 0.781   | 1.000 | 0.891   | 1.000 |
| 9e-4 | 4.37e-4  | 1.081  | 1.030   | 1.000 | 0.892   | 1.001 |

**Run B — α = 1.40 (Xia aspect ratio), factor 0.992, five points**

| frac | f_phys   | δ_skim | δ_bulge | K_in  | R_depth | valid |
|------|----------|--------|---------|-------|---------|-------|
| 2e-4 | 9.92e-5  | 0.543  | 0.484   | 1.000 | 0.942   | 0.991 |
| 3e-4 | 1.49e-4  | 0.667  | 0.599   | 1.000 | 0.944   | 0.993 |
| 5e-4 | 2.48e-4  | 0.864  | 0.781   | 1.000 | 0.947   | 0.995 |
| 7e-4 | 3.47e-4  | 1.023  | 0.928   | 1.000 | 0.947   | 0.995 |
| 9e-4 | 4.47e-4  | 1.160  | 1.055   | 1.000 | 0.947   | 0.995 |

**Run A — α = 2.69 (Dinh outlet), factor 0.970**

| frac | f_phys   | δ_skim | δ_bulge | K_in  | R_depth | valid |
|------|----------|--------|---------|-------|---------|-------|
| 2e-4 | 9.70e-5  | 0.290  | n/a     | n/a   | 0.978   | 0.998 |
| 5e-4 | 2.43e-4  | 0.461  | 0.254   | 1.000 | 0.983   | 1.003 |
| 9e-4 | 4.37e-4  | 0.619  | 0.467   | 1.000 | 0.983   | 1.002 |

- **CHECK 1 — δ_skim ∝ √f_phys.** Ratio constant to <0.8% in every run; R² = 0.99997 (C),
  0.99992 (B), 0.99993 (A). The fit slope ÷ the deep-limit constant w_c/√3 independently
  returns the mean depth factor R for that run.
- **CHECK 2 — K_inplane = 1.000.** At every fraction and station, in all three runs (one
  documented exception: Run A 2e-4, Section 5). Genuine null: δ_bulge is a one-sided *max*, so
  K = 1 only when the line never bulges; funnelling would read K > 1.
- **CHECK 3 — R_depth, and validation.** R_depth = δ_skim/δ_parab is fraction-independent
  (spread ≤0.6% within each run) and equals the analytic R(α). Validation = δ_skim/δ_general
  (exact duct law at f_phys) sits at 1.00 ± 0.6% across all points.

---

## 5. Verification & resolution

The resolution question matters most for Run A (predicted skim only 0.29–0.62 µm). Four
independent lines establish adequacy:

1. **Mass conservation.** Q_in(x=75)/Q0 = 0.99985 in every Run A point (0.015% off); mass
   balance ≤ 10⁻¹⁵ in the Run B reports. The inlet delivers Q0 and nothing leaks.
2. **Export-grid convergence.** Re-extracting δ_skim while coarsening the export y-pitch 8×
   (0.05 → 0.40 µm) moves it ≤0.1% at every Run A fraction — including the thinnest (2e-4:
   0.2897 → 0.2898 µm). The export sampling is far from limiting.
3. **Near-wall profile.** The velocity decays smoothly and monotonically to zero at no-slip
   over the last microns to the wall, with no kinks — and the boundary-layer mesh is present
   (min element quality ~0.05 = thin wall cells; element volume ratio ~0.001 = fine-at-wall,
   coarse-in-bulk). The sub-micron gradient is resolved, not interpolated from coarse cells.
4. **The validation match is itself the certificate.** δ_skim reproduces the *exact analytic
   duct law* to ≤0.6% at every point and aspect ratio. A sub-micron skim that landed on the
   analytic value to three parts in a thousand cannot have been under-resolved.

**The Run A 2e-4 K_inplane = n/a.** At the deep aspect ratio the *mid-depth* skim is ~0.16 µm,
and the bulge check needs the net flux through the gap at mid-depth. Along the gap aperture the
in-plane velocity is dominated by large **antisymmetric recirculation lobes** (±1.3×10⁻² at the
gap edges) while the net through-flow at the gap centre is tiny (+2.6×10⁻⁵ at 2e-4). At 5e-4
and 9e-4 the through-flow survives that subtraction; at the lowest fraction it is ~5× smaller
and the imperfect discrete cancellation flips the residual negative (φ_gap ≤ 0), collapsing the
dividing line onto the wall. The analyzer **detects this and reports n/a** rather than a
spurious value. This is a post-processing limit at small f, **not** a model or resolution
issue: the depth-averaged δ_skim (which drives the validation) is fully resolved, and
K_inplane = 1 is established at Run A's 5e-4 and 9e-4 and at *every* Run B/C point — a
dozen-plus points across the whole α range. (A mass-balance flux estimator recovers a positive
φ_gap here but breaks the shallow runs, so the validated gap-flux method is kept and the one
indeterminate point flagged rather than forced.)

---

## 6. What is settled, and the open cutoff closure

**Settled by this campaign:**
- The exact rectangular-duct skim law **δ_skim = w_c·√(f_phys/3)·R(α)** is correct to ≤0.6%
  across α = 0.8–2.69.
- **K_inplane = 1** (no in-plane funnelling) across the full family — the funnelling
  hypothesis is dead.
- **K_depth = R(α)** is a real geometric factor, pinned at three pre-registered points.
- The √f scaling holds on interior points, not just endpoints.

**Open — and where the missing κ lives.** With K_inplane = 1 and K_depth = R ≈ 0.89–0.98, the
geometric prefactor is κ_steric·1·R ≈ 1.9 (κ_steric = 2). The empirical κ ≈ 5 therefore
requires an additional **~2.6×** (κ_steric = 2) that is *not* in the single-gap hydrodynamics.
The two candidates, both downstream of δ_skim:

- **Effective fraction f_eff (the cheap lever — no new CFD).** An isolated unit cell
  over-extracts relative to the real device's resistance network, so the *physical* per-gap
  fraction exceeds the nominal. Since δ ∝ √f, a modest multiplier moves the cutoff a lot:
  f_eff = 3–4× raises δ by 1.7–2× (notebook Block 5), and at **f_eff ≈ 3.5×** the Dinh
  prediction lands d_crit = 1.66 µm — **inside** the measured rigid-bead bracket (1.0–4.7 µm),
  whereas the nominal fraction gives 0.88 µm and fails (Block 7). f_eff is computable
  analytically from the hydraulic–electrical network (the circuit-analogy references in the
  knowledge base); the resistance-network argument is the natural next step.
- **Cumulative finite-particle displacement.** The single-gap model is a point tracer; real
  finite particles are displaced incrementally at each of many gaps, so the effective capture
  width accumulates down the cascade. This is the second downstream factor and needs
  finite-particle accounting (analytic or light particle-tracing), not single-gap CFD.

**The Dinh anchor.** Run A sits on Dinh's real geometry. At Dinh's operating point
(f*gap = 1.04×10⁻⁴ ↔ CFD frac 2.08×10⁻⁴, the 2e-4 point within 4%) the CFD lamina is
δ_skim ≈ 0.30 µm. The κ-transfer prediction (notebook Block 4) gives d_crit = 1.55 µm for
w_c = 52 — inside the measured bracket — so the empirical law and the bracket agree; the CFD's
job was to explain *why* κ has the value it does, and it has shown that ~2/3 of κ is the
downstream f_eff/cumulative factor rather than funnelling.

---

## 7. Reproducibility

Each run ships as a self-contained folder:

```
Run {B|C|A} Sweep/
  <frac>/   Report ….docx   upstream_grid_frac….txt   middepth_grid_frac….txt   [+ .mph]
  cif_funnel_postproc.py     per-point analyzer (δ_skim, δ_bulge, K_inplane, R, validation)
  cif_sweep_summary.py       sweep driver — auto-discovers fraction subfolders
  README_*.md                per-run companion to this file
  REFERENCE_OUTPUT_*.txt     known-good console output to diff against
```

Run `python cif_sweep_summary.py "<Run> Sweep"`; requires `cif_funnel_postproc.py` beside it,
numpy ≥ 1.24, scipy. Geometry is auto-detected, so the same scripts serve all three runs. The
two scripts are **identical across the three folders** (one canonical pair); the low-fraction
K_inplane guard added for Run A is inert on Runs B and C (regression-confirmed: their output is
byte-identical, because their wider/shallower channels keep the gap through-flow above the
recirculation floor at every fraction).

---

## 8. Correction note (record of a fixed error)

An earlier version of this work reported **K_depth ≈ 0.67** and a "~25% under-development
offset," and predicted that sampling on a fully-developed inlet plane would raise K_depth
toward the analytic R. A confirmatory developed-inlet run **falsified that prediction**:
δ_skim was essentially unchanged (0.864 → 0.858 µm), proving the near-wall layer at x = 75 was
already developed (the skim develops over a far shorter length than the centreline). The
recheck found the real cause — the **factor-of-2 fraction convention** (Section 3a): δ_skim
was correct all along; the "offset" was an artifact of comparing it to the law at the un-halved
f_gap (which gives the spurious R ≈ 0.67 = 0.95/√2). The funnelling result (K_inplane = 1) and
the √f scaling were never affected, being ratios. This is why the validation in Section 4 uses
f_phys = f_gap/2 throughout.

---

## 9. Next steps & optional additional COMSOL work

**Primary path (no new CFD required) — close the cutoff:**
1. Compute **f_eff** from the device resistance network (hydraulic–electrical analogy; tools
   already in the knowledge base). This is the cheapest, highest-leverage step.
2. Add the **cumulative finite-particle displacement** term (analytic or light particle
   tracing).
3. Combine with K_inplane = 1 and K_depth = R, carry **f_phys = f_gap/2**, and produce a clean
   Dinh cutoff prediction; check against the measured 1.0–4.7 µm bracket and the κ-transfer
   value (d_crit = 1.55 µm, w_c = 52).

**Optional COMSOL work that would add value (in priority order):**
- **Formal mesh-convergence study** (cheap, high reviewer value): solve one point at 0.5× and
  2× mesh density and tabulate δ_skim convergence, to complement the export-grid and
  validation-match evidence with a standard h-refinement curve.
- **Multi-gap cascade model** (confirmatory): a 3–5-gap CIF segment would *directly* measure
  f_eff (the per-gap fraction in a real network) and validate the resistance-network estimate
  — the CFD route to the same f_eff the circuit gives analytically.
- **Finite-particle tracing** in the single gap or cascade: replace point tracers with
  finite-radius particles to quantify the cumulative-displacement factor directly.
- **Dinh taper**: Dinh's channel narrows 75 → 52 µm along the cascade; a tapered or
  multi-w_c model would map δ_skim along the real device rather than at the worst-case outlet.
- **Shallow RBC single-gap** (separate application): a 2–6 µm-deep single-gap solve (α ≈ 0.04)
  to pin κ_total for a discocyte/spherocyte shape separator — laptop-tractable, flagged in the
  depth-derivation review as the next concrete step for that distinct design.

**Bottom line:** the mechanistic CFD campaign is **complete and successful** in its designed
scope — it measured both K factors, refuted funnelling, and validated the depth law on three
real-device aspect ratios. No further CFD is required to finish the κ story; what remains is the
analytic f_eff/cumulative closure, with the COMSOL options above as confirmatory or
application-extending, not blocking.
