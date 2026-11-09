# Run B f-Sweep — generation procedure (COMSOL 6.3 UI)

**Procedure of record** for how the five fraction solves in this archive were produced.
Geometry is the deep Run B family (w_c = 100, depth = 140 µm, α = 1.4); only the per-gap
boundary fraction changes between points.

**Goal of the sweep:** vary the per-gap fraction to (1) validate the skim-extraction method
(δ_skim ∝ √f), (2) show K_inplane ≈ 1 holds at every operating point (no in-plane
funnelling), and (3) confirm δ_skim reproduces the exact developed-duct skim law.

**Factor of 2 (carry this through):** COMSOL reports `f_gap = intop_sout/intop_in` on the
quarter cell. The physical one-wall/full-duct fraction the analytic law uses is
**f_physical = f_gap/2** (the lamina is 2× the quarter lamina by z-symmetry; the duct is 4×
the quarter inlet by y- and z-symmetry). δ_skim, K_inplane, and the √f scaling are
unaffected; only absolute fraction/cutoff predictions need the ÷2.

**Time:** changing the fraction does **not** trigger a re-mesh — it only changes a boundary
velocity — so each point is a ~90 s re-solve on the existing mesh. Five points + exports is
well under an hour.

---

## Step 0 — Make a separate sweep file (keep canonical Run B frozen)

1. **File → Open** → your canonical model (currently frac = 5e-4).
2. **File → Save As** → a sweep working copy. **All steps below happen in this copy.**
3. **Re-use the 5e-4 point.** The 5e-4 solve already *is* the f ≈ 4.96e-4 point; you only
   need to compute the other fractions.

**Points to run:** `2e-4`, `3e-4`, `5e-4` (re-used), `7e-4`, `9e-4` — five points give a
clean √f line, a K_inplane(f) trend, and a fraction-independent depth-factor check.

---

## Step 1 — Set up the two Export nodes once (then re-fire per fraction)

Two planes are exported per fraction. Set both Export nodes up now so each point is just
"change the target folder, click Export."

### 1a. Upstream plane (feeds δ_skim)
1. **Results → right-click Export → Data** → node "Data 1".
2. **Settings:**
   - **Dataset:** `Study 1/Solution 1`
   - **Points to evaluate in:** `Grid`
   - **x:** `75`  **y:** `range(0,0.1,50)`  **z:** `range(0,1,70)`
   - **Expressions:** add three rows — `u`, `v`, `w`
   - **Output → Filename:** point it into the fraction's folder (any name; the post-processor
     classifies the plane automatically).

### 1b. Mid-depth plane (feeds δ_bulge / K_inplane)
1. **Results → right-click Export → Data** → node "Data 2".
2. **Settings:**
   - **Dataset:** `Study 1/Solution 1`
   - **Points to evaluate in:** `Grid`
   - **x:** `range(60,0.5,140)`  **y:** `range(0,0.1,75)`  **z:** `70`
   - **Expressions:** `u`, `v`, `w`

> The post-processor accepts whichever layout COMSOL writes (column or grid-block), and it
> identifies upstream vs mid-depth by which coordinate is fixed — so file names are free.

---

## Step 2 — The per-fraction loop

Do this block once for **each** fraction:

### 2a. Change the fraction
1. **Global Definitions → Parameters 1.**
2. Find the **`frac`** row, type the new value (e.g. `2e-4`), press **Enter**.
3. **Confirm it propagated:** the **`U_side`** Value should change (it is `frac*Q0/A_side`);
   **`U_in`** should stay **0.095238 m/s**. Spot checks:
   - frac = 2e-4 → U_side ≈ **1.587e-5 m/s**
   - frac = 9e-4 → U_side ≈ **7.143e-5 m/s**

### 2b. Solve
1. **Study 1 → Compute** (F8).
2. Healthy = converges in **~4 iterations**, no warnings, DOF stays **292 917** (same mesh),
   in-core (~8–9 GB). If it takes 20+ iterations or memory spikes, recheck Step 2a — a
   mis-typed fraction or a parameter that lost its unit is the usual cause.

### 2c. Verify (record the achieved fraction)
1. **Results → Derived Values → Global Evaluation 1 → Evaluate (=).**
2. Read the **mass balance** `(intop_in−intop_cout−intop_sout)/intop_in` (should be
   ~1e-15 to 1e-5) and the **f_gap** `intop_sout/intop_in` (the achieved fraction, ~0.8%
   below the nominal `frac`). The report `.docx` records both; the summary script reads
   f_gap from it automatically.

### 2d. Export both planes
1. **Export → Data 1** (upstream) → into this fraction's folder → **Export**.
2. **Export → Data 2** (mid-depth) → into this fraction's folder → **Export**.

Repeat 2a–2d for each fraction. When done, each fraction's folder holds two `.txt` exports;
add the `.mph` (Save As per fraction, or one model holding the last solve) and the report.

---

## Step 3 — Post-process (one command)

Arrange the archive as below, then run the summary from the root. It discovers each
fraction subfolder, classifies the two exports, reads f_gap from each report, and prints the
table plus the three checks and the robustness block.

```
python cif_sweep_summary.py            # from inside the Run B Sweep folder
```

No file lists to edit. Requires `cif_funnel_postproc.py` beside it; `numpy >= 1.24`, `scipy`.
The output should match `REFERENCE_OUTPUT.txt`.

**Expected result:** K_inplane = 1.000 at every f; δ_skim ∝ √f (R² = 0.99992); R_depth ≈
0.95 and validation ≈ 0.99 (δ_skim equals the exact duct law at f_physical = f_gap/2).

---

## Step 4 — Archive layout

```
Run B Sweep/
├── README_sweep.md
├── cif_funnel_postproc.py
├── cif_sweep_summary.py
├── REFERENCE_OUTPUT.txt
├── CIF_RunB_Sweep_StepByStep.md   (this file)
└── 2e-4/ 3e-4/ 5e-4/ 7e-4/ 9e-4/  (each: .mph model, 2 .txt exports, .docx report)
```

The `5e-4` folder is the canonical single-point Run B; the sweep is the robustness layer.

---

## Optional — extend to other geometries

The same procedure (Save As → edit the geometry parameters → re-mesh → same pre-flight
gates → the same two exports → the same post-processor) extends to other aspect ratios as
confirmatory cross-checks of K_inplane ≈ 1. They are not required for the conclusions here.

### Alternative — COMSOL Parametric Sweep (more self-documenting)
**Study 1 → Parametric Sweep**, add `frac` with the value list (`2e-4 3e-4 5e-4 7e-4 9e-4`),
**Compute**. Point each Export node's **Dataset** at `Study 1/Parametric Solutions 1` and set
its **Parameter selection (frac)** per export. Tidier in the `.mph`, identical data.
