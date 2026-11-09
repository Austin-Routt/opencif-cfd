# cfd/ — Python CFD

An independent finite-element solve of the two COMSOL geometries of the dissertation, written with open-source tools so
that anyone can repeat it without a COMSOL licence. "Python CFD" is the name the dissertation and this archive give to
this solver. It is the author's own code, not a published package: the scripts in `src/`, built on scikit-fem
(finite-element assembly), PyAMG (algebraic multigrid) and SciPy (MINRES). It solves:

- the **single-gap quarter cell** of the June 2026 campaign (Runs A, B and C; Section 4.9, C4.9-P06), and
- the **five-gap outlet slab** at Dinh et al.'s V1 and V3 outlet widths (Section 4.9, C4.9-P07; Section 5.5).

Both solve steady Stokes flow with P2-P1 (Taylor-Hood) tetrahedra on graded tensor-product meshes, so every face of the
geometry is a mesh plane and the geometry is exact. The saddle-point system is solved by MINRES with a block
preconditioner (algebraic multigrid on the velocity block, the lumped pressure mass matrix on the pressure block).

The code was first written during the October 2026 reviews. This version (9 November 2026) is the one the dissertation
uses. Its changes are listed at the end.

## Install

Python 3.10 or later. The runs in `runs/` used Python 3.13.16, numpy 2.5.3, scipy 1.18.1, scikit-fem 12.0.2,
pyamg 5.4.1 and matplotlib 3.11.2 on Linux x86_64 (each run's `meta.json` records its versions).

    python3 -m pip install -r requirements.txt

## Check the installation (about one minute)

    python3 src/smoke_test.py

It solves Run B (developed inlet) and the V1 slab on coarse meshes and checks the skim width against the exact duct
law (within 3%), F against the converged value (within 0.003), the flow balance (1e-4) and the solver residual (1e-6).
It prints `SMOKE TEST PASSED` and exits with status 0.

## Repeat every run (about an hour on two cores, at most 5 GB of memory)

From the root of the archive (the folder that holds `comsol/` and `cfd/`):

    bash cfd/run_all.sh comsol
    python3 cfd/src/compare.py .
    python3 cfd/src/figures.py .

`run_all.sh` solves each case, post-processes it against COMSOL's exports of the same case and deletes the solution
file (`solution.npz`). A run whose folder already holds `metrics.json` is skipped, so the script can be restarted; to
repeat the archived runs from scratch, rename `cfd/runs` first. `compare.py` reads only the archived files and writes
`cfd/results/comparison.json`; `report.py` writes `cfd/results/COMPARISON.md` from it, and `figures.py` draws
`cfd/figures/`.

| Run folder | Case | Inlet | Fraction Q_side/Q_in | Mesh scales |
|---|---|---|---|---|
| `B_uni_s*` | Run B (w_c 100 µm, depth 140 µm) | uniform, as COMSOL Run B | 4.962120761577921 × 10⁻⁴ | 2, 1.414, 1 |
| `B_dev_s*` | Run B | developed duct profile | 4.8751210307083063 × 10⁻⁴ | 2, 1.414, 1 |
| `A_dev_s*` | Run A (52 µm, 140 µm) | developed, as COMSOL Run A | 4.8518 × 10⁻⁴ | 2, 1.414, 1 |
| `C_dev_s*` | Run C (100 µm, 80 µm) | developed, as COMSOL Run C | 4.8550 × 10⁻⁴ | 2, 1.414, 1 |
| `V1_s*`, `V3_s*` | five-gap slab | 100 Pa, as COMSOL | — | 3, 2, 1.7 |

The fractions are COMSOL's achieved values (the re-evaluated and re-solved Run B tables at full precision, and the
Run A and C reports), so the two codes are compared at equal fraction. Mesh scale s sets the first cell at the
gap-side wall (0.1 s µm in the single gap), the cells at the gap edges (s µm) and the bulk cells (8 s µm single gap,
5 s µm five gaps).

## Files

| File | What it does |
|---|---|
| `src/stokes_gap.py` | single-gap quarter cell: mesh, assembly, boundary conditions, MINRES; writes `solution.npz` and `meta.json` |
| `src/stokes_5gap.py` | five-gap slab (half width, half depth by symmetry); the stress form, so the pressure inlet and outlets carry COMSOL's traction condition |
| `src/postproc_gap.py` | probes the solution on the points of COMSOL's upstream and mid-depth exports; skim widths, the duct law, the dividing streamline, the gap-mouth lobes, for both codes with the same routines |
| `src/postproc_5gap.py` | flows, F, the ten probe pressures, the flow through each gap on COMSOL's gap-plane points, and the mid-depth plane |
| `src/skim.py` | the exact rectangular-duct law, the skim width (native grid and the campaign's resampling), the streamline tracer |
| `src/fivegap.py` | the slab geometry and the per-gap flow integral, shared by both codes' post-processing |
| `src/comsol_io.py` | readers for COMSOL's table and plane exports (both layouts) |
| `src/compare.py` | the comparison and the ASME V&V 20 grid estimates |
| `src/figures.py`, `src/pubstyle.py` | the figures, in the dissertation's style |
| `src/smoke_test.py` | the installation check |
| `run_all.sh` | every run, in order |
| `runs/<run>/` | `meta.json` (mesh, solver, versions), `metrics.json` (every result), the solve and post-processing logs, the traced dividing streamlines and, for the finest mesh of each case, the probed planes (`*.txt.gz`) |
| `runs/run_all.log`, `runs/peak_memory.txt` | the order and times of the runs, and the peak memory of each step |
| `results/comparison.json` | every number of the comparison |
| `figures/` | PDF (vector) and PNG (600 dpi) |

## Boundary conditions and what differs from COMSOL

- Single gap: inlet velocity (uniform or the developed duct profile, unit mean), center outlet at zero pressure, side
  outlet at a uniform outflow velocity that sets the fraction, symmetry at the channel midline (y = 0) and at mid-depth,
  no slip elsewhere. The weak form is the vector Laplacian; its natural outlet condition is zero pseudo-traction, which
  developed outflow satisfies, whereas COMSOL applies zero traction. The outlet lies 60 µm past the gap and 105 µm past
  the sampling plane. With a uniform inlet Python CFD imposes the inlet velocity up to the wall edges, so its
  inflow equals U_in A_in exactly, whereas COMSOL's falls 1.75% short; the comparison is made at each code's achieved
  fraction, which removes this difference.
- Five gaps: the stress form 2μ ε(u):ε(v), a pressure load on the inlet with zero tangential velocity (COMSOL's
  "normal flow"), traction-free outlets at zero pressure, symmetry at the midline and at mid-depth (COMSOL modeled the
  full depth; the mid-depth plane is a symmetry plane of the problem), no slip elsewhere.

## Changes from the October 2026 packages

1. Runs A and C use the developed inlet that COMSOL's Runs A and C use; the October runs used a uniform inlet.
2. Every single-gap run is made at COMSOL's achieved fraction instead of the nominal 5 × 10⁻⁴.
3. The post-processing no longer imports the campaign's script from a fixed path. `skim.py` implements the skim
   width from its definition and is checked against the campaign's reference output (Run B 0.864, Run A 0.461,
   Run C 0.804 µm with the campaign's resampling); both codes are measured with the same routines on the same points.
4. COMSOL's values are read from its exported tables and planes, not typed into the scripts.
5. The flow through each gap and the mid-depth plane are computed on COMSOL's export grids.
6. The observed order of the five-gap mesh study uses V&V 20's representative cell size, h proportional to N^(−1/3);
   the October package used the mesh-scale setting.
7. `--tol`, `--inlet` with checked choices, the library versions in every record, a smoke test, this README.
8. The name: the October packages and the first November drafts called this solver "the second code"; it is now
   "Python CFD" throughout. The figures take the colors of the dissertation's Chapters 1 and 2 (V1 blue, V3 vermillion;
   runs B blue, A vermillion, C green).
