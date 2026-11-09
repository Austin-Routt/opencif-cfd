# opencif-cfd: the COMSOL models and Python CFD of the CIF dissertation

[![DOI](https://zenodo.org/badge/1451755952.svg)](https://doi.org/10.5281/zenodo.24091946)

This repository holds every computational-fluid-dynamics file behind Austin Routt's dissertation on Controlled
Incremental Filtration (CIF) design: the COMSOL Multiphysics models, reports and exports, an independent open-source
solve of the same geometries, and the comparison of the two. It closes the CFD work of Sections 4.9 and 5.5, and
`DISSERTATION_MAP.md` says which file stands behind each figure, paragraph and number of the dissertation that uses
CFD. Its content is the CFD archive of 9 November 2026 (`OpenCIF_CFD_Archive_2026-11-09`), file for file, with the
files that make it a repository added (this README, `DISSERTATION_MAP.md`, `CITATION.cff`, the licences,
`checks/verify_sha256.py`, `.gitattributes`, `.gitignore`) and the dissertation's CFD figures
(`dissertation_figures/`).

"Python CFD" names the independent solver throughout: the author's finite-element code in `cfd/src/`, built on
scikit-fem and PyAMG. It is not a published package, and the October 2026 packages called it "the second code".

Read in this order:

1. this file (what is here, the results, what was checked, the limits);
2. `comsol/README.md` (each COMSOL model: settings, versions, what was re-evaluated or re-solved);
3. `cfd/README.md` (Python CFD: how to install, check and repeat it);
4. `cfd/results/COMPARISON.md` (every comparison, with the numbers);
5. `DISSERTATION_MAP.md` (where the dissertation uses each file);
6. `MANIFEST.md` (every file with its size and SHA-256).

## Related repositories

- [`cif-reference-implementation`](https://github.com/Austin-Routt/cif-reference-implementation): the preprint's
  reference implementation (Routt, 2026, engrXiv, https://doi.org/10.31224/7152). Its `COMSOL/` folder holds the
  five-gap models as the preprint published them (COMSOL 6.3, tag v3.20) and, since 9 November 2026, copies of the
  re-evaluated models that `comsol/five_gap/` holds here. The SHA-256 in `SHA256SUMS.txt` identify the files, so
  either copy can be checked against the other. The single-gap campaign, the Run B re-evaluations and Python CFD are
  only here.
- [`opencif`](https://github.com/Austin-Routt/opencif): the targets-to-photomask generator (release 1.0.0,
  https://doi.org/10.5281/zenodo.23620100). It does not run CFD; the dissertation compares its one-dimensional results
  with the files here.

## What the CFD is for

Two geometries are solved, each in steady Stokes flow (creeping flow; Reynolds numbers are below 10):

- **The single-gap quarter cell** (Runs A, B and C of June 2026). One gap of a CIF array, one quarter of the cross
  section, with a side outlet that draws a set fraction of the center flow through the gap. It checks two analytic
  results of Chapter 3: the depth factor of the skim width (the width of the near-wall layer that a gap withdraws)
  and the absence of in-plane funneling toward the gap.
- **The five-gap outlet slab** at Dinh et al.'s V1 and V3 outlet widths. Five gaps, both outlets at one pressure. Its
  filtrate fraction F is compared in Section 5.5 with Dinh et al.'s measured F.

## Results

**Five-gap slab, filtrate fraction F.** COMSOL's numerical uncertainty is ASME V&V 20-2009's two-grid estimate (factor
of safety 3); Python CFD's is its three-grid estimate (factor of safety 1.25).

| | V1 | V3 |
|---|---|---|
| COMSOL, coarse mesh (115,904 / 104,705 elements) | 0.50702 | 0.62799 (the stored May solution: 0.62808) |
| COMSOL, extra-coarse mesh (24,783 / 21,448 elements) | 0.50558 | 0.62552 |
| COMSOL u_num | 2.1 × 10⁻³ | 3.4 × 10⁻³ |
| COMSOL extrapolated (order 2 assumed) | 0.50783 | 0.62931 |
| Python CFD, finest of three meshes (149,904 / 156,600 tetrahedra) | 0.50814 | 0.62980 |
| Python CFD, observed order; extrapolated | 2.1; 0.50785 | 2.6; 0.62958 |
| Python CFD u_num | 3 × 10⁻⁴ | 2 × 10⁻⁴ |
| Dinh et al., measured (mean ± SD) | 0.505 ± 0.033 | 0.649 ± 0.036 |

The pressure difference across each gap agrees within 3.6% (about COMSOL's own change between its two meshes) and the
flow through each gap within 1.5%. The total inflow converges from opposite sides: Python CFD's finest mesh
carries 0.75% and 0.89% more than COMSOL's coarse mesh and 0.20% and 0.23% more than COMSOL's extrapolation.

**Single gap, skim width at x = 75 µm** (fraction 5 × 10⁻⁴; both codes measured with the same routine on the same
points, at COMSOL's achieved fraction and with COMSOL's inlet).

| | Run B, uniform inlet | Run B, developed inlet | Run A, developed | Run C, developed |
|---|---|---|---|---|
| COMSOL, one mesh (µm) | 0.8646 | 0.8596 | 0.4614 | 0.8058 |
| Python CFD, finest of three meshes (µm) | 0.8625 | 0.8575 | 0.4597 | 0.8012 |
| Python CFD / COMSOL | 0.9976 | 0.9976 | 0.9962 | 0.9943 |
| spread of Python CFD over its three meshes | 1.5 × 10⁻⁴ | 1.0 × 10⁻⁵ | 4.5 × 10⁻⁵ | 3.7 × 10⁻⁵ |
| Python CFD / exact duct law | 0.9938 | 0.9968 | 1.0000 | 0.9969 |
| COMSOL / exact duct law | 0.9962 | 0.9992 | 1.0038 | 1.0026 |

The developed inlet changes Run B's skim width at equal fraction by +0.30% in both codes. Python CFD is converged
to 1.5 × 10⁻⁴, so the 0.24% to 0.57% between the codes most likely sits in COMSOL's single mesh at the wall; it is
smaller than the 0.53% between the simulated and the analytic depth factors that Chapter 4 reports.

## What was checked, and what was found

**COMSOL (the files of 9 November 2026, against the closeout plan and the May and June reports).**

- Every table holds one row at 17 significant figures, and in every model the inflow equals the two outflows to
  10⁻¹⁵.
- The element counts are those of the May reports: 115,904, 24,783 and 21,448. The extra-coarse models now give their
  own F (0.50558 and 0.62552), which replaces the tables that cif-reference-implementation's v3.20 models had copied
  from the coarse models.
- V1 coarse reproduces its May tables to every printed digit (F = 0.5070248).
- V3 coarse: its stored May solution gives F = 0.6280814, against 0.6280649 in the tables it carried (the preprint's
  production solve on 104,651 elements). Its re-solve rebuilt the mesh at unchanged settings (104,705 elements,
  550,173 DOF, 8 min 42 s, 13.1 GB) and gives 0.6279926. The stored solution therefore matches the model's
  parameters within the spread of mesh realizations, 9 × 10⁻⁵ in F.
- The flow through the five gaps, integrated over the gap mid-plane export, equals the side-outlet flow within 0.03%
  in every model.
- On the mid-depth plane, a symmetry plane of the problem, the out-of-plane velocity is at most 0.2% of the peak speed
  on the coarse meshes and 0.4% on the extra-coarse ones, as a discretization error should behave. NaN appears only
  inside the wall.
- The Run B passage export matches the June mid-depth export at the 5,617 points they share to 10⁻¹⁶ m/s, so it reads
  the canonical stored solution, whose achieved fraction is 4.962120761577921 × 10⁻⁴.
- Run B with a developed inlet: inlet set as in Runs A and C, the same 56,526-element mesh, four Newton iterations to
  1.7 × 10⁻¹⁵; its achieved fraction, 4.8751210 × 10⁻⁴, and its skim width with the campaign's resampling,
  0.858 µm, are those of the run that the campaign README describes.
- A separate check, written from scratch without the scripts of `cfd/` (`checks/independent_recheck/`), recomputed
  every COMSOL number above from the raw exports and agreed to every printed digit: F, u_num, the gap flows, the
  half-lives, the validation distances and the four skim widths (its duct law, from a different Fourier expansion,
  agrees to 10⁻¹⁰ µm). It also found that the ten probe pressures of each model's Tables 4 and 5 equal the pressures in
  its mid-depth export.
- To note: the files were saved by COMSOL 6.4 (build 429), so they open in 6.4 or later; the V3 "before" folder has no
  Word report (the May report covers that solution); `v1 coarse mesh.png` was a byte-identical copy of
  `V1_coarse_mesh_xy.png` and is left out; the V3 model's description of `Ly_half` reads "89 µm" while its value is
  102 µm (a label only).

**Python CFD (`cfd/`).**

- The October packages ran Runs A and C with a uniform inlet, whereas COMSOL's Runs A and C use a developed inlet, so
  their agreement figures compared two different inlets. Every single-gap run is now made with COMSOL's inlet and at
  COMSOL's achieved fraction.
- The post-processing no longer imports the campaign's script from a fixed path; `skim.py` implements the skim width
  from its definition and reproduces the campaign's reference output (0.864, 0.461 and 0.804 µm for Runs B, A and C).
- Re-run with the final code, the five-gap meshes reproduce the October records to 2 × 10⁻⁸ (relative), so the code
  is deterministic across library versions.
- `smoke_test.py` passes in about one minute.

## Limits

- COMSOL's single-gap runs have one mesh each; Python CFD's three meshes per run stand in for a COMSOL refinement
  study, and the remaining difference between the codes, a few tenths of a percent, most likely sits in COMSOL's
  near-wall mesh.
- COMSOL's five-gap uncertainty rests on two mesh levels from presets that are not geometrically similar, which is
  why V&V 20's factor of safety of 3 applies. A finer COMSOL mesh does not fit in 16 GB of memory.
- Python CFD's meshes are graded tensor grids refined by a scale factor: structured but not geometrically
  similar, so its observed order depends on how the cell size is defined. This archive uses V&V 20's representative
  size h = (V/N)^(1/3).
- The traced mid-depth dividing streamline and the inflow and outflow lobes at the gap mouth sit at a no-slip wall and
  a corner; they differ between the codes by several percent and move with the mesh. Use them qualitatively.
- Neither geometry is Dinh et al.'s device (the slab has a different gap and no outlet path), so the slab's agreement
  with their measurements is consistency, not a validation of opencif.

## Layout

| Folder | Contents |
|---|---|
| `comsol/five_gap/` | the five-gap models: V1 and V3 at two meshes, and V3's coarse model before and after its re-solve |
| `comsol/single_gap/campaign/` | the June 2026 single-gap campaign (Runs A, B, C), as the author archived it |
| `comsol/single_gap/runB_5e-4_passage/` | the canonical Run B re-evaluated, with a fine export of the gap passage |
| `comsol/single_gap/runB_5e-4_developed/` | Run B re-solved with a fully developed inlet |
| `cfd/` | Python CFD, its runs, the comparison and its figures |
| `checks/independent_recheck/` | a from-scratch recomputation of the COMSOL numbers, which agrees to every printed digit |
| `checks/verify_sha256.py` | checks every file against `SHA256SUMS.txt` (any platform, Python 3) |
| `DISSERTATION_MAP.md` | each figure, paragraph and number of the dissertation that uses CFD, with its files here |
| `dissertation_figures/` | the dissertation's CFD figures (F3-5, F4-1, F4-2, F4-3, F4-5, F4-6, F5-1) as it prints them |

## Software

COMSOL Multiphysics 6.3 (build 420) built every model; the November 2026 files were saved by COMSOL 6.4 (build 429)
and open in 6.4 or later. Python CFD needs Python 3.10 or later with numpy, scipy, scikit-fem, pyamg and
matplotlib (`cfd/requirements.txt`).

## Integrity

`SHA256SUMS.txt` holds the SHA-256 of every file except itself and `MANIFEST.md`. From the repository's root, run
`python3 checks/verify_sha256.py` (any platform) or `sha256sum -c SHA256SUMS.txt` (Linux, macOS, Git Bash).
`.gitattributes` turns off line-ending conversion, so a checkout on Windows keeps every byte and every digest.

## Citing

`CITATION.cff` gives the citation; GitHub shows it under "Cite this repository". Each release is archived on Zenodo.
Cite the version you used: version 1.0.0, the files of 9 November 2026 that the dissertation uses, is
https://doi.org/10.5281/zenodo.24091947. The DOI https://doi.org/10.5281/zenodo.24091946 (the badge above) stands
for all versions and always resolves to the newest.

> Routt, A. H. (2026). *opencif-cfd: The COMSOL models and Python CFD of the CIF dissertation* (Version 1.0.0)
> [Computer software]. Zenodo. https://doi.org/10.5281/zenodo.24091947

## Author and licences

Austin Routt (ORCID 0000-0002-7286-2924). The code (`cfd/src/`, `cfd/run_all.sh`, `checks/`) is under the MIT licence
(`LICENSE`, also `cfd/LICENSE`). The COMSOL models, reports, exports, run records, results and figures are under the
Creative Commons Attribution 4.0 International licence (`LICENSE-DATA.md`); `CITATION.cff`, and so the Zenodo
record, names that licence, because Zenodo takes one. COMSOL Multiphysics is needed to open the
`.mph` files; the exports and reports can be read without it.
