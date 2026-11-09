# Where the dissertation uses these files

Every figure, table, paragraph and number of the dissertation that rests on CFD, with the files of this repository that
stand behind it. Paragraph IDs are those of the dissertation's outline (for example C4.9-P07 is Section 4.9, paragraph
7); number IDs are those of its numbers register. The dissertation always uses the files here, which are the updated
values: where the preprint (Routt, 2026) printed an earlier value, the preprint is updated, not the dissertation.
The dissertation cites this repository as version 1.0.0, https://doi.org/10.5281/zenodo.24091947 (reference key
`opencifcfd100` in its outline).

## Figures and tables

| Item | What it shows | Files here |
|---|---|---|
| F3-5 | the depth factor against aspect ratio, with the single-gap solves | `comsol/single_gap/campaign/Run {A,B,C} Sweep/REFERENCE_OUTPUT*.txt` (the fraction-averaged simulated depth factors) |
| T3-3 | the depth factor at the simulated aspect ratios and at the calibration widths | the same reference outputs |
| F4-1 | Run B's mid-depth plane, with the dividing streamline | `comsol/single_gap/campaign/Run B Sweep/5e-4/middepth_grid_frac5e-4.txt` |
| F4-2 | Run B's cross-section 25 µm upstream of the gap, against the duct field | `comsol/single_gap/campaign/Run B Sweep/5e-4/upstream_grid_frac5e-4 v3.txt` |
| F4-3 | the single-gap model of Runs C, B and A (a drawing) | the parameters of `comsol/single_gap/campaign/Run B Sweep/5e-4/Run B - Report.docx` and `comsol/single_gap/campaign/CIF_Funnelling_CFD_Unified_README.md` |
| F4-5 | the eleven single-gap solves against the duct law | the campaign's reference outputs and sweep summaries (`comsol/single_gap/campaign/`) |
| F4-6 | Python CFD against COMSOL, single gap and five-gap slab | `cfd/results/comparison.json` (the dissertation's copy, `kb/cfd_comparison_2026-11-09.json`, is byte-identical: SHA-256 `b57d690dd533964f192c…`, the full digest in `SHA256SUMS.txt`) |
| F5-1 | the five-gap slab and its center-to-side pressure differences | `comsol/five_gap/V1_coarse/table1.txt` to `table5.txt` and `comsol/five_gap/V3_coarse_resolved/Table1_after.txt` to `Table5_after.txt` |

The figures themselves, as the dissertation prints them, are in `dissertation_figures/`. The archive's own figures
(`cfd/figures/`: F4-1c, F4-5d, F5-1c and `two_code_five_gap`) are drawn by `cfd/src/figures.py`
from the same files; the dissertation's F4-6 combines F4-5d and F5-1c.

## Paragraphs

| Paragraph | What it states | Files here |
|---|---|---|
| C3.6-P02, C4.9-P06 | the analytic depth factor matches the single-gap solves to 0.53%; no in-plane funneling | `comsol/single_gap/campaign/` (Runs A, B, C: models, reports, exports, reference outputs); Run B with a developed inlet: `comsol/single_gap/runB_5e-4_developed/`; Python CFD's runs `cfd/runs/{B_uni,B_dev,A_dev,C_dev}_s*` |
| C4.9-P07 | the five-gap slab and its numerical uncertainty from two meshes | `comsol/five_gap/` (all five folders, with `comsol/README.md`); the two-grid estimate in `cfd/results/comparison.json` (`five_gap.grid.comsol`) |
| C4.9-P08 | Run B's dividing streamline and upstream field | as F4-1 and F4-2; the fine passage export `comsol/single_gap/runB_5e-4_passage/` |
| C4.9-P09 | Python CFD agrees with both COMSOL solves | `cfd/` (code in `cfd/src/`, runs in `cfd/runs/`, every comparison in `cfd/results/COMPARISON.md`); the independent recomputation `checks/independent_recheck/` |
| C4.11-P01, Table T4-4 (row 3) | the slab's numerical uncertainty as its own category | as C4.9-P07 |
| C5.5-P05 | the slab against Dinh et al.'s measured F; the pressure half-lives | the tables of `comsol/five_gap/V1_coarse/` and `comsol/five_gap/V3_coarse_resolved/`; `cfd/results/comparison.json` (`five_gap.compare`, `five_gap.half_lives`) |
| C6.7-P01, C6.8-P04 | the CFD as part of the research artifact; RQ4's answer | this repository as a whole |
| Appendix (V&V 20 subsection), erratum D109 | u_num of the slab from two meshes, the validation macros | `cfd/results/comparison.json` (`five_gap.grid.comsol`, `five_gap.compare.*.validation_comsol`) |

## Numbers

| ID | Value | What | Files here |
|---|---|---|---|
| N3-19 | 0.53% | the largest difference of the analytic and simulated depth factors | the campaign's reference outputs |
| N4-01 to N4-04 | 1; 0.891, 0.945, 0.981 | the in-plane funneling factor; the simulated depth factors of Runs C, B, A | the campaign's reference outputs |
| N4-05, N4-06 | 115,904; 104,705 | elements of the V1 and V3 slabs, finer mesh (V3 re-solved) | `comsol/five_gap/V1_coarse/report.docx`, `comsol/five_gap/V3_coarse_resolved/report_after.docx` |
| N4-17, N4-18 | 24,783; 21,448 | elements, coarser mesh | the reports in `comsol/five_gap/V1_extra_coarse/` and `V3_extra_coarse/` |
| N4-19, N4-20 | 2.1 × 10⁻³; 3.4 × 10⁻³ | u_num of F, V1 and V3 | `cfd/results/comparison.json`, `five_gap.grid.comsol` |
| N4-21, N4-22 | 1.4 × 10⁻³; 2.5 × 10⁻³ | change of F between the two meshes | the same |
| N4-23 | 0.30% | the developed inlet's change of Run B's skim width | `comsol/single_gap/runB_5e-4_developed/`; `comparison.json`, `single_gap.comsol.inlet_effect_B_native` |
| N4-24 to N4-27 | 2.1, 2.6; 0.50785, 0.62958 | Python CFD's observed order and extrapolated F | `comparison.json`, `five_gap.grid.V1_python`, `V3_python` |
| N4-28, N4-29 | 0.24% to 0.57%; 1.5 × 10⁻⁴ | Python CFD's skim widths below COMSOL's; their largest change between meshes | `comparison.json`, `single_gap.runs` |
| N5-01, N5-02 | 0.5070; 0.6280 | F of the slab, V1 and V3 (V3 re-solved) | the tables of `comsol/five_gap/V1_coarse/` and `V3_coarse_resolved/` (F = 1 − Q_c/Q_in) |
| N5-07, N5-08 | 0.944; 1.042 | the pressure half-lives, in gaps | the same tables (Tables 4 and 5, the probe pressures); `comparison.json`, `five_gap.half_lives` |

N5-03 to N5-06 (Dinh et al.'s measured F and its standard deviations), N5-20 and N3-22 to N3-23 (the one-dimensional
outlet identity) come from the papers and from opencif, not from CFD; the network half-lives (0.966 and 1.090 gaps)
come from the preprint.

## Earlier values

The preprint printed V3's F as 0.6281, its half-life as 1.041 gaps and their agreement with the network as 4.7%, from
the May 2026 production solve (104,651 elements). The re-solved model here gives 0.6280, 1.042 gaps and 4.6%. The three
V3 meshes at one setting differ by 9 × 10⁻⁵ in F, about 1/40 of u_num (`comsol/README.md`).
