# COMSOL files

Every file here was made by Austin Routt in COMSOL Multiphysics on a 4-core, 16 GB Windows 11 computer. The models were
built in COMSOL 6.3 (build 420). The November 2026 files were opened, re-evaluated or re-solved and saved in COMSOL 6.4
(build 429), so they open in 6.4 or later. The 6.3 originals of the five-gap models stay in the public repository
`cif-reference-implementation` (folder `COMSOL/`, tag v3.20), which since 9 November 2026 also holds copies of the
re-evaluated models in `five_gap/`; their SHA-256 in `../SHA256SUMS.txt` identify the same files. File names inside each folder are the ones COMSOL wrote
on the author's computer; `../MANIFEST.md` lists every file with its size and SHA-256.

## five_gap/ — the outlet slab of Section 4.9 (C4.9-P07) and Section 5.5

The model holds half of the width (symmetry plane at the channel midline) and the full depth. Center channel
52 µm (26 µm modeled), side channel 44 µm (V1) or 57 µm (V3), wall 19 µm across the flow, five gap openings 22 µm
along the flow, posts 27 µm, pitch 49 µm, depth 140 µm, length 272 µm; side channel closed upstream. Creeping Flow,
P2 + P1 elements, PARDISO, fully coupled, four Newton iterations. Inlet (center, x = 0) 100 Pa with normal flow; both
outlets 0 Pa with backflow suppressed; viscosity 1.2 mPa·s.

| Folder | Mesh (physics-controlled presets) | Elements | DOF | Solution | Newton error |
|---|---|---|---|---|---|
| `V1_coarse` | global Extra coarse, walls Coarse, 1 boundary layer | 115,904 | 606,563 | stored, solved 7 May 2026 | 4.8 × 10⁻¹² |
| `V1_extra_coarse` | global Extremely coarse, walls Extra coarse | 24,783 | 138,474 | stored, solved 7 May 2026 | 7.6 × 10⁻¹¹ |
| `V3_coarse_resolved` | as V1 coarse | 104,705 | 550,173 | **re-solved 9 Nov 2026** (8 min 42 s, 13.1 GB) | 6.8 × 10⁻¹² |
| `V3_coarse_stored_May7` | as V1 coarse | 104,206 | 547,783 | stored, solved 7 May 2026 | 4.9 × 10⁻¹² |
| `V3_extra_coarse` | global Extremely coarse, walls Extra coarse | 21,448 | 121,788 | stored, solved 7 May 2026 | 3.9 × 10⁻¹⁰ |

Each folder holds:

- the model (`*_reeval*.mph`) and, except `V3_coarse_stored_May7`, its complete Word report;
- five table exports at 17 significant figures (Results → Derived Values → Clear and Evaluate All): Table 1 Q_in,
  Table 2 Q_c_out (center outlet), Table 3 Q_s_out_top (side outlet), Table 4 the center pressure at
  x = 38, 87, 136, 185, 234 µm (y = 0, z = 70 µm), Table 5 the side pressure at the same x (y = 67 µm in V1, 73.5 µm
  in V3). Each table holds one row. Flows are for the modeled half width at full depth, in m³/s;
- the mid-depth plane z = 70 µm (`*_middepth_z70*.txt`; x 0 to 272 µm every 0.5 µm, y 0 to 89 or 102 µm every
  0.5 µm; u, v, w in m/s and p in Pa; NaN inside the wall);
- the gap mid-plane y = 35.5 µm (`*_gapplane_y35p5*.txt`; x 0 to 272 µm every 0.25 µm, z 0 to 140 µm every 1 µm);
- snapshots: the mesh seen from above (`*_mesh_xy*.png`) and the speed on the mid-depth plane, with and without a
  title (`*_middepth_speed*.png`), 3000 × 2250 px at 300 dpi.

**V3 coarse.** Its model file's global definitions were edited 19 s after its 7 May solve ended, so it was re-solved.
The re-solve rebuilt the mesh from the unchanged settings, to 104,705 elements against the stored solution's 104,206.
`V3_coarse_stored_May7` holds the file and exports as they were before the re-solve; `V3_coarse_resolved` holds the
re-solved model. The three V3 coarse solutions known (the preprint's production solve on 104,651 elements, whose
tables the 7 May file still carried; the stored 104,206; the re-solved 104,705) give F = 0.6280649, 0.6280814 and
0.6279926: a spread of 9 × 10⁻⁵ from mesh realization at fixed settings, about 1/40 of the two-grid uncertainty
below.

### Results (from the table exports; `../cfd/results/comparison.json` holds every digit)

| | V1 coarse | V1 extra coarse | V3 coarse (re-solved) | V3 coarse (stored) | V3 extra coarse |
|---|---|---|---|---|---|
| Q_in (m³/s) | 3.213970 × 10⁻¹⁰ | 3.182220 × 10⁻¹⁰ | 3.644964 × 10⁻¹⁰ | 3.645205 × 10⁻¹⁰ | 3.599678 × 10⁻¹⁰ |
| F = 1 − Q_c/Q_in | 0.5070248 | 0.5055835 | 0.6279926 | 0.6280814 | 0.6255172 |
| (Q_in − Q_c − Q_s)/Q_in | 7 × 10⁻¹⁶ | −8 × 10⁻¹⁷ | −1 × 10⁻¹⁶ | 1 × 10⁻¹⁵ | −3 × 10⁻¹⁶ |
| sum of the five gap flows / Q_s (gap plane) | 0.99973 | 0.99985 | 0.99967 | 0.99968 | 0.99990 |

Numerical uncertainty of F (ASME V&V 20-2009, two grids, Fs = 3, p = 2, refinement ratio from the element counts,
u_num = GCI/1.15): V1 r = 1.672, ΔF = 1.441 × 10⁻³, **u_num = 2.09 × 10⁻³**; V3 r = 1.696, ΔF = 2.475 × 10⁻³,
**u_num = 3.44 × 10⁻³** (3.58 × 10⁻³ with the stored solution). The presets' meshes are not geometrically similar,
which V&V 20 counts against a two-grid estimate; Fs = 3 is its factor for that case. Both F and Q_in rise with
refinement (Q_in by 1.0% and 1.3%).

## single_gap/

`campaign/` is the June 2026 single-gap campaign as the author archived it (Runs A, B and C; Run B at five fractions,
A and C at three), with lock files and Python byte-code removed: each fraction's model, report and two plane exports
(upstream plane x = 75 µm; mid-depth plane), the READMEs, reference outputs and the campaign's own post-processing
scripts. Run B uses a uniform inlet velocity; Runs A and C use a fully developed inlet.

`runB_5e-4_passage/` is the canonical Run B (fraction 5 × 10⁻⁴, uniform inlet; the stored solution of 13 June 2026,
56,526 elements, 292,917 DOF) re-opened and re-evaluated: its four tables at full precision (the achieved fraction
Q_sout/Q_in = 4.962120761577921 × 10⁻⁴), the gap passage at mid-depth on a fine grid (x 96 to 124 µm, y 46 to 58 µm,
every 0.05 µm, z = 70 µm; u, v, w; NaN inside the posts), speed and mesh snapshots, and its report.

`runB_5e-4_developed/` is the same model with the inlet set to Fully developed flow (average velocity U_in, each
disjoint selection separately, as in Runs A and C), re-solved on 9 November 2026 on the same mesh (56,526 elements,
292,918 DOF, four Newton iterations to 1.7 × 10⁻¹⁵, 86 s, 9.3 GB): four tables (achieved fraction
4.8751210307083063 × 10⁻⁴, the value of the developed-inlet run that the campaign README describes), the campaign's
two planes on the campaign's grids, and its report. At equal fraction the developed inlet raises the skim width at
x = 75 µm by 0.30% (0.27% with the campaign's resampling).
