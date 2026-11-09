# Figures recolored to the palette of Chapters 1 and 2 (9 November 2026)

In Chapters 1 and 2 a series takes blue, vermillion or green (Okabe-Ito), and black marks a single curve, an analytic
reference or a frame. Two entity colors of the later chapters broke that rule and are changed in `tools/pubstyle.py`,
which every figure script shares, so each entity keeps one color across the dissertation:

| entity | before | now |
|---|---|---|
| Dinh et al.'s V1 (VCOL) | black | blue #0072B2 |
| CFD run C, α = 0.80 (ACOL) | black | green #009E73 |
| V3, run B, run A | vermillion, blue, vermillion | unchanged |

Line styles and markers are unchanged. All three pairs of blue, vermillion and green pass the color-blind check
(worst ΔE 11.0 deutan, 18.7 normal vision).

- `dissertation_figures/`: the four figures that draw V1 or run C, rebuilt by `tools/figures_kb.py` and
  `tools/figures_cfd_campaign.py` with matplotlib 3.10.9 (the version of the original build): F3-5, F3-7, F4-5, F5-1,
  and the proof sheet. Rebuilt the same way, every other figure is byte-identical to its previous file, and both
  check files are unchanged (24 checks, 0 mismatches). No caption names V1's or run C's color, so no caption changes.
- `cfd_panels/`: the CFD panels of the archive (F4-1c, F4-5d, F5-1c and the two-code five-gap figure), in the same
  colors, labeled "Python CFD". In F4-1c COMSOL's dividing streamline is vermillion, as in F4-1(a, b), and Python
  CFD's is blue, dashed, with a thin white edge so that it shows on the dark end of the color map.
- `tools/`: the two changed scripts of the outline (`pubstyle.py`: VCOL and ACOL; `figures_kb.py`: F3-7 and F5-1 now
  read those maps instead of hard-coded colors).
