# Independent recheck of the COMSOL numbers

Scripts written from scratch, without the code of `cfd/`, from a brief that gave the file formats and the definitions
but not the expected values. They recompute, from the raw COMSOL exports in `comsol/`:

- `fivegap.py comsol/five_gap` — Q_in, Q_c, Q_s, F, the flow balance, the per-gap pressure differences and their
  half-lives, the two-grid V&V 20 estimate, the flow through each gap and the out-of-plane velocity at mid-depth
  (`fivegap_out.txt`); `fivegap_xcheck.py` checks the probe pressures of Tables 4 and 5 against the mid-depth export;
- `singlegap.py LABEL FILE W_C H F_Q [LABEL FILE W_C H F_Q ...]` — the skim width at x = 75 µm and the exact duct law
  for each upstream export given (`singlegap_out.txt` holds Run B uniform and developed, Run A and Run C, with the
  fractions of `cfd/README.md`); `duct_alt.py` solves the duct law with a second Fourier expansion;
- `table4_check.py FILE ...` — the achieved fractions in the two Run B models' Table 4; `inspect_plane.py` — the
  export layouts.

Every number agrees with `cfd/results/comparison.json` to the digits printed there.
