#!/bin/bash
# Every solve of Python CFD, in priority order, each followed by its post-processing.
# usage: bash run_all.sh COMSOL_ROOT   (COMSOL_ROOT holds single_gap/campaign/Run {A,B,C} Sweep and single_gap/runB_5e-4_developed)
# The fractions are COMSOL's achieved Q_sout/Q_in, so that the two codes are compared at equal fraction:
#   Run B uniform 4.962120761577921e-4 and developed 4.8751210307083063e-4 (Table 4 of the re-evaluated and re-solved
#   models, full precision); Runs A and C 4.8518e-4 and 4.8550e-4 (their reports, five digits).
set -u
R="$1"; HERE="$(cd "$(dirname "$0")" && pwd)"; SRC="$HERE/src"; RUNS="$HERE/runs"; mkdir -p "$RUNS"
CAMP="$R/single_gap/campaign"
B_UP="$CAMP/Run B Sweep/5e-4/upstream_grid_frac5e-4 v3.txt"; B_MID="$CAMP/Run B Sweep/5e-4/middepth_grid_frac5e-4.txt"
BD_UP="$R/single_gap/runB_5e-4_developed/upstream_grid_frac5e-4 v3.txt"; BD_MID="$R/single_gap/runB_5e-4_developed/middepth_grid_frac5e-4.txt"
A_UP="$CAMP/Run A Sweep/5e-4/upstream_official_5e-4.txt"; A_MID="$CAMP/Run A Sweep/5e-4/middepth_official_5e-4.txt"
C_UP="$CAMP/Run C Sweep/5e-4/upstream_official_5e-4.txt"; C_MID="$CAMP/Run C Sweep/5e-4/middepth_official_5e-4.txt"
FB=4.962120761577921e-4; FBD=4.8751210307083063e-4; FA=4.8518e-4; FC=4.8550e-4

gap() {   # name run inlet frac scale up mid
  local d="$RUNS/$1"
  if [ -f "$d/metrics.json" ]; then echo "skip $1"; return; fi
  echo "=== $1 $(date +%T)"
  python3 -I "$SRC/stokes_gap.py" --run "$2" --inlet "$3" --frac "$4" --scale "$5" --out "$d" > "$d.solve.log" 2>&1 || { echo "FAILED solve $1"; return; }
  python3 -I "$SRC/postproc_gap.py" "$d" --comsol_upstream "$6" --comsol_middepth "$7" --comsol_fq "$4" > "$d.post.log" 2>&1 || { echo "FAILED post $1"; return; }
  rm -f "$d/solution.npz"; tail -3 "$d.post.log"
}
five() {  # name variant scale
  local d="$RUNS/$1"
  if [ -f "$d/metrics.json" ]; then echo "skip $1"; return; fi
  echo "=== $1 $(date +%T)"
  python3 -I "$SRC/stokes_5gap.py" --variant "$2" --scale "$3" --out "$d" > "$d.solve.log" 2>&1 || { echo "FAILED solve $1"; return; }
  python3 -I "$SRC/postproc_5gap.py" "$d" > "$d.post.log" 2>&1 || { echo "FAILED post $1"; return; }
  rm -f "$d/solution.npz"; grep -o '"F": [0-9.]*' "$d.post.log" | head -1
}

gap A_dev_s1.414 A developed $FA 1.414 "$A_UP" "$A_MID"
gap C_dev_s1.414 C developed $FC 1.414 "$C_UP" "$C_MID"
gap B_dev_s1.414 B developed $FBD 1.414 "$BD_UP" "$BD_MID"
gap B_uni_s1.414 B uniform $FB 1.414 "$B_UP" "$B_MID"
gap B_uni_s2 B uniform $FB 2 "$B_UP" "$B_MID"
gap A_dev_s2 A developed $FA 2 "$A_UP" "$A_MID"
gap C_dev_s2 C developed $FC 2 "$C_UP" "$C_MID"
five V1_s3 V1 3
five V3_s3 V3 3
five V1_s2 V1 2
five V3_s2 V3 2
five V1_s1.7 V1 1.7
five V3_s1.7 V3 1.7
gap B_uni_s1 B uniform $FB 1 "$B_UP" "$B_MID"
gap B_dev_s2 B developed $FBD 2 "$BD_UP" "$BD_MID"
gap B_dev_s1 B developed $FBD 1 "$BD_UP" "$BD_MID"
gap A_dev_s1 A developed $FA 1 "$A_UP" "$A_MID"
gap C_dev_s1 C developed $FC 1 "$C_UP" "$C_MID"
echo "ALL DONE $(date +%T)"
