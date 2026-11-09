"""Independent five-gap checker. Usage: python3 -I fivegap.py <five_gap_dir>"""
import sys, os, glob, math
import numpy as np

ROOT = sys.argv[1]
MODELS = ["V1_coarse", "V1_extra_coarse", "V3_coarse_resolved", "V3_coarse_stored_May7", "V3_extra_coarse"]


def read_table(path):
    rows = []
    with open(path, encoding="utf-8", errors="replace") as f:
        for line in f:
            if line.startswith("%") or not line.strip():
                continue
            rows.append([float(t) for t in line.split()])
    return rows


def find_table(folder, n):
    c = [p for p in glob.glob(os.path.join(folder, "*")) if os.path.basename(p).lower().startswith(f"table{n}")]
    assert len(c) == 1, (folder, n, c)
    return c[0]


def load_points(fn):
    a = np.loadtxt(fn, comments="%", encoding="utf-8")
    return a


def plane_grid(a, axes):
    """Put point data on a regular 2D grid. axes = (i,j) column indices of the varying coords."""
    i, j = axes
    xi = np.unique(a[:, i]); xj = np.unique(a[:, j])
    ii = np.searchsorted(xi, a[:, i]); jj = np.searchsorted(xj, a[:, j])
    nv = a.shape[1] - 3
    G = np.full((nv, len(xi), len(xj)), np.nan)
    cnt = np.zeros((len(xi), len(xj)), int)
    np.add.at(cnt, (ii, jj), 1)
    assert cnt.max() == 1 and cnt.min() == 1, "grid not complete/unique"
    for k in range(nv):
        G[k, ii, jj] = a[:, 3 + k]
    return xi, xj, G


def halflife(dp, n):
    k = np.arange(n); y = np.log(np.asarray(dp[:n]))
    s = np.polyfit(k, y, 1)[0]
    return math.log(2) / (-s), s


out = {}
for m in MODELS:
    d = os.path.join(ROOT, m)
    T = {}
    for n in range(1, 6):
        p = find_table(d, n)
        rows = read_table(p)
        if len(rows) != 1:
            print(f"WARNING {m} table{n}: {len(rows)} data rows")
        T[n] = rows[0]
    Qin, Qc, Qs = T[1][0], T[2][0], T[3][0]
    F = 1 - Qc / Qin
    imb = (Qin - Qc - Qs) / Qin
    dp = np.array(T[4]) - np.array(T[5])
    hl5, s5 = halflife(dp, 5); hl4, s4 = halflife(dp, 4)
    r = dict(Qin=Qin, Qc=Qc, Qs=Qs, F=F, imb=imb, pc=T[4], ps=T[5], dp=dp, hl5=hl5, hl4=hl4, s5=s5, s4=s4)

    # gap plane y = 35.5
    gp = glob.glob(os.path.join(d, "*gapplane_y35p5*.txt"))
    assert len(gp) == 1
    a = load_points(gp[0])
    assert np.allclose(a[:, 1], a[0, 1])
    xs, zs, G = plane_grid(a, (0, 2))  # G[var, ix, iz], var: u v w p
    v = np.nan_to_num(G[1], nan=0.0)
    gaps = []
    nan_inside = []
    for k in range(5):
        x0, x1 = 27 + 49 * k, 49 + 49 * k
        sel = (xs >= x0 - 1e-9) & (xs <= x1 + 1e-9)
        xg = xs[sel]
        assert abs(xg[0] - x0) < 1e-9 and abs(xg[-1] - x1) < 1e-9
        vg = v[sel]
        qz = np.trapezoid(vg, zs, axis=1)
        q = np.trapezoid(qz, xg) * 1e-12  # um^2 -> m^2
        gaps.append(q)
        # NaN audit: interior points (not on gap edges, not on z walls) that are NaN
        raw = G[1][sel]
        interior = raw[1:-1, 1:-1]
        nan_inside.append(int(np.isnan(interior).sum()))
        # edge values
    r["gapQ"] = np.array(gaps)
    r["gap_nan_interior"] = nan_inside
    # wall checks on gap plane: v on z=0, z=140, x-edges
    finite = ~np.isnan(G[1])
    r["gap_wall_max_v"] = {
        "z=0": np.nanmax(np.abs(G[1][:, 0])), "z=140": np.nanmax(np.abs(G[1][:, -1])),
    }
    edge_vals = []
    for k in range(5):
        for xe in (27 + 49 * k, 49 + 49 * k):
            ix = np.argmin(abs(xs - xe))
            edge_vals.append(np.nanmax(np.abs(G[1][ix])) if np.any(~np.isnan(G[1][ix])) else np.nan)
    r["gap_edge_max_v"] = np.nanmax(edge_vals) if np.any(~np.isnan(edge_vals)) else np.nan
    r["gap_edge_nan_cols"] = int(np.sum(np.isnan(edge_vals)))
    r["vmax_gap"] = np.nanmax(np.abs(G[1]))
    # non-NaN outside gaps?
    ingap = np.zeros(len(xs), bool)
    for k in range(5):
        ingap |= (xs >= 27 + 49 * k - 1e-9) & (xs <= 49 + 49 * k + 1e-9)
    r["finite_outside_gaps"] = int(finite[~ingap].sum())
    r["nan_inside_gaps_total"] = int((~finite[ingap]).sum())

    # mid-depth plane
    mp = glob.glob(os.path.join(d, "*middepth_z70*.txt"))
    assert len(mp) == 1
    b = load_points(mp[0])
    assert np.allclose(b[:, 2], 70)
    ok = ~np.isnan(b[:, 3:6]).any(axis=1)
    U = b[ok, 3:6]
    spd = np.sqrt((U ** 2).sum(1))
    r["wmax"] = np.abs(U[:, 2]).max(); r["spdmax"] = spd.max(); r["wratio"] = r["wmax"] / r["spdmax"]
    r["mid_nan"] = int((~ok).sum()); r["mid_n"] = len(b)
    # mid-depth u at x=0 inlet (for info)
    out[m] = r

np.set_printoptions(precision=9)
for m, r in out.items():
    print(f"\n### {m}")
    print(f"Q_in={r['Qin']:.10e} Q_c={r['Qc']:.10e} Q_s={r['Qs']:.10e}")
    print(f"F={r['F']:.10f} (7dp {r['F']:.7f})  imbalance={r['imb']:.6e}")
    print("p_c =", np.array(r['pc'])); print("p_s =", np.array(r['ps']))
    print("dp  =", r['dp'])
    print(f"half-life 5 gaps={r['hl5']:.8f} (slope {r['s5']:.8f}); 4 gaps={r['hl4']:.8f} (slope {r['s4']:.8f})")
    print("gap Q =", r['gapQ'])
    print(f"sum gap Q = {r['gapQ'].sum():.10e}; sum/Q_s = {r['gapQ'].sum()/r['Qs']:.8f}; sum/(Qin-Qc) = {r['gapQ'].sum()/(r['Qin']-r['Qc']):.8f}")
    print("gap Q / Q_in =", r['gapQ'] / r['Qin'])
    print("NaN interior per gap:", r['gap_nan_interior'], " NaN in gap columns total:", r['nan_inside_gaps_total'],
          " finite outside gaps:", r['finite_outside_gaps'])
    print("wall |v| max:", r['gap_wall_max_v'], " edge-column |v| max:", r['gap_edge_max_v'], " edge cols all-NaN:", r['gap_edge_nan_cols'], " max|v| gap:", r['vmax_gap'])
    print(f"middepth: max|w|={r['wmax']:.6e} max speed={r['spdmax']:.8e} ratio={r['wratio']:.6e}  NaN rows {r['mid_nan']}/{r['mid_n']}")

# GCI
N = {"V1_coarse": 115904, "V1_extra_coarse": 24783, "V3_coarse_resolved": 104705, "V3_coarse_stored_May7": 104206, "V3_extra_coarse": 21448}
print("\n### GCI")
gci = {}
for fine, coarse in [("V1_coarse", "V1_extra_coarse"), ("V3_coarse_resolved", "V3_extra_coarse"), ("V3_coarse_stored_May7", "V3_extra_coarse")]:
    Ff, Fc = out[fine]["F"], out[coarse]["F"]
    r = (N[fine] / N[coarse]) ** (1 / 3)
    G = 3 * abs(Ff - Fc) / (r ** 2 - 1)
    un = G / 1.15
    Fext = Ff + (Ff - Fc) / (r ** 2 - 1)
    gci[fine] = un
    print(f"{fine} vs {coarse}: Nf={N[fine]} Nc={N[coarse]} r={r:.8f} Ff={Ff:.10f} Fc={Fc:.10f} dF={Ff-Fc:.8e} GCI={G:.8e} u_num={un:.8e} F_ext={Fext:.10f}")

print("\n### Validation")
for fine, mean, sd, fac in [("V1_coarse", 0.505, 0.033, 0.031), ("V3_coarse_resolved", 0.649, 0.036, 0.015)]:
    F = out[fine]["F"]; uD = sd; uin = fac * uD; un = gci[fine]
    E = F - mean
    uval = math.sqrt(un ** 2 + uin ** 2 + uD ** 2)
    print(f"{fine}: F={F:.10f} E={E:.8e} E/uD={E/uD:.8f} u_num={un:.8e} u_input={uin:.8e} u_val={uval:.10e} E/u_val={E/uval:.8f} (u_val/uD-1)%={100*(uval/uD-1):.8f} uD/u_num={uD/un:.8f}")
