"""Readers for COMSOL's text exports (tables and plane data), so that Python CFD compares itself with the COMSOL
files as exported, not with values typed into a script.

  read_table(path)  -> list of floats (the one data row of a Table export), and the header lines
  read_plane(path)  -> dict(axes=(xs, ys, zs) or None, x, y, z, u, v, w[, p]) for both layouts COMSOL writes:
                       (a) 'Spreadsheet': one row per point, columns x y z expr...;
                       (b) 'Grid': three coordinate lines, then each expression's values; the singleton axis is the
                           plane's fixed coordinate and the first of the other two axes varies fastest.
The 'Grid' layout is checked when read: the no-slip walls of the exported plane must come out at zero velocity
(see check_grid_order in the single-gap post-processor)."""
import numpy as np


def _data_lines(path):
    with open(path, encoding='utf-8', errors='replace') as fh:
        return [l for l in fh if l.strip() and not l.lstrip().startswith('%')]


def header(path):
    with open(path, encoding='utf-8', errors='replace') as fh:
        return [l.rstrip('\n') for l in fh if l.startswith('%')]


def read_table(path):
    rows = _data_lines(path)
    if len(rows) != 1:
        raise ValueError('%s: expected one data row, found %d (was Evaluate All used instead of Clear and Evaluate '
                         'All?)' % (path, len(rows)))
    return [float(t) for t in rows[0].split()]


def read_plane(path):
    rows = _data_lines(path)
    hdr = header(path)
    is_grid = any(l.strip() == '% Grid' for l in hdr)
    nexpr = None
    for l in hdr:
        if l.startswith('% Expressions:'):
            nexpr = int(l.split(':')[1])
    if not is_grid:
        M = np.array([[float(t) for t in l.split()] for l in rows])
        out = dict(axes=None, x=M[:, 0], y=M[:, 1], z=M[:, 2])
        names = ['u', 'v', 'w', 'p'][:M.shape[1] - 3]
        for k, nm in enumerate(names):
            out[nm] = M[:, 3 + k]
        return out
    a = np.array([float(t) for t in rows[0].split()])
    b = np.array([float(t) for t in rows[1].split()])
    c = np.array([float(t) for t in rows[2].split()])
    data = np.array([float(t) for l in rows[3:] for t in l.split()])
    axes = [a, b, c]
    n = [len(a), len(b), len(c)]
    sing = int(np.argmin(n))
    others = [i for i in range(3) if i != sing]
    A1, A2 = axes[others[0]], axes[others[1]]
    if nexpr is None:
        nexpr = data.size // (len(A1) * len(A2))
    if data.size != nexpr * len(A1) * len(A2):
        raise ValueError('%s: %d values for %d expressions on %d x %d points' % (path, data.size, nexpr, len(A1),
                                                                                 len(A2)))
    blk = data.reshape(nexpr, len(A2), len(A1))          # first other axis fastest
    G2, G1 = np.meshgrid(A2, A1, indexing='ij')
    coords = {others[0]: G1.ravel(), others[1]: G2.ravel(), sing: np.full(G1.size, axes[sing][0])}
    out = dict(axes=(a, b, c), x=coords[0], y=coords[1], z=coords[2])
    for k, nm in enumerate(['u', 'v', 'w', 'p'][:nexpr]):
        out[nm] = blk[k].ravel()
    return out


def as_grid(P, i_fast, i_slow):
    """Reshape a plane (dict from read_plane, scattered or grid) onto its tensor grid: returns (A_fast, A_slow,
    fields as arrays [n_fast, n_slow])."""
    keys = ['x', 'y', 'z']
    cf, cs = P[keys[i_fast]], P[keys[i_slow]]
    af, as_ = np.unique(np.round(cf, 9)), np.unique(np.round(cs, 9))
    jf = np.searchsorted(af, np.round(cf, 9)); js = np.searchsorted(as_, np.round(cs, 9))
    out = {}
    for nm in ('u', 'v', 'w', 'p'):
        if nm in P:
            G = np.full((len(af), len(as_)), np.nan)
            G[jf, js] = P[nm]
            out[nm] = G
    return af, as_, out
