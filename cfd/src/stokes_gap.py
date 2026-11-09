#!/usr/bin/env python3
"""Independent single-gap Stokes solve of the CIF quarter cell (open-source: scikit-fem + pyamg + scipy).

Geometry and boundary conditions are read off the COMSOL model files of the author's campaign
(cif_funnel_run{A,B,C}*.mph: smodel.json parameters, dmodel.xml blocks and box selections):
  center half-channel  x in [0, L_x], y in [0, w_c/2], z in [0, h/2]
  gap aperture         x in [L_entry, L_entry + L_gap_x], y in [w_c/2, w_c/2 + gap]
  side channel         x in [0, L_x], y in [w_c/2 + gap, w_c/2 + gap + w_s]
  inlet (x = 0, center): normal inflow velocity U_in = Q0 / A_in (uniform), or the developed duct profile
  center outlet (x = L_x): pressure 0 (natural)
  side outlet (x = L_x, side): normal outflow velocity U_side = frac Q0 / A_side (uniform)
  symmetry at y = 0 (v = 0) and z = h/2 (w = 0); no slip elsewhere (incl. the side channel's upstream end)
Discretization: P2-P1 Taylor-Hood tetrahedra on a graded tensor mesh (every box face is a mesh plane, so the
geometry is exact), MINRES on the symmetric saddle-point system with a block preconditioner (smoothed-aggregation
AMG on the velocity Laplacian, lumped pressure mass matrix).  Lengths in um, velocities in units of U_in, mu = 1
(Stokes flow is linear, so every reported width depends on the fraction alone).
"""
import sys, time, json, argparse
import numpy as np
import scipy.sparse as sp
import scipy.sparse.linalg as spla
from skfem import MeshHex, Basis, ElementTetP2, ElementTetP1, ElementVector, asm, condense
from skfem.models.poisson import vector_laplace
from skfem.models.general import divergence
import pyamg


def graded(a, b, h_a, h_b, hmax, g):
    """Nodes from a to b with first cell h_a at a (None: coarse), h_b at b (None: coarse), cells growing by g up
    to hmax; the interior is filled with cells <= hmax and the whole array rescaled to end exactly at b."""
    L = b - a
    left = []
    if h_a is not None:
        h = h_a
        while sum(left) + h < L / 2 and h < hmax:
            left.append(h); h *= g
    right = []
    if h_b is not None:
        h = h_b
        while sum(right) + h < L / 2 and h < hmax:
            right.append(h); h *= g
    rem = L - sum(left) - sum(right)
    n_mid = max(1, int(np.ceil(rem / hmax)))
    cells = left + [rem / n_mid] * n_mid + right[::-1]
    cells = np.array(cells) * (L / np.sum(cells))
    return a + np.concatenate([[0.0], np.cumsum(cells)])


def build_mesh(G, s):
    """Graded tensor mesh of the quarter cell at scale s (s = 1 is the reference; sizes scale with s)."""
    wc2, gap, ws, h2 = G['w_c'] / 2, G['gap'], G['w_s'], G['depth'] / 2
    x0, x1, Lx = G['L_entry'], G['L_entry'] + G['L_gap_x'], G['L_x']
    hw = G['h_wall'] * s          # first cell at the gap-side wall (COMSOL: boundary layer 0.1 um)
    hmax = G['hmax'] * s          # bulk (COMSOL: 8 um)
    hgap = G['h_gap'] * s         # near the gap edges (COMSOL: 3 um on the wall patch)
    g = G['growth']
    xs = np.concatenate([graded(0, x0, None, hgap, hmax, g)[:-1],
                         graded(x0, x1, hgap, hgap, min(hmax, hgap * 2), g)[:-1],
                         graded(x1, Lx, hgap, None, hmax, g)])
    ys = np.concatenate([graded(0, wc2, None, hw, hmax, g)[:-1],
                         graded(wc2, wc2 + gap, hw * 4, hw * 4, min(hmax, 2.5 * s), g)[:-1],
                         graded(wc2 + gap, wc2 + gap + ws, hw * 4, None, hmax, g)])
    zs = graded(0, h2, G['h_floor'] * s, None, hmax, g)
    m = MeshHex.init_tensor(xs, ys, zs)
    c = m.p[:, m.t].mean(axis=1)                       # element centers
    wall = (c[1] > wc2) & (c[1] < wc2 + gap) & ~((c[0] > x0) & (c[0] < x1))
    m = m.remove_elements(np.nonzero(wall)[0])
    m = m.to_meshtet()
    tol = 1e-6
    m = m.with_boundaries({
        'inlet': lambda x: (x[0] < tol) & (x[1] < wc2 + tol),
        'outlet_c': lambda x: (x[0] > Lx - tol) & (x[1] < wc2 + tol),
        'outlet_s': lambda x: (x[0] > Lx - tol) & (x[1] > wc2 + gap - tol),
        'sym_y': lambda x: x[1] < tol,
        'sym_z': lambda x: x[2] > h2 - tol,
    })
    return m, (xs, ys, zs)


def duct_profile(y, z, a, b, N=60):
    """Developed rectangular-duct velocity (unit mean) on |y| <= a, 0 <= z <= 2b (floor at z = 0)."""
    zz = z - b
    u = np.zeros_like(y, dtype=float)
    for n in range(N):
        k = (2 * n + 1) * np.pi / (2 * b)
        u += ((-1) ** n / (2 * n + 1) ** 3) * (1 - np.cosh(k * y) / np.cosh(k * a)) * np.cos(k * zz)
    # mean over the section
    ny, nz = 201, 201
    yg = np.linspace(-a, a, ny); zg = np.linspace(0, 2 * b, nz)
    YY, ZZ = np.meshgrid(yg, zg, indexing='ij')
    um = np.zeros_like(YY)
    for n in range(N):
        k = (2 * n + 1) * np.pi / (2 * b)
        um += ((-1) ** n / (2 * n + 1) ** 3) * (1 - np.cosh(k * YY) / np.cosh(k * a)) * np.cos(k * (ZZ - b))
    mean = np.trapezoid(np.trapezoid(um, zg, axis=1), yg) / (2 * a * 2 * b)
    return u / mean


def solve(G, s, inlet='uniform', tol=1e-9, log=print):
    t0 = time.time()
    m, grids = build_mesh(G, s)
    log('mesh: %d tets, %d vertices (scale %.3g); %d x %d x %d grid' % (m.t.shape[1], m.p.shape[1], s,
                                                                      len(grids[0]), len(grids[1]), len(grids[2])))
    bu = Basis(m, ElementVector(ElementTetP2()), intorder=4)
    bp = bu.with_element(ElementTetP1())
    # assemble in element chunks (the COO triplets of 184k P2 tets at once exceed this container's memory)
    nchunk = max(1, int(np.ceil(m.t.shape[1] / 30000)))
    A = None; B = None
    for idx in np.array_split(np.arange(m.t.shape[1]), nchunk):
        bu_c = Basis(m, ElementVector(ElementTetP2()), intorder=4, elements=idx)
        bp_c = bu_c.with_element(ElementTetP1())
        Ac = asm(vector_laplace, bu_c).tocsr(); Bc = asm(divergence, bu_c, bp_c).tocsr()
        A = Ac if A is None else (A + Ac); B = Bc if B is None else (B + Bc)
        del bu_c, bp_c, Ac, Bc
    A = A.tocsr(); B = B.tocsr()                    # B u = int q div u
    nu, npre = A.shape[0], B.shape[0]
    log('assembled: %d velocity dofs, %d pressure dofs, nnz(A) %d (%.0f s)' % (nu, npre, A.nnz, time.time() - t0))
    wc2, gap, ws, h2 = G['w_c'] / 2, G['gap'], G['w_s'], G['depth'] / 2
    A_in = wc2 * h2; A_side = ws * h2
    U_side = G['frac'] * A_in / A_side              # U_in = 1
    # Dirichlet data: walls first, then symmetry, then the inlet and side outlet (later features override earlier
    # ones on shared edges, as COMSOL's feature order does)
    x = np.zeros(nu)
    D = set()
    wall_facets = np.setdiff1d(m.boundary_facets(), np.concatenate([m.boundaries[k] for k in
                               ('inlet', 'outlet_c', 'outlet_s', 'sym_y', 'sym_z')]))
    dofs_wall = bu.get_dofs(facets=wall_facets).all()
    D.update(dofs_wall.tolist()); x[dofs_wall] = 0.0
    d = bu.get_dofs(facets=m.boundaries['sym_y']).all(['u^2']); D.update(d.tolist()); x[d] = 0.0
    d = bu.get_dofs(facets=m.boundaries['sym_z']).all(['u^3']); D.update(d.tolist()); x[d] = 0.0
    din = bu.get_dofs(facets=m.boundaries['inlet'])
    # the inlet facets' nodes that lie on walls keep no slip (edges y = wc2 and z = 0); COMSOL's later Inlet feature
    # overrides the wall there, so we do the same except on the symmetry planes, where the normal component stays 0
    d1 = din.all(['u^1']); d23 = din.all(['u^2', 'u^3'])
    D.update(d1.tolist()); D.update(d23.tolist()); x[d23] = 0.0
    if inlet == 'uniform':
        x[d1] = 1.0
    elif inlet == 'developed':
        pts = bu.doflocs[:, d1]
        x[d1] = duct_profile(pts[1], pts[2], wc2, h2)
    else:
        raise ValueError('inlet must be uniform or developed')
    dso = bu.get_dofs(facets=m.boundaries['outlet_s'])
    d1s = dso.all(['u^1']); d23s = dso.all(['u^2', 'u^3'])
    D.update(d1s.tolist()); D.update(d23s.tolist()); x[d1s] = U_side; x[d23s] = 0.0
    # re-impose no slip on wall dofs that the inlet/outlet patches touched (shared edges): COMSOL order is
    # Wall 1 < Inlet 1 < Outlet 2, so the inlet/outlet win; keep that (comment out to let the wall win)
    D = np.array(sorted(D))
    K = sp.bmat([[A, B.T], [B, None]], format='csr')
    xfull = np.concatenate([x, np.zeros(npre)])
    Kc, fc, xc, I = condense(K, np.zeros(nu + npre), x=xfull, D=D)
    log('condensed: %d unknowns (%.0f s)' % (Kc.shape[0], time.time() - t0))
    Iu = I[I < nu]; Ip = I[I >= nu]
    nfu = len(Iu)
    Auu = Kc[:nfu, :nfu].tocsr()
    ml = pyamg.ruge_stuben_solver(Auu, strength=('classical', {'theta': 0.5}), max_coarse=500)
    Pamg = ml.aspreconditioner(cycle='V')
    from skfem.models.poisson import mass
    Mp = asm(mass, bp).diagonal()[Ip - nu]          # lumped pressure mass (mu = 1)
    def prec(v):
        out = np.empty_like(v)
        out[:nfu] = Pamg.matvec(v[:nfu])
        out[nfu:] = v[nfu:] / Mp
        return out
    Mop = spla.LinearOperator(Kc.shape, matvec=prec, dtype=float)
    log('AMG levels %d, operator complexity %.2f (%.0f s)' % (len(ml.levels), ml.operator_complexity(), time.time() - t0))
    it = [0]
    def cb(xk):
        it[0] += 1
        if it[0] % 50 == 0:
            log('  minres it %d, residual %.2e' % (it[0], np.linalg.norm(Kc @ xk - fc) / np.linalg.norm(fc)))
    sol, info = spla.minres(Kc, fc, M=Mop, rtol=tol, maxiter=4000, callback=cb)
    res = np.linalg.norm(Kc @ sol - fc) / np.linalg.norm(fc)
    log('minres info %d, %d iterations, relative residual %.2e (%.0f s)' % (info, it[0], res, time.time() - t0))
    xfull[I] = sol
    u = xfull[:nu]; p = xfull[nu:]
    return dict(mesh=m, bu=bu, bp=bp, u=u, p=p, grids=grids, A_in=A_in, U_side=U_side,
                n_tets=m.t.shape[1], n_vert=m.p.shape[1], n_dofs=Kc.shape[0], minres_it=it[0], residual=res,
                seconds=time.time() - t0)


def probe(bu, u, pts):
    """Velocity (3, n) at the points pts (3, n); points outside the mesh give nan."""
    v = bu.probes(pts) @ u
    return v.reshape(3, -1)


def versions():
    import skfem, scipy, platform
    return dict(python=platform.python_version(), numpy=np.__version__, scipy=scipy.__version__,
                scikit_fem=skfem.__version__, pyamg=pyamg.__version__, machine=platform.machine())


GEOM = {'A': dict(w_c=52., w_s=44., gap=19., depth=140.), 'B': dict(w_c=100., w_s=60., gap=20., depth=140.),
        'C': dict(w_c=100., w_s=60., gap=20., depth=80.)}


if __name__ == '__main__':
    ap = argparse.ArgumentParser(description='Single-gap quarter-cell Stokes solve (Runs A, B, C of the COMSOL campaign)')
    ap.add_argument('--run', default='B', choices=sorted(GEOM))
    ap.add_argument('--frac', type=float, default=5e-4,
                    help='quarter-cell fraction Q_side/Q_in; pass COMSOL\'s achieved value to compare at equal fraction')
    ap.add_argument('--scale', type=float, default=1.0, help='mesh scale s (first wall cell 0.1 s um, bulk 8 s um)')
    ap.add_argument('--inlet', default='uniform', choices=['uniform', 'developed'],
                    help='uniform (COMSOL Run B) or the developed duct profile (COMSOL Runs A and C)')
    ap.add_argument('--tol', type=float, default=1e-9, help='MINRES relative residual')
    ap.add_argument('--out', default='out')
    a = ap.parse_args()
    G = dict(GEOM[a.run], L_entry=100., L_gap_x=20., L_exit=60., L_x=180., frac=a.frac,
             h_wall=0.1, h_floor=0.5, h_gap=1.0, hmax=8.0, growth=1.3)
    import os
    os.makedirs(a.out, exist_ok=True)
    r = solve(G, a.scale, a.inlet, tol=a.tol)
    meta = {k: r[k] for k in ('n_tets', 'n_vert', 'n_dofs', 'minres_it', 'residual', 'seconds')}
    meta.update(run=a.run, frac=a.frac, scale=a.scale, inlet=a.inlet, tol=a.tol, G=G, versions=versions())
    np.savez_compressed(os.path.join(a.out, 'solution.npz'), p=r['mesh'].p, t=r['mesh'].t, u=r['u'], pr=r['p'],
                        xs=r['grids'][0], ys=r['grids'][1], zs=r['grids'][2])
    json.dump(meta, open(os.path.join(a.out, 'meta.json'), 'w'), indent=1)
    print(json.dumps(meta))
