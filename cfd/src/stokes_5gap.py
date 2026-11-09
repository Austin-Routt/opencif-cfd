#!/usr/bin/env python3
"""Independent Stokes solve of the preprint's five-gap outlet slab (Dinh 2024 V1 and V3 outlet widths), as modeled in
the COMSOL files 'v1 5gap coarse.mph' / 'v3_5gap_coarse.mph' (cif-reference-implementation, COMSOL/):

  center half-channel  x in [0, L_total], y in [0, w_c/2], z in [0, depth]        (symmetry at y = 0)
  side channel         x in [0, L_total], y in [w_c/2 + gap, w_c/2 + gap + w_s]
  five gap passages    x in [L_post + k pitch, L_post + k pitch + L_gap_x], y in [w_c/2, w_c/2 + gap], k = 0..4
  w_c 52, gap 19 (the wall across the flow), L_post 27, L_gap_x 22 (the opening along the flow), pitch 49,
  L_total 272, depth 140; w_s 44 (V1) or 57 (V3); mu 1.2 mPa s
  inlet (x = 0, center): static pressure p_in = 100 Pa, normal flow (tangential velocity zero)
  outlets (x = L_total, center and side): static pressure 0, full traction condition
  everything else no slip (including the side channel's upstream end)
The boundary conditions are uniform in z, so the solution is symmetric about mid-depth: this solve models the half
depth z in [0, depth/2] with a symmetry plane at z = depth/2 (COMSOL modeled the full depth).

Discretization: P2-P1 tetrahedra on a graded tensor mesh (exact geometry), the stress form 2 mu eps(u):eps(v) so that the
pressure boundary condition is the one COMSOL applies ((-pI + 2 mu eps) n = -p0 n), MINRES with a block preconditioner
(classical AMG on the vector Laplacian, lumped pressure mass).  Solved with mu = 1 and p_in = 1; velocities scale as
p_in/mu (Stokes flow is linear), so u_phys [um/s] = u * (100 Pa / 1.2e-3 Pa s) and p_phys = p * 100 Pa.
"""
import sys, os, time, json, argparse
import numpy as np
import scipy.sparse as sp
import scipy.sparse.linalg as spla
from skfem import MeshHex, Basis, FacetBasis, ElementTetP2, ElementTetP1, ElementVector, asm, condense, BilinearForm, LinearForm
from skfem.helpers import sym_grad, ddot, div, dot
from skfem.models.poisson import vector_laplace, mass
import pyamg

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from stokes_gap import graded                                                 # noqa: E402

GEOM = {'V1': dict(w_c=52., w_s=44., gap=19., L_post=27., L_gap_x=22., depth=140., n_gaps=5),
        'V3': dict(w_c=52., w_s=57., gap=19., L_post=27., L_gap_x=22., depth=140., n_gaps=5)}


def build_mesh(G, s):
    wc2, gap, ws, h2 = G['w_c'] / 2, G['gap'], G['w_s'], G['depth'] / 2
    pitch = G['L_post'] + G['L_gap_x']; n = G['n_gaps']
    Lx = (n + 1) * G['L_post'] + n * G['L_gap_x']
    edges = [G['L_post'] + k * pitch for k in range(n)]                        # gap upstream edges
    hg, hb, hw, hf, g = G['h_gap'] * s, G['h_bulk'] * s, G['h_wall'] * s, G['h_floor'] * s, G['growth']
    xs = [graded(0, edges[0], None, hg, hb, g)[:-1]]
    for k in range(n):
        x0 = edges[k]; x1 = x0 + G['L_gap_x']
        xs.append(graded(x0, x1, hg, hg, hb, g)[:-1])
        xend = edges[k + 1] if k + 1 < n else Lx
        xs.append(graded(x1, xend, hg, hg if k + 1 < n else None, hb, g)[:-1])
    xs = np.concatenate(xs + [[Lx]])
    ys = np.concatenate([graded(0, wc2, None, hw, hb, g)[:-1],
                         graded(wc2, wc2 + gap, hw, hw, hb, g)[:-1],
                         graded(wc2 + gap, wc2 + gap + ws, hw, hw, hb, g)])
    zs = graded(0, h2, hf, None, hb, g)
    m = MeshHex.init_tensor(xs, ys, zs)
    c = m.p[:, m.t].mean(axis=1)
    ingap = np.zeros(c.shape[1], dtype=bool)
    for x0 in edges:
        ingap |= (c[0] > x0) & (c[0] < x0 + G['L_gap_x'])
    solid = (c[1] > wc2) & (c[1] < wc2 + gap) & ~ingap
    m = m.remove_elements(np.nonzero(solid)[0]).to_meshtet()
    tol = 1e-6
    m = m.with_boundaries({
        'inlet': lambda x: (x[0] < tol) & (x[1] < wc2 + tol),
        'outlet_c': lambda x: (x[0] > Lx - tol) & (x[1] < wc2 + tol),
        'outlet_s': lambda x: (x[0] > Lx - tol) & (x[1] > wc2 + gap - tol),
        'sym_y': lambda x: x[1] < tol,
        'sym_z': lambda x: x[2] > h2 - tol,
    })
    return m, (xs, ys, zs), Lx


@BilinearForm
def stress_form(u, v, w):
    return 2.0 * ddot(sym_grad(u), sym_grad(v))


@BilinearForm
def div_form(u, q, w):
    return div(u) * q


def solve(G, s, tol=1e-9, log=print):
    t0 = time.time()
    m, grids, Lx = build_mesh(G, s)
    log('mesh: %d tets, %d vertices (scale %.3g); grid %d x %d x %d' % (m.t.shape[1], m.p.shape[1], s,
                                                                      len(grids[0]), len(grids[1]), len(grids[2])))
    bu = Basis(m, ElementVector(ElementTetP2()), intorder=4)
    bp = bu.with_element(ElementTetP1())
    nchunk = max(1, int(np.ceil(m.t.shape[1] / 25000)))
    A = L = B = None
    for idx in np.array_split(np.arange(m.t.shape[1]), nchunk):
        bu_c = Basis(m, ElementVector(ElementTetP2()), intorder=4, elements=idx)
        bp_c = bu_c.with_element(ElementTetP1())
        Ac = asm(stress_form, bu_c).tocsr(); Lc = asm(vector_laplace, bu_c).tocsr(); Bc = asm(div_form, bu_c, bp_c).tocsr()
        A = Ac if A is None else A + Ac; L = Lc if L is None else L + Lc; B = Bc if B is None else B + Bc
        del bu_c, bp_c, Ac, Lc, Bc
    nu, npre = A.shape[0], B.shape[0]
    log('assembled: %d velocity dofs, %d pressure dofs, nnz(A) %d (%.0f s)' % (nu, npre, A.nnz, time.time() - t0))
    # pressure boundary load: -p0 int n.v on the inlet (p0 = 1) and 0 on the outlets
    fb_in = FacetBasis(m, ElementVector(ElementTetP2()), facets=m.boundaries['inlet'], intorder=4)

    @LinearForm
    def pload(v, w):
        return -1.0 * dot(w.n, v)
    f = asm(pload, fb_in)
    f = np.concatenate([f, np.zeros(npre)])
    # Dirichlet: no slip on walls; symmetry; inlet tangential (normal flow)
    wc2, gap, ws, h2 = G['w_c'] / 2, G['gap'], G['w_s'], G['depth'] / 2
    x = np.zeros(nu + npre)
    D = set()
    named = np.concatenate([m.boundaries[k] for k in ('inlet', 'outlet_c', 'outlet_s', 'sym_y', 'sym_z')])
    wall_facets = np.setdiff1d(m.boundary_facets(), named)
    D.update(bu.get_dofs(facets=wall_facets).all().tolist())
    D.update(bu.get_dofs(facets=m.boundaries['sym_y']).all(['u^2']).tolist())
    D.update(bu.get_dofs(facets=m.boundaries['sym_z']).all(['u^3']).tolist())
    D.update(bu.get_dofs(facets=m.boundaries['inlet']).all(['u^2', 'u^3']).tolist())
    D = np.array(sorted(D))
    K = sp.bmat([[A, -B.T], [-B, None]], format='csr')              # symmetric indefinite
    f[nu:] *= -1.0
    Kc, fc, xc, I = condense(K, f, x=x, D=D)
    Iu = I[I < nu]; nfu = len(Iu); Ip = I[I >= nu]
    Lc_free = L[Iu][:, Iu].tocsr()
    ml = pyamg.ruge_stuben_solver(Lc_free, strength=('classical', {'theta': 0.5}), max_coarse=500)
    Pamg = ml.aspreconditioner(cycle='V')
    Mp = asm(mass, bp).diagonal()[Ip - nu]
    def prec(v):
        out = np.empty_like(v); out[:nfu] = Pamg.matvec(v[:nfu]); out[nfu:] = v[nfu:] / Mp; return out
    Mop = spla.LinearOperator(Kc.shape, matvec=prec, dtype=float)
    log('condensed %d unknowns; AMG levels %d, complexity %.2f (%.0f s)' % (Kc.shape[0], len(ml.levels), ml.operator_complexity(), time.time() - t0))
    it = [0]
    def cb(xk):
        it[0] += 1
        if it[0] % 100 == 0:
            log('  minres it %d, residual %.2e' % (it[0], np.linalg.norm(Kc @ xk - fc) / np.linalg.norm(fc)))
    sol, info = spla.minres(Kc, fc, M=Mop, rtol=tol, maxiter=6000, callback=cb)
    res = np.linalg.norm(Kc @ sol - fc) / np.linalg.norm(fc)
    log('minres info %d, %d iterations, relative residual %.2e (%.0f s)' % (info, it[0], res, time.time() - t0))
    x[I] = sol
    return dict(mesh=m, u=x[:nu], p=x[nu:], grids=grids, Lx=Lx, n_tets=m.t.shape[1], n_vert=m.p.shape[1],
                n_dofs=Kc.shape[0], minres_it=it[0], residual=res, seconds=time.time() - t0)


if __name__ == '__main__':
    ap = argparse.ArgumentParser()
    ap.add_argument('--variant', default='V1', choices=sorted(GEOM))
    ap.add_argument('--scale', type=float, default=1.0)
    ap.add_argument('--tol', type=float, default=1e-9, help='MINRES relative residual')
    ap.add_argument('--out', default='runs/v1')
    a = ap.parse_args()
    G = dict(GEOM[a.variant], h_gap=1.0, h_bulk=5.0, h_wall=1.0, h_floor=1.0, growth=1.3)
    os.makedirs(a.out, exist_ok=True)
    r = solve(G, a.scale, tol=a.tol)
    meta = {k: r[k] for k in ('n_tets', 'n_vert', 'n_dofs', 'minres_it', 'residual', 'seconds', 'Lx')}
    from stokes_gap import versions
    meta.update(variant=a.variant, scale=a.scale, tol=a.tol, G=G, mu_Pa_s=1.2e-3, p_in_Pa=100.0, versions=versions())
    np.savez_compressed(os.path.join(a.out, 'solution.npz'), p=r['mesh'].p, t=r['mesh'].t, u=r['u'], pr=r['p'],
                        xs=r['grids'][0], ys=r['grids'][1], zs=r['grids'][2])
    json.dump(meta, open(os.path.join(a.out, 'meta.json'), 'w'), indent=1)
    print(json.dumps(meta))
