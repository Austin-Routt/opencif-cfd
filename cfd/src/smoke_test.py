#!/usr/bin/env python3
"""A two-minute check that an installation reproduces Python CFD's results on coarse meshes.

  1. Run B, developed inlet, mesh scale 4: the skim width on the plane x = 75 um must lie within 3% of the exact
     rectangular-duct law (the converged runs lie within 0.1%), the flow balance must close to 1e-4 and MINRES must reach
     a true relative residual below 1e-6 (the runs stop at 5e-8 to 3e-7; a residual of 7e-11 changes the skim width by
     less than 2e-7, relative).
  2. The V1 five-gap slab, mesh scale 4: F must lie within 0.003 of 0.5079 (the converged runs' extrapolation) and
     the flow balance must close to 1e-4.
Exit status 0 when both pass.  usage: python3 smoke_test.py
"""
import sys, os, time
import numpy as np
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import stokes_gap as SG                                                       # noqa: E402
import stokes_5gap as S5                                                      # noqa: E402
import skim                                                                   # noqa: E402
from postproc_gap import tensor_finder                                        # noqa: E402
from skfem import Basis, ElementTetP2, ElementVector                          # noqa: E402


def probe_field(mesh, grids, u, pts):
    fnd = tensor_finder(mesh, *grids)
    mesh.element_finder = lambda mapping=None: fnd
    bu = Basis(mesh, ElementVector(ElementTetP2()), intorder=4)
    return (bu.probes(pts) @ u).reshape(3, -1)


def main():
    ok = True
    t0 = time.time()
    G = dict(SG.GEOM['B'], L_entry=100., L_gap_x=20., L_exit=60., L_x=180., frac=5e-4,
             h_wall=0.1, h_floor=0.5, h_gap=1.0, hmax=8.0, growth=1.3)
    r = SG.solve(G, 4.0, 'developed', log=lambda *a: None)
    wc2, h2, eps = 50.0, 70.0, 1e-7
    ys = np.arange(0, wc2 + 1e-9, 0.05); zs = np.arange(0, h2 + 1e-9, 1.0)
    Y, Z = np.meshgrid(ys, zs, indexing='ij')
    V = probe_field(r['mesh'], r['grids'], r['u'],
                    np.vstack([np.full(Y.size, 75.0), np.minimum(Y.ravel(), wc2 - eps), np.clip(Z.ravel(), eps, h2 - eps)]))
    U = V[0].reshape(Y.shape); U[-1] = 0.0; U[:, 0] = 0.0
    d, _ = skim.skim_native(ys, zs, U, wc2, 5e-4)
    law = skim.duct_delta(2.5e-4, 100.0, 140.0)
    rel = d / law - 1
    print('single gap, scale 4: %d tets, MINRES residual %.1e; skim %.4f um, duct law %.4f um (%+.2f%%)  [%s]'
          % (r['n_tets'], r['residual'], d, law, 100 * rel, 'pass' if abs(rel) < 0.03 else 'FAIL'))
    ok &= abs(rel) < 0.03 and r['residual'] < 1e-6
    G5 = dict(S5.GEOM['V1'], h_gap=1.0, h_bulk=5.0, h_wall=1.0, h_floor=1.0, growth=1.3)
    r5 = S5.solve(G5, 4.0, log=lambda *a: None)
    Lx = r5['Lx']
    def flux(x, ylo, yhi):
        yy = np.linspace(ylo + eps, yhi - eps, 161); zz = np.linspace(eps, 70 - eps, 141)
        YY, ZZ = np.meshgrid(yy, zz, indexing='ij')
        V = probe_field(r5['mesh'], r5['grids'], r5['u'], np.vstack([np.full(YY.size, x), YY.ravel(), ZZ.ravel()]))
        return np.trapezoid(np.trapezoid(V[0].reshape(YY.shape), zz, axis=1), yy)
    Q_in, Q_c, Q_s = flux(eps, 0, 26), flux(Lx - eps, 0, 26), flux(Lx - eps, 45, 89)
    F = 1 - Q_c / Q_in; imb = (Q_in - Q_c - Q_s) / Q_in
    print('five gaps V1, scale 4: %d tets, MINRES residual %.1e; F %.4f (expected 0.5079 +- 0.003), balance %.1e  [%s]'
          % (r5['n_tets'], r5['residual'], F, imb, 'pass' if abs(F - 0.5079) < 0.003 and abs(imb) < 1e-4 else 'FAIL'))
    ok &= abs(F - 0.5079) < 0.003 and abs(imb) < 1e-4 and r5['residual'] < 1e-6
    print('%s (%.0f s)' % ('SMOKE TEST PASSED' if ok else 'SMOKE TEST FAILED', time.time() - t0))
    sys.exit(0 if ok else 1)


if __name__ == '__main__':
    main()
