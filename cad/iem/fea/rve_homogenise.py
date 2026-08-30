"""JOB 1a -- direct 3D FEA of the real graded gyroid, one wall thickness at a time.

The wing lattice cannot be volume-meshed whole: a 12 x 14 x 7 mm blade at the
0.2 mm wall / 1.2 mm cell the generator actually emits needs ~0.05 mm elements,
which is O(10^7) tets -- far past what a Python solver will carry.  So the
lattice is characterised here on a multi-cell coupon (direct micro-FE on the
exact ``generate.gyroid`` field, no homogenisation theory in the geometry), and
the resulting effective modulus is carried into the macro wing model.

Test: uniaxial compression of an N-cell coupon.  Bottom face rollered
(u_z = 0, tangentially free), top face given u_z = -delta (tangentially free),
lateral faces traction-free.  That is a free-lateral-expansion uniaxial-stress
test, so E_eff = (F/A_gross) / (delta/L) directly.

Run:  .venv/bin/python fea/rve_homogenise.py
"""

import json
import os
import sys
import time

import numpy as np

sys.stdout.reconfigure(line_buffering=True) if hasattr(sys.stdout, 'reconfigure') else None

import _common as K
import generate as gen


def coupon(wall, cell, ncell=(2, 2, 3), h=None):
    """Effective Young's modulus of the sheet gyroid at a given wall thickness."""
    h = h or wall / 4.0
    Lx, Ly, Lz = (n * cell for n in ncell)
    # snap the grid so an integer number of voxels spans the coupon
    nx, ny, nz = (max(1, int(round(L / h))) for L in (Lx, Ly, Lz))
    hx, hy, hz = Lx / nx, Ly / ny, Lz / nz
    h = min(hx, hy, hz)
    nx, ny, nz = (int(round(L / h)) for L in (Lx, Ly, Lz))

    def fn(X, Y, Z, C):
        return gen.gyroid(X, Y, Z, cell, wall)

    m, cent = K.voxel_hex(fn, ((0, nx * h), (0, ny * h), (0, nz * h)), h)
    nel = m.t.shape[1]
    rel_rho = nel * h ** 3 / (Lx * Ly * Lz)

    basis = K.make_basis(m)
    Kmat = K.stiffness(basis, K.E_TI)

    p = m.p
    tol = 0.25 * h
    bot = np.where(p[2] < tol)[0]
    top = np.where(p[2] > nz * h - tol)[0]
    if bot.size == 0 or top.size == 0:
        raise RuntimeError("coupon lost a platen face")

    delta = 1.0e-3 * Lz          # 0.1% nominal strain, linear so the value is arbitrary
    dofs = lambda nodes, c: 3 * nodes + c          # noqa: E731  (x,y,z interleaved)
    D = np.concatenate([dofs(bot, 2), dofs(top, 2)])
    xp = np.zeros(Kmat.shape[0])
    xp[dofs(top, 2)] = -delta
    # kill the two remaining rigid-body modes (x, y translation + z rotation)
    anchor = bot[np.lexsort((p[1, bot], p[0, bot]))]
    D = np.concatenate([D, dofs(anchor[:1], 0), dofs(anchor[:1], 1),
                        dofs(anchor[-1:], 0)])

    u, r = K.solve_fixed(Kmat, basis, D, xp, verbose=True)
    F = -float(r[dofs(top, 2)].sum())              # compressive reaction, N
    A = Lx * Ly
    E_eff = (F / A) / (delta / Lz)
    vm = K.von_mises(basis, u, K.E_TI)
    # stress concentration factor: peak wall stress per unit *apparent* stress
    kt = float(np.percentile(vm, 99.5)) / (F / A)
    return dict(wall=wall, cell=cell, h=h, nelem=nel, ndof=int(Kmat.shape[0]),
                rel_rho=rel_rho, E_eff=E_eff, E_ratio=E_eff / K.E_TI, kt=kt,
                ncell=list(ncell))


def main():
    g, P = K.geom()
    cell = P["gyroid_cell"]
    walls = [P["wall_face"], 0.25, 0.30, 0.35, P["wall_root"]]

    print(f"Ti-6Al-4V  E = {K.E_TI/1e3:.0f} GPa  nu = {K.NU_TI}")
    print(f"sheet gyroid, cell = {cell} mm, coupon 2 x 2 x 3 cells, "
          f"free lateral faces\n")
    print(f"{'wall':>5} {'h':>6} {'elems':>8} {'dof':>9} {'rho*/rho':>9} "
          f"{'E_eff GPa':>10} {'E*/Es':>8} {'Kt':>6}")
    rows = []
    for w in walls:
        t0 = time.time()
        r = coupon(w, cell)
        r['secs'] = time.time() - t0
        rows.append(r)
        print(f"{r['wall']:5.2f} {r['h']:6.3f} {r['nelem']:8d} {r['ndof']:9d} "
              f"{r['rel_rho']:9.4f} {r['E_eff']/1e3:10.2f} {r['E_ratio']:8.4f} "
              f"{r['kt']:6.2f}")

    # ---- mesh convergence on the governing (thinnest) wall -----------------
    print("\nmesh convergence, wall = %.2f mm:" % walls[0])
    conv = []
    for div in (2, 3, 4, 5):
        r = coupon(walls[0], cell, ncell=(1, 1, 2), h=walls[0] / div)
        conv.append(r)
        print(f"  h = wall/{div} = {r['h']:.4f} mm  {r['nelem']:7d} elems  "
              f"rho* {r['rel_rho']:.4f}  E_eff {r['E_eff']/1e3:6.2f} GPa")

    # ---- coupon-size convergence ------------------------------------------
    print("\ncoupon-size convergence, wall = %.2f mm, h = wall/4:" % walls[0])
    size = []
    for nc in [(1, 1, 2), (2, 2, 3), (3, 3, 3)]:
        r = coupon(walls[0], cell, ncell=nc)
        size.append(r)
        print(f"  {nc[0]}x{nc[1]}x{nc[2]} cells  {r['nelem']:7d} elems  "
              f"rho* {r['rel_rho']:.4f}  E_eff {r['E_eff']/1e3:6.2f} GPa")

    out = os.path.join(os.path.dirname(os.path.abspath(__file__)), "rve_results.json")
    with open(out, "w") as fh:
        json.dump(dict(walls=rows, mesh_convergence=conv, size_convergence=size),
                  fh, indent=2)
    print(f"\nwrote {out}")


if __name__ == "__main__":
    main()
