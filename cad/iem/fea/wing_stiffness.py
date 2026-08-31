"""JOB 1b -- compression of the graded-gyroid wing and of the jacket skin.

Model
-----
Body      the free-standing wing span (everything outboard of y = y_root), taken
          from ``generate.wing_envelope_mesh``'s SDF -- the same macro envelope
          the generator prints.
Material  Ti-6Al-4V, graded: every element gets the effective modulus that
          ``rve_homogenise.py`` measured for the gyroid wall thickness the
          generator's own grading law puts at that element's root distance.
          Elements inside ``solid_root`` are fully dense Ti.
Root      the y = y_root cut face is encastre (all 3 dofs).
Platen    rigid frictionless flat plane advancing along -n_hat, where n_hat is
          the in-plane normal to the root->tip chord -- the blade's weak bending
          axis, i.e. the lowest-stiffness transverse direction and therefore the
          most favourable case for producing a soft plateau.
Contact   active-set node-to-plane.  The mesh is pre-rotated so n_hat is +z, so
          each contact constraint is a single Dirichlet dof (frictionless: the
          two tangential dofs stay free).  Nodes whose reaction goes tensile are
          released and the step is re-solved.  Contact area therefore grows with
          the platen stroke -- that nonlinearity IS in the model.
Not in    geometric nonlinearity (see the justification the script prints), and
the model plasticity.  Both are argued in docs/MECH_VALIDATION.md.

The jacket patch (``run_jacket``) is done the other way round -- its mesh resolves
the real gyroid walls and the real skin perforations, so every element there is
SOLID Ti.  Applying the homogenised modulus to an already-porous mesh would knock
the stiffness down twice for the same porosity.

Run:  .venv/bin/python fea/wing_stiffness.py
"""

import json
import math
import os
import sys
import time

import numpy as np

import _common as K
import generate as gen


# ---------------------------------------------------------------------------
# effective modulus of the graded lattice, from the RVE study
# ---------------------------------------------------------------------------

def rve_table():
    p = os.path.join(os.path.dirname(os.path.abspath(__file__)), "rve_results.json")
    if not os.path.exists(p):
        sys.exit("run rve_homogenise.py first (rve_results.json missing)")
    d = json.load(open(p))["walls"]
    w = np.array([r["wall"] for r in d])
    E = np.array([r["E_eff"] for r in d])
    rho = np.array([r["rel_rho"] for r in d])
    o = np.argsort(w)
    return w[o], E[o], rho[o]


def graded_E(cent, P, y_root, wtab, Etab):
    """Element modulus from the generator's own grading law."""
    x, y, z = cent
    dz = np.maximum(0.0, -z)
    dy = np.maximum(0.0, y - y_root)
    root_dist = np.sqrt(dz * dz + dy * dy)
    t = np.clip(root_dist / P["grade_len"], 0.0, 1.0)
    wall = P["wall_root"] + (P["wall_face"] - P["wall_root"]) * t
    E = np.interp(wall, wtab, Etab)
    E[root_dist < P["solid_root"]] = K.E_TI       # solid collar
    return E, wall, root_dist


# ---------------------------------------------------------------------------
# active-set frictionless rigid-platen contact on a linear-elastic body
# ---------------------------------------------------------------------------

def press(Kmat, basis, p_rot, fixed_nodes, delta, s0, surf=None, max_it=10,
          verbose=False):
    """Advance a rigid plane (normal +z in the rotated frame) by `delta`.

    Only *surface* nodes can touch the platen, so the candidate set is the
    mesh boundary above the platen plane.  Returns (F, u, active) with F the
    total normal force on the platen (N).
    """
    nd = basis.nodal_dofs
    if surf is None:
        surf = basis.mesh.boundary_nodes()
    s_plat = s0 - delta
    cand = surf[p_rot[2, surf] > s_plat - 1e-12]
    cand = np.setdiff1d(cand, fixed_nodes)
    active = cand.copy()
    n = Kmat.shape[0]
    for it in range(max_it):
        xp = np.zeros(n)
        D = [nd[c, fixed_nodes] for c in range(3)]
        if active.size:
            D.append(nd[2, active])
            xp[nd[2, active]] = s_plat - p_rot[2, active]
        D = np.concatenate(D)
        u, r = K.solve_fixed(Kmat, basis, D, xp)
        if active.size == 0:
            return 0.0, u, active
        rz = r[nd[2, active]]
        # The platen advances in -z, so at a genuinely-contacting node the platen
        # must PUSH the body down: the reaction r = (K u) at that dof is negative.
        # A node with r_z > 0 is being pulled onto the platen -- release it.
        bad = active[rz > 0]
        if bad.size == 0:
            return float(-rz.sum()), u, active
        if verbose:
            print(f"      contact it{it}: releasing {bad.size} of {active.size}")
        active = np.setdiff1d(active, bad)
    return float(-r[nd[2, active]].sum()), u, active


def rot_to_z(nhat):
    """Rotation matrix taking unit vector nhat to +z."""
    nhat = np.asarray(nhat, float)
    nhat = nhat / np.linalg.norm(nhat)
    z = np.array([0.0, 0.0, 1.0])
    v = np.cross(nhat, z)
    c = float(np.dot(nhat, z))
    if np.linalg.norm(v) < 1e-12:
        return np.eye(3) if c > 0 else -np.eye(3)
    vx = np.array([[0, -v[2], v[1]], [v[2], 0, -v[0]], [-v[1], v[0], 0]])
    return np.eye(3) + vx + vx @ vx / (1 + c)


# ---------------------------------------------------------------------------
# the wing
# ---------------------------------------------------------------------------

def wing_model(h=0.15):
    g, P = K.geom()
    pts = gen._bezier_pts(g.wing_p0, g.wing_p1, g.wing_p2)

    def fn(X, Y, Z, C):
        d2 = gen._bezier_dist2d(C, pts) - 0.5 * P["wing_thick"]
        rate = math.tan(math.radians(P["wing_taper_deg"]))
        shrink = rate * np.maximum(0.0, (-Z) - (P["wing_width"] - P["wing_taper"]))
        wing = np.maximum(d2 + shrink,
                          np.maximum(Z - P["wing_z_top"], -Z - P["wing_width"]))
        return gen.I(wing, g.y_root - Y)      # free span only, root plane at y_root

    ymax = max(p[1] for p in pts) + 1.0
    xmin = min(p[0] for p in pts) - 3.5
    xmax = max(p[0] for p in pts) + 3.5
    b = ((xmin, xmax), (g.y_root, ymax), (-P["wing_width"] - 1.0, 0.5))
    # snap the root plane onto a node plane so the encastre face is exact
    b = ((xmin, xmax), (g.y_root, g.y_root + h * math.ceil((ymax - g.y_root) / h)),
         (-P["wing_width"] - 1.0, 0.5))
    m, cent = K.voxel_hex(fn, b, h)
    return g, P, m, cent, pts


def run_wing(h=0.15, steps=None, verbose=True):
    wtab, Etab, rhotab = rve_table()
    g, P, m, cent, pts = wing_model(h)
    nel = m.t.shape[1]
    vol = nel * h ** 3
    E, wall, rd = graded_E(cent, P, g.y_root, wtab, Etab)

    # ---- platen direction: in-plane normal to the root->tip chord ----------
    root_xy = np.array([g.core_cx + P["wing_root_dx"], g.y_root])
    tip_xy = np.array(g.wing_p2)
    chord = tip_xy - root_xy
    nhat = np.array([chord[1], -chord[0], 0.0])
    nhat /= np.linalg.norm(nhat)
    if nhat[0] < 0:
        nhat = -nhat                              # point outboard (+x side)

    fixed = np.where(m.p[1] < g.y_root + 0.25 * h)[0]
    # Rotate the MESH, not just the coordinates used to pick contact nodes.
    # Assembling K on the unrotated mesh while constraining dof 2 presses along
    # GLOBAL z, whatever n_hat says -- the bug this line fixes.
    R = rot_to_z(nhat)
    m = K.rotate_mesh(m, R)
    basis = K.make_basis(m)
    Kmat = K.stiffness(basis, E)
    p_rot = m.p
    s0 = float(p_rot[2].max())

    if verbose:
        print(f"  mesh h={h:.3f} mm  {nel} hexes  {Kmat.shape[0]} dof  "
              f"vol {vol:.1f} mm^3  root face {fixed.size} nodes")
        print(f"  platen normal ({nhat[0]:+.3f}, {nhat[1]:+.3f}, {nhat[2]:+.3f})  "
              f"free span {np.linalg.norm(chord):.2f} mm")
        print(f"  graded wall {wall.min():.3f}-{wall.max():.3f} mm -> "
              f"E {E.min()/1e3:.1f}-{E.max()/1e3:.1f} GPa  "
              f"(solid collar: {(rd < P['solid_root']).sum()} elems)")

    steps = steps if steps is not None else [0.25 * i for i in range(1, 11)]
    surf = m.boundary_nodes()
    rows = []
    for d in steps:
        t0 = time.time()
        F, u, act = press(Kmat, basis, p_rot, fixed, d, s0, surf=surf)
        vm = K.von_mises(basis, u, E)
        rows.append(dict(delta=d, F=F, k=F / d, n_contact=int(act.size),
                         vm_max=float(vm.max()), vm_p999=float(np.percentile(vm, 99.9)),
                         secs=time.time() - t0))
        if verbose:
            print(f"    d={d*1e3:8.1f} um  F={F:12.4g} N  k={F/d:11.4g} N/mm  "
                  f"contact nodes {act.size:5d}  vM max {vm.max():10.4g} MPa  "
                  f"({rows[-1]['secs']:.0f}s)")
    return dict(h=h, nelem=nel, ndof=int(Kmat.shape[0]), vol=vol,
                nhat=nhat.tolist(), span=float(np.linalg.norm(chord)),
                E_min=float(E.min()), E_max=float(E.max()), steps=rows)


# ---------------------------------------------------------------------------
# the jacket skin, through-thickness compression
# ---------------------------------------------------------------------------

def run_jacket(h=0.06, half=1.5, verbose=True):
    """Through-thickness compression of a real patch of the jacket wall.

    A square patch under the core's bottom pole -- the flattest part of the
    shell, and clear of every magnet pocket and locating pin -- taken from the
    full ``part_jacket_wing`` SDF, so it carries the real perforated 0.6 mm skin
    plus the graded gyroid behind it.

    The shell is a curved offset of the core, so face selection is done on the
    core's own signed distance (``core_outer``), not on a z-plane:

      inner face  env ~ clearance                 -> encastre (the core behind it
                                                     is 1.2 mm solid Ti, ~70x
                                                     stiffer than this shell)
      outer face  env ~ clearance + jacket_thick  -> uniform u_z = -delta,
                                                     tangentially free
                                                     (frictionless flat platen)
    """
    wtab, Etab, rhotab = rve_table()
    g, P = K.geom()
    fn_full, _ = gen.part_jacket_wing(g)
    cx = g.core_cx
    z_in = -(P["core_rz"] + P["clearance"])
    z_out = z_in - P["jacket_thick"]
    b = ((cx - half, cx + half), (-half, half), (z_out - 0.5, z_in + 0.5))

    m, cent = K.voxel_hex(fn_full, b, h)
    nel = m.t.shape[1]

    # This mesh resolves the ACTUAL lattice walls and the actual perforated skin,
    # so every element is solid Ti.  Using the homogenised E from rve_results here
    # would knock the modulus down for a porosity the mesh already represents.
    E = np.full(nel, K.E_TI)

    # distance from the core surface, at nodes (the shell is a curved offset of
    # the core, so faces cannot be picked off a z-plane)
    def env_at(pts):
        C = gen.Ctx(pts[0], pts[1], pts[2])
        return g.core_outer(pts[0], pts[1], pts[2], C)

    d_n = env_at(m.p) - P["clearance"]              # 0 at the inner face

    basis = K.make_basis(m)
    Kmat = K.stiffness(basis, E)
    nd = basis.nodal_dofs

    inner = np.where(d_n < 0.75 * h)[0]
    outer = np.where(d_n > P["jacket_thick"] - 0.75 * h)[0]
    outer = np.setdiff1d(outer, inner)
    if inner.size == 0 or outer.size == 0:
        raise RuntimeError("jacket patch lost a face")

    delta = 1.0e-4                                # 0.1 um, deep in the linear range
    n = Kmat.shape[0]
    xp = np.zeros(n)
    xp[nd[2, outer]] = -delta
    D = np.concatenate([nd[0, inner], nd[1, inner], nd[2, inner], nd[2, outer]])
    u, r = K.solve_fixed(Kmat, basis, D, xp)
    F = float(-r[nd[2, outer]].sum())
    k = F / delta
    A = (2 * half) ** 2
    fill = nel * h ** 3 / (A * P["jacket_thick"])
    if verbose:
        print(f"  patch {2*half:.1f} x {2*half:.1f} mm, h={h:.3f} mm, {nel} hexes, "
              f"{n} dof   ({inner.size} inner / {outer.size} outer nodes)")
        print(f"  solid volume {nel*h**3:.2f} mm^3 of a {A*P['jacket_thick']:.2f} mm^3 "
              f"envelope -> {fill:.3f} fill")
        print(f"  through-thickness k = {k:.4g} N/mm over {A:.1f} mm^2 "
              f"-> {k/A:.4g} N/mm per mm^2")
    return dict(h=h, half=half, nelem=nel, ndof=int(n), area=A, fill=fill,
                k=k, k_per_area=k / A, F=F, delta=delta)


# ---------------------------------------------------------------------------

def gibson_ashby(rho):
    """Sheet-lattice collapse estimates, Gibson & Ashby scaling."""
    return dict(rel_rho=rho,
                sigma_plastic=0.3 * K.SY_TI * rho ** 1.5,   # plastic collapse plateau
                sigma_elastic=0.05 * K.E_TI * rho ** 3)     # elastic cell-wall buckling


def main():
    sys.stdout.reconfigure(line_buffering=True)
    out = {}

    print("=" * 78)
    print("WING -- graded gyroid, transverse compression on the weak bending axis")
    print("=" * 78)
    fine = [1e-5, 3e-5, 1e-4, 3e-4, 1e-3]
    coarse = [0.25 * i for i in range(1, 11)]
    out["wing"] = run_wing(h=0.18, steps=fine + coarse)

    print("\nmesh convergence (single 0.25 mm step):")
    conv = []
    for h in (0.30, 0.24, 0.18, 0.15):
        r = run_wing(h=h, steps=[0.25], verbose=False)
        conv.append(r)
        st = r["steps"][0]
        print(f"  h={h:.2f} mm  {r['nelem']:7d} hexes  vol {r['vol']:7.2f} mm^3  "
              f"k = {st['k']:.4g} N/mm")
    out["wing_convergence"] = conv

    print("\nmesh convergence in the ELASTIC regime (0.1 um step) -- this is the "
          "number\nthe verdict rests on:")
    econv = []
    for h in (0.24, 0.18, 0.15):
        r = run_wing(h=h, steps=[1e-4], verbose=False)
        econv.append(r)
        st = r["steps"][0]
        print(f"  h={h:.2f} mm  {r['nelem']:7d} hexes  k = {st['k']:.0f} N/mm  "
              f"({st['n_contact']} contact nodes)")
    out["wing_elastic_convergence"] = econv

    print("\n" + "=" * 78)
    print("JACKET SKIN -- through-thickness compression onto the core")
    print("=" * 78)
    out["jacket"] = [run_jacket(h=hh) for hh in (0.09, 0.07, 0.06)]

    # ---- lattice collapse estimates --------------------------------------
    wtab, Etab, rhotab = rve_table()
    print("\n" + "=" * 78)
    print("IS THERE A BUCKLING / SOFTENING PLATEAU?  (Gibson-Ashby, sheet lattice)")
    print("=" * 78)
    ga = []
    for w, rho in zip(wtab, rhotab):
        d = gibson_ashby(rho)
        d["wall"] = w
        ga.append(d)
        gov = "plastic" if d["sigma_plastic"] < d["sigma_elastic"] else "elastic buckling"
        print(f"  wall {w:.2f} mm  rho* {rho:.3f}  plastic collapse "
              f"{d['sigma_plastic']:7.1f} MPa   elastic buckling "
              f"{d['sigma_elastic']:8.1f} MPa   -> {gov} governs")
    out["gibson_ashby"] = ga

    # ---- analytic cross-check --------------------------------------------
    g, P = K.geom()
    Ic = P["wing_width"] * P["wing_thick"] ** 3 / 12.0
    A = P["wing_width"] * P["wing_thick"]
    L = out["wing"]["span"]
    Eb = 20.0e3            # GPa-ish mid-grade effective modulus
    Gb = Eb / (2 * (1 + K.NU_TI))
    k_eb = 3 * Eb * Ic / L ** 3
    k_tim = 1.0 / (L ** 3 / (3 * Eb * Ic) + L / ((5.0 / 6.0) * Gb * A))
    print("\n" + "=" * 78)
    print("ANALYTIC CROSS-CHECK -- tip-loaded cantilever, %.1fx%.1f mm section"
          % (P["wing_thick"], P["wing_width"]))
    print("=" * 78)
    print(f"  L = {L:.2f} mm, I = {Ic:.2f} mm^4, E_eff = {Eb/1e3:.0f} GPa")
    print(f"  Euler-Bernoulli   k = {k_eb:8.1f} N/mm")
    print(f"  Timoshenko        k = {k_tim:8.1f} N/mm  (shear adds "
          f"{100*(k_eb/k_tim-1):.0f}% compliance at L/t = {L/P['wing_thick']:.1f})")
    print(f"  FEA (flat platen) k = {out['wing']['steps'][0]['k']:8.1f} N/mm "
          f"at the smallest step -- higher than the beam value because a flat "
          f"platen\n    contacts a patch, not a point, so local indentation "
          f"stiffness adds to bending.")
    out["analytic"] = dict(L=L, I=Ic, E=Eb, k_euler=k_eb, k_timoshenko=k_tim)

    dest = os.path.join(os.path.dirname(os.path.abspath(__file__)), "wing_results.json")
    with open(dest, "w") as fh:
        json.dump(out, fh, indent=2)
    print(f"\nwrote {dest}")


if __name__ == "__main__":
    main()
