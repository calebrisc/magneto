#!/usr/bin/env python3
"""
stability.py -- quasi-static retention-under-load for a seated IEM.

The contact contract checks that each part TOUCHES what it should.  It says
nothing about whether the assembly STAYS PUT when something pulls on it, which is
a different failure and the one a wearer actually notices.  This is a rigid-body
force/moment balance -- screw theory with Coulomb friction cones -- not FEA.

LOADS TO RESIST
    (a) skirt preload reaction   the compressed sealing skirt pushes the shell
                                 back out.  Modelled as its contact normals with
                                 a fixed total of SKIRT_PRELOAD N, whose resultant
                                 IS the outward axial push -- so it appears once,
                                 as both the destabilising force and the source of
                                 friction at the land.  Always on, never scaled.
    (b) cable tug                CABLE_TUG N at the cable exit, worst direction in
                                 a downward-backward cone.
    (c) inertial                 G_LOAD x g on ASSEMBLY_MASS, worst direction over
                                 the sphere.

RESISTANCE
    Each contact i can push the shell along the ear's outward surface normal n_i
    (flesh pushes, never pulls) with N_i >= 0, plus friction in the tangent plane
    bounded by |t_i| <= MU * N_i.  The friction cone is linearised into an
    8-sided pyramid, which is conservative (inscribed).

    Normal-force budgets are set by what actually presses each contact:
      skirt land   total normal fixed at SKIRT_PRELOAD (set by its compression)
      wing pad     total normal fixed at k_wing x interference (a leaf spring)
      jacket face  free reaction, >= 0 -- it pushes only while something presses
                   the shell onto the concha floor, which equilibrium enforces
      plungers     PLUNGER_FORCE each, if the build has any (this one does not)

    Geometric interlock needs no special term: it falls out of the contact normals
    opposing the load direction.

SCORE
    One LP per sampled load direction, with the load scale `s` as a variable:
    maximise s subject to equilibrium, the friction pyramids and the normal
    budgets.  s* >= 1 means the spec loads are resisted, and s* is the margin.
    The reported verdict is the WORST s* over all sampled directions.

    stable    worst s* >= 1.5
    marginal  worst s* >= 1.0
    FAIL      worst s* <  1.0, and the failing direction is reported.

ASSUMPTIONS worth arguing with:
  * MU = 0.4 for silicone/Ti against skin.  Published skin friction against
    smooth elastomers spans roughly 0.3-1.0 depending on moisture, so 0.4 is a
    deliberately dry/conservative pick.  Raise it and everything gets easier.
  * The cable exit is NOT modelled in the build; its position is proxied by the
    2-pin socket pocket centre and flagged.
  * Rigid body, rigid ear, no skin deformation and no time dependence.
"""

from __future__ import annotations

import numpy as np
from scipy.optimize import linprog

MU = 0.40
SKIRT_PRELOAD = 0.31          # N, total normal at the sealing land
CABLE_TUG = 0.50              # N
G_LOAD = 3.0
ASSEMBLY_MASS = 0.008         # kg
PLUNGER_FORCE = (0.18, 0.49)  # N each, when the build has plungers
K_WING = 0.294                # N/mm
CONTACT_BAND = 1.0            # mm; |signed distance| under this counts as contact
PYRAMID_SIDES = 8
CONE_HALF_DEG = 45.0
JACKET_CAP = 20.0             # N, effectively a free reaction


def _basis(n):
    a = np.array([0.0, 0.0, 1.0])
    if abs(n @ a) > 0.9:
        a = np.array([1.0, 0.0, 0.0])
    u = np.cross(n, a)
    u /= np.linalg.norm(u)
    return u, np.cross(n, u)


def _generators(pts, nrms, com):
    """Friction-pyramid wrench generators for a set of contacts.

    Returns (W, owner) where W is 6 x m and owner[k] is the contact index of
    generator k, so per-contact normal budgets can be summed.
    """
    W, owner = [], []
    th = np.linspace(0, 2 * np.pi, PYRAMID_SIDES, endpoint=False)
    for i, (p, n) in enumerate(zip(pts, nrms)):
        n = n / max(np.linalg.norm(n), 1e-12)
        u, v = _basis(n)
        r = p - com
        for t in th:
            g = n + MU * (np.cos(t) * u + np.sin(t) * v)
            W.append(np.concatenate([g, np.cross(r, g)]))
            owner.append(i)
    if not W:
        return np.zeros((6, 0)), np.zeros(0, int)
    return np.array(W).T, np.array(owner)


def _cone_dirs(axis, half_deg, n):
    axis = axis / np.linalg.norm(axis)
    u, v = _basis(axis)
    out = [axis]
    for ring in (0.5, 1.0):
        for t in np.linspace(0, 2 * np.pi, 6, endpoint=False):
            a = np.radians(half_deg) * ring
            out.append(np.cos(a) * axis
                       + np.sin(a) * (np.cos(t) * u + np.sin(t) * v))
    return [d / np.linalg.norm(d) for d in out]


def _sphere_dirs(n=14):
    g = (1 + 5 ** 0.5) / 2
    i = np.arange(n)
    z = 1 - 2 * (i + 0.5) / n
    r = np.sqrt(np.clip(1 - z * z, 0, 1))
    ph = 2 * np.pi * i / g
    return [np.array([r[k] * np.cos(ph[k]), r[k] * np.sin(ph[k]), z[k]])
            for k in range(n)]


def _max_scale(W, owner, groups, w_always, w_scaled):
    """Largest s with  W.lam + w_always + s.w_scaled = 0  under the budgets.

    Normal-force budgets are CAPS, not equalities.  A preload bounds how hard a
    contact can push, and therefore how much friction it offers; it does not
    force that contact to be loaded.  Writing them as equalities over-constrains
    the balance and reports every pose as infeasible -- which is exactly what a
    first pass of this did, returning margin 0.00 on a pose that is clearly not
    falling out of the ear.
    """
    m = W.shape[1]
    if m == 0:
        return 0.0
    A_eq = np.hstack([W, w_scaled.reshape(6, 1)])
    b_eq = -w_always
    rows_ub, rhs_ub = [], []
    for idx, (_kind, cap) in groups.items():
        sel = np.zeros(m + 1)
        sel[:m] = (owner == idx).astype(float)
        rows_ub.append(sel)
        rhs_ub.append(cap)
    c = np.zeros(m + 1)
    c[-1] = -1.0                                           # maximise s
    bounds = [(0, None)] * m + [(0, 20.0)]
    res = linprog(c,
                  A_ub=np.array(rows_ub) if rows_ub else None,
                  b_ub=np.array(rhs_ub) if rows_ub else None,
                  A_eq=A_eq, b_eq=b_eq, bounds=bounds, method="highs")
    return float(res.x[-1]) if res.success else 0.0


def stability_check(rec, P, field, transform_fn, contacts=None,
                    plungers=0, cable_point=None, com=None):
    """Worst-case load scale the seated pose can resist.  See module docstring."""
    M = np.array(rec["transform"], float)

    # ---- gather contacts, grouped by the part that presses them ------------- #
    groups, pts, nrms, owner_of = {}, [], [], []

    def add(key, keys, kind, cap):
        got = []
        for k in keys:
            if k not in P or not len(P[k]):
                continue
            w = transform_fn(P[k], M)
            d = field.query(w)
            for j in np.where(np.abs(d) <= CONTACT_BAND)[0]:
                got.append(w[j])
        if not got:
            return 0
        # decimate to at most 12 representative points per group
        got = np.array(got)
        if len(got) > 12:
            got = got[np.linspace(0, len(got) - 1, 12).astype(int)]
        gi = len(groups)
        for p in got:
            _, idx = field.tree.query(p[None, :])
            pts.append(p); nrms.append(field.nrm[int(idx[0])]); owner_of.append(gi)
        groups[gi] = (kind, cap)
        return len(got)

    n_skirt = add("skirt", ["rim", "soft"], "cap", SKIRT_PRELOAD)
    tip = float(np.median(field.query(transform_fn(P["wing_tip"], M))))
    f_wing = K_WING * max(-tip, 0.0)
    n_wing = add("wing", ["wing_tip", "wing_mid"], "cap", f_wing) if f_wing > 0 else 0
    n_jack = add("jacket", ["jacket"], "cap", JACKET_CAP)
    if plungers:
        add("plungers", ["plunger"], "cap", plungers * PLUNGER_FORCE[0])

    if not pts:
        return dict(verdict="FAIL", margin=0.0, reason="no contacts at all",
                    n_skirt=0, n_wing=0, n_jacket=0, f_wing=f_wing, mu=MU)

    pts = np.array(pts); nrms = np.array(nrms)
    owner = np.array(owner_of)
    if com is None:
        com = pts.mean(axis=0)

    W, own2 = [], []
    th = np.linspace(0, 2 * np.pi, PYRAMID_SIDES, endpoint=False)
    for i, (p, n) in enumerate(zip(pts, nrms)):
        n = n / max(np.linalg.norm(n), 1e-12)
        u, v = _basis(n)
        r = p - com
        for t in th:
            g = n + MU * (np.cos(t) * u + np.sin(t) * v)
            W.append(np.concatenate([g, np.cross(r, g)]))
            own2.append(owner[i])
    W = np.array(W).T
    own2 = np.array(own2)

    # ---- disturbance wrenches ---------------------------------------------- #
    if cable_point is None:
        cable_point = com
    f_inert = G_LOAD * ASSEMBLY_MASS * 9.81
    # downward-backward cone: -z inferior, -x posterior (scan frame)
    cable_dirs = _cone_dirs(np.array([-1.0, 0.0, -1.0]), CONE_HALF_DEG, 13)
    inert_dirs = _sphere_dirs(14)

    # (a) skirt preload reaction: always on, never scaled.  Acts outward along
    # the nozzle axis, applied at the rim centroid.
    from earfit import NOZZLE_AXIS
    nax = M[:3, :3] @ NOZZLE_AXIS
    nax = nax / np.linalg.norm(nax)
    rim_c = transform_fn(P["rim"], M).mean(axis=0)
    F_sk = -SKIRT_PRELOAD * nax               # outward = away from the canal
    w_always = np.concatenate([F_sk, np.cross(rim_c - com, F_sk)])

    # Diagnostics that make a FAIL actionable: how much straight pull-out the
    # contacts can hold at all, and how many of them actually INTERLOCK (normal
    # opposing the escape direction) rather than pushing the shell out.
    pull = -nax
    w_pull = np.concatenate([pull, np.cross(np.zeros(3), pull)])
    pullout_capacity = _max_scale(W, own2, groups, np.zeros(6), w_pull)
    comp = nrms @ pull
    interlock = int((comp < 0).sum())
    demand = SKIRT_PRELOAD + CABLE_TUG + f_inert

    worst, worst_dirs = 1e9, None
    for cd in cable_dirs:
        for idd in inert_dirs:
            F = CABLE_TUG * cd + f_inert * idd
            Mo = np.cross(cable_point - com, CABLE_TUG * cd)
            w = np.concatenate([F, Mo])
            sc = _max_scale(W, own2, groups, w_always, w)
            if sc < worst:
                worst, worst_dirs = sc, (cd, idd)

    verdict = "stable" if worst >= 1.5 else "marginal" if worst >= 1.0 else "FAIL"
    return dict(verdict=verdict, margin=worst,
                pullout_capacity=pullout_capacity, demand=demand,
                interlock=interlock, n_contacts=len(pts),
                friction_budget=MU * (SKIRT_PRELOAD + f_wing),
                cable_dir=None if worst_dirs is None else worst_dirs[0].tolist(),
                inert_dir=None if worst_dirs is None else worst_dirs[1].tolist(),
                n_skirt=n_skirt, n_wing=n_wing, n_jacket=n_jack,
                f_wing=f_wing, tip=tip, mu=MU,
                f_inert=f_inert, plungers=plungers)
