"""Mag-plunger wing: sizing study for the magnet-sprung contact rail.

Competing concept to the Ti macro-gyroid sheet (docs/MECH_VALIDATION.md §3-4).
Instead of a Ti spring, a rigid contact rail (~14 x 6 mm, soft pad on top) rides
on guide pins from the jacket and is pushed toward the antihelix by repelling
bonded-NdFeB pairs.  A cam preset takes up coarse per-ear engagement, so the
magnets only have to cover +/- 1.5 mm of DYNAMIC travel about a rest point.

Requirement: total force (all pairs summed) inside 0.15-0.5 N across the whole
+/- 1.5 mm window; lowest max/min ratio wins.

Method matches docs/MAGFLOAT_MAGNETS.md exactly so the two studies are
comparable: magpylib 5.0.1 + magpylib-force 0.3.1 (Maxwell stress / meshed
dipole, a direct force solve, not a dU/dz finite difference).  Rings are
CylinderSegment over phi = 0-360, discs are Cylinder; axial magnetisation
M = Br/mu0; source +z, target -z so the facing poles repel; force is
getFT(source, target, anchor=target_position)[0][2].  200-cell meshing.

Run:  .venv/bin/python fea/mag_plunger.py
"""

import itertools
import json
import math
import os
import sys

import numpy as np

import magpylib as magpy
from magpylib_force import getFT

MU0 = 4e-7 * math.pi

MATERIALS = {
    "bonded NdFeB": dict(Br=0.65, rho=6.0),
    "N35 sintered": dict(Br=1.20, rho=7.5),
}

# (label, OD, ID, thickness) in mm -- ID = 0 is a solid disc
CANDIDATES = [
    ("3x1",     3.0, 0.0, 1.0),
    ("4x1",     4.0, 0.0, 1.0),
    ("4x2x1",   4.0, 2.0, 1.0),
    ("5x2.5x1", 5.0, 2.5, 1.0),
]

TRAVEL = 1.5        # mm, +/- about the rest point
F_LO, F_HI = 0.15, 0.50
GAP_MIN = 0.30      # mm, closest approach allowed (encapsulation + clearance)


MM = 1e-3        # magpylib 5 is SI: lengths in METRES, polarisation in tesla.
                 # Feeding it millimetres silently scales every force by 1e6
                 # (force goes as length^2) -- caught by validate() below.


def magnet(od, idd, t, Br, z, sign):
    """One axially magnetised puck centred at z (mm), magnetisation sign * Br."""
    pol = (0, 0, sign * Br)
    if idd <= 0:
        m = magpy.magnet.Cylinder(polarization=pol,
                                  dimension=(od * MM, t * MM),
                                  position=(0, 0, z * MM))
    else:
        m = magpy.magnet.CylinderSegment(
            polarization=pol,
            dimension=(idd / 2 * MM, od / 2 * MM, t * MM, 0, 360),
            position=(0, 0, z * MM))
    m.meshing = 200
    return m


def force(od, idd, t, Br, gap):
    """Axial repulsion (N) between two identical coaxial pucks at a face gap."""
    src = magnet(od, idd, t, Br, -(gap + t) / 2.0, +1)
    tgt = magnet(od, idd, t, Br, +(gap + t) / 2.0, -1)
    F = getFT(src, tgt, anchor=tgt.position)[0]
    return float(F[2])


def mass_g(od, idd, t, rho):
    return math.pi / 4.0 * (od ** 2 - idd ** 2) * t * rho * 1e-3


def stray_mT(od, idd, t, Br, dist):
    """On-axis |B| at `dist` mm beyond the back face of one puck."""
    m = magnet(od, idd, t, Br, 0.0, +1)
    B = m.getB([0, 0, (t / 2 + dist) * MM])
    return float(np.linalg.norm(B)) * 1e3


def validate():
    """Reproduce a published row of docs/MAGFLOAT_MAGNETS.md as a method check."""
    want = {0.5: 5.625, 1.0: 3.328, 2.0: 1.490, 3.0: 0.821}   # 7x3x1.5, N52 Br=1.45
    print("method check vs docs/MAGFLOAT_MAGNETS.md, 7x3x1.5 mm N52 (Br 1.45 T):")
    ok = True
    for g, ref in want.items():
        got = force(7.0, 3.0, 1.5, 1.45, g)
        err = 100 * (got / ref - 1)
        ok &= abs(err) < 2.0
        print(f"    gap {g:.1f} mm   this run {got:6.3f} N   published {ref:6.3f} N"
              f"   {err:+.2f}%")
    print(f"  -> {'MATCH' if ok else 'MISMATCH'}\n")
    return ok


def best_window(od, idd, t, Br, npairs):
    """Best rest gap for a given puck and pair count, or None if nothing fits."""
    best = None
    for rest in np.arange(GAP_MIN + TRAVEL, 4.001, 0.05):
        gaps = (rest - TRAVEL, rest, rest + TRAVEL)
        f = [npairs * force(od, idd, t, Br, g) for g in gaps]
        if min(f) < F_LO or max(f) > F_HI:
            continue
        ratio = max(f) / min(f)
        if best is None or ratio < best["ratio"]:
            best = dict(rest=float(rest), f_near=f[0], f_rest=f[1], f_far=f[2],
                        ratio=float(ratio))
    return best


DEPTH_BUDGET = 5.0   # mm of jacket depth available for the plunger stack


def stack_depth(t, rest):
    """Axial depth the plunger needs: fixed puck + rest gap + moving puck."""
    return 2 * t + rest


def search(cands, mats, pairs, rest_hi=9.0, depth_cap=None, travel=None):
    """Every (puck, material, npairs) that holds the band over +/- TRAVEL."""
    tr = TRAVEL if travel is None else travel
    hits = []
    for (lab, od, idd, t), mat, n in itertools.product(cands, mats, pairs):
        Br = MATERIALS[mat]["Br"]
        for rest in np.arange(GAP_MIN + tr, rest_hi + 1e-9, 0.05):
            if depth_cap is not None and stack_depth(t, rest) > depth_cap:
                continue
            f = [n * force(od, idd, t, Br, g)
                 for g in (rest - tr, rest, rest + tr)]
            if min(f) < F_LO or max(f) > F_HI:
                continue
            hits.append(dict(puck=lab, od=od, id=idd, t=t, material=mat, Br=Br,
                             npairs=n, travel=tr, rest=float(rest), f_near=f[0],
                             f_rest=f[1], f_far=f[2], ratio=max(f) / min(f),
                             depth=stack_depth(t, rest),
                             mass_mg=1e3 * 2 * n * mass_g(od, idd, t,
                                                          MATERIALS[mat]["rho"])))
            break        # smallest feasible rest gap for this combo
    return hits


def closest_miss(od, idd, t, Br, n, rest_hi=9.0, depth_cap=None):
    """Best (lowest) ratio reachable, and what stops it qualifying."""
    best = None
    for rest in np.arange(GAP_MIN + TRAVEL, rest_hi + 1e-9, 0.05):
        if depth_cap is not None and stack_depth(t, rest) > depth_cap:
            continue
        f = [n * force(od, idd, t, Br, g)
             for g in (rest - TRAVEL, rest, rest + TRAVEL)]
        r = max(f) / min(f)
        if best is None or r < best["ratio"]:
            binds = []
            if min(f) < F_LO:
                binds.append(f"F_far {min(f):.3f} < {F_LO} N")
            if max(f) > F_HI:
                binds.append(f"F_near {max(f):.3f} > {F_HI} N")
            if r > F_HI / F_LO:
                binds.append(f"ratio {r:.2f} > {F_HI/F_LO:.2f} allowed by the band")
            best = dict(rest=float(rest), ratio=r, f_near=f[0], f_far=f[2],
                        depth=stack_depth(t, rest), binds="; ".join(binds))
    return best


def main():
    sys.stdout.reconfigure(line_buffering=True)
    out = {"validation_ok": validate(), "sweep": [], "best": None}

    print(f"window: rest +/- {TRAVEL} mm ({2*TRAVEL} mm total), total force in "
          f"[{F_LO}, {F_HI}] N -> the band itself allows a ratio of only "
          f"{F_HI/F_LO:.2f}\n")

    listed = [c for c in CANDIDATES]
    print("STEP 1 -- the four listed pucks, inside the %.0f mm depth budget:"
          % DEPTH_BUDGET)
    hits = search(listed, MATERIALS, (1, 2), depth_cap=DEPTH_BUDGET)
    out["listed_in_budget"] = hits
    if hits:
        for h in sorted(hits, key=lambda r: r["ratio"]):
            print(f"    {h['puck']:>9} {h['material']:>13} {h['npairs']} pair "
                  f"rest {h['rest']:.2f}  ratio {h['ratio']:.2f}")
    else:
        print("    NONE FIT.  Closest miss per puck (2 pairs, best material):")
        for lab, od, idd, t in listed:
            for mat in MATERIALS:
                cm = closest_miss(od, idd, t, MATERIALS[mat]["Br"], 2,
                                  depth_cap=DEPTH_BUDGET)
                print(f"      {lab:>9} {mat:>13}  best rest {cm['rest']:.2f} mm "
                      f"(depth {cm['depth']:.1f}) ratio {cm['ratio']:6.2f}  "
                      f"F {cm['f_near']:.3f} -> {cm['f_far']:.3f} N")
                out.setdefault("closest_miss", []).append(dict(puck=lab, material=mat, **cm))

    print("\nSTEP 2 -- why: the ratio is set by geometry, not by magnet strength.")
    print("    Over +/- 1.5 mm the gap swings from (rest-1.5) to (rest+1.5), so a")
    print("    small rest gap is a huge relative swing.  Ratio vs rest gap and OD")
    print("    (N35, ratio must be <= %.2f):" % (F_HI / F_LO))
    big = [("5x2.5x1", 5, 2.5, 1.0), ("8x4x1.5", 8, 4, 1.5), ("10x5x2", 10, 5, 2.0),
           ("14x7x3", 14, 7, 3.0)]
    rests = (2.0, 3.0, 4.0, 5.0, 6.0, 7.0)
    print("      %-10s" % "puck" + "".join("%8s" % ("r=%.0f" % r) for r in rests))
    grid = []
    for lab, od, idd, t in big:
        row = []
        for r in rests:
            f = [force(od, idd, t, 1.20, g) for g in (r - TRAVEL, r + TRAVEL)]
            row.append(f[0] / f[1])
        grid.append(dict(puck=lab, rests=list(rests), ratios=row))
        print("      %-10s" % lab + "".join("%8.2f" % v for v in row))
    out["ratio_grid"] = grid
    print("\n    Depth needed = 2t + rest gap.  Budget is %.0f mm, so:" % DEPTH_BUDGET)
    for lab, od, idd, t in big:
        print(f"      {lab:>9}: max in-budget rest gap {DEPTH_BUDGET - 2*t:.1f} mm")

    print("\nSTEP 3 -- what WOULD work if the depth budget is relaxed:")
    ext = [("5x2.5x1", 5, 2.5, 1.0), ("6x3x1.5", 6, 3, 1.5), ("8x4x1.5", 8, 4, 1.5),
           ("10x5x2", 10, 5, 2.0), ("12x6x2", 12, 6, 2.0), ("14x7x3", 14, 7, 3.0)]
    feas = search(ext, MATERIALS, (1, 2), rest_hi=9.0)
    out["feasible_relaxed"] = feas
    if not feas:
        print("    nothing in the extended set fits either.")
        json.dump(out, open(os.path.join(os.path.dirname(os.path.abspath(__file__)),
                                         "mag_plunger_results.json"), "w"), indent=2)
        return out
    feas.sort(key=lambda r: (r["depth"], r["ratio"]))
    print(f"    {'puck':>9} {'material':>13} {'pairs':>5} {'rest':>6} {'depth':>6} "
          f"{'F near':>7} {'F rest':>7} {'F far':>7} {'ratio':>6} {'mass mg':>8}")
    for h in feas[:8]:
        print(f"    {h['puck']:>9} {h['material']:>13} {h['npairs']:>5} "
              f"{h['rest']:6.2f} {h['depth']:6.2f} {h['f_near']:7.3f} "
              f"{h['f_rest']:7.3f} {h['f_far']:7.3f} {h['ratio']:6.2f} "
              f"{h['mass_mg']:8.1f}")
    mf = feas[0]
    out["best_relaxed"] = mf
    print(f"\n    MINIMUM-DEPTH FEASIBLE: {mf['puck']} {mf['material']}, "
          f"{mf['npairs']} pair(s), rest {mf['rest']:.2f} mm -> needs "
          f"{mf['depth']:.2f} mm of depth vs the {DEPTH_BUDGET:.0f} mm budget "
          f"({mf['depth']/DEPTH_BUDGET:.1f}x over).")

    print("\nSTEP 4 -- the fallback: keep +/- 0.75 mm dynamic travel instead.")
    f075 = search(listed, MATERIALS, (1, 2), depth_cap=DEPTH_BUDGET, travel=0.75)
    out["travel_0p75"] = f075
    if f075:
        f075.sort(key=lambda r: r["ratio"])
        for h in f075[:4]:
            print(f"    {h['puck']:>9} {h['material']:>13} {h['npairs']} pair "
                  f"rest {h['rest']:.2f} depth {h['depth']:.2f}  F "
                  f"{h['f_near']:.3f}->{h['f_far']:.3f} N  ratio {h['ratio']:.2f}"
                  f"  {h['mass_mg']:.0f} mg")
    else:
        print("    still nothing inside the depth budget.")

    # profile the best design that actually FITS the depth budget
    b = (sorted(f075, key=lambda r: r["ratio"])[0] if f075 else mf)
    out["profiled"] = b
    print(f"\n  profiling the best IN-BUDGET design: {b['puck']} {b['material']}, "
          f"{b['npairs']} pair(s), +/-{b['travel']} mm travel, "
          f"{b['depth']:.2f} mm deep")

    # ---- mechanical envelope of the winner --------------------------------
    print("\n" + "=" * 78)
    print("MECHANICAL ENVELOPE OF THE WINNER")
    print("=" * 78)
    od, idd, t, Br = b["od"], b["id"], b["t"], b["Br"]
    n = b["npairs"]

    # axial magnetic stiffness at the rest point (central difference)
    dh = 0.05
    kmag = n * (force(od, idd, t, Br, b["rest"] - dh)
                - force(od, idd, t, Br, b["rest"] + dh)) / (2 * dh)
    print(f"  axial magnetic stiffness at rest   {kmag:.4f} N/mm "
          f"(the Ti sheet's k is 0.15-0.35 by comparison -- but note a magnet "
          f"spring\n                                     is SOFTER for the same "
          f"force, which is the whole point)")

    # guide pins: L/D >= 2 for no cocking, inside a 5 mm jacket depth budget
    L_guide = 5.0
    D_max = L_guide / 2.0
    print(f"  guide pin: {L_guide:.1f} mm engaged length in the 5 mm depth budget "
          f"-> D <= {D_max:.1f} mm for L/D >= 2")
    print(f"             two pins on a 14 mm rail; pin OD 2.0 mm gives L/D = "
          f"{L_guide/2.0:.1f}")

    # rocking restoring torque for a two-pair layout at the rail ends
    if n == 2:
        s = 10.0        # mm between the two plungers on a 14 mm rail
        k_tors = 0.5 * (kmag / n) * s ** 2       # N.mm per rad
        print(f"  two plungers {s:.0f} mm apart -> torsional restoring stiffness "
              f"{k_tors:.2f} N.mm/rad")
        print(f"             = {k_tors*math.radians(1):.3f} N.mm per degree of "
              f"rail rock; the rail self-levels, the pins only take residual")
    else:
        print("  ONE pair: no restoring torque about the rail axis -- the guide "
              "pins carry\n             all the rocking moment.")

    # masses
    m_mag = b["mass_mg"]
    rail = dict(l=14.0, w=6.0, t=0.8, rho_ti=4.43)
    m_rail = rail["l"] * rail["w"] * rail["t"] * rail["rho_ti"] * 1e-3 * 1e3
    m_pins = 2 * math.pi / 4 * 2.0 ** 2 * 5.0 * 4.43e-3 * 1e3
    m_pad = 14.0 * 6.0 * 1.0 * 1.1e-3 * 1e3
    print(f"\n  added mass: magnets {m_mag:.0f} mg + rail "
          f"({rail['l']:.0f}x{rail['w']:.0f}x{rail['t']:.1f} Ti) {m_rail:.0f} mg "
          f"+ 2 pins {m_pins:.0f} mg + pad {m_pad:.0f} mg")
    print(f"              TOTAL {m_mag + m_rail + m_pins + m_pad:.0f} mg "
          f"per side")
    out["mass_mg"] = dict(magnets=m_mag, rail=m_rail, pins=m_pins, pad=m_pad,
                          total=m_mag + m_rail + m_pins + m_pad)
    out["k_axial_N_per_mm"] = kmag

    # stray field at the bone-conduction sensor
    print("\n  stray field at the bone sensor, 8 mm from the puck face:")
    for lab, od_, id_, t_ in CANDIDATES:
        s8 = stray_mT(od_, id_, t_, MATERIALS[b["material"]]["Br"], 8.0)
        mark = "  <- winner" if lab == b["puck"] else ""
        print(f"    {lab:>9}  {s8:7.3f} mT{mark}")
    out["stray_mT_at_8mm"] = stray_mT(od, idd, t, Br, 8.0)

    dest = os.path.join(os.path.dirname(os.path.abspath(__file__)),
                        "mag_plunger_results.json")
    with open(dest, "w") as fh:
        json.dump(out, fh, indent=2)
    print(f"\nwrote {dest}")
    return out


if __name__ == "__main__":
    main()
