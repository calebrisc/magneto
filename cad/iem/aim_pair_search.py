#!/usr/bin/env python3
"""
aim_pair_search.py -- joint search over the two plunger aims.

Every aim so far was a guess checked afterwards.  This searches the pair space
directly against the four objectives, and reports what the best pair cannot do as
well as what it can.

OBJECTIVES
  1 coverage    both stacks find cartilage within stroke, at floating
                equilibrium, on as many of the 13 ears as possible
  2 reaction    the reaction sum on the body points INTO the concha floor --
                the one surface every ear presents and the jacket already rests
                on.  Floor-inward is (+0.18, +0.39, -0.90), measured.  A near
                cancelling pair (small |sum|) also satisfies this, since then the
                skirt and jacket carry only a couple.
  3 no press    no ear bottoms a stack above 2 N
  4 stability   >= 1.00x at mu 0.6 / over-ear routing on BOTH a small and a
                large ear

STAGES
  A  per-ear surface distance for every candidate aim at the nominal pose --
     one batched ray cast per ear, cheap, prunes aims that reach nobody
  B  analytic pair score on objective 2 -- free
  C  floating equilibrium (3-DOF translation) for the surviving pairs
  D  stability matrix for the winner

Usage:
    python aim_pair_search.py --grid 15 --top 12
"""

from __future__ import annotations

import argparse
import glob
import json
import os
import sys

import numpy as np
import trimesh
from scipy.optimize import least_squares

from earfit import ALIGNED, EarField, iem_points, transform

FLOOR_IN = np.array([0.18, 0.39, -0.90])
FLOOR_IN /= np.linalg.norm(FLOOR_IN)
F_PL = 0.33
K_SKIRT, SKIRT_MAX, K_JAC, K_RIG = 0.124, 0.31, 0.50, 20.0
BOTTOM_LIMIT = 2.0


def gen_bits():
    import generate
    P = generate.PARAMS
    g = generate.G(P)
    return (np.array(g.core_c, float), np.array(g.core_r, float),
            P["clearance"] + P["jacket_thick"], P)


def base_for(cc, cr, off, a):
    return cc + (cr ** 2 * a) / np.linalg.norm(cr * a) + a * off


def sphere_grid(step_deg):
    out = []
    for pol in np.arange(step_deg, 180.0, step_deg):
        n_az = max(1, int(round(360.0 / step_deg * np.sin(np.radians(pol)))))
        for az in np.linspace(0, 2 * np.pi, n_az, endpoint=False):
            t = np.radians(pol)
            out.append(np.array([np.sin(t) * np.cos(az), np.cos(t),
                                 np.sin(t) * np.sin(az)]))
    out.append(np.array([0.0, 1.0, 0.0]))
    return [v / np.linalg.norm(v) for v in out]


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--json-dir", default=None)
    ap.add_argument("--grid", type=float, default=15.0)
    ap.add_argument("--top", type=int, default=12)
    ap.add_argument("--compacted", type=float, default=4.5)
    ap.add_argument("--stroke", type=float, default=6.0)
    a = ap.parse_args()

    js = sorted(glob.glob(os.path.join(a.json_dir or ALIGNED, "*.json")))
    cc, cr, off, GP = gen_bits()
    P = iem_points()
    L_c, S = a.compacted, a.stroke

    ears = []
    for p in js:
        rec = json.load(open(p))
        patch = trimesh.load(os.path.join(ALIGNED, rec["patch"]), force="mesh")
        ears.append((rec, patch, EarField(patch, seed=0)))
    n = len(ears)

    aims = sphere_grid(a.grid)
    A = np.array(aims)
    print(f"n = {n} ears | {len(aims)} candidate aims on a {a.grid:.0f} deg grid")
    print(f"stack: compacted {L_c:.2f} mm, stroke {S:.2f} mm "
          f"-> window [{L_c:.2f}, {L_c + S:.2f}] mm\n")

    # ---- stage A: reach of every single aim, nominal pose ------------------ #
    D = np.full((n, len(aims)), np.inf)
    for i, (rec, patch, _f) in enumerate(ears):
        M = np.array(rec["transform"], float); R = M[:3, :3]
        org = np.array([transform(base_for(cc, cr, off, v)[None, :], M)[0] for v in A])
        dirs = (A @ R.T)
        dirs /= np.linalg.norm(dirs, axis=1, keepdims=True)
        loc, ir, _ = patch.ray.intersects_location(org, dirs, multiple_hits=False)
        for k, j in enumerate(ir):
            d = np.linalg.norm(loc[k] - org[j])
            if d < D[i, j]:
                D[i, j] = d
    inwin = (D >= L_c) & (D <= L_c + S)
    reach = inwin.sum(axis=0)
    keep = np.where(reach >= 3)[0]
    print(f"stage A: {len(keep)} aims reach >= 3 ears inside the stroke window "
          f"(best single aim {reach.max()}/{n})")

    # ---- stage B: analytic pair score on the reaction direction ------------ #
    pairs = []
    for x in range(len(keep)):
        for y in range(x + 1, len(keep)):
            j1, j2 = keep[x], keep[y]
            ssum = A[j1] + A[j2]
            mag = np.linalg.norm(ssum)
            react = -ssum
            align = (react @ FLOOR_IN) / mag if mag > 1e-6 else 1.0
            # good = reaction into the floor, OR a near-cancelling couple
            obj2 = max(align, 1.0 - mag)
            both = int((inwin[:, j1] & inwin[:, j2]).sum())
            pairs.append((both, obj2, mag, align, j1, j2))
    pairs.sort(key=lambda t: (t[0] + 2.0 * t[1]), reverse=True)
    print(f"stage B: {len(pairs)} pairs; best analytic obj2 "
          f"{max(p[1] for p in pairs):+.2f} "
          f"(1.0 = reaction straight into the floor)")
    best_align = max(pairs, key=lambda t: t[3])
    print(f"         best achievable floor-alignment over ALL pairs: "
          f"{best_align[3]:+.2f} at |sum| {best_align[2]:.2f}")
    cand = pairs[:a.top]

    # ---- stage C: floating equilibrium ------------------------------------- #
    print(f"\nstage C: floating equilibrium for the top {len(cand)} pairs\n")
    print(f"{'aim 1':>22} {'aim 2':>22} {'both':>5} {'obj2':>6} {'maxN':>7}")
    results = []
    for both0, obj2, mag, align, j1, j2 in cand:
        a1, a2 = A[j1], A[j2]
        nb, worstN = 0, 0.0
        det = []
        for rec, patch, field in ears:
            M = np.array(rec["transform"], float); R = M[:3, :3]
            bs = []
            for v in (a1, a2):
                d = R @ v; d /= np.linalg.norm(d)
                bs.append((d, transform(base_for(cc, cr, off, v)[None, :], M)[0]))
            ring = transform(__import__("seal_compliance").rim_from_mesh(n=72), M)
            jac = transform(P["jacket"], M)
            rig = transform(np.vstack([P["core_s"], P["face_s"]]), M)

            def gp(o, d):
                loc, ir, _ = patch.ray.intersects_location(
                    o[None, :], d[None, :], multiple_hits=False)
                return float(np.linalg.norm(loc[0] - o)) if len(ir) else np.inf

            def resid(t):
                F = np.zeros(3)
                for d, b in bs:
                    dd = gp(b + t, d)
                    if L_c <= dd <= L_c + S:
                        F += F_PL * (-d)
                    elif np.isfinite(dd) and dd < L_c:
                        F += K_RIG * (L_c - dd) * (-d)
                for pts, k, cap in ((ring, K_SKIRT, SKIRT_MAX),
                                    (jac, K_JAC, None), (rig, K_RIG, None)):
                    q = pts + t
                    sd = field.query(q)
                    pen = np.clip(-sd, 0.0, None)
                    if not pen.any():
                        continue
                    _, idx = field.tree.query(q[pen > 0])
                    f = (k * pen[pen > 0])[:, None] * field.nrm[idx]
                    tot = f.sum(axis=0)
                    if cap is not None and np.linalg.norm(tot) > cap:
                        tot = tot / np.linalg.norm(tot) * cap
                    F += tot
                return F

            r = least_squares(resid, np.zeros(3), bounds=(-6, 6),
                              xtol=1e-2, ftol=1e-2, diff_step=0.3)
            eng, bn = [], 0.0
            for d, b in bs:
                dd = gp(b + r.x, d)
                eng.append(L_c <= dd <= L_c + S)
                if np.isfinite(dd) and dd < L_c:
                    bn = max(bn, K_RIG * (L_c - dd))
            if all(eng):
                nb += 1
            worstN = max(worstN, bn)
            det.append((rec["ear_id"], eng, bn))
        results.append(dict(a1=a1, a2=a2, both=nb, obj2=obj2, mag=mag,
                            align=align, worstN=worstN, det=det))
        print(f"{np.round(a1,2)!s:>22} {np.round(a2,2)!s:>22} "
              f"{nb:5d} {obj2:6.2f} {worstN:7.2f}")

    results.sort(key=lambda r: (r["both"], r["obj2"], -r["worstN"]), reverse=True)
    b = results[0]
    print(f"\n=== BEST PAIR ===")
    print(f"  aim 1 {np.round(b['a1'],3)}   aim 2 {np.round(b['a2'],3)}")
    print(f"  both stacks engaged: {b['both']}/{n}")
    print(f"  reaction sum |{b['mag']:.2f}|, floor alignment {b['align']:+.2f} "
          f"(1.0 = straight into the floor)")
    print(f"  worst bottomed strut {b['worstN']:.2f} N "
          f"({'OK' if b['worstN'] <= BOTTOM_LIMIT else 'EXCEEDS 2 N'})")
    print(f"\n  per-ear: {'ear':<12}{'aim1':>6}{'aim2':>6}{'bottom N':>10}")
    for eid, eng, bn in b["det"]:
        print(f"           {eid:<12}{'yes' if eng[0] else 'no':>6}"
              f"{'yes' if eng[1] else 'no':>6}{bn:10.2f}")
    json.dump({k: (v.tolist() if isinstance(v, np.ndarray) else
                   (v if k != "det" else [[d[0], list(map(bool, d[1])), d[2]]
                                          for d in v]))
               for k, v in b.items()},
              open(os.path.join(ALIGNED, "aim_pair_best.json"), "w"), indent=1)
    return 0


if __name__ == "__main__":
    sys.exit(main())
