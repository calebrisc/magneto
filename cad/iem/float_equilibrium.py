#!/usr/bin/env python3
"""
float_equilibrium.py -- coupled body-float equilibrium for expanding plungers.

`plunger_stroke_spec.py` measured the gap at a RIGID seat: body pinned, distance
from boss to flesh.  That is the wrong physics for a self-filling plunger.  The
body FLOATS: each stack pushes until skirt, jacket and the other stack balance
it, so the stroke a site actually consumes is an equilibrium result, not a
measurement of the static gap.

MODEL (3-DOF translation, which the brief allows as the minimum)
    Unknown: body translation `t` from the seated pose.  Forces on the body:

    plungers   each stack delivers a FLAT force F_p along -aim (reaction of
               pushing the ear) whenever the surface distance from the jacket
               outer face along its aim lies inside [L_c, L_c + S] -- compacted
               length to full extension.  Outside that it is either bottomed
               (surface closer than compacted, handled as rigid interference) or
               out of reach, contributing nothing.
    skirt      penalty spring over the 24-point sealing ring: each penetrating
               point pushes back along the ear normal, total capped at its
               preload.  Its compliance window is the same deformation budget
               the seal metric uses.
    jacket     same, on the ear-facing gyroid skin.
    rigid      core/faceplate get a stiff penalty so the solve cannot resolve
               the balance by burying the shell.

    Solved by least squares on the residual force.  Rotation is not solved; the
    brief permits translation-only, and the honest consequence is stated in the
    output rather than hidden.

Usage:
    python float_equilibrium.py --json-dir DIR [--stroke 6.0] [--compacted 4.5]
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

from earfit import ALIGNED, NOZZLE_AXIS, EarField, iem_points, transform

F_PLUNGER = 0.33      # N, mid of the 0.18-0.49 band (flat-force approximation)
SKIRT_MAX = 0.31      # N total
K_SKIRT = 0.124       # N/mm per unit penetration, ~0.31 N over a 2.5 mm window
K_JACKET = 0.50
K_RIGID = 20.0


def ear_bits(rec):
    patch = trimesh.load(os.path.join(ALIGNED, rec["patch"]), force="mesh")
    return patch, EarField(patch, seed=0)


def site_geometry():
    import generate
    P = generate.PARAMS
    g = generate.G(P)
    cc, cr = np.array(g.core_c, float), np.array(g.core_r, float)
    off = P["clearance"] + P["jacket_thick"]
    out = []
    for pl in g.plungers:
        a = np.asarray(pl["aim"], float); a /= np.linalg.norm(a)
        surf = cc + (cr ** 2 * a) / np.linalg.norm(cr * a)
        out.append((pl["name"], a, surf + a * off))     # base = jacket outer face
    return out


def solve_ear(rec, P, patch, field, sites, L_c, S):
    M = np.array(rec["transform"], float)
    R = M[:3, :3]
    ring = transform(__import__("seal_compliance").rim_from_mesh(n=72), M)
    jac = transform(P["jacket"], M)
    rig = transform(np.vstack([P["core_s"], P["face_s"]]), M)
    bases = [(nm, R @ a / np.linalg.norm(R @ a), transform(b[None, :], M)[0])
             for nm, a, b in sites]

    def gap(org, d):
        loc, ir, _ = patch.ray.intersects_location(org[None, :], d[None, :],
                                                   multiple_hits=False)
        return float(np.linalg.norm(loc[0] - org)) if len(ir) else np.inf

    def resid(t):
        F = np.zeros(3)
        for _nm, d, b in bases:
            D = gap(b + t, d)
            if L_c <= D <= L_c + S:
                F += F_PLUNGER * (-d)
        for pts, k, cap in ((ring, K_SKIRT, SKIRT_MAX), (jac, K_JACKET, None),
                            (rig, K_RIGID, None)):
            q = pts + t
            dd = field.query(q)
            pen = np.clip(-dd, 0.0, None)
            if not pen.any():
                continue
            _, idx = field.tree.query(q[pen > 0])
            n = field.nrm[idx]
            f = (k * pen[pen > 0])[:, None] * n
            tot = f.sum(axis=0)
            if cap is not None and np.linalg.norm(tot) > cap:
                tot = tot / np.linalg.norm(tot) * cap
            F += tot
        return F

    r = least_squares(resid, np.zeros(3), bounds=(-6.0, 6.0), xtol=1e-3,
                      ftol=1e-3, diff_step=0.2)
    t = r.x
    out = {}
    for nm, d, b in bases:
        D = gap(b + t, d)
        eng = L_c <= D <= L_c + S
        out[nm] = dict(D=D, consumed=(D - L_c) if np.isfinite(D) else np.inf,
                       engaged=bool(eng))
    out["_t"] = t.tolist()
    out["_resid"] = float(np.linalg.norm(r.fun))
    return out


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--json-dir", default=None)
    ap.add_argument("--compacted", type=float, default=4.5)
    ap.add_argument("--stroke", type=float, default=6.0)
    a = ap.parse_args()
    js = sorted(glob.glob(os.path.join(a.json_dir or ALIGNED, "*.json")))
    if not js:
        sys.exit("no seatings")
    P = iem_points()
    sites = site_geometry()
    L_c, S = a.compacted, a.stroke
    print(f"\nn = {len(js)} ears | compacted L_c = {L_c:.2f} mm, stroke S = "
          f"{S:.2f} mm -> working window on surface distance "
          f"[{L_c:.2f}, {L_c + S:.2f}] mm")
    print(f"flat plunger force {F_PLUNGER:.2f} N per stack "
          f"(band 0.18-0.49); 3-DOF translation only\n")

    per = {nm: [] for nm, _, _ in sites}
    eng = {nm: 0 for nm, _, _ in sites}
    shifts = []
    print(f"{'ear':<22}{'shift mm':>9} " +
          " ".join(f"{nm[:12]:>14}" for nm, _, _ in sites))
    for p in js:
        rec = json.load(open(p))
        patch, field = ear_bits(rec)
        r = solve_ear(rec, P, patch, field, sites, L_c, S)
        t = np.array(r["_t"]); shifts.append(np.linalg.norm(t))
        cells = []
        for nm, _, _ in sites:
            d = r[nm]
            per[nm].append(d["consumed"])
            eng[nm] += int(d["engaged"])
            cells.append(("--" if not np.isfinite(d["D"]) else f"{d['consumed']:+.2f}")
                         + ("*" if d["engaged"] else " "))
        print(f"{rec['ear_id']:<22}{np.linalg.norm(t):9.2f} " +
              " ".join(f"{c:>14}" for c in cells))
    print("\n(* = stack engaged at equilibrium; value = stroke consumed from "
          "compacted)\n")
    print(f"body shift at equilibrium: median {np.median(shifts):.2f} mm, "
          f"max {max(shifts):.2f} mm\n")
    print("| site | engaged | stroke consumed min | median | p90 | max |")
    print("|---|---|---|---|---|---|")
    for nm, _, _ in sites:
        v = np.array(per[nm]); v = v[np.isfinite(v)]
        if not len(v):
            print(f"| {nm} | {eng[nm]}/{len(js)} | -- | -- | -- | -- |")
            continue
        print(f"| {nm} | {eng[nm]}/{len(js)} | {v.min():.2f} | "
              f"{np.median(v):.2f} | {np.percentile(v,90):.2f} | {v.max():.2f} |")
    return 0


if __name__ == "__main__":
    sys.exit(main())
