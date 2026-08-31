#!/usr/bin/env python3
"""
viz_equilibrium.py -- seated scenes at the FLOATING-EQUILIBRIUM pose.

Exports one GLB per ear with the body at the pose the coupled equilibrium found
(not the seating-optimiser pose) and each plunger stack drawn at the extension it
actually consumed on that ear.

The stacks are drawn as primitives -- a green cylinder for the magnet column, a
red disc for the pad -- because the winning aims come out of `aim_pair_search.py`
and no STL exists for them yet.  The primitives are dimensionally honest (they
start at the jacket outer face and end at the measured contact distance) but they
are a schematic of a stack, not a manufacturable part.

Usage:
    python viz_equilibrium.py --a1 x,y,z --a2 x,y,z --ears small=ID,median=ID,large=ID
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

import stability as stab
from earfit import ALIGNED, HERE, NOZZLE_T, EarField, iem_points, transform
from viz_scene import cluster_decimate

VIZ = os.path.join(HERE, "viz")
F_PL, K_SKIRT, SKIRT_MAX, K_JAC, K_RIG = 0.33, 0.124, 0.31, 0.50, 20.0
COL = {"ear": (205, 175, 135, 255), "core": (190, 192, 200, 255),
       "faceplate": (205, 207, 214, 255), "jacket_wing": (62, 108, 198, 255),
       "carrier": (232, 140, 48, 255), "nozzle_insert_short": (128, 128, 130, 255),
       "stack": (110, 200, 110, 255), "pad": (220, 90, 90, 255)}


def gen_bits():
    import generate
    P = generate.PARAMS
    g = generate.G(P)
    return (np.array(g.core_c, float), np.array(g.core_r, float),
            P["clearance"] + P["jacket_thick"], P)


def frame_for(a):
    a = a / np.linalg.norm(a)
    r = np.array([0.0, 0.0, 1.0])
    if abs(r @ a) > 0.9:
        r = np.array([1.0, 0.0, 0.0])
    u = np.cross(a, r); u /= np.linalg.norm(u)
    return a, u, np.cross(a, u)


def solve(rec, patch, field, P, bases, L_c, S):
    M = np.array(rec["transform"], float)
    ring = transform(__import__("seal_compliance").rim_from_mesh(n=72), M)
    jac = transform(P["jacket"], M)
    rig = transform(np.vstack([P["core_s"], P["face_s"]]), M)

    def gp(o, d):
        loc, ir, _ = patch.ray.intersects_location(o[None, :], d[None, :],
                                                   multiple_hits=False)
        return float(np.linalg.norm(loc[0] - o)) if len(ir) else np.inf

    def resid(t):
        F = np.zeros(3)
        for d, b in bases:
            dd = gp(b + t, d)
            if L_c <= dd <= L_c + S:
                F += F_PL * (-d)
            elif np.isfinite(dd) and dd < L_c:
                F += K_RIG * (L_c - dd) * (-d)
        for pts, k, cap in ((ring, K_SKIRT, SKIRT_MAX), (jac, K_JAC, None),
                            (rig, K_RIG, None)):
            q = pts + t
            pen = np.clip(-field.query(q), 0.0, None)
            if not pen.any():
                continue
            _, idx = field.tree.query(q[pen > 0])
            tot = ((k * pen[pen > 0])[:, None] * field.nrm[idx]).sum(axis=0)
            if cap is not None and np.linalg.norm(tot) > cap:
                tot = tot / np.linalg.norm(tot) * cap
            F += tot
        return F

    r = least_squares(resid, np.zeros(3), bounds=(-6, 6), xtol=1e-3, ftol=1e-3,
                      diff_step=0.2)
    out = []
    for d, b in bases:
        dd = gp(b + r.x, d)
        out.append(dict(D=dd, consumed=(dd - L_c) if np.isfinite(dd) else None,
                        engaged=bool(L_c <= dd <= L_c + S)))
    return r.x, out


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--json-dir", default=None)
    ap.add_argument("--a1", required=True)
    ap.add_argument("--a2", required=True)
    ap.add_argument("--ears", required=True, help="small=ID,median=ID,large=ID")
    ap.add_argument("--compacted", type=float, default=4.5)
    ap.add_argument("--stroke", type=float, default=6.0)
    a = ap.parse_args()
    aims = [np.array([float(x) for x in s.split(",")]) for s in (a.a1, a.a2)]
    aims = [v / np.linalg.norm(v) for v in aims]
    want = dict(kv.split("=") for kv in a.ears.split(","))
    cc, cr, off, GP = gen_bits()
    P = iem_points()
    os.makedirs(VIZ, exist_ok=True)
    sdir = os.path.join(HERE, "stl", "right")

    for label, eid in want.items():
        fs = glob.glob(os.path.join(a.json_dir or ALIGNED, f"*_{eid}_right.json"))
        if not fs:
            print(f"{label}: {eid} not found"); continue
        rec = json.load(open(fs[0]))
        patch = trimesh.load(os.path.join(ALIGNED, rec["patch"]), force="mesh")
        field = EarField(patch, seed=0)
        M = np.array(rec["transform"], float); R = M[:3, :3]
        bases = []
        for v in aims:
            d = R @ v; d /= np.linalg.norm(d)
            b = transform((cc + (cr ** 2 * v) / np.linalg.norm(cr * v)
                           + v * off)[None, :], M)[0]
            bases.append((d, b))
        t, sites = solve(rec, patch, field, P, bases, a.compacted, a.stroke)

        scene = trimesh.Scene()
        ear = cluster_decimate(patch, 55000)
        ear.visual.face_colors = COL["ear"]
        scene.add_geometry(ear, node_name="ear", geom_name="ear")
        T = np.eye(4); T[:3, 3] = t                     # equilibrium float
        for nm in ("core", "faceplate", "jacket_wing", "carrier",
                   "nozzle_insert_short"):
            m = trimesh.load(os.path.join(sdir, f"{nm}.stl"), force="mesh")
            m = cluster_decimate(m, 15000)
            if nm in ("carrier", "nozzle_insert_short"):
                m.apply_transform(NOZZLE_T)
            m.apply_transform(M); m.apply_transform(T)
            m.visual.face_colors = COL[nm]
            scene.add_geometry(m, node_name=nm, geom_name=nm)

        for k, ((d, b), s) in enumerate(zip(bases, sites)):
            nm = ["site1", "site2"][k]
            if not np.isfinite(s["D"]):
                continue
            L = max(s["D"] - 1.0, 0.5)                  # stack column, pad on top
            col = trimesh.creation.cylinder(radius=2.6, height=L)
            col.apply_transform(trimesh.geometry.align_vectors([0, 0, 1], d))
            col.apply_translation(b + t + d * (L / 2.0))
            col.visual.face_colors = COL["stack"]
            scene.add_geometry(col, node_name=f"stack_{nm}", geom_name=f"stack_{nm}")
            pad = trimesh.creation.cylinder(radius=3.0, height=1.0)
            pad.apply_transform(trimesh.geometry.align_vectors([0, 0, 1], d))
            pad.apply_translation(b + t + d * (s["D"] - 0.5))
            pad.visual.face_colors = COL["pad"]
            scene.add_geometry(pad, node_name=f"pad_{nm}", geom_name=f"pad_{nm}")

        out = os.path.join(VIZ, f"seated_{label}.glb")
        scene.export(out)
        mb = os.path.getsize(out) / 1e6

        # Stability at this equilibrium.  Float the POSE rather than the points:
        # M_f = T(t) @ M keeps skirt, jacket and pads consistently displaced,
        # which point-by-point translation of only the pads would not.
        M_f = M.copy(); M_f[:3, 3] = M[:3, 3] + t
        rec_f = dict(rec); rec_f["transform"] = M_f.tolist()
        Pl = dict(P)
        pads, meta = [], []
        for k, (v, s_) in enumerate(zip(aims, sites)):
            if not np.isfinite(s_["D"]):
                continue
            base_d = cc + (cr ** 2 * v) / np.linalg.norm(cr * v) + v * off
            tip_d = base_d + v * s_["D"]
            _, u, w = frame_for(v)
            pads.append(tip_d)
            pads.extend(tip_d + 2.7 * (np.cos(q) * u + np.sin(q) * w)
                        for q in np.linspace(0, 2 * np.pi, 8, endpoint=False))
            meta.append(dict(name=f"site{k+1}", aim=v.tolist()))
        cells = {}
        if pads:
            Pl["plunger"] = np.array(pads)
            Pl["_plunger"] = meta
            cp = transform(P["_cable_exit"][None, :], M_f)[0]
            com = transform(P["_com"][None, :], M_f)[0]
            for mu in (0.4, 0.6, 0.8, 1.0):
                for tug in (0.5, 0.2):
                    r = stab.stability_check(rec_f, Pl, field, transform,
                                             cable_point=cp, com=com,
                                             mu=mu, cable_tug=tug)
                    cells[f"mu{mu}_tug{tug}"] = round(r["margin"], 3)
        mj = dict(label=label, ear=f"{rec['dataset']}/{rec['ear_id']}",
                  aim_pair=[aims[0].tolist(), aims[1].tolist()],
                  equilibrium_translation_mm=t.tolist(),
                  sites=[dict(aim=aims[i].tolist(), surface_distance_mm=sites[i]["D"],
                              stroke_consumed_mm=sites[i]["consumed"],
                              engaged=sites[i]["engaged"]) for i in range(2)],
                  stability_margin=cells,
                  note=("Body at the coupled floating-equilibrium pose. Plunger "
                        "stacks are SCHEMATIC primitives (green column, red pad) "
                        "drawn at the extension consumed on this ear -- the "
                        "searched aims have no STL yet. Ear/shell meshes are "
                        "decimated for viewing only."))
        json.dump(mj, open(os.path.join(VIZ, f"seated_{label}_meta.json"), "w"),
                  indent=1)
        print(f"{label:<7} {rec['ear_id']:<8} shift {np.linalg.norm(t):.2f} mm  "
              f"D {sites[0]['D']:.2f}/{sites[1]['D']:.2f}  "
              f"eng {sites[0]['engaged']}/{sites[1]['engaged']}  "
              f"{mb:.2f} MB  mu0.6/0.2N {cells.get('mu0.6_tug0.2','--')}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
