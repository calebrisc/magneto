#!/usr/bin/env python3
"""
hook_config.py -- lip-hook retention as the PRIMARY design, not an upper bound.

Reframe: the magnet stacks are DEPLOYMENT (they carry the pad out to the lip);
retention is interlock plus friction.  Each pad is a hook whose contact patch
seats on the UNDERSIDE of a cartilage overhang -- the cymba lip and the antihelix
undercut this report already located.  Once hooked, the contact normal is set by
the cartilage, not by the pad's aim, so the reaction-direction impossibility that
killed the preload approach does not apply.

MODEL
  1 cast from the jacket outer face along the site aim to find the lip;
  2 collect ear samples near that hit whose normal OPPOSES pull-out
    (n . pull < UNDERCUT) -- that is the lip's underside;
  3 the hook's contact patch is placed on those samples, so it inherits their
    interlocking normals;
  4 the stack must be able to deliver the pad there: the along-aim distance to
    the patch must lie inside [compacted, compacted + stroke];
  5 usable WRAP DEPTH is the extent of that undercut patch along the pull-out
    direction -- how far under the lip a pad could actually reach.  A hook that
    assumes 1.5-2 mm of wrap needs that much lip to wrap around.

Stability is then run over ALL pull directions including a straight outward yank
(`cable_mode="sphere"`), which is the case a hook must survive and a
friction-only fit cannot.

Usage:
    python hook_config.py --json-dir DIR --ears small=ID,median=ID,large=ID
"""

from __future__ import annotations

import argparse
import glob
import json
import os
import sys

import numpy as np
import trimesh

import stability as stab
from earfit import ALIGNED, HERE, NOZZLE_AXIS, NOZZLE_T, EarField, iem_points, transform
from viz_scene import cluster_decimate

VIZ = os.path.join(HERE, "viz")
UNDERCUT = -0.20
WRAP_WANT = (1.5, 2.0)
SEARCH_R = 9.0
MIN_PTS = 4
COL = {"ear": (205, 175, 135, 255), "core": (190, 192, 200, 255),
       "faceplate": (205, 207, 214, 255), "jacket_wing": (62, 108, 198, 255),
       "carrier": (232, 140, 48, 255), "nozzle_insert_short": (128, 128, 130, 255),
       "stack": (110, 200, 110, 255), "hook": (220, 90, 90, 255)}


def gen_bits():
    import generate
    P = generate.PARAMS
    g = generate.G(P)
    return (np.array(g.core_c, float), np.array(g.core_r, float),
            P["clearance"] + P["jacket_thick"], P)


def hook_at(rec, patch, field, aim_d, cc, cr, off, L_c, S):
    """Find the lip underside this aim can reach, and whether a stack reaches it."""
    M = np.array(rec["transform"], float); R = M[:3, :3]
    pull = -(R @ NOZZLE_AXIS); pull /= np.linalg.norm(pull)
    base_d = cc + (cr ** 2 * aim_d) / np.linalg.norm(cr * aim_d) + aim_d * off
    b = transform(base_d[None, :], M)[0]
    d = R @ aim_d; d /= np.linalg.norm(d)
    loc, ir, _ = patch.ray.intersects_location(b[None, :], d[None, :],
                                               multiple_hits=False)
    res = dict(reach=np.inf, seated=False, n_pts=0, wrap=0.0, patch=None,
               hit=None)
    if not len(ir):
        return res
    H = loc[0]
    res["hit"] = H
    near = field.tree.query_ball_point(H, SEARCH_R)
    if not near:
        return res
    k = field.nrm[near] @ pull
    sel = [near[i] for i in np.where(k < UNDERCUT)[0]]
    if len(sel) < MIN_PTS:
        return res
    pts = field.pts[sel]
    # keep the undercut cluster closest to where the aim arrives
    c0 = pts[np.argmin(np.linalg.norm(pts - H, axis=1))]
    keep = sel_i = [s for s, p in zip(sel, pts)
                    if np.linalg.norm(p - c0) <= 4.0]
    if len(keep) < MIN_PTS:
        return res
    pp = field.pts[keep]
    C = pp.mean(axis=0)
    res["reach"] = float((C - b) @ d)
    res["seated"] = bool(L_c <= res["reach"] <= L_c + S)
    res["n_pts"] = len(keep)
    # usable wrap = extent of the undercut patch along the pull-out direction
    proj = pp @ pull
    res["wrap"] = float(proj.max() - proj.min())
    idx = np.argsort(np.linalg.norm(pp - C, axis=1))[:9]
    res["patch"] = pp[idx]
    res["patch_design"] = (pp[idx] - M[:3, 3]) @ R
    return res


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--json-dir", default=None)
    ap.add_argument("--ears", required=True)
    ap.add_argument("--compacted", type=float, default=4.5)
    ap.add_argument("--stroke", type=float, default=6.0)
    ap.add_argument("--aims", default=None, help="two 'x,y,z' aims, ';'-separated")
    a = ap.parse_args()
    cc, cr, off, GP = gen_bits()
    P = iem_points()
    if a.aims:
        aims = [np.array([float(x) for x in s.split(",")]) for s in a.aims.split(";")]
        names = ["site1", "site2"]
    else:
        import generate
        g = generate.G(generate.PARAMS)
        aims = [np.asarray(pl["aim"], float) for pl in g.plungers]
        names = [pl["name"] for pl in g.plungers]
    aims = [v / np.linalg.norm(v) for v in aims]
    want = dict(kv.split("=") for kv in a.ears.split(","))
    js = sorted(glob.glob(os.path.join(a.json_dir or ALIGNED, "*.json")))
    os.makedirs(VIZ, exist_ok=True)
    sdir = os.path.join(HERE, "stl", "right")

    print(f"aims: " + "  ".join(f"{n} {np.round(v,2)}" for n, v in zip(names, aims)))
    print(f"stack window [{a.compacted:.2f}, {a.compacted + a.stroke:.2f}] mm; "
          f"hook wants {WRAP_WANT[0]}-{WRAP_WANT[1]} mm of wrap\n")
    print(f"{'ear':<12}" + "".join(f"{n[:11]+' reach':>18}" for n in names)
          + f"{'both seated':>13}{'min wrap':>10}")
    rows = []
    for p in js:
        rec = json.load(open(p))
        if not isinstance(rec, dict) or "patch" not in rec:
            continue                      # sidecar JSONs live in the same dir
        patch = trimesh.load(os.path.join(ALIGNED, rec["patch"]), force="mesh")
        field = EarField(patch, seed=0)
        hs = [hook_at(rec, patch, field, v, cc, cr, off, a.compacted, a.stroke)
              for v in aims]
        both = all(h["seated"] for h in hs)
        wraps = [h["wrap"] for h in hs if h["patch"] is not None]
        rows.append((rec, patch, field, hs, both))
        cells = "".join(
            (f"{'--':>18}" if not np.isfinite(h['reach'])
             else f"{h['reach']:>13.2f}{'*' if h['seated'] else ' ':>2}   ")
            for h in hs)
        print(f"{rec['ear_id']:<12}{cells}{str(both):>13}"
              f"{(min(wraps) if wraps else 0):>10.2f}")
    nb = sum(1 for r in rows if r[4])
    print(f"\nboth hooks seat within stroke on {nb}/{len(rows)} ears")
    allw = [h["wrap"] for _r, _p, _f, hs, _b in rows for h in hs
            if h["patch"] is not None]
    if allw:
        allw = np.array(allw)
        print(f"usable wrap depth: min {allw.min():.2f}  median {np.median(allw):.2f}"
              f"  p90 {np.percentile(allw,90):.2f} mm  "
              f"(hook needs {WRAP_WANT[0]}-{WRAP_WANT[1]})")
        print(f"  sites with < {WRAP_WANT[0]} mm of lip to wrap: "
              f"{int((allw < WRAP_WANT[0]).sum())}/{len(allw)}")

    # ---- stability + scenes for the three named ears ----------------------- #
    print("\n=== stability, ALL pull directions incl. straight outward yank ===")
    for label, eid in want.items():
        m = [r for r in rows if r[0]["ear_id"] == eid]
        if not m:
            print(f"{label}: {eid} not in the set"); continue
        rec, patch, field, hs, both = m[0]
        M = np.array(rec["transform"], float)
        Pl = dict(P)
        pd = [h["patch_design"] for h in hs if h["patch"] is not None]
        Pl["plunger"] = np.vstack(pd) if pd else np.zeros((0, 3))
        Pl["_plunger"] = [dict(name=n) for n, h in zip(names, hs)
                          if h["patch"] is not None]
        cp = transform(P["_cable_exit"][None, :], M)[0]
        com = transform(P["_com"][None, :], M)[0]
        cells = {}
        print(f"\n  {label.upper():<7} {eid}  both seated: {both}  "
              f"wrap {min([h['wrap'] for h in hs if h['patch'] is not None] or [0]):.2f} mm")
        print("  | mu | 0.5 N | 0.2 N |"); print("  |---|---|---|")
        for mu in (0.4, 0.6, 0.8):
            r = []
            for tug in (0.5, 0.2):
                x = stab.stability_check(rec, Pl, field, transform,
                                         cable_point=cp, com=com, mu=mu,
                                         cable_tug=tug, cable_mode="sphere")
                cells[f"mu{mu}_tug{tug}"] = round(x["margin"], 3)
                r.append(f"{x['margin']:.2f}x{'P' if x['margin'] >= 1 else 'f'}")
            print(f"  | {mu:.1f} | {r[0]} | {r[1]} |")

        scene = trimesh.Scene()
        ear = cluster_decimate(patch, 55000); ear.visual.face_colors = COL["ear"]
        scene.add_geometry(ear, node_name="ear", geom_name="ear")
        for nm in ("core", "faceplate", "jacket_wing", "carrier",
                   "nozzle_insert_short"):
            mm = trimesh.load(os.path.join(sdir, f"{nm}.stl"), force="mesh")
            mm = cluster_decimate(mm, 14000)
            if nm in ("carrier", "nozzle_insert_short"):
                mm.apply_transform(NOZZLE_T)
            mm.apply_transform(M); mm.visual.face_colors = COL[nm]
            scene.add_geometry(mm, node_name=nm, geom_name=nm)
        R = M[:3, :3]
        for nmn, v, h in zip(names, aims, hs):
            if h["patch"] is None or not np.isfinite(h["reach"]):
                continue
            d = R @ v; d /= np.linalg.norm(d)
            b = transform((cc + (cr ** 2 * v) / np.linalg.norm(cr * v)
                           + v * off)[None, :], M)[0]
            L = max(h["reach"] - 1.2, 0.5)
            col = trimesh.creation.cylinder(radius=2.4, height=L)
            col.apply_transform(trimesh.geometry.align_vectors([0, 0, 1], d))
            col.apply_translation(b + d * (L / 2.0))
            col.visual.face_colors = COL["stack"]
            scene.add_geometry(col, node_name=f"stack_{nmn}", geom_name=f"stack_{nmn}")
            # hook sketched as a wrap-thick pad on the lip underside
            hook = trimesh.creation.box(extents=(5.0, 5.0, max(h["wrap"], 1.0)))
            pull = -(R @ NOZZLE_AXIS); pull /= np.linalg.norm(pull)
            hook.apply_transform(trimesh.geometry.align_vectors([0, 0, 1], pull))
            hook.apply_translation(h["patch"].mean(axis=0))
            hook.visual.face_colors = COL["hook"]
            scene.add_geometry(hook, node_name=f"hook_{nmn}", geom_name=f"hook_{nmn}")
        out = os.path.join(VIZ, f"hook_{label}.glb")
        scene.export(out)
        meta = dict(label=label, ear=f"{rec['dataset']}/{rec['ear_id']}",
                    aims={n: v.tolist() for n, v in zip(names, aims)},
                    both_hooks_seated=bool(both),
                    sites=[dict(name=n, reach_mm=(None if not np.isfinite(h["reach"])
                                                  else h["reach"]),
                                seated=h["seated"], undercut_pts=h["n_pts"],
                                usable_wrap_mm=h["wrap"])
                           for n, h in zip(names, hs)],
                    wrap_required_mm=list(WRAP_WANT),
                    stability_margin_all_directions=cells,
                    note=("Hook-seated config: pad contact patches placed on the "
                          "UNDERSIDE of the local overhang, so normals are the "
                          "cartilage's. Stacks and hooks are schematic primitives. "
                          "Stability swept over all pull directions including a "
                          "straight outward yank."))
        json.dump(meta, open(os.path.join(VIZ, f"hook_{label}_meta.json"), "w"),
                  indent=1)
        print(f"  -> {out}  {os.path.getsize(out)/1e6:.2f} MB")
    return 0


if __name__ == "__main__":
    sys.exit(main())
