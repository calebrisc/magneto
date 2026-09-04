#!/usr/bin/env python3
"""
crus_bands.py -- the 9/1 cymba-fold clamp, sized S/M/L by fit.

The clamp: body in the cavum under the crus helicis, one magnetic pad expanding
into the cymba so the crus ridge is pinched between pad and body.  The pad's
contact is the crus overhang -- the surface whose normal opposes pull-out --
so the retention is interlock, and the sizing variable is how far the pad has
to travel from the jacket face to that overhang.

Per ear (only ears where gyro_arm_variance found a crus, n~41):
  crus_w   the crus point (datum frame -> scan frame)
  patch    undercut samples within R of the crus point, n.pull < UNDERCUT
  travel   distance from the jacket outer face (at the cymba aim) to the patch
           centroid, projected on the aim -- the pad extension length
  stability  pad on the patch, skirt + jacket as usual, all pull directions,
           mu 0.6, tug 0.5 / 0.2 N

Bands: 3 contiguous bands over travel minimising the widest one, then the
stack must fit: (n_gap+1) x 1 mm magnets + 1 mm pad + 0.8 mm plate <=
compacted (+ recess, if the pad pocket is sunk into the jacket).

Usage:  python crus_bands.py            (runs everything, ~2 min)
"""
from __future__ import annotations
import glob, json, os, sys, time
import numpy as np, trimesh
import stability as stab
from earfit import ALIGNED, HERE, NOZZLE_AXIS, EarField, iem_points, transform
from hook_config import gen_bits
from gyro_arm_variance import datum
from size_bands import ears, bands_from, sizing, AIM, MU, TUGS, MARGIN, GAP_STROKE, MAG_T, PAD_T, PLATE_T

R_SEARCH = 6.0; UNDERCUT = -0.20; MIN_PTS = 4
CACHE = os.path.join(ALIGNED, "crus_bands.json")
OUT_MD = os.path.join(HERE, "..", "..", "docs", "CLAMP_SIZE_BANDS.md")


def run():
    cc, cr, off, _ = gen_bits(); P = iem_points()
    gav = {e["id"]: e for e in json.load(open(os.path.join(ALIGNED, "gyro_arm_variance.json")))}
    rows = []; t0 = time.time()
    for rec in ears():
        g = gav.get(rec["ear_id"]); 
        if not g or g.get("crus") is None: continue
        patch = trimesh.load(os.path.join(ALIGNED, rec["patch"]), force="mesh")
        field = EarField(patch, seed=0)
        M = np.array(rec["transform"], float); Rm = M[:3, :3]
        ap, B = datum(rec)
        crus_w = ap + B @ np.array(g["crus"], float)
        pull = -(Rm @ NOZZLE_AXIS); pull /= np.linalg.norm(pull)
        near = field.tree.query_ball_point(crus_w, R_SEARCH)
        row = dict(ear=rec["ear_id"], travel=None, wrap=None, n_pts=0)
        if near:
            k = field.nrm[near] @ pull
            sel = [near[i] for i in np.where(k < UNDERCUT)[0]]
            if len(sel) >= MIN_PTS:
                pp = field.pts[sel]; C = pp.mean(axis=0)
                # travel is measured toward the crus patch itself, not along the
                # cymba aim: direction from the core centre to the patch centroid
                # (design frame), base point on the ellipsoid surface + jacket.
                C_d = (C - M[:3, 3]) @ Rm
                a_d = C_d - cc; a_d /= np.linalg.norm(a_d)
                base_d = cc + (cr ** 2 * a_d) / np.linalg.norm(cr * a_d) + a_d * off
                b = transform(base_d[None, :], M)[0]; d = Rm @ a_d; d /= np.linalg.norm(d)
                row["aim_design"] = a_d.tolist()
                row["travel"] = float((C - b) @ d)
                row["lateral"] = float(np.linalg.norm((C - b) - ((C - b) @ d) * d))
                proj = pp @ pull; row["wrap"] = float(proj.max() - proj.min()); row["n_pts"] = len(sel)
                idx = np.argsort(np.linalg.norm(pp - C, axis=1))[:9]
                Pl = dict(P); Pl["plunger"] = (pp[idx] - M[:3, 3]) @ Rm; Pl["_plunger"] = [dict(name="crus")]
                cp = transform(P["_cable_exit"][None, :], M)[0]; com = transform(P["_com"][None, :], M)[0]
                for tug in TUGS:
                    x = stab.stability_check(rec, Pl, field, transform, cable_point=cp, com=com,
                                             mu=MU, cable_tug=tug, cable_mode="sphere")
                    row[f"tug{tug}"] = round(float(x["margin"]), 3)
                    row[f"interlock{tug}"] = int(x.get("interlock", 0))
                    row[f"pullout{tug}"] = round(float(x.get("pullout_capacity", 0)), 3)
        rows.append(row)
        print(f"{rec['ear_id']:<10} travel {row['travel'] if row['travel'] is None else round(row['travel'],2)!s:<8} "
              f"lat {row.get('lateral', float('nan')):.1f} wrap {row['wrap'] if row['wrap'] is None else round(row['wrap'],2)!s:<6} "
              f"{row.get('tug0.5','-')}x {row.get('tug0.2','-')}x  [{time.time()-t0:.0f}s]", flush=True)
    json.dump(rows, open(CACHE, "w"), indent=1)
    return rows


def report(rows):
    ok = [r for r in rows if r["travel"] is not None]
    tr = np.array([r["travel"] for r in ok])
    L = ["\n\n# Crus-pinch clamp (the 9/1 design), sized S/M/L by fit\n",
         f"{len(rows)} ears with a detected crus (of 102); undercut found at the crus on {len(ok)}. "
         "`travel` = distance from the jacket face (on the line core-centre -> patch) to the crus overhang patch. Script: `cad/iem/crus_bands.py`.\n",
         f"travel: min {tr.min():.2f}  p10 {np.percentile(tr,10):.2f}  median {np.median(tr):.2f}  "
         f"p90 {np.percentile(tr,90):.2f}  max {tr.max():.2f} mm (span {tr.max()-tr.min():.2f})\n",
         "| size | band (travel, mm) | compacted | stroke | gaps | stack length | fits (recess 0 / 3 mm) | ears | pass @0.5 N | pass @0.2 N | median margin 0.5 N |",
         "|---|---|---|---|---|---|---|---|---|---|---|"]
    b = bands_from(tr)
    for name, (lo, hi) in zip("SML", b):
        sz = sizing(lo, hi); mem = [r for r in ok if lo - 1e-9 <= r["travel"] <= hi + 1e-9]
        a = np.array([r["tug0.5"] for r in mem]); c = np.array([r["tug0.2"] for r in mem])
        L.append(f"| {name} | {lo:.2f}–{hi:.2f} | {sz['compacted']:.2f} | {sz['stroke']:.2f} | {sz['n_gap']} | {sz['stack_len']:.2f} | "
                 f"{'yes' if sz['fits'] else 'no'} / {'yes' if sz['stack_len'] <= sz['compacted'] + 3 else 'no'} | {len(mem)} | "
                 f"{(a>=1).sum()} ({100*(a>=1).mean():.0f}%) | {(c>=1).sum()} ({100*(c>=1).mean():.0f}%) | {np.median(a):.2f}x |")
    a = np.array([r["tug0.5"] for r in ok]); c = np.array([r["tug0.2"] for r in ok])
    il = np.array([r["interlock0.5"] for r in ok]); po = np.array([r["pullout0.5"] for r in ok])
    L.append(f"\nAll {len(ok)} ears, pad seated by construction: pass {(a>=1).sum()} @0.5 N, {(c>=1).sum()} @0.2 N; "
             f"median straight pull-out capacity {np.median(po):.2f} N vs demand ~1.05 N (0.31 skirt + 0.5 tug + 0.24 inertial); "
             f"median interlocking contacts {np.median(il):.0f}.\n")
    L.append("Universal-vs-banded is not the question here: with the pad placed on the overhang by construction every ear is "
             "'seated', so the pass rate is the ceiling for a single crus pad. Sizing only decides whether a physical stack "
             "can deliver the pad to that overhang (the fits column).\n")
    md = "\n".join(L); open(OUT_MD, "a").write(md); print(md)


if __name__ == "__main__":
    rows = run() if not os.path.exists(CACHE) or "--rerun" in sys.argv else json.load(open(CACHE))
    report(rows)
