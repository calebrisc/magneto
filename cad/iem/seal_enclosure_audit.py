#!/usr/bin/env python3
"""
seal_enclosure_audit.py -- does the sealed rim actually ENCLOSE the canal?

Raised by human review of viz/seated_scene.glb: the nozzle looked like it was
sealing against flesh *near* the canal aperture rather than *around* it.  If so,
that is a hole in the scoring, not a cosmetic issue -- `seal_compliance.py` asks
only whether the rim forms a continuous contact loop against flesh.  A rim laid
flat on the concha floor, nowhere near the canal, satisfies that perfectly while
sealing nothing at all.  An acoustic seal requires the contact loop to *surround*
the canal entrance.

Three checks per ear, all in the ear's scan frame:

  (a) AIM       angle between the nozzle axis and the direction from the rim
                centre to the detected canal aperture.  0 deg = the bore points
                straight at the canal.  Also reported against the inward concha
                normal for context (anatomy puts the canal axis 30-60 deg off it,
                so that one is not expected to be zero).

  (b) OFFSET    in-plane distance from the rim-circle centre to the aperture,
                measured in the rim plane.  This is the crisp "is the bore over
                the canal" number: compare it with the rim radius.

  (c) ENCLOSURE the aperture must project INSIDE the rim ring *and* the contact
                loop must be continuous around it.  Only then does the sealed
                loop surround the canal entrance.  Reported at the conservative
                2.5 mm budget, with the same travel search seal_compliance uses.

`axial` is the aperture's position along the nozzle axis relative to the rim
plane: positive means the canal mouth lies ahead of the rim, which is what a
bore opening into the canal looks like.  Strongly negative means the rim has
been driven past the aperture into flesh.

Usage:
    python seal_enclosure_audit.py --json-dir /path/to/seatings
    python seal_enclosure_audit.py --ear sonicom/P0023
"""

from __future__ import annotations

import argparse
import glob
import json
import os
import sys

import numpy as np
import trimesh

import seal_compliance as sc
from earfit import ALIGNED, NOZZLE_AXIS, EarField, transform

BUDGET = 2.5


def audit_ear(rec, patch, rim_design):
    field = EarField(patch, seed=0)
    M = np.array(rec["transform"], float)
    A = np.array(rec["aperture"], float)
    n_out = np.array(rec["floor_normal"], float)

    n_ax = M[:3, :3] @ NOZZLE_AXIS
    n_ax /= np.linalg.norm(n_ax)
    rim = transform(rim_design, M)
    C = rim.mean(axis=0)
    R = float(np.median(np.linalg.norm(rim - C, axis=1)))

    d = A - C
    axial = float(d @ n_ax)
    lateral = float(np.linalg.norm(d - axial * n_ax))
    nd = np.linalg.norm(d)
    aim = float(np.degrees(np.arccos(np.clip((d / nd) @ n_ax, -1, 1)))) if nd > 1e-9 else 0.0
    aim_norm = float(np.degrees(np.arccos(np.clip(n_ax @ (-n_out), -1, 1))))

    # continuity of the contact loop at the conservative budget
    best = None
    for t in np.linspace(0.0, sc.TRAVEL_MM, sc.TRAVEL_STEPS):
        dist = field.query(rim + t * n_ax)
        sealed = dist <= BUDGET
        glen, _ = sc.longest_false_run(sealed)
        gap = 360.0 * glen / len(sealed)
        key = (float(sealed.mean()), -gap)
        if best is None or key > best[0]:
            best = (key, float(sealed.mean()), gap)
    _, cover, gap = best
    continuous = (cover >= sc.COVER_MIN) and (gap <= sc.GAP_MAX_DEG)

    inside = lateral < R
    return dict(dataset=rec["dataset"], ear_id=rec["ear_id"],
                aim=aim, aim_norm=aim_norm, lateral=lateral, axial=axial,
                rim_r=R, inside=bool(inside), continuous=bool(continuous),
                encloses=bool(inside and continuous), cover=cover, gap=gap)


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--json-dir", default=None)
    ap.add_argument("--stl-dir", default=None)
    ap.add_argument("--ear", default=None, help="highlight one ear, e.g. sonicom/P0023")
    a = ap.parse_args()

    js = sorted(glob.glob(os.path.join(a.json_dir or ALIGNED, "*.json")))
    if not js:
        sys.exit("no seatings found")
    rim_design = sc.rim_from_mesh(a.stl_dir)

    rows = []
    for p in js:
        rec = json.load(open(p))
        patch = trimesh.load(os.path.join(ALIGNED, rec["patch"]), force="mesh")
        rows.append(audit_ear(rec, patch, rim_design))

    n = len(rows)
    arr = lambda k: np.array([r[k] for r in rows], float)   # noqa: E731
    print(f"\nn = {n} ears; rim radius {rows[0]['rim_r']:.2f} mm\n")

    print("| check | median | p10 | p90 | max |")
    print("|---|---|---|---|---|")
    for k, lab in (("aim", "(a) aim: nozzle axis vs rim-centre→aperture (deg)"),
                   ("aim_norm", "     nozzle axis vs inward concha normal (deg)"),
                   ("lateral", "(b) lateral offset rim centre→aperture (mm)"),
                   ("axial", "     aperture axial position vs rim plane (mm)")):
        v = arr(k)
        print(f"| {lab} | {np.median(v):.1f} | {np.percentile(v,10):.1f} | "
              f"{np.percentile(v,90):.1f} | {v.max():.1f} |")

    ins = sum(r["inside"] for r in rows)
    con = sum(r["continuous"] for r in rows)
    enc = sum(r["encloses"] for r in rows)
    print(f"\n(c) ENCLOSURE at the {BUDGET} mm budget")
    print(f"    aperture projects inside the rim ring : {ins}/{n} ({100*ins/n:.0f} %)")
    print(f"    contact loop continuous               : {con}/{n} ({100*con/n:.0f} %)")
    print(f"    BOTH (loop actually surrounds canal)  : {enc}/{n} ({100*enc/n:.0f} %)")
    print(f"    sealed-but-NOT-enclosing (false pass) : {con - enc}/{n}")

    if a.ear:
        ds, eid = a.ear.split("/")
        r = next((x for x in rows if x["dataset"] == ds and x["ear_id"] == eid), None)
        if r:
            print(f"\n--- {a.ear} ---")
            for k in ("aim", "aim_norm", "lateral", "axial", "rim_r", "cover",
                      "gap", "inside", "continuous", "encloses"):
                print(f"    {k:<11} {r[k]}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
