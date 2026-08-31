#!/usr/bin/env python3
"""
retention_analysis.py -- why does the wing fail to retain, and what would fix it?

Retention became the #2 failure axis at v3 and has drifted the wrong way at every
revision as the body moved under it.  `tryon.py` only reports the number
(`wing_tip`, a median signed distance); this asks *why* each failing ear fails,
and sweeps the two geometric knobs to see what would fix the most ears.

THE PHYSICS.  The wing is a compliant leaf spring at ~0.25 N/mm, so it is meant
to be a light interference fit: contact with roughly 0 to 1.5 mm of overlap is
absorbed by deflection and is exactly what retains.  A positive gap is no
retention at all; more than ~1.5 mm of overlap means the spring is bottoming and
levering the shell back out.  So the target band used here is

    tip signed distance in [-1.5, 0.0] mm      ("retained")

which is the spring's working range, and is narrower than tryon.py's grading band
because it asks for real retention rather than a passing grade.

TAXONOMY.  Each failing ear is put in exactly one bucket, in priority order:

  BLOCKED     the wing is already overlapping >= 1.5 mm somewhere along its span
              (mid or tip) -- it bottoms out before it can seat.  Making the wing
              longer makes this worse.
  MISDIRECTED the wing has a gap at the tip, but its reach direction is nearly
              parallel to the ear surface it is supposed to press (angle between
              the wing's +Y growth axis and the local inward surface normal
              >= 60 deg).  Extending it slides it along the antihelix instead of
              into it, so length will not fix these.
  SHORT       the wing has a gap and points at the surface (angle < 60 deg).
              Pure reach problem; extra length converts directly into press.
  OVERPRESSED overlap beyond the spring's range with nothing blocking -- the wing
              is simply too long for this ear.

THE SWEEP.  At each ear's frozen seated pose, the tip sample points are extended
by L mm along the design +Y (the wing's growth direction) and/or the whole wing
is splayed by theta degrees about the design +X axis through the wing root, and
the tip distance re-queried.  Frozen pose, so this is a lower bound: a wing that
retains better would also let the shell sit slightly differently.

Usage:
    python retention_analysis.py
    python retention_analysis.py --json-dir /path/to/saved/seatings
"""

from __future__ import annotations

import argparse
import glob
import json
import os
import sys

import numpy as np
import trimesh

from earfit import ALIGNED, HERE, EarField, iem_points, transform

RETAIN_LO, RETAIN_HI = -1.5, 0.0      # spring working range
BLOCK_OVERLAP = -1.5                  # beyond this the spring bottoms out
MISDIRECT_DEG = 60.0


def wing_root(stl_dir=None):
    """Centroid of the wing's root band (the jacket rim it cantilevers from)."""
    sdir = stl_dir or os.path.join(HERE, "stl", "right")
    jw = trimesh.load(os.path.join(sdir, "jacket_wing.stl"), force="mesh")
    v = np.asarray(jw.vertices)
    band = v[(v[:, 1] > 7.5) & (v[:, 1] < 8.5)]
    return band.mean(axis=0) if len(band) else np.array([0.0, 8.0, -3.0])


def morph(pts, root, ext, theta_deg):
    """Extend the wing by `ext` mm along +Y and splay it `theta_deg` about +X
    through the root."""
    q = pts.copy()
    if theta_deg:
        t = np.radians(theta_deg)
        c, s = np.cos(t), np.sin(t)
        d = q - root
        y, z = d[:, 1].copy(), d[:, 2].copy()
        d[:, 1], d[:, 2] = c * y - s * z, s * y + c * z
        q = root + d
    if ext:
        # extend along the (possibly splayed) growth direction
        g = np.array([0.0, 1.0, 0.0])
        if theta_deg:
            t = np.radians(theta_deg)
            g = np.array([0.0, np.cos(t), np.sin(t)])
        q = q + ext * g
    return q


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--json-dir", default=None)
    ap.add_argument("--tryon-csv", default=None,
                    help="tryon.csv, to split out the ears that FAIL the graded "
                         "retention axis from the wider spring-band population")
    ap.add_argument("--stl-dir", default=None)
    a = ap.parse_args()

    jdir = a.json_dir or ALIGNED
    js = sorted(glob.glob(os.path.join(jdir, "*.json")))
    if not js:
        sys.exit("no seatings found")

    graded = {}
    if a.tryon_csv:
        import csv as _csv
        for r in _csv.DictReader(open(a.tryon_csv)):
            graded[(r["dataset"], r["ear_id"])] = r["g_retention"]

    P = iem_points(stl_dir=a.stl_dir)
    root = wing_root(a.stl_dir)
    tip, mid = P["wing_tip"], P["wing_mid"]

    EXT = [-3.0, -2.5, -2.0, -1.5, -1.0, -0.5, 0.0, 0.5, 1.0, 1.5]
    THETA = [-15.0, -10.0, -5.0, 0.0, 5.0, 10.0, 15.0]

    rows = []
    grid = np.zeros((len(EXT), len(THETA)), int)
    for p in js:
        rec = json.load(open(p))
        patch = trimesh.load(os.path.join(ALIGNED, rec["patch"]), force="mesh")
        field = EarField(patch, seed=0)
        M = np.array(rec["transform"], float)

        d_tip = float(np.median(field.query(transform(tip, M))))
        d_mid = float(np.median(field.query(transform(mid, M))))
        worst = min(float(field.query(transform(tip, M)).min()),
                    float(field.query(transform(mid, M)).min()))

        # local surface normal where the tip is aiming
        tw = transform(tip, M)
        _, idx = field.tree.query(tw.mean(axis=0)[None, :])
        n_ear = field.nrm[int(idx[0])]
        growth = M[:3, :3] @ np.array([0.0, 1.0, 0.0])
        growth /= np.linalg.norm(growth)
        # angle between the wing's reach direction and pressing straight in
        ang = float(np.degrees(np.arccos(np.clip(growth @ (-n_ear), -1.0, 1.0))))

        if worst <= BLOCK_OVERLAP and d_tip > RETAIN_HI:
            bucket = "blocked"
        elif d_tip > RETAIN_HI:
            bucket = "misdirected" if ang >= MISDIRECT_DEG else "short"
        elif d_tip < RETAIN_LO:
            bucket = "overpressed"
        else:
            bucket = "retained"

        rows.append(dict(dataset=rec["dataset"], ear_id=rec["ear_id"],
                         tip=d_tip, mid=d_mid, worst=worst, angle=ang,
                         bucket=bucket,
                         graded=graded.get((rec["dataset"], rec["ear_id"]), "?")))

        for i, e in enumerate(EXT):
            for j, th in enumerate(THETA):
                q = transform(morph(tip, root, e, th), M)
                dv = float(np.median(field.query(q)))
                if RETAIN_LO <= dv <= RETAIN_HI:
                    grid[i, j] += 1

    n = len(rows)
    print(f"\nn = {n} ears, frozen at their seated poses\n")

    print("### Retention taxonomy (spring band: tip in [-1.5, 0.0] mm)\n")
    print("| outcome | ears | share |")
    print("|---|---|---|")
    for b in ("retained", "short", "misdirected", "blocked", "overpressed"):
        k = sum(1 for r in rows if r["bucket"] == b)
        print(f"| {b} | {k} | {100*k/n:.0f} % |")

    for b in ("short", "misdirected", "blocked", "overpressed"):
        sub = [r for r in rows if r["bucket"] == b]
        if not sub:
            continue
        g = np.array([r["tip"] for r in sub])
        an = np.array([r["angle"] for r in sub])
        print(f"\n{b}: n={len(sub)}  tip median {np.median(g):+.2f} mm "
              f"(p90 {np.percentile(g,90):+.2f})  aim angle median {np.median(an):.0f} deg")

    fails = [r for r in rows if r["graded"] == "fail"]
    if fails:
        print(f"\n### The {len(fails)} ears that FAIL the graded retention axis\n")
        print("| ear | tip mm | worst overlap mm | aim angle | diagnosis |")
        print("|---|---|---|---|---|")
        for r in sorted(fails, key=lambda r: r["tip"]):
            print(f"| {r['dataset']} {r['ear_id']} | {r['tip']:+.2f} | {r['worst']:+.2f} "
                  f"| {r['angle']:.0f}° | {r['bucket']} |")

    print("\n### Ears retained after a wing change (extension x splay)\n")
    print("| ext \\\\ splay | " + " | ".join(f"{t:+.0f}°" for t in THETA) + " |")
    print("|" + "---|" * (len(THETA) + 1))
    for i, e in enumerate(EXT):
        print(f"| **{e:+.1f} mm** | " + " | ".join(str(grid[i, j]) for j in range(len(THETA))) + " |")

    i, j = np.unravel_index(int(np.argmax(grid)), grid.shape)
    print(f"\nbest: {EXT[i]:+.1f} mm extension, {THETA[j]:+.0f}° splay -> "
          f"{grid[i,j]}/{n} retained (from {grid[EXT.index(0.0), THETA.index(0.0)]}/{n} as built)")
    return 0


if __name__ == "__main__":
    sys.exit(main())
