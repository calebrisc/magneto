#!/usr/bin/env python3
"""
seal_compliance.py -- compliance-aware seal scoring on already-seated ears.

WHY THIS EXISTS
---------------
`tryon.py`'s seal axis treats the skirt as RIGID: it scores the fraction of the
Ø19 rim circle lying within 1.5 mm of flesh, and calls anything else a gap.  The
physical skirt is nothing like that.  It is a 0.35 mm Shore-A-10/15 silicone
flare that folds and drapes over 2-4 mm, and it is carried on a mag-float
carrier with **1.5 mm of axial travel** to seat it.  A rigid-rim score therefore
counts as "leak" a great many gaps the real part simply closes by deforming.

It also asks the wrong QUESTION.  Coverage fraction is not what seals: a rim
touching over 90 % of its perimeter with one continuous 36 deg hole leaks, while
one touching 88 % broken into a dozen 3 deg specks does not.  A seal is a
CLOSED LOOP.  So this module scores angular CONTINUITY, not coverage.

WHAT IT DOES -- no new seating optimisation
-------------------------------------------
Reads the seatings `align_ear.py` already produced (ears/aligned/*.json) and
re-scores only the seal, three ways:

    budget 1.5 mm   the rigid skirt, for comparison with the existing metric
    budget 2.5 mm   CONSERVATIVE compliance -- drape only, little help from travel
    budget 4.0 mm   OPTIMISTIC compliance -- full drape plus seated carrier

For each ear and budget:

  1. Sample the rim at 360 points (1 deg; an 18 deg gap threshold needs finer
     resolution than the 72 points tryon.py uses).
  2. Model the carrier's axial travel explicitly: translate the rim along the
     nozzle axis, into the ear, over t in [0, 1.5] mm, and keep the best t.
     This is a 1-D seating-depth search inside the seal metric -- NOT a re-run of
     the 6-DOF pose optimiser.  The pose is exactly as `align_ear.py` left it.
  3. A rim sample is SEALED where its signed distance <= budget: negative means
     the lip is already pressed into flesh, positive means the silicone has to
     span that gap and can do so within its deformation budget.
  4. The rim PASSES if it is sealed over >= 95 % of the perimeter AND its largest
     single unsealed arc is <= 18 deg.  Both conditions, because either alone is
     gameable in the way described above.

INTERTRAGIC NOTCH
-----------------
The notch is the soft gap between tragus and antitragus at the inferior-anterior
margin of the concha -- the known risk spot, because there is no cartilage wall
there for the skirt to land on.  Each ear's largest unsealed arc is classified by
where its midpoint falls: the notch sector is +-45 deg about the
inferior-anterior direction, projected into the rim plane (+x anterior, +z
superior in the scan frame; right ears).

Usage:
    python seal_compliance.py                    # all aligned ears
    python seal_compliance.py --md out.md --csv out.csv
"""

from __future__ import annotations

import argparse
import csv
import glob
import json
import os
import sys

import numpy as np
import trimesh

import os as _os

from earfit import (ALIGNED, HERE, NOZZLE_AXIS, NOZZLE_T, EarField, _rim_circle,
                    transform)

BUDGETS = (1.5, 2.5, 4.0)
BUDGET_LABEL = {1.5: "rigid", 2.5: "conservative", 4.0: "optimistic"}
N_RIM = 360                  # 1 deg resolution
COVER_MIN = 0.95             # >= 95 % of the perimeter must be sealed
GAP_MAX_DEG = 18.0           # any single unsealed arc above this leaks
TRAVEL_MM = 1.5              # mag-float carrier axial travel
TRAVEL_STEPS = 16
NOTCH_HALF_DEG = 45.0        # half-width of the intertragic-notch sector


def rim_from_mesh(stl_dir=None, n=N_RIM):
    """Sample the AS-BUILT skirt rim off carrier.stl, one point per azimuth.

    The analytic circle in earfit is a perfect Ø19 ring, which was right until the
    intertragic-notch sector landed: the rim now flares ~1 mm further out over the
    inferior 90 deg, and a circle samples straight past the very feature the
    sector was added to provide.  Taking the max-radius locus per azimuth bin
    picks up whatever shape the generator actually produced, including future
    changes.

    carrier.stl is written in the NOZZLE-LOCAL frame, so the contour is built
    there and then canted into the assembly frame with NOZZLE_T.
    """
    sdir = stl_dir or _os.path.join(HERE, "stl", "right")
    m = trimesh.load(_os.path.join(sdir, "carrier.stl"), force="mesh")
    v = np.asarray(m.vertices, float)
    r = np.hypot(v[:, 1], v[:, 2])
    th = np.degrees(np.arctan2(v[:, 2], v[:, 1])) % 360.0
    # only the lip band, so a wide flange further down the cone cannot win
    lip = v[:, 0] > v[:, 0].max() - 0.60
    v, r, th = v[lip], r[lip], th[lip]
    edges = np.linspace(0.0, 360.0, n + 1)
    idx = np.clip(np.digitize(th, edges) - 1, 0, n - 1)
    pts = np.zeros((n, 3))
    for b in range(n):
        sel = np.where(idx == b)[0]
        if len(sel) == 0:                     # empty bin: fall back to the circle
            a = np.radians(0.5 * (edges[b] + edges[b + 1]))
            from earfit import SKIRT_RIM_R, SKIRT_RIM_X
            pts[b] = [SKIRT_RIM_X, SKIRT_RIM_R * np.cos(a), SKIRT_RIM_R * np.sin(a)]
        else:
            pts[b] = v[sel[np.argmax(r[sel])]]
    return pts @ NOZZLE_T[:3, :3].T + NOZZLE_T[:3, 3]


def longest_false_run(mask):
    """(length, start index) of the longest circular run of False."""
    n = len(mask)
    if mask.all():
        return 0, -1
    if not mask.any():
        return n, 0
    best_len = best_start = 0
    run_len = 0
    run_start = 0
    # walk twice to handle wrap-around
    for i in range(2 * n):
        if not mask[i % n]:
            if run_len == 0:
                run_start = i % n
            run_len += 1
            if run_len > best_len:
                best_len, best_start = run_len, run_start
        else:
            run_len = 0
    return min(best_len, n), best_start


def notch_reference(M):
    """Unit vector for the intertragic notch direction, in the rim plane.

    Inferior and slightly anterior: the notch sits below the tragus at the
    front-bottom of the concha.  Built in the SCAN frame (+x anterior,
    +z superior) and projected into the plane of the rim.
    """
    axis = M[:3, :3] @ NOZZLE_AXIS
    axis = axis / np.linalg.norm(axis)
    ref = np.array([0.45, 0.0, -1.0])            # anterior-ish, inferior
    ref = ref - (ref @ axis) * axis              # into the rim plane
    n = np.linalg.norm(ref)
    return ref / n if n > 1e-9 else None


def score_ear(rec, patch, rim_design):
    field = EarField(patch, seed=0)
    M = np.array(rec["transform"], float)
    axis = M[:3, :3] @ NOZZLE_AXIS
    axis = axis / np.linalg.norm(axis)
    rim0 = transform(rim_design, M)
    centre = rim0.mean(axis=0)

    # angular coordinate of each rim sample, measured about the nozzle axis from
    # the notch reference direction
    ref = notch_reference(M)
    perp = np.cross(axis, ref)
    off = rim0 - centre
    ang = np.degrees(np.arctan2(np.einsum('ij,j->i', off, perp),
                                np.einsum('ij,j->i', off, ref))) % 360.0
    order = np.argsort(ang)
    ang_sorted = ang[order]

    # the signed distances depend only on the travel t, not on the budget, so
    # query once per t and evaluate all three budgets against the same field
    dists = []
    for t in np.linspace(0.0, TRAVEL_MM, TRAVEL_STEPS):
        dists.append((float(t), field.query(rim0 + t * axis)[order]))

    out = {}
    for B in BUDGETS:
        best = None
        for t, d in dists:
            sealed = d <= B
            cover = float(sealed.mean())
            glen, gstart = longest_false_run(sealed)
            gap_deg = 360.0 * glen / len(sealed)
            key = (cover, -gap_deg)
            if best is None or key > best[0]:
                best = (key, t, cover, gap_deg, gstart, glen)
        _, t, cover, gap_deg, gstart, glen = best
        ok = (cover >= COVER_MIN) and (gap_deg <= GAP_MAX_DEG)
        if glen > 0:
            mid = ang_sorted[(gstart + glen // 2) % len(ang_sorted)]
            dev = abs((mid + 180.0) % 360.0 - 180.0)      # deg from notch dir
            in_notch = dev <= NOTCH_HALF_DEG
        else:
            mid, dev, in_notch = float("nan"), float("nan"), False
        # drape only, with the carrier left where the seating put it -- isolates
        # how much of the gain is silicone and how much is the 1.5 mm of travel
        d0 = dists[0][1]
        s0 = d0 <= B
        g0len, _ = longest_false_run(s0)
        g0 = 360.0 * g0len / len(s0)
        ok0 = (float(s0.mean()) >= COVER_MIN) and (g0 <= GAP_MAX_DEG)

        out[B] = dict(travel=t, cover=cover, gap_deg=gap_deg,
                      pass_=bool(ok), gap_mid_deg=float(mid),
                      gap_dev_deg=float(dev), notch=bool(in_notch),
                      cover_notravel=float(s0.mean()), gap_notravel=g0,
                      pass_notravel=bool(ok0))
    return out


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--rim", choices=("mesh", "circle"), default="mesh",
                    help="'mesh' samples the as-built rim off carrier.stl "
                         "(captures the intertragic-notch flare); 'circle' uses "
                         "the analytic Ø19 ring, for comparison with older runs")
    ap.add_argument("--stl-dir", default=None)
    ap.add_argument("--json-dir", default=None,
                    help="directory of per-ear seating JSONs (default ears/aligned); "
                         "patches are always read from ears/aligned")
    ap.add_argument("--csv", default=os.path.join(ALIGNED, "seal_compliance.csv"))
    ap.add_argument("--md", default=os.path.join(ALIGNED, "seal_compliance.md"))
    a = ap.parse_args()

    js = sorted(glob.glob(os.path.join(a.json_dir or ALIGNED, "*.json")))
    if not js:
        sys.exit("nothing in ears/aligned -- run align_ear.py first")

    rim_design = (rim_from_mesh(a.stl_dir, N_RIM) if a.rim == "mesh"
                  else _rim_circle(n=N_RIM))
    print(f"rim model: {a.rim}", flush=True)
    rows = []
    for p in js:
        rec = json.load(open(p))
        patch = trimesh.load(os.path.join(ALIGNED, rec["patch"]), force="mesh")
        s = score_ear(rec, patch, rim_design)
        row = dict(dataset=rec["dataset"], ear_id=rec["ear_id"])
        for B in BUDGETS:
            tag = f"b{B}"
            for k, v in s[B].items():
                row[f"{tag}_{k}"] = v
        rows.append(row)
        print(f"{row['dataset']:>9}/{row['ear_id']:<7} "
              + "  ".join(f"{B}mm cov {s[B]['cover']:.2f} gap {s[B]['gap_deg']:5.1f}deg "
                          f"{'PASS' if s[B]['pass_'] else 'leak'}"
                          + ("(notch)" if s[B]['notch'] and not s[B]['pass_'] else "")
                          for B in BUDGETS), flush=True)

    with open(a.csv, "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=list(rows[0].keys()))
        w.writeheader()
        w.writerows(rows)

    L = ["| skirt budget | model | seal pass | fail | pass rate |",
         "|---|---|---|---|---|"]
    for B in BUDGETS:
        n = sum(1 for r in rows if r[f"b{B}_pass_"])
        L.append(f"| {B:.1f} mm | {BUDGET_LABEL[B]} | **{n}** | {len(rows)-n} | "
                 f"{100*n/len(rows):.0f} % |")
    L += ["", "| skirt budget | failing ears | gap at intertragic notch | elsewhere | notch share |",
          "|---|---|---|---|---|"]
    for B in BUDGETS:
        bad = [r for r in rows if not r[f"b{B}_pass_"]]
        nn = sum(1 for r in bad if r[f"b{B}_notch"])
        share = f"{100*nn/len(bad):.0f} %" if bad else "--"
        L.append(f"| {B:.1f} mm | {len(bad)} | **{nn}** | {len(bad)-nn} | {share} |")
    L += ["", "| skirt budget | pass, drape only | pass, drape + 1.5 mm travel | travel contributes |",
          "|---|---|---|---|"]
    for B in BUDGETS:
        n0 = sum(1 for r in rows if r[f"b{B}_pass_notravel"])
        n1 = sum(1 for r in rows if r[f"b{B}_pass_"])
        L.append(f"| {B:.1f} mm | {n0} | **{n1}** | +{n1-n0} |")
    L += ["", "| skirt budget | median sealed arc | median worst gap | median travel used |",
          "|---|---|---|---|"]
    for B in BUDGETS:
        cov = np.array([r[f"b{B}_cover"] for r in rows])
        gap = np.array([r[f"b{B}_gap_deg"] for r in rows])
        tr = np.array([r[f"b{B}_travel"] for r in rows])
        L.append(f"| {B:.1f} mm | {100*np.median(cov):.0f} % | {np.median(gap):.0f} deg | "
                 f"{np.median(tr):.2f} mm |")
    md = "\n".join(L)
    with open(a.md, "w") as f:
        f.write(md + "\n")
    print("\n" + md)
    print(f"\n{len(rows)} ears -> {a.csv}, {a.md}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
