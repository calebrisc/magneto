#!/usr/bin/env python3
"""
shrink_estimate.py -- what would protrusion be if the shell were smaller?

A first-order estimate for a driver-downsizing decision, WITHOUT rebuilding the
geometry.  Protrusion is

    max over faceplate points p of  (M @ p - tragus) . floor_normal

so displacing each faceplate point by a small delta in the DESIGN frame changes
it by (M[:3,:3] @ delta) . floor_normal, exactly, at a fixed pose.  That makes a
shell shrink cheap to evaluate: offset the points and re-take the max.

Default models a shell 2 mm smaller in Y and Z -- roughly an 8 mm driver in place
of the current one:

    Y (width, symmetric about y = 0): 1 mm inward per side, delta_y = -sign(y)
    Z (stack height): the whole 2 mm comes off the +Z / faceplate side, because
      the -Z jacket face still has to bed on the concha floor.

WHAT THIS IS NOT.  The pose is frozen at the seating `align_ear.py` found for the
full-size shell.  A genuinely smaller shell would also re-seat -- generally
deeper and flatter in the bowl -- so this is the CONSERVATIVE half of the answer:
it captures the geometry that disappears, not the better fit that follows.  Treat
it as a lower bound on the benefit, and rebuild before committing.

Usage:
    python shrink_estimate.py                    # 2 mm in Y and Z
    python shrink_estimate.py --dy 2 --dz 2 --dy 3
"""

from __future__ import annotations

import argparse
import glob
import json
import os
import sys

import numpy as np

from earfit import ALIGNED, iem_points, transform
from tryon import PROT_MARGINAL, PROT_PASS


def protrusion(P, recs, dy=0.0, dz=0.0):
    """Protrusion per ear with the faceplate shrunk by dy (total, symmetric in Y)
    and dz (total, off the +Z side)."""
    fp = P["faceplate"].copy()
    fp[:, 1] -= np.sign(fp[:, 1]) * (0.5 * dy)
    fp[:, 2] -= dz
    out = []
    for rec in recs:
        M = np.array(rec["transform"], float)
        n = np.array(rec["floor_normal"], float)
        trg = np.array(rec["tragus"], float)
        w = transform(fp, M)
        out.append(float(np.max(np.einsum("ij,j->i", w - trg, n))))
    return np.array(out)


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--json-dir", default=None,
                    help="directory of per-ear seating JSONs (default ears/aligned). "
                         "Point this at a saved copy when another run may be "
                         "rewriting ears/aligned underneath you.")
    ap.add_argument("--case", action="append", default=[],
                    help="a 'dy,dz' shrink in mm; repeatable (default 0,0 and 2,2)")
    a = ap.parse_args()
    cases = [tuple(float(x) for x in c.split(",")) for c in a.case] or [(0.0, 0.0),
                                                                       (2.0, 2.0)]
    jdir = a.json_dir or ALIGNED
    recs = [json.load(open(p)) for p in sorted(glob.glob(os.path.join(jdir, "*.json")))]
    if not recs:
        sys.exit("nothing in ears/aligned")
    P = iem_points()

    print(f"n = {len(recs)} ears, pose frozen at the as-seated full-size shell\n")
    print("| shell shrink | protrusion median | p90 | max | pass ≤10 mm | ≤14 mm | fail >14 mm |")
    print("|---|---|---|---|---|---|---|")
    for dy, dz in cases:
        v = protrusion(P, recs, dy, dz)
        lab = "none (as built)" if (dy == 0 and dz == 0) else f"−{dy:g} mm Y, −{dz:g} mm Z"
        print(f"| {lab} | {np.median(v):.2f} mm | {np.percentile(v,90):.2f} | "
              f"{v.max():.2f} | {100*np.mean(v<=PROT_PASS):.0f} % | "
              f"{100*np.mean(v<=PROT_MARGINAL):.0f} % | "
              f"{100*np.mean(v>PROT_MARGINAL):.0f} % |")
    return 0


if __name__ == "__main__":
    sys.exit(main())
