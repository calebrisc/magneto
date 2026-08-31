#!/usr/bin/env python3
"""
plunger_stroke_spec.py -- size a REPELLING MAGNET STACK per plunger site.

A single repelling pair gives one bounded force over a short stroke, so the cam
has to preset stand-off per wearer.  Put N pairs in series on the guide pin and
the strokes ADD at the same bounded force: the plunger self-fills whatever gap
that ear presents, and the cam disappears.  This sizes the stroke each site needs
from the ears themselves.

MEASUREMENT
    For each ear and each site aim, cast a ray from the boss MOUNT FACE along the
    aim to the ear surface.  That distance -- `D_mount` -- is exactly what the
    plunger stack has to span, since the stack occupies mount -> pad tip.

SIZING RULE (1 mm margin each end, as specified)
    compacted length   L_c  <=  min(D_mount) - 1      fits the closest ear
    expanded reach     L_e  >=  max(D_mount) + 1      reaches the farthest
    stroke             S    =   L_e - L_c  =  range(D_mount) + 2

LIP SEATING (cymba)
    Reaching the first surface is not the same as seating UNDER the overhang.
    For the cymba site we also measure how far along the aim the genuinely
    undercut surface sits (ear normal opposing pull-out), so the stroke can be
    sized to hook rather than merely touch.

PER-GAP CURVE
    MECH_VALIDATION 5 characterises a 5x2.5x1 N35 repelling pair as usable over
    +-0.75 mm about a 2.75 mm rest gap, i.e. 1.5 mm of usable stroke per gap.
    Stacking may shift the optimum per-gap geometry; anything proposed here is
    FLAGGED for the magnet agent to verify rather than asserted.

Usage:
    python plunger_stroke_spec.py --json-dir DIR
"""

from __future__ import annotations

import argparse
import glob
import json
import os
import sys

import numpy as np
import trimesh

from earfit import ALIGNED, NOZZLE_AXIS, EarField, transform

GAP_STROKE = 1.5      # mm usable per repelling gap (+-0.75 about 2.75 rest)
GAP_REST = 2.75       # mm rest gap
MAG_T = 1.0           # mm magnet thickness
UNDERCUT = -0.2       # n . pull_out below this counts as genuine overhang
MARGIN = 1.0          # mm at each end


def geom():
    import generate
    P = generate.PARAMS
    g = generate.G(P)
    return (np.array(g.core_c, float), np.array(g.core_r, float),
            P["clearance"] + P["jacket_thick"], P["plunger_boss_h"],
            g.pl_pad1 + P["plunger_rocker"], P["plunger_pad_t"],
            P["plunger_plate"], g.plungers)


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--json-dir", default=None)
    a = ap.parse_args()
    js = sorted(glob.glob(os.path.join(a.json_dir or ALIGNED, "*.json")))
    if not js:
        sys.exit("no seatings")
    core_c, core_r, off, boss_h, pad_reach, pad_t, plate_t, plungers = geom()

    sites = [(pl["name"], np.asarray(pl["aim"], float)) for pl in plungers]
    print(f"\nn = {len(js)} ears, poses as already seated")
    print(f"boss height {boss_h:.2f} mm; as-built stack mount->pad tip "
          f"{pad_reach:.2f} mm\n")

    results = {}
    for name, aim in sites:
        aim = aim / np.linalg.norm(aim)
        Dm, Dlip = [], []
        for p in js:
            rec = json.load(open(p))
            patch = trimesh.load(os.path.join(ALIGNED, rec["patch"]), force="mesh")
            M = np.array(rec["transform"], float)
            R = M[:3, :3]
            pull = -(R @ NOZZLE_AXIS); pull /= np.linalg.norm(pull)
            surf = core_c + (core_r ** 2 * aim) / np.linalg.norm(core_r * aim)
            mount = surf + aim * (off + boss_h)
            org = transform(mount[None, :], M)[0]
            d_w = R @ aim; d_w /= np.linalg.norm(d_w)
            loc, ir, _ = patch.ray.intersects_location(
                org[None, :], d_w[None, :], multiple_hits=False)
            if not len(ir):
                Dm.append(np.nan); Dlip.append(np.nan)
                continue
            hit = loc[0]
            Dm.append(float(np.linalg.norm(hit - org)))
            # how far along the aim does genuine undercut surface sit?
            field = EarField(patch, seed=0)
            near = field.tree.query_ball_point(hit, 8.0)
            best = np.nan
            if near:
                k = field.nrm[near] @ pull
                sel = [near[i] for i in np.where(k < UNDERCUT)[0]]
                if sel:
                    s = (field.pts[sel] - org) @ d_w
                    s = s[s > 0]
                    if len(s):
                        best = float(np.median(s))
            Dlip.append(best)
        results[name] = (np.array(Dm, float), np.array(Dlip, float))

    for name, (Dm, Dlip) in results.items():
        v = Dm[np.isfinite(Dm)]
        miss = int((~np.isfinite(Dm)).sum())
        print(f"=== {name} ===")
        if not len(v):
            print("  no ear reached\n")
            continue
        print(f"  pad-to-surface engagement D_mount over {len(v)} ears "
              f"({miss} no-hit)")
        print(f"    min {v.min():.2f}  median {np.median(v):.2f}  "
              f"p90 {np.percentile(v,90):.2f}  max {v.max():.2f} mm")
        lo, hi = v.min() - MARGIN, v.max() + MARGIN
        S = hi - lo
        n_gap = int(np.ceil(S / GAP_STROKE))
        per = S / n_gap
        compact = (n_gap + 1) * MAG_T + pad_t + plate_t
        print(f"    sizing: compacted <= {lo:.2f}, expanded >= {hi:.2f} "
              f"-> STROKE {S:.2f} mm")
        print(f"    -> N = {n_gap} gaps x {per:.2f} mm "
              f"(vs {GAP_STROKE:.2f} mm/gap as characterised)")
        print(f"    compacted length {(n_gap+1)}x{MAG_T:.0f} mm magnets + pad "
              f"{pad_t:.1f} + plate {plate_t:.1f} = {compact:.2f} mm "
              f"vs budget {lo:.2f} mm -> "
              f"{'FITS' if compact <= lo else 'DOES NOT FIT'}")
        if compact > lo:
            n_fit = int(np.floor((lo - pad_t - plate_t) / MAG_T)) - 1
            print(f"       at most N = {max(n_fit,0)} gaps fit the closest ear; "
                  f"that is {max(n_fit,0)*GAP_STROKE:.2f} mm of stroke against "
                  f"{S:.2f} mm needed")
        w = Dlip[np.isfinite(Dlip)]
        if len(w):
            print(f"  undercut (lip) depth along the aim, {len(w)} ears: "
                  f"median {np.median(w):.2f}  p90 {np.percentile(w,90):.2f}  "
                  f"max {w.max():.2f} mm")
            hl = w.max() + MARGIN
            print(f"    to SEAT UNDER THE LIP on every ear: expanded reach "
                  f">= {hl:.2f} mm -> stroke {hl - lo:.2f} mm "
                  f"(N = {int(np.ceil((hl-lo)/GAP_STROKE))} gaps)")
        print()
    return 0


if __name__ == "__main__":
    sys.exit(main())
