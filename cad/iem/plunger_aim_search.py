#!/usr/bin/env python3
"""
plunger_aim_search.py -- derive a plunger aim from the ears, not from intuition.

The antitragus leg of the tripod misses P0023 by 14.6 mm (see TRYON_REPORT.md),
so its aim (0.20, -0.96, -0.20) was picked without checking where cartilage
actually is.  This sweeps candidate aims and reports where the surface really
sits, across the ears already aligned and seated.  GEOMETRY QUERIES ONLY -- no
re-seating; each ear keeps the pose it already has.

METHOD
    For a candidate unit aim `a` in the design frame, the generator places the
    plunger boss by walking out of the core along `a`:

        surf  = core_c + (r^2 * a) / |r * a|          (core ellipsoid surface)
        base  = surf + a * (clearance + jacket_thick)  (jacket outer surface)
        mount = base + a * plunger_boss_h
        pad   = mount + a * (pl_pad1 + plunger_rocker)

    So the leg's reach from the jacket's outer surface to the pad tip is
    STACK = boss_h + pl_pad1 + rocker.  We cast a ray from `base` along `a` into
    the ear and record the distance to the first hit.  That distance is what the
    leg has to span.

WINDOWS
    as-built   the existing stack can serve surface distances in
               [STACK - cam, STACK + travel] -- the cam presets stand-off over
               its range and the spring adds travel.
    target     the coordinator's stated 1.5-6.0 mm working range, measured from
               the same reference.  It is much shorter than the as-built stack,
               so it describes a shorter leg; both are reported rather than
               guessing which is meant.

Usage:
    python plunger_aim_search.py --json-dir DIR [--grid 10] [--extra x,y,z]
"""

from __future__ import annotations

import argparse
import glob
import json
import os
import sys

import numpy as np
import trimesh

from earfit import ALIGNED, EarField, transform

TARGET_LO, TARGET_HI = 1.5, 6.0


def geom():
    import generate
    P = generate.PARAMS
    g = generate.G(P)
    stack = (P["plunger_boss_h"] + g.pl_pad1 + P["plunger_rocker"])
    return (np.array(g.core_c, float), np.array(g.core_r, float),
            P["clearance"] + P["jacket_thick"], stack,
            P["plunger_cam_range"], P["plunger_travel"], P["plunger_aims"])


def base_point(core_c, core_r, off, a):
    surf = core_c + (core_r ** 2 * a) / np.linalg.norm(core_r * a)
    return surf + a * off


def hemisphere(step_deg, axis=np.array([0.0, -1.0, 0.0]), max_off=90.0):
    """Unit directions on a `step_deg` grid within max_off of `axis`."""
    out = []
    for pol in np.arange(0.0, max_off + 1e-9, step_deg):
        if pol == 0.0:
            out.append(axis / np.linalg.norm(axis))
            continue
        n_az = max(1, int(round(360.0 / step_deg * np.sin(np.radians(pol)))))
        ref = np.array([1.0, 0.0, 0.0])
        u = np.cross(axis, ref); u /= np.linalg.norm(u)
        v = np.cross(axis, u)
        for az in np.linspace(0, 2 * np.pi, n_az, endpoint=False):
            d = (np.cos(np.radians(pol)) * axis / np.linalg.norm(axis)
                 + np.sin(np.radians(pol)) * (np.cos(az) * u + np.sin(az) * v))
            out.append(d / np.linalg.norm(d))
    return out


def probe(recs, aims, core_c, core_r, off):
    """Distance to the ear along each aim, and whether that surface INTERLOCKS.

    D[i, j]  distance from the boss base along aims[j] to ear i's surface
    K[i, j]  n_ear . pull_out at the hit.  Negative means the surface faces back
             against the escape direction -- geometric interlock, which resists
             pull-out far beyond friction.  Positive means the surface would help
             eject the shell, so a pad there only ever contributes friction.
    """
    D = np.full((len(recs), len(aims)), np.inf)
    K = np.full((len(recs), len(aims)), np.nan)
    A = np.array(aims)
    from earfit import NOZZLE_AXIS
    for i, (rec, patch) in enumerate(recs):
        M = np.array(rec["transform"], float)
        R = M[:3, :3]
        pull = -(R @ NOZZLE_AXIS)
        pull /= np.linalg.norm(pull)
        field = EarField(patch, seed=0)
        org = np.array([transform(base_point(core_c, core_r, off, a)[None, :], M)[0]
                        for a in A])
        dirs = np.array([R @ a for a in A])
        dirs /= np.linalg.norm(dirs, axis=1, keepdims=True)
        loc, ir, _ = patch.ray.intersects_location(org, dirs, multiple_hits=False)
        if len(ir):
            d = np.linalg.norm(loc - org[ir], axis=1)
            for k, j in enumerate(ir):
                if d[k] < D[i, j]:
                    D[i, j] = d[k]
                    _, ni = field.tree.query(loc[k][None, :])
                    K[i, j] = float(field.nrm[int(ni[0])] @ pull)
    return D, K


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--json-dir", default=None)
    ap.add_argument("--grid", type=float, default=10.0)
    ap.add_argument("--extra", action="append", default=[],
                    help="additional aim 'x,y,z' to evaluate explicitly")
    a = ap.parse_args()

    js = sorted(glob.glob(os.path.join(a.json_dir or ALIGNED, "*.json")))
    if not js:
        sys.exit("no seatings found")
    recs = []
    for p in js:
        rec = json.load(open(p))
        recs.append((rec, trimesh.load(os.path.join(ALIGNED, rec["patch"]),
                                       force="mesh")))
    core_c, core_r, off, stack, cam, travel, builtin = geom()
    lo_built, hi_built = stack - cam, stack + travel
    n = len(recs)
    print(f"\nn = {n} ears, poses as already seated (no re-seating)")
    print(f"leg stack, jacket outer surface -> pad tip: {stack:.2f} mm")
    print(f"as-built window [{lo_built:.2f}, {hi_built:.2f}] mm "
          f"(cam {cam:.1f} + travel {travel:.2f})")
    print(f"target window   [{TARGET_LO:.2f}, {TARGET_HI:.2f}] mm (stated)\n")

    aims = hemisphere(a.grid)
    labels = [None] * len(aims)
    for nm, v in builtin:
        aims.append(np.array(v, float) / np.linalg.norm(v)); labels.append(nm)
    for e in a.extra:
        v = np.array([float(x) for x in e.split(",")], float)
        aims.append(v / np.linalg.norm(v)); labels.append("extra " + e)

    D, K = probe(recs, aims, core_c, core_r, off)
    A = np.array(aims)

    def score(lo, hi):
        return ((D >= lo) & (D <= hi)).sum(axis=0)

    def lock(lo, hi):
        """ears both in range AND landing on an interlocking surface."""
        return ((D >= lo) & (D <= hi) & (K < 0)).sum(axis=0)

    # A leg aimed at the concha floor scores well on distance and is useless for
    # retention -- it presses where the jacket already presses.  Rank inferior-leg
    # candidates only among directions with a real inferior component.
    infer = A[:, 1] <= -0.30
    for lo, hi, tag in ((TARGET_LO, TARGET_HI, "target 1.5-6.0"),
                        (lo_built, hi_built, f"as-built {lo_built:.1f}-{hi_built:.1f}")):
        c, kk = score(lo, hi), lock(lo, hi)
        for sub, sub_tag in ((np.ones(len(aims), bool), "all directions"),
                             (infer, "inferior only (a_y <= -0.30)")):
            order = [j for j in np.argsort(-(c + 0.5 * kk)) if sub[j]]
            print(f"--- best aims, {tag} mm window, {sub_tag} ---")
            print("| aim (x,y,z) | in range | interlocking | median dist | p10-p90 |")
            print("|---|---|---|---|---|")
            seen = 0
            for j in order:
                if seen >= 5:
                    break
                fin = D[:, j][np.isfinite(D[:, j])]
                if not len(fin):
                    continue
                lab = f" `{labels[j]}`" if labels[j] else ""
                print(f"| ({A[j][0]:+.2f}, {A[j][1]:+.2f}, {A[j][2]:+.2f}){lab} "
                      f"| **{c[j]}**/{n} | {kk[j]}/{n} | {np.median(fin):.1f} mm "
                      f"| {np.percentile(fin,10):.1f}-{np.percentile(fin,90):.1f} |")
                seen += 1
            print()

    print("--- named aims (built-in + extra) ---")
    print("| aim | in target | in as-built | median | min | max | misses |")
    print("|---|---|---|---|---|---|---|")
    print("| aim | in target | in as-built | interlock | median | min | max | misses |")
    for j, lab in enumerate(labels):
        if lab is None:
            continue
        v = D[:, j]
        fin = v[np.isfinite(v)]
        it = int(((v >= TARGET_LO) & (v <= TARGET_HI)).sum())
        ib = int(((v >= lo_built) & (v <= hi_built)).sum())
        kl = int(np.nansum(K[:, j] < 0))
        med = f"{np.median(fin):.1f}" if len(fin) else "--"
        mn = f"{fin.min():.1f}" if len(fin) else "--"
        mx = f"{fin.max():.1f}" if len(fin) else "--"
        print(f"| {lab} ({A[j][0]:+.2f},{A[j][1]:+.2f},{A[j][2]:+.2f}) "
              f"| {it}/{n} | {ib}/{n} | {kl}/{n} | {med} | {mn} | {mx} "
              f"| {int((~np.isfinite(v)).sum())} |")

    # per-ear spread for the best target-window aim
    c, kk = score(TARGET_LO, TARGET_HI), lock(TARGET_LO, TARGET_HI)
    cand = [j for j in range(len(aims)) if infer[j]]
    j = max(cand, key=lambda k: (c[k] + 0.5 * kk[k])) if cand else int(np.argmax(c))
    print(f"\n--- per-ear distance for the winning INFERIOR target aim "
          f"({aims[j][0]:+.2f}, {aims[j][1]:+.2f}, {aims[j][2]:+.2f}) ---")
    print("| ear | distance | in target window | interlock |")
    print("|---|---|---|---|")
    for i, (rec, _) in enumerate(recs):
        d = D[i, j]
        ok = TARGET_LO <= d <= TARGET_HI
        ds = "no hit" if not np.isfinite(d) else f"{d:.2f} mm"
        kv = K[i, j]
        ks = "--" if not np.isfinite(kv) else ("yes" if kv < 0 else "no")
        print(f"| {rec['dataset']}/{rec['ear_id']} | {ds} | "
              f"{'yes' if ok else 'NO'} | {ks} |")
    fin = D[:, j][np.isfinite(D[:, j])]
    if len(fin):
        print(f"\nspread: median {np.median(fin):.2f}, p10 {np.percentile(fin,10):.2f}, "
              f"p90 {np.percentile(fin,90):.2f}, min {fin.min():.2f}, max {fin.max():.2f} mm; "
              f"{int((~np.isfinite(D[:, j])).sum())} ears with no hit")
    return 0


if __name__ == "__main__":
    sys.exit(main())
