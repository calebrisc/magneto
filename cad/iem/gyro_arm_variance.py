#!/usr/bin/env python3
"""
gyro_arm_variance.py -- how much do the gyro-arm target regions actually vary?

Feeds the staged-deployment concept: an arm reaches the cymba, then a distal
segment continues under the antihelix/crus when the proximal stage feels
resistance.  The question this answers is how many stages and degrees of freedom
the measured anatomy actually demands -- as opposed to how many seem prudent.

SKIRT-DATUM FRAME (per ear, so every ear is measured in the same place)
    origin  the detected aperture centroid
    +Z      the canal axis, CONSTRUCTED as the inward concha-floor normal raked
            40 deg anteriorly (the scans have no resolvable canal; see
            viz_marked.py).  Points medially, into the head.
    +Y      superior, orthogonalised against +Z
    +X      +Y x +Z, anterior for a right ear
    Cylindrical (r, theta, z): theta measured from +X (anterior) toward +Y
    (superior), so 90 deg is straight up, 180 deg posterior.

FEATURES
    rim path      for each azimuth, the most LATERAL surface point at mid radius
                  -- the concha rim crest.  Resampled on a fixed azimuth grid so
                  ears are point-wise comparable.
    cymba pocket  deepest point in the superior sector, its depth below the local
                  rim, and the centroid of the OVERHANGING lip above it
    antihelix     the undercut lip in the superior-posterior sector, with a
                  circle fitted in the datum XY plane: centre, radius, arc extent
    crus helicis  the anterior ridge between cymba and cavum -- the most lateral
                  point in the anterior sector at mid radius

An UNDERCUT is surface whose outward normal has a medial component (n . +Z >
UNDERCUT), i.e. the underside of an overhang: the only geometry that resists
pull-out without relying on friction.

SCALE VS SHAPE
    Every statistic is computed twice: raw, and after dividing each ear's
    lengths by its own size scale (median rim radius).  If the SD collapses
    under normalisation the population differs only in SIZE, and one
    spring-loaded arc with a single scale degree of freedom covers it.  If it
    does not, the SHAPE differs and the arm needs articulation.

Usage:
    python gyro_arm_variance.py --json-dir DIR
"""

from __future__ import annotations

import argparse
import glob
import json
import os
import sys

import numpy as np
import trimesh

from earfit import ALIGNED

UNDERCUT = 0.20
RIM_R = (4.0, 15.0)      # mm; the concha bowl, NOT the pinna
Z_BAND = 12.0            # mm along the canal axis either side of the aperture
AZ_STEP = 10.0
SEC = dict(crus=(-30.0, 40.0), cymba=(55.0, 125.0), antihelix=(110.0, 215.0))
CANAL_RAKE = 40.0


def datum(rec):
    ap = np.array(rec["aperture"], float)
    n_in = -np.array(rec["floor_normal"], float)
    n_in /= np.linalg.norm(n_in)
    ant = np.array([1.0, 0.0, 0.0])
    t = ant - (ant @ n_in) * n_in
    t /= np.linalg.norm(t)
    r = np.radians(CANAL_RAKE)
    ez = np.cos(r) * n_in + np.sin(r) * t
    ez /= np.linalg.norm(ez)
    up = np.array([0.0, 0.0, 1.0])
    ey = up - (up @ ez) * ez
    ey /= np.linalg.norm(ey)
    ex = np.cross(ey, ez)
    ex /= np.linalg.norm(ex)
    return ap, np.column_stack([ex, ey, ez])


def to_datum(P, ap, B):
    return (P - ap) @ B


def sector(th, lo, hi):
    t = (th + 180.0) % 360.0 - 180.0
    return (t >= lo) & (t <= hi)


def fit_circle(xy):
    A = np.column_stack([2 * xy[:, 0], 2 * xy[:, 1], np.ones(len(xy))])
    b = (xy ** 2).sum(axis=1)
    sol, *_ = np.linalg.lstsq(A, b, rcond=None)
    c = sol[:2]
    R = float(np.sqrt(max(sol[2] + c @ c, 1e-9)))
    ang = np.degrees(np.arctan2(xy[:, 1] - c[1], xy[:, 0] - c[0]))
    ang = np.sort((ang + 360.0) % 360.0)
    gaps = np.diff(np.concatenate([ang, ang[:1] + 360.0]))
    extent = 360.0 - gaps.max()
    return c, R, float(extent)


def measure(rec, patch):
    ap, B = datum(rec)
    V = to_datum(np.asarray(patch.vertices), ap, B)
    C = to_datum(patch.triangles_center, ap, B)
    N = np.asarray(patch.face_normals) @ B
    r = np.hypot(C[:, 0], C[:, 1])
    th = np.degrees(np.arctan2(C[:, 1], C[:, 0]))
    out = {}

    # --- rim path: most lateral (min z) point per azimuth ------------------- #
    az = np.arange(-180.0, 180.0, AZ_STEP)
    # The cached patch is the whole 85 mm ear window -- pinna, helix and a slab
    # of scalp.  Without this constraint the "rim" detector happily returns the
    # helix or the skull behind it, which is what a first pass did: it reported
    # 30-48 mm concha depths against an anatomical 9-17 mm.
    band = (r > RIM_R[0]) & (r < RIM_R[1]) & (np.abs(C[:, 2]) < Z_BAND)
    rim = np.full((len(az), 3), np.nan)
    for i, a0 in enumerate(az):
        m = band & sector(th, a0, a0 + AZ_STEP)
        if m.sum() >= 3:
            j = np.argmin(C[m][:, 2])
            rim[i] = C[m][j]
    out["rim"] = rim
    out["rim_radius"] = float(np.nanmedian(np.hypot(rim[:, 0], rim[:, 1])))

    # --- cymba pocket ------------------------------------------------------- #
    m = band & sector(th, *SEC["cymba"])
    if m.sum() >= 5:
        floor = C[m][np.argmax(C[m][:, 2])]
        out["cymba_floor"] = floor
        i = int(np.argmin(np.abs(az - np.degrees(np.arctan2(floor[1], floor[0])))))
        out["cymba_depth"] = float(floor[2] - rim[i, 2]) if np.isfinite(rim[i, 2]) \
            else np.nan
        u = m & (N[:, 2] > UNDERCUT)
        out["cymba_lip"] = C[u].mean(axis=0) if u.sum() >= 3 else np.full(3, np.nan)
        out["cymba_lip_area"] = float(patch.area_faces[u].sum()) if u.sum() else 0.0
    # --- antihelix undercut ------------------------------------------------- #
    m = band & sector(th, *SEC["antihelix"]) & (N[:, 2] > UNDERCUT)
    if m.sum() >= 6:
        P = C[m]
        out["ah_closest"] = P[np.argmin(np.linalg.norm(P, axis=1))]
        c, R, ext = fit_circle(P[:, :2])
        out["ah_c"] = c; out["ah_R"] = R; out["ah_extent"] = ext
        out["ah_area"] = float(patch.area_faces[m].sum())
    # --- crus helicis ------------------------------------------------------- #
    m = band & sector(th, *SEC["crus"]) & (r > 5.0) & (r < 14.0)
    if m.sum() >= 5:
        out["crus"] = C[m][np.argmin(C[m][:, 2])]
    return out


def stats(vals, lab, unit="mm"):
    v = np.array([x for x in vals if x is not None and np.isfinite(x)], float)
    if not len(v):
        return f"| {lab} | -- | -- | -- | -- |"
    return (f"| {lab} | {v.mean():.2f} | {v.std(ddof=1):.2f} | "
            f"{v.min():.2f}–{v.max():.2f} | {len(v)} |")


def main():
    ap_ = argparse.ArgumentParser(description=__doc__)
    ap_.add_argument("--json-dir", default=None)
    a = ap_.parse_args()
    js = sorted(glob.glob(os.path.join(a.json_dir or ALIGNED, "*.json")))
    if not js:
        sys.exit("no seatings")
    ALL = []
    for p in js:
        rec = json.load(open(p))
        patch = trimesh.load(os.path.join(ALIGNED, rec["patch"]), force="mesh")
        m = measure(rec, patch)
        m["id"] = rec["ear_id"]
        # INDEPENDENT size proxy.  Normalising by the rim radius this script
        # measures would divide the features by a quantity carrying the same
        # detector noise, inflating CV rather than removing scale.  basin_
        # inscribed comes from align_ear's depth map and knows nothing about
        # these detectors.
        m["scale"] = float(rec["basin_inscribed"])
        m["synthetic"] = rec["dataset"] == "synthetic"
        ALL.append(m)
    # Synthetic corners are parametric bowls, not anatomy.  Mixing them into a
    # variance study of real ears inflates exactly the number being measured.
    M = [m for m in ALL if not m["synthetic"]]
    n = len(M)
    print(f"\n{len(ALL)} ears measured; {len(ALL)-n} synthetic corners excluded "
          f"from the statistics (parametric bowls, not anatomy)")
    print(f"\nn = {n} ears, skirt-datum frame (origin = aperture centroid, "
          f"+Z = canal axis)")
    print(f"NOTE: this is the 13-ear short list; the same script runs over all "
          f"~103 aligned ears unchanged.\n")

    scales = np.array([m["scale"] for m in M])
    rr = np.array([m["rim_radius"] for m in M])
    print(f"measured rim radius: mean {rr.mean():.2f}, SD {rr.std(ddof=1):.2f} mm")
    print(f"size proxy (basin inscribed radius, independent): mean {scales.mean():.2f}, "
          f"SD {scales.std(ddof=1):.2f}, range {scales.min():.2f}–{scales.max():.2f} mm "
          f"({100*scales.std(ddof=1)/scales.mean():.0f} % CV)\n")

    def col(key, idx=None, norm=False):
        out = []
        for m, s in zip(M, scales):
            v = m.get(key)
            if v is None:
                out.append(np.nan); continue
            x = v if idx is None else np.asarray(v)[idx]
            x = float(x)
            out.append(x / s if norm else x)
        return out

    def dist(key, norm=False):
        out = []
        for m, s in zip(M, scales):
            v = m.get(key)
            if v is None or not np.all(np.isfinite(np.asarray(v, float))):
                out.append(np.nan); continue
            d = float(np.linalg.norm(np.asarray(v, float)))
            out.append(d / s if norm else d)
        return out

    for tag, norm in (("RAW (mm)", False), ("SCALE-NORMALISED (x / ear scale)", True)):
        print(f"--- {tag} ---")
        print("| feature | mean | SD | range | n |")
        print("|---|---|---|---|---|")
        print(stats(dist("cymba_floor", norm), "cymba floor, distance from origin"))
        print(stats(col("cymba_depth", None, norm), "cymba pocket depth"))
        print(stats(dist("cymba_lip", norm), "cymba lip centroid, distance"))
        print(stats(dist("ah_closest", norm), "antihelix lip, closest point"))
        print(stats(col("ah_R", None, norm), "antihelix arc radius"))
        print(stats(col("ah_extent"), "antihelix arc extent (deg)"))
        print(stats(dist("crus", norm), "crus helicis, distance from origin"))
        print()

    print("--- does normalising by size remove the variance? (CV = SD/mean) ---")
    print("| feature | CV raw | CV scaled | verdict |")
    print("|---|---|---|---|")
    feats = [("cymba floor dist", lambda nz: dist("cymba_floor", nz)),
             ("cymba pocket depth", lambda nz: col("cymba_depth", None, nz)),
             ("cymba lip dist", lambda nz: dist("cymba_lip", nz)),
             ("antihelix closest", lambda nz: dist("ah_closest", nz)),
             ("antihelix arc radius", lambda nz: col("ah_R", None, nz)),
             ("crus dist", lambda nz: dist("crus", nz))]
    for lab, fn in feats:
        a0 = np.array([x for x in fn(False) if np.isfinite(x)])
        a1 = np.array([x for x in fn(True) if np.isfinite(x)])
        if len(a0) < 3 or len(a1) < 3:
            continue
        c0, c1 = a0.std(ddof=1) / a0.mean(), a1.std(ddof=1) / a1.mean()
        v = "SCALE" if c1 < 0.75 * c0 else ("shape" if c1 > 1.1 * c0 else "mixed")
        print(f"| {lab} | {c0:.2f} | {c1:.2f} | {v} |")
    print()

    # --- rim path point-wise variance --------------------------------------- #
    az = np.arange(-180.0, 180.0, AZ_STEP)
    R = np.stack([m["rim"] for m in M])                     # ears x az x 3
    for tag, norm in (("RAW", False), ("SCALE-NORMALISED", True)):
        A = R / scales[:, None, None] if norm else R
        sd = np.sqrt(np.nanmean(np.nanvar(A, axis=0, ddof=1), axis=1))
        ok = np.isfinite(sd)
        print(f"--- concha rim path, point-wise across-ear SD ({tag}) ---")
        print(f"  overall: mean {np.nanmean(sd):.2f}, min {np.nanmin(sd):.2f}, "
              f"max {np.nanmax(sd):.2f}"
              + ("" if norm else " mm"))
        lo = [(az[i], sd[i]) for i in range(len(az)) if ok[i] and sd[i] < 1.0]
        print(f"  azimuths with SD < 1.0: {len(lo)}/{int(ok.sum())}"
              + ("" if not lo else "  -> " +
                 ", ".join(f"{int(t)}deg:{v:.2f}" for t, v in lo[:10])))
        band_lbl = [("anterior/crus", -30, 40), ("cymba", 55, 125),
                    ("antihelix", 110, 215), ("posterior", 150, 215)]
        for nm, l, h in band_lbl:
            m2 = np.array([ok[i] and (l <= ((az[i] + 180) % 360 - 180) <= h)
                           for i in range(len(az))])
            if m2.any():
                print(f"    {nm:<15} mean SD {np.nanmean(sd[m2]):.2f}")
        print()

    json.dump([{k: (v.tolist() if isinstance(v, np.ndarray) else v)
                for k, v in m.items()} for m in M],
              open(os.path.join(ALIGNED, "gyro_arm_variance.json"), "w"), indent=1)
    return 0


if __name__ == "__main__":
    sys.exit(main())
