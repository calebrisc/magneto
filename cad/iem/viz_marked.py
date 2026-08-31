#!/usr/bin/env python3
"""
viz_marked.py -- one ear, landmarks painted, assembly seated BY CONSTRUCTION.

No optimiser anywhere in this file.  The pose is built from the detected
landmarks by direct construction, so the render shows what the landmark data
actually says rather than what a cost function negotiated:

    nozzle axis  ->  the canal direction probed from the mesh
    skirt rim centre -> the detected aperture centroid
    faceplate (+Z)   -> forced to face out of the head

MARKERS
    magenta ring   the aperture boundary, found by casting rays radially out
                   from the aperture centroid in the plane perpendicular to the
                   canal axis and taking the first surface hit -- so it traces
                   the real funnel wall, it is not a drawn circle of assumed size
    magenta sphere the aperture centroid itself
    tinted floor   ear faces within CAVUM_R of the aperture, lightly recoloured

Usage:
    python viz_marked.py --json-dir DIR --ear P0023
"""

from __future__ import annotations

import argparse
import glob
import json
import os
import sys

import numpy as np
import trimesh

from align_ear import cone_dirs
from earfit import ALIGNED, HERE, NOZZLE_AXIS, NOZZLE_T, RIM_CENTRE, transform
from viz_scene import cluster_decimate

VIZ = os.path.join(HERE, "viz")
CAVUM_R = 11.0
CANAL_RAKE = 40.0        # deg, canal axis vs concha-floor normal (anatomical)
COL = {"ear": (205, 175, 135, 255), "cavum": (232, 205, 170, 255),
       "core": (190, 192, 200, 255), "faceplate": (205, 207, 214, 255),
       "jacket_wing": (62, 108, 198, 255), "carrier": (232, 140, 48, 255),
       "nozzle_insert_short": (128, 128, 130, 255),
       "mark": (255, 0, 190, 255)}


def build_pose(ap, canal_dir, outward):
    """Rotation+translation placing the rim on the aperture, nozzle down the canal."""
    a1 = NOZZLE_AXIS / np.linalg.norm(NOZZLE_AXIS)
    up_d = np.array([0.0, 1.0, 0.0])                    # design superior
    a2 = up_d - (up_d @ a1) * a1
    a2 /= np.linalg.norm(a2)
    a3 = np.cross(a1, a2)

    b1 = canal_dir / np.linalg.norm(canal_dir)
    up_w = np.array([0.0, 0.0, 1.0])                    # world superior
    b2 = up_w - (up_w @ b1) * b1
    if np.linalg.norm(b2) < 1e-6:
        b2 = np.cross(b1, outward)
    b2 /= np.linalg.norm(b2)
    b3 = np.cross(b1, b2)

    R = np.column_stack([b1, b2, b3]) @ np.column_stack([a1, a2, a3]).T
    if (R @ np.array([0.0, 0.0, 1.0])) @ outward < 0:   # faceplate must face out
        R = np.column_stack([b1, -b2, -b3]) @ np.column_stack([a1, a2, a3]).T
    M = np.eye(4)
    M[:3, :3] = R
    M[:3, 3] = ap - R @ RIM_CENTRE
    return M


def canal_axis(patch, ap, inward, half=55.0, n=400, cap=25.0):
    """Deepest penetrating direction within a cone about the inward floor normal.

    align_ear's unconstrained canal_probe is unusable here: in these scans the
    canal is a dimple, so the longest unobstructed run heads ALONG the concha
    wall rather than into the ear -- on P0023 it came out 98 deg from the floor
    normal, i.e. skimming the surface.  Constraining the search to a cone about
    the inward normal keeps it a canal axis rather than a wall-grazing ray.
    """
    dirs = cone_dirs(inward, half, n)
    org = np.tile(ap + inward * -0.8, (len(dirs), 1))       # start just outside
    loc, ir, _ = patch.ray.intersects_location(org, dirs, multiple_hits=False)
    best, bd = inward, -1.0
    runs = np.full(len(dirs), cap)
    for k, j in enumerate(ir):
        runs[j] = min(runs[j], np.linalg.norm(loc[k] - org[j]))
    j = int(np.argmax(runs))
    return float(runs[j]), dirs[j]


def aperture_ring(patch, ap, axis, outward, n=96, reach=14.0):
    """Trace the funnel wall around the aperture: radial rays, first hit each."""
    r = np.array([0.0, 0.0, 1.0])
    if abs(r @ axis) > 0.9:
        r = np.array([1.0, 0.0, 0.0])
    u = np.cross(axis, r); u /= np.linalg.norm(u)
    v = np.cross(axis, u)
    # lift the origin off the surface, or every ray hits its own start point
    start = ap + outward * 0.8
    org, dirs = [], []
    for t in np.linspace(0, 2 * np.pi, n, endpoint=False):
        org.append(start)
        dirs.append(np.cos(t) * u + np.sin(t) * v)
    loc, ir, _ = patch.ray.intersects_location(np.array(org), np.array(dirs),
                                               multiple_hits=False)
    pts = {}
    for k, j in enumerate(ir):
        d = np.linalg.norm(loc[k] - start)
        if 0.8 <= d <= reach and (j not in pts or d < pts[j][1]):
            pts[j] = (loc[k], d)
    return np.array([p for p, _d in (pts[j] for j in sorted(pts))]) if pts \
        else np.zeros((0, 3))


def main():
    ap_ = argparse.ArgumentParser(description=__doc__)
    ap_.add_argument("--json-dir", default=None)
    ap_.add_argument("--ear", default="P0023")
    a = ap_.parse_args()
    fs = glob.glob(os.path.join(a.json_dir or ALIGNED, f"*_{a.ear}_right.json"))
    if not fs:
        sys.exit(f"{a.ear} not found")
    rec = json.load(open(fs[0]))
    patch = trimesh.load(os.path.join(ALIGNED, rec["patch"]), force="mesh")
    ap = np.array(rec["aperture"], float)
    outward = np.array(rec["outward"], float)
    medial = -outward

    n_in = -np.array(rec["floor_normal"], float)
    n_in /= np.linalg.norm(n_in)
    run, probed = canal_axis(patch, ap, n_in)
    print(f"{a.ear}: constrained canal probe -> run {run:.2f} mm, "
          f"dir {np.round(probed, 3)}")

    # The probe is reported but NOT used.  1-2 mm of run means the scan has no
    # canal to find -- structured light cannot see past the first bend, as this
    # report documents (median canal_run 1.8 mm across all 107 ears).  So the
    # axis is CONSTRUCTED from anatomy instead of pretending to measure it:
    # the inward concha-floor normal raked anteriorly by CANAL_RAKE, the angle
    # docs/EAR_ANTHROPOMETRY.md gives and make_synthetic_ear.py builds to.
    ant = np.array([1.0, 0.0, 0.0])
    t_ant = ant - (ant @ n_in) * n_in
    t_ant /= np.linalg.norm(t_ant)
    r = np.radians(CANAL_RAKE)
    cdir = np.cos(r) * n_in + np.sin(r) * t_ant
    cdir /= np.linalg.norm(cdir)
    off_norm = float(np.degrees(np.arccos(np.clip(cdir @ n_in, -1, 1))))
    print(f"  axis CONSTRUCTED: inward floor normal raked {CANAL_RAKE:.0f} deg "
          f"anteriorly -> {np.round(cdir, 3)} ({off_norm:.0f} deg off normal)")
    M = build_pose(ap, cdir, outward)
    print(f"  constructed pose: nozzle -> canal, rim centre -> aperture, "
          f"faceplate outward")
    chk = transform(np.array([RIM_CENTRE]), M)[0]
    print(f"  rim centre lands {np.linalg.norm(chk - ap):.3f} mm from the aperture")

    scene = trimesh.Scene()
    ear = cluster_decimate(patch, 60000)
    fc = np.tile(np.array(COL["ear"], np.uint8), (len(ear.faces), 1))
    near = np.linalg.norm(ear.triangles_center - ap, axis=1) < CAVUM_R
    fc[near] = COL["cavum"]
    ear.visual.face_colors = fc
    scene.add_geometry(ear, node_name="ear", geom_name="ear")
    print(f"  cavum tint: {int(near.sum())} faces within {CAVUM_R:.0f} mm")

    ring = aperture_ring(patch, ap, cdir, outward)
    print(f"  aperture ring: {len(ring)} boundary points, "
          f"radius {np.linalg.norm(ring - ap, axis=1).mean():.2f} mm mean"
          if len(ring) else "  aperture ring: none found")
    for i, p in enumerate(ring):
        b = trimesh.creation.uv_sphere(radius=0.42, count=[8, 8])
        b.apply_translation(p)
        b.visual.face_colors = COL["mark"]
        scene.add_geometry(b, node_name=f"ring{i}", geom_name=f"ring{i}")
    cen = trimesh.creation.uv_sphere(radius=1.3, count=[16, 16])
    cen.apply_translation(ap)
    cen.visual.face_colors = COL["mark"]
    scene.add_geometry(cen, node_name="aperture", geom_name="aperture")

    sdir = os.path.join(HERE, "stl", "right")
    for nm in ("core", "faceplate", "jacket_wing", "carrier",
               "nozzle_insert_short"):
        m = trimesh.load(os.path.join(sdir, f"{nm}.stl"), force="mesh")
        m = cluster_decimate(m, 15000)
        if nm in ("carrier", "nozzle_insert_short"):
            m.apply_transform(NOZZLE_T)
        m.apply_transform(M)
        m.visual.face_colors = COL[nm]
        scene.add_geometry(m, node_name=nm, geom_name=nm)

    os.makedirs(VIZ, exist_ok=True)
    out = os.path.join(VIZ, f"marked_{a.ear}.glb")
    scene.export(out)
    print(f"  -> {out}  {os.path.getsize(out)/1e6:.2f} MB")
    json.dump(dict(ear=f"{rec['dataset']}/{rec['ear_id']}",
                   seating="BY CONSTRUCTION - no optimiser",
                   aperture=ap.tolist(), canal_dir=cdir.tolist(),
                   canal_dir_source="CONSTRUCTED from anatomy, not probed",
                   canal_rake_deg=CANAL_RAKE,
                   probe_run_mm=run, probe_dir=list(map(float, probed)),
                   canal_off_floor_normal_deg=off_norm, transform=M.tolist(),
                   markers=dict(ring_points=len(ring), ring_colour="magenta",
                                centroid_sphere="magenta",
                                cavum_tint_radius_mm=CAVUM_R),
                   note=("Pose built directly from landmarks, no optimiser. Skirt "
                         "rim centre on the detected aperture centroid; faceplate "
                         "forced outward. The canal axis is CONSTRUCTED as the "
                         "inward concha-floor normal raked 40 deg anteriorly -- "
                         "the scan has no resolvable canal (probe run 1.5 mm), so "
                         "probing it would be fiction. The magenta ring traces "
                         "the real funnel wall by radial ray casting, not a drawn "
                         "circle of assumed size.")),
              open(os.path.join(VIZ, f"marked_{a.ear}_meta.json"), "w"), indent=1)
    return 0


if __name__ == "__main__":
    sys.exit(main())
