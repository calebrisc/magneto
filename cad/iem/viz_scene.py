#!/usr/bin/env python3
"""
viz_scene.py -- export one seated ear + IEM as a GLB for human eyeball checking.

Every number in docs/TRYON_REPORT.md comes from a signed-distance field and a
pose matrix.  Those can be self-consistently wrong -- the v2 nozzle-frame bug and
the faceplate-into-the-skull roll were both invisible in the metrics until
someone reasoned about the geometry.  This writes the seated scene out so a human
can look at it.

Picks a representative ear (nearest the median protrusion among ears that pass
the seal at the conservative 2.5 mm budget), applies the seating transform to
every part, and writes:

    viz/seated_scene.glb        ear + IEM in the seated pose, parts coloured
    viz/seated_scene_meta.json  ear id, its scores, and the 4x4 seating transform

Parts are placed exactly as the try-on scores them, including the nozzle-local
-> assembly cant for the carrier and nozzle insert -- so if the render looks
wrong, the metrics are wrong too, which is the point.

DECIMATION.  The raw parts are ~1.25 M triangles, far past a 4 MB budget, and no
quadric simplifier is installed.  `cluster_decimate` does vertex clustering:
snap to a grid, weld, drop degenerate and duplicate faces, with the grid size
binary-searched to hit a face target.  Visual quality only -- never use these
meshes for measurement.

Usage:
    python viz_scene.py
    python viz_scene.py --ear sonicom/P0007
"""

from __future__ import annotations

import argparse
import csv
import json
import os
import sys

import numpy as np
import trimesh

from earfit import ALIGNED, HERE, NOZZLE_T, transform

VIZ = os.path.join(HERE, "viz")
NOZZLE_FRAME_PARTS = {"carrier", "nozzle_insert_short"}
COLORS = {
    "ear":                 (205, 175, 135, 255),   # tan
    "core":                (190, 192, 200, 255),   # silver
    "faceplate":           (205, 207, 214, 255),   # silver, brighter
    "jacket_wing":         ( 62, 108, 198, 255),   # blue
    "carrier":             (232, 140,  48, 255),   # orange
    "nozzle_insert_short": (128, 128, 130, 255),   # grey
}


def cluster_decimate(mesh, target_faces):
    """Vertex-clustering decimation to roughly `target_faces`."""
    if len(mesh.faces) <= target_faces:
        return mesh
    lo, hi = 1e-3, float(max(mesh.extents)) / 4.0
    best = mesh
    for _ in range(18):
        g = 0.5 * (lo + hi)
        key = np.floor(mesh.vertices / g).astype(np.int64)
        _, inv = np.unique(key, axis=0, return_inverse=True)
        nv = inv.max() + 1
        verts = np.zeros((nv, 3))
        cnt = np.bincount(inv, minlength=nv).astype(float)
        for k in range(3):
            verts[:, k] = np.bincount(inv, weights=mesh.vertices[:, k], minlength=nv) / cnt
        f = inv[mesh.faces]
        f = f[(f[:, 0] != f[:, 1]) & (f[:, 1] != f[:, 2]) & (f[:, 0] != f[:, 2])]
        f = np.unique(np.sort(f, axis=1), axis=0)
        if len(f) == 0:
            hi = g
            continue
        m = trimesh.Trimesh(vertices=verts, faces=f, process=False)
        if len(f) > target_faces:
            lo = g
        else:
            best, hi = m, g
        if abs(len(f) - target_faces) < 0.05 * target_faces:
            return m
    return best


def pick_ear(rows, seal, want=None):
    ok = {(r["dataset"], r["ear_id"]) for r in seal if r["b2.5_pass_"] == "True"}
    if want:
        ds, eid = want.split("/")
        return next(r for r in rows if r["dataset"] == ds and r["ear_id"] == eid)
    cand = [r for r in rows if (r["dataset"], r["ear_id"]) in ok]
    if not cand:
        cand = rows
    med = np.median([float(r["protrusion"]) for r in rows])
    return min(cand, key=lambda r: abs(float(r["protrusion"]) - med))


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--ear", default=None, help="dataset/ear_id, e.g. sonicom/P0007")
    ap.add_argument("--ear-faces", type=int, default=60000)
    ap.add_argument("--part-faces", type=int, default=18000)
    a = ap.parse_args()

    rows = list(csv.DictReader(open(os.path.join(ALIGNED, "tryon.csv"))))
    seal = list(csv.DictReader(open(os.path.join(ALIGNED, "seal_compliance.csv"))))
    row = pick_ear(rows, seal, a.ear)
    ds, eid = row["dataset"], row["ear_id"]
    rec = json.load(open(os.path.join(ALIGNED, f"{ds}_{eid}_right.json")))
    M = np.array(rec["transform"], float)
    print(f"ear {ds}/{eid}: protrusion {float(row['protrusion']):.2f} mm, "
          f"grade {row['grade']}")

    os.makedirs(VIZ, exist_ok=True)
    scene = trimesh.Scene()

    ear = trimesh.load(os.path.join(ALIGNED, rec["patch"]), force="mesh")
    ear = cluster_decimate(ear, a.ear_faces)
    ear.visual.face_colors = COLORS["ear"]
    scene.add_geometry(ear, node_name="ear", geom_name="ear")

    sdir = os.path.join(HERE, "stl", "right")
    for name in ("core", "faceplate", "jacket_wing", "carrier", "nozzle_insert_short"):
        m = trimesh.load(os.path.join(sdir, f"{name}.stl"), force="mesh")
        m = cluster_decimate(m, a.part_faces)
        if name in NOZZLE_FRAME_PARTS:          # written in the nozzle-local frame
            m.apply_transform(NOZZLE_T)
        m.apply_transform(M)                    # into the ear's scan frame
        m.visual.face_colors = COLORS[name]
        scene.add_geometry(m, node_name=name, geom_name=name)
        print(f"  {name:<20} {len(m.faces):6d} faces")

    glb = os.path.join(VIZ, "seated_scene.glb")
    scene.export(glb)
    mb = os.path.getsize(glb) / 1e6
    print(f"\n{glb}  {mb:.2f} MB")

    meta = dict(
        ear=f"{ds}/{eid}", dataset=ds, ear_id=eid, side=rec["side"],
        source=rec["source"],
        note=("Seated pose exactly as scored in docs/TRYON_REPORT.md. Carrier and "
              "nozzle insert are modelled in the nozzle-local frame and carry "
              "NOZZLE_T before the seating transform. Meshes are decimated by "
              "vertex clustering for viewing only -- do not measure them."),
        scores={k: row[k] for k in (
            "grade", "worst_metric", "rim_cover", "rim_gap", "rim_press",
            "wing_tip", "wing_mid", "jacket_mean", "protrusion", "hard_min",
            "g_seal", "g_retention", "g_clearance", "g_protrusion")},
        seal_compliance={r["b2.5_pass_"]: None for r in [] } or {},
        seating_transform=np.array(rec["transform"], float).tolist(),
        nozzle_local_to_assembly=NOZZLE_T.tolist(),
        landmarks=dict(aperture=rec["aperture"], tragus=rec["tragus"],
                       floor_normal=rec["floor_normal"], outward=rec["outward"]),
        colors={k: list(v) for k, v in COLORS.items()},
    )
    srow = next((r for r in seal if r["dataset"] == ds and r["ear_id"] == eid), None)
    if srow:
        meta["seal_compliance"] = {k: srow[k] for k in srow
                                   if k.startswith(("b1.5", "b2.5", "b4.0"))}
    mp = os.path.join(VIZ, "seated_scene_meta.json")
    with open(mp, "w") as f:
        json.dump(meta, f, indent=1)
    print(mp)
    return 0


if __name__ == "__main__":
    sys.exit(main())
