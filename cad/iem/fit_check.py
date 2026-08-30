#!/usr/bin/env python3
"""
fit_check.py -- measure the Magneto IEM against a scanned ear mesh.

Loads an ear-surface STL (SONICOM or HUTUBS export, or any watertight/open ear
mesh in mm) plus the generated jacket_wing.stl and carrier.stl, and reports how
far the two contact features sit from the ear surface:

  * wing tip   -- the single most distal point of the wing (max +Y), plus the
                  20 vertices nearest to it, so a single spike cannot dominate.
  * skirt rim  -- 72 points sampled around the sealing lip of the mag-float
                  carrier's skirt (the ring of maximum radius about the X axis).

For each it prints min / mean / max distance to the nearest ear surface and, if
the ear mesh is watertight, the signed distance (negative = the part is INSIDE
the ear surface, i.e. interference).

ALIGNMENT IS THE HARD PART AND IS NOT AUTOMATED.  The IEM lives in the design
frame (origin = nozzle base, +X into the ear, +Z out of the head, +Y superior);
an ear scan lives in whatever frame the capture rig used.  Supply the rigid
transform that maps design frame -> scan frame with --transform, either as
6 numbers "tx,ty,tz,rx,ry,rz" (mm and degrees, XYZ intrinsic Euler, rotation
applied first) or as a path to a text file holding a 4x4 row-major matrix.

To obtain that transform in practice:
  1. In MeshLab/Blender pick three landmarks on the scan -- the centre of the
     canal aperture, the deepest point of the cavum concha, and the antihelix
     crest directly above the aperture.
  2. The aperture centre is the design origin; aperture-centre -> concha-floor
     normal is +X; aperture-centre -> antihelix crest is +Y.
  3. Build the transform from that frame and pass it here.

Usage:
    python fit_check.py --ear /path/to/subject_042_ear.stl
    python fit_check.py --ear ear.stl --transform "0,0,0,0,-25,10"
    python fit_check.py --ear ear.stl --transform align.txt --stl-dir stl/right
"""

from __future__ import annotations

import argparse
import math
import os
import sys

import numpy as np

try:
    import trimesh
    from trimesh.proximity import ProximityQuery
except ImportError:  # pragma: no cover
    sys.exit("trimesh missing -- activate cad/iem/.venv")


def parse_transform(spec):
    if spec is None:
        return np.eye(4)
    if os.path.exists(spec):
        M = np.loadtxt(spec)
        if M.shape != (4, 4):
            sys.exit(f"{spec}: expected a 4x4 matrix, got {M.shape}")
        return M
    vals = [float(v) for v in spec.replace(";", ",").split(",")]
    if len(vals) != 6:
        sys.exit("--transform needs 6 numbers 'tx,ty,tz,rx,ry,rz' or a 4x4 matrix file")
    tx, ty, tz, rx, ry, rz = vals
    M = trimesh.transformations.euler_matrix(
        math.radians(rx), math.radians(ry), math.radians(rz), "rxyz")
    M[:3, 3] = (tx, ty, tz)
    return M


def wing_tip_points(mesh, n=20):
    """The most distal wing vertices (largest +Y)."""
    v = mesh.vertices
    idx = np.argsort(v[:, 1])[-n:]
    return v[idx]


def skirt_rim_points(mesh, n=72, band=0.35):
    """Points around the sealing lip: max radius about the X axis."""
    v = mesh.vertices
    r = np.hypot(v[:, 1], v[:, 2])
    rmax = r.max()
    sel = v[r > rmax - band]
    if len(sel) < n:
        return sel
    th = np.arctan2(sel[:, 2], sel[:, 1])
    out = []
    edges = np.linspace(-np.pi, np.pi, n + 1)
    for i in range(n):
        m = (th >= edges[i]) & (th < edges[i + 1])
        if m.any():
            out.append(sel[m][np.argmax(np.hypot(sel[m][:, 1], sel[m][:, 2]))])
    return np.asarray(out)


def report(label, pts, ear, pq):
    if len(pts) == 0:
        print(f"  {label:12s} -- no points found")
        return
    _, dist, _ = pq.on_surface(pts)
    line = (f"  {label:12s} n={len(pts):3d}  nearest-surface distance  "
            f"min {dist.min():6.2f}  mean {dist.mean():6.2f}  max {dist.max():6.2f} mm")
    if ear.is_watertight:
        sd = pq.signed_distance(pts)          # >0 inside the closed surface
        pen = sd[sd > 0]
        line += (f"  |  inside-surface points {len(pen):3d}"
                 + (f", worst interference {pen.max():.2f} mm" if len(pen) else ""))
    print(line)


def main():
    ap = argparse.ArgumentParser(description="Magneto IEM fit check against an ear mesh")
    ap.add_argument("--ear", required=True, help="ear surface STL/PLY/OBJ (mm)")
    ap.add_argument("--stl-dir", default=None, help="directory holding the IEM STLs")
    ap.add_argument("--transform", default=None,
                    help="'tx,ty,tz,rx,ry,rz' (mm, deg) or a 4x4 matrix file")
    ap.add_argument("--scale", type=float, default=1.0,
                    help="scale factor applied to the ear mesh (use 1000 for metres)")
    args = ap.parse_args()

    here = os.path.dirname(os.path.abspath(__file__))
    sdir = args.stl_dir or os.path.join(here, "stl", "right")

    ear = trimesh.load(args.ear, force="mesh")
    if args.scale != 1.0:
        ear.apply_scale(args.scale)
    print(f"ear mesh   {args.ear}")
    print(f"           {len(ear.faces)} faces, watertight={ear.is_watertight}, "
          f"bbox {np.round(ear.extents, 1)} mm")

    M = parse_transform(args.transform)
    if args.transform is None:
        print("           NOTE: no --transform given; the IEM is being compared in the "
              "raw design frame.  Numbers are meaningless unless the scan happens to "
              "share that frame.")

    pq = ProximityQuery(ear)

    for name, extract, label in (("jacket_wing.stl", wing_tip_points, "wing tip"),
                                 ("carrier.stl", skirt_rim_points, "skirt rim")):
        path = os.path.join(sdir, name)
        if not os.path.exists(path):
            print(f"  missing {path} -- run generate.py --all first")
            continue
        m = trimesh.load(path, force="mesh")
        pts = extract(m)
        pts = trimesh.transform_points(pts, M)
        report(label, pts, ear, pq)

    print("\nInterpretation: a well-fitting wing tip sits 0.0-0.5 mm off the antihelix "
          "(light contact, no preload) and the skirt rim wants 0.2-1.0 mm of "
          "interference against the concha floor so the mag-float preload has "
          "something to push into.")
    return 0


if __name__ == "__main__":
    sys.exit(main())
