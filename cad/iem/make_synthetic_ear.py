#!/usr/bin/env python3
"""
make_synthetic_ear.py -- parametric ears at the corners of the design envelope.

The scan datasets tell us what the middle of the population looks like; they do
not reliably contain the tails.  `docs/EAR_ANTHROPOMETRY.md` gives the envelope
the one-size IEM has to survive:

    cavum concha length  11-23 mm      aperture height  7-18 mm
    concha width         10-24 mm      aperture width   4.5-14 mm
    concha depth          8-18 mm

so we build the four corners of the two dimensions that actually load the part
-- **aperture size** (what the Ø19 skirt has to seal against) and **concha
depth** (what the 13.65 mm protrusion has to fit inside) -- and let the concha
plan dimensions track the depth, because shallow conchas are also small ones.

    xs_shallow   aperture  7.0 x 4.5   concha 11 x 10 x  8   (5th pct-ish)
    xs_deep      aperture  7.0 x 4.5   concha 23 x 24 x 18
    xl_shallow   aperture 18.0 x 14.0  concha 11 x 10 x  8
    xl_deep      aperture 18.0 x 14.0  concha 23 x 24 x 18   (95th pct-ish)

Each ear is a signed-distance field polygonised with marching cubes, written in
the same head frame the scan datasets use (+x anterior, +y left, +z superior,
mm), so align_ear.py and tryon.py treat them exactly like a real scan.  The
model is:

  * a 78 mm sphere of "head" whose lateral pole sits just outboard of the
    concha rim (the right ear's lateral direction is -y);
  * a **concha bowl** -- half-ellipsoid of L x W x D carved into it;
  * a **canal funnel** running anteromedially out of the bowl floor at 40
    degrees to the floor normal (the real angle; this is the thing the IEM's
    orthogonal design frame does not have), tapering over 3 mm from twice the
    aperture down to the aperture ellipse and then running 12 mm at constant
    section;
  * an **antihelix ridge** -- a raised torus segment along the posterior and
    superior rim of the bowl, which is what the wing hooks under;
  * a **tragus** lobe anterior of the aperture, which is what the faceplate
    protrusion is measured against;
  * a **helix** rim outboard of all of it, so the hull-relative depth map in
    align_ear.py sees the same kind of picture it sees on a real ear.

Ground-truth landmarks go into ears/synthetic/<name>.landmarks.json, so
align_ear.py's detector can be scored against them.

Usage:
    python make_synthetic_ear.py            # all four, into ears/synthetic/
    python make_synthetic_ear.py --voxel 0.3
"""

from __future__ import annotations

import argparse
import json
import math
import os
import sys

import numpy as np

try:
    from skimage import measure
except ImportError:  # pragma: no cover
    sys.exit("scikit-image missing -- activate cad/iem/.venv")

import trimesh

from earfit import EARS

CORNERS = {
    "xs_shallow": dict(ap=(7.0, 4.5), concha=(11.0, 10.0, 8.0)),
    "xs_deep":    dict(ap=(7.0, 4.5), concha=(23.0, 24.0, 18.0)),
    "xl_shallow": dict(ap=(18.0, 14.0), concha=(11.0, 10.0, 8.0)),
    "xl_deep":    dict(ap=(18.0, 14.0), concha=(23.0, 24.0, 18.0)),
}

CANAL_RAKE_DEG = 40.0      # canal axis vs concha-floor normal
# Only ~4 mm of canal, on purpose.  A structured-light scanner cannot see round
# the first bend, so on SONICOM and HUTUBS the canal is a funnel plus a short
# stub (our canal probe measures 0.6-3 mm of free run on real ears).  Modelling
# a full 25 mm canal here would make the synthetic ears *easier* to landmark in
# a way no scan is, and would put the depth map's deepest point 10 mm down the
# bore instead of at the mouth where the skirt has to seal.
CANAL_LEN = 4.0            # mm of constant-section canal past the funnel


# --------------------------------------------------------------------------- #
# SDF primitives -- local frame u=anterior, v=superior, w=lateral(outward)
# --------------------------------------------------------------------------- #

def _ellipsoid(U, V, W, c, r):
    return (np.sqrt(((U - c[0]) / r[0]) ** 2 + ((V - c[1]) / r[1]) ** 2
                    + ((W - c[2]) / r[2]) ** 2) - 1.0) * min(r)


def _smax(a, b, k):
    return 0.5 * (a + b + np.sqrt((a - b) ** 2 + k * k))


def _smin(a, b, k):
    return 0.5 * (a + b - np.sqrt((a - b) ** 2 + k * k))


def _elliptic_tube(U, V, W, origin, axis, ra, rb, half_up, t0, t1, flare):
    """Elliptic cylinder along `axis`, radii (ra, rb) scaled by `flare` at t0."""
    ax = axis / np.linalg.norm(axis)
    up = half_up - (half_up @ ax) * ax
    up /= np.linalg.norm(up)
    side = np.cross(ax, up)
    du, dv, dw = U - origin[0], V - origin[1], W - origin[2]
    t = du * ax[0] + dv * ax[1] + dw * ax[2]
    a = du * side[0] + dv * side[1] + dw * side[2]
    b = du * up[0] + dv * up[1] + dw * up[2]
    s = np.clip((t - t0) / max(t1 - t0, 1e-6), 0.0, 1.0)
    k = flare + (1.0 - flare) * np.clip(s * (t1 - t0) / 3.0, 0.0, 1.0)
    rad = np.sqrt((a / (ra * k)) ** 2 + (b / (rb * k)) ** 2) - 1.0
    rad = rad * min(ra, rb)
    cap = np.maximum(t0 - t, t - t1)
    return _smax(rad, cap, 0.4)


def _torus_seg(U, V, W, c, R, r, keep):
    """Torus in the (u, v) plane at w = c[2]; `keep` masks the wanted arc."""
    du, dv, dw = U - c[0], V - c[1], W - c[2]
    q = np.sqrt(du ** 2 + dv ** 2) - R
    d = np.sqrt(q ** 2 + dw ** 2) - r
    return np.where(keep, d, 1e3)


def build_field(cfg, grid):
    U, V, W = grid
    L, Wd, D = cfg["concha"]
    aph, apw = cfg["ap"]

    # head: a 78 mm sphere whose lateral pole sits at w = +2, so the ear
    # window / convex-hull logic in align_ear.py sees the same kind of skull
    # curvature it sees on a real scan
    head = np.sqrt(U ** 2 + V ** 2 + (W + 76.0) ** 2) - 78.0

    # concha bowl carved into it, its floor sloping down to the canal mouth --
    # the deepest point of a real cavum *is* the canal mouth, and a symmetric
    # bowl with an off-centre hole would hand the detector an error that no
    # real ear has
    bowl_u = 0.18 * L
    bowl = _ellipsoid(U, V, W, (bowl_u, 0.0, 0.0), (0.5 * L, 0.5 * Wd, D))
    solid = _smax(head, -bowl, 1.2)

    # antihelix ridge: raised arc along the posterior + superior rim
    ang = np.arctan2(V, U - bowl_u)
    keep = (ang > math.radians(35.0)) & (ang < math.radians(215.0))
    ridge = _torus_seg(U, V, W, (bowl_u, 0.0, 1.0),
                       0.5 * max(L, Wd) + 1.8, 3.2, keep)
    solid = _smin(solid, ridge, 1.4)

    # helix: a bigger outboard arc so the depth map sees a pinna outline
    keep_h = (ang > math.radians(15.0)) & (ang < math.radians(250.0))
    helix = _torus_seg(U, V, W, (bowl_u, 0.0, 2.0),
                       0.5 * max(L, Wd) + 9.0, 4.0, keep_h)
    solid = _smin(solid, helix, 1.6)

    # tragus: a lobe just anterior of the aperture
    trg = _ellipsoid(U, V, W, (0.5 * L + bowl_u + 2.0, -2.0, 0.5), (3.2, 5.0, 4.5))
    solid = _smin(solid, trg, 1.0)

    # canal: bored out of the bowl floor, raked anteromedially
    th = math.radians(CANAL_RAKE_DEG)
    axis = np.array([math.sin(th), 0.0, -math.cos(th)])     # into the head
    apc = np.array([bowl_u, 0.0, -D + 0.6])                 # bowl floor = canal mouth
    canal = _elliptic_tube(U, V, W, apc, axis,
                           0.5 * apw, 0.5 * aph, np.array([0.0, 1.0, 0.0]),
                           -2.0, CANAL_LEN, flare=2.0)
    solid = _smax(solid, -canal, 0.5)
    return solid, apc, axis


def polygonise(field, origin, spacing):
    v, f, _, _ = measure.marching_cubes(field, level=0.0, spacing=(spacing,) * 3)
    return trimesh.Trimesh(vertices=v + origin, faces=f, process=True)


def make(name, cfg, voxel):
    L, Wd, D = cfg["concha"]
    pad = 38.0
    bu = 0.18 * L
    lo = np.array([-0.5 * L - pad + bu, -0.5 * Wd - pad, -D - CANAL_LEN - 14.0])
    hi = np.array([0.5 * L + pad + bu, 0.5 * Wd + pad, 12.0])
    n = np.ceil((hi - lo) / voxel).astype(int) + 1
    ax = [lo[i] + voxel * np.arange(n[i]) for i in range(3)]
    U = ax[0][:, None, None]; V = ax[1][None, :, None]; W = ax[2][None, None, :]
    fld, apc, axis = build_field(cfg, (U, V, W))
    m = polygonise(np.ascontiguousarray(fld), lo, voxel)
    m.fix_normals()

    # local (u,v,w) -> head frame (x anterior, y left, z superior) for a RIGHT
    # ear, whose lateral direction is -y:  x=u, y=-w, z=v
    M = np.array([[1.0, 0.0, 0.0, 0.0],
                  [0.0, 0.0, -1.0, 0.0],
                  [0.0, 1.0, 0.0, 0.0],
                  [0.0, 0.0, 0.0, 1.0]])
    m.apply_transform(M)
    ap_head = (M[:3, :3] @ apc).tolist()
    ax_head = (M[:3, :3] @ axis).tolist()
    return m, ap_head, ax_head


def main():
    p = argparse.ArgumentParser(description=__doc__,
                                formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--voxel", type=float, default=0.40)
    p.add_argument("--out", default=os.path.join(EARS, "synthetic"))
    a = p.parse_args()
    os.makedirs(a.out, exist_ok=True)
    truth = {}
    for name, cfg in CORNERS.items():
        m, ap, ax = make(name, cfg, a.voxel)
        path = os.path.join(a.out, f"{name}.stl")
        m.export(path)
        truth[name] = dict(aperture=ap, canal_axis=ax,
                           aperture_hw=cfg["ap"], concha=cfg["concha"])
        print(f"{name:12s} {len(m.faces):7d} faces  watertight={m.is_watertight}  "
              f"bbox {np.round(m.extents, 1)}  -> {path}")
    with open(os.path.join(a.out, "landmarks_truth.json"), "w") as f:
        json.dump(truth, f, indent=1)
    print(f"ground truth -> {os.path.join(a.out, 'landmarks_truth.json')}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
