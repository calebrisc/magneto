#!/usr/bin/env python3
"""
earfit.py -- shared machinery for the Magneto IEM virtual try-on.

Holds the three things align_ear.py and tryon.py both need:

  1. Dataset conventions and mesh loading (HUTUBS is in metres, SONICOM in mm;
     both use the same head frame).
  2. Ear-region extraction + the hull-relative depth map that everything else
     keys off.
  3. A fast signed-distance oracle for the ear surface, and the tagged sample
     points taken off the generated IEM STLs.

HEAD FRAME (both datasets, verified by inspection -- nose at max +x, vertex at
max +z, ears at the |y| extremes):

    +x anterior (nose)      +y toward the LEFT ear      +z superior

so the RIGHT ear -- the one the master STLs are cut for -- lives at -y and its
lateral (outward) direction is -y.

IEM DESIGN FRAME (from generate.py / README.md):

    origin  nozzle base = core face
    +X      nozzle axis, toward the canal
    +Y      superior, toward the antihelix
    +Z      lateral, toward the faceplate
    right ear is the master; skirt rim sits at x = +13.65, Ø19.0
"""

from __future__ import annotations

import os

import numpy as np
import trimesh
from scipy.spatial import cKDTree

HERE = os.path.dirname(os.path.abspath(__file__))
EARS = os.path.join(HERE, "ears")
ALIGNED = os.path.join(EARS, "aligned")

# ---- IEM geometry constants, mirrored from generate.py PARAMS -------------- #
SKIRT_RIM_X = 13.65        # mm, design X of the skirt rim / seal plane
SKIRT_RIM_R = 9.325        # mm, rim centreline radius (19.0/2 - 0.35/2)
CARRIER_X0 = 4.65          # mm, carrier proximal face

DATASETS = {
    "hutubs": dict(scale=1000.0, glob="3D head meshes/*.ply",
                   id_from=lambda p: os.path.basename(p).split("_")[0]),
    "sonicom": dict(scale=1.0, glob="*.stl",
                    id_from=lambda p: os.path.basename(p).split("_")[0]),
    "synthetic": dict(scale=1.0, glob="*.stl",
                      id_from=lambda p: os.path.basename(p).rsplit(".", 1)[0]),
}


# --------------------------------------------------------------------------- #
# ear extraction
# --------------------------------------------------------------------------- #

def load_head(path, scale):
    m = trimesh.load(path, force="mesh")
    if scale != 1.0:
        m.apply_scale(scale)
    return m


def ear_window(head, side="right", pad_post=45.0, pad_ant=40.0,
               pad_inf=45.0, pad_sup=40.0, depth=50.0):
    """Axis-aligned box around one ear.

    Anchored on the single most-lateral vertex of that half of the head, which
    on every scan we have seen is a point on the helix rim.  Failure mode: a
    scan that still has a shoulder or a hand wider than the pinna -- guarded by
    only considering vertices in the upper half of the mesh.
    """
    v = head.vertices
    lat = -1.0 if side == "right" else 1.0            # sign of y that is lateral
    zmid = 0.5 * (head.bounds[0][2] + head.bounds[1][2])
    cand = v[v[:, 2] > zmid - 40.0]
    p = cand[np.argmin(lat * -cand[:, 1])] if lat < 0 else cand[np.argmax(cand[:, 1])]
    if lat < 0:
        p = cand[np.argmin(cand[:, 1])]
    x0, x1 = p[0] - pad_post, p[0] + pad_ant
    z0, z1 = p[2] - pad_inf, p[2] + pad_sup
    if lat < 0:
        y0, y1 = p[1] - 5.0, p[1] + depth
    else:
        y0, y1 = p[1] - depth, p[1] + 5.0
    return dict(anchor=p, x=(x0, x1), y=(y0, y1), z=(z0, z1), lat=lat)


def crop(head, win):
    v = head.vertices
    m = ((v[:, 0] > win["x"][0]) & (v[:, 0] < win["x"][1]) &
         (v[:, 1] > win["y"][0]) & (v[:, 1] < win["y"][1]) &
         (v[:, 2] > win["z"][0]) & (v[:, 2] < win["z"][1]))
    f = np.where(m[head.faces].all(axis=1))[0]
    if len(f) == 0:
        raise RuntimeError("ear window is empty")
    return head.submesh([f], append=True)


def depth_map(patch, win, step=0.5):
    """Hull-relative depth of the ear, rasterised on the (x, z) plane.

    Rays are cast along the medial direction; for each cell we record the first
    hit on the ear and the first hit on the patch's own convex hull, and return
    the difference.  Subtracting the hull is what makes the concha -- and not
    the curvature of the skull behind the ear -- the deepest thing in frame.
    """
    lat = win["lat"]
    x0, x1 = win["x"]
    z0, z1 = win["z"]
    xs = np.arange(x0, x1, step)
    zs = np.arange(z0, z1, step)
    XX, ZZ = np.meshgrid(xs, zs, indexing="ij")
    ystart = win["anchor"][1] + lat * 20.0
    org = np.stack([XX.ravel(), np.full(XX.size, ystart), ZZ.ravel()], axis=1)
    dirs = np.tile([0.0, -lat, 0.0], (len(org), 1))

    def first_hit(mesh):
        loc, ir, _ = mesh.ray.intersects_location(org, dirs, multiple_hits=False)
        d = np.full(XX.size, np.nan)
        d[ir] = loc[:, 1]
        return d.reshape(XX.shape)

    surf = first_hit(patch)
    hull = first_hit(patch.convex_hull)
    rel = (surf - hull) * (-lat)          # positive = deeper than the hull
    return dict(xs=xs, zs=zs, rel=rel, surf=surf, ystart=ystart, lat=lat)


# --------------------------------------------------------------------------- #
# signed-distance oracle
# --------------------------------------------------------------------------- #

class EarField:
    """Nearest-surface signed distance from a dense point sampling.

    We sample the ear patch densely (~0.25 mm) and keep each sample's face
    normal.  distance(p) = |p - q|, sign = sign((p - q) . n_q), so positive is
    outside the flesh and negative is inside it (interference).

    This is ~1000x faster than trimesh's exact signed_distance and is accurate
    to the sampling density away from thin ridges; on a helix rim the sign can
    flip within ~0.3 mm of the surface, which does not matter because every
    decision we make has a threshold well outside that band.
    """

    def __init__(self, patch, n=250_000):
        pts, fid = trimesh.sample.sample_surface(patch, n)
        self.pts = np.asarray(pts)
        self.nrm = patch.face_normals[fid]
        self.tree = cKDTree(self.pts)

    def query(self, p):
        d, i = self.tree.query(p, workers=-1)
        s = np.einsum("ij,ij->i", p - self.pts[i], self.nrm[i])
        return np.where(s < 0, -d, d)

    def unsigned(self, p):
        d, _ = self.tree.query(p, workers=-1)
        return d


# --------------------------------------------------------------------------- #
# IEM sample points
# --------------------------------------------------------------------------- #

def _rim_circle(n=72):
    th = np.linspace(0, 2 * np.pi, n, endpoint=False)
    return np.stack([np.full(n, SKIRT_RIM_X),
                     SKIRT_RIM_R * np.cos(th),
                     SKIRT_RIM_R * np.sin(th)], axis=1)


def iem_points(stl_dir=None, seed=0):
    """Tagged sample points on the IEM, in the design frame.

    rim        72 points on the Ø19 skirt-rim circle (the seal band)
    wing_tip   the 40 most-distal wing vertices
    wing_mid   40 vertices at mid-span of the wing
    jacket     the ear-facing (-Z) skin of the jacket
    rigid      everything printed in Ti -- core, faceplate, jacket + wing.
               This is the set that cannot deform, so it drives the hard
               interference test.
    shell      rigid minus the wing: the surfaces that are not supposed to
               press into anything, which is what the seating search minimises
    soft       the silicone carrier body and skirt
    faceplate  the outermost +Z shell points, for the protrusion measurement
    """
    sdir = stl_dir or os.path.join(HERE, "stl", "right")
    rng = np.random.default_rng(seed)
    out = {}

    def load(n):
        return trimesh.load(os.path.join(sdir, n), force="mesh")

    core = load("core.stl")
    face = load("faceplate.stl")
    jw = load("jacket_wing.stl")
    car = load("carrier.stl")

    out["rim"] = _rim_circle()

    jv = jw.vertices
    tipidx = np.argsort(jv[:, 1])[-40:]
    out["wing_tip"] = jv[tipidx]
    ymid = 0.5 * (jv[:, 1].max() + 8.0)
    band = np.where(np.abs(jv[:, 1] - ymid) < 0.6)[0]
    out["wing_mid"] = jv[rng.choice(band, min(40, len(band)), replace=False)]

    # jacket ear face: the -Z skin, excluding the wing span
    body = jv[(jv[:, 1] < 8.0) & (jv[:, 2] < -3.0)]
    out["jacket"] = body[rng.choice(len(body), min(140, len(body)), replace=False)]

    s_core = trimesh.sample.sample_surface(core, 260)[0]
    s_face = trimesh.sample.sample_surface(face, 140)[0]
    s_jw = trimesh.sample.sample_surface(jw, 700)[0]
    out["rigid"] = np.vstack([s_core, s_face, s_jw])
    # `shell` is `rigid` minus the wing.  The wing is *meant* to press into the
    # antihelix, so counting it as interference would make every seated pose
    # look like a crash; the shell is not meant to press anywhere.
    out["shell"] = np.vstack([s_core, s_face, s_jw[s_jw[:, 1] < 8.0]])

    cv = car.vertices
    r = np.hypot(cv[:, 1], cv[:, 2])
    out["soft"] = cv[rng.choice(len(cv), 400, replace=False)]
    out["skirt"] = cv[r > 6.0]

    fv = face.vertices
    out["faceplate"] = fv[fv[:, 2] > fv[:, 2].max() - 1.5]

    out["_bbox"] = np.vstack([core.bounds, face.bounds, jw.bounds, car.bounds])
    return out


def transform(pts, M):
    with np.errstate(over="ignore", invalid="ignore", divide="ignore"):
        return pts @ M[:3, :3].T + M[:3, 3]


def rt_matrix(t, rot_deg, base):
    """base @ (rotation about the design origin, then translation), all in the
    ear frame: M = base @ T(t) @ R."""
    t = np.clip(np.asarray(t, float), -1e3, 1e3)
    rx, ry, rz = np.radians(np.clip(np.asarray(rot_deg, float), -720, 720))
    R = trimesh.transformations.euler_matrix(rx, ry, rz, "rxyz")
    R[:3, 3] = t
    return base @ R
