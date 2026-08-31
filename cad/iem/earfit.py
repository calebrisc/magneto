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
def _nozzle_frame(cant=None):
    """Nozzle-local -> assembly frame, read straight out of generate.py.

    `cant` overrides `nozzle_cant_deg` (degrees) so a single checkout can score
    the canted default and the collinear baseline side by side.

    Since `nozzle_cant_deg`, the carrier/skirt and the nozzle inserts are modelled
    about the nozzle axis and their STLs are written in that local frame; only the
    assembly applies `nozzle_T`.  Loading carrier.stl and treating its coordinates
    as assembly-frame coordinates would put the seal 45 deg away from where it
    actually is, so the transform is imported rather than copied -- a hardcoded
    cant here would silently rot the moment the generator's default changes.

    Returns (nozzle_T, rim_x, rim_r) with rim_x/rim_r in the NOZZLE-LOCAL frame.
    """
    import generate                                   # same directory
    P = generate.PARAMS
    if cant is not None:
        P = dict(P, nozzle_cant_deg=float(cant))
    g = generate.G(P)
    # generate's skirt_rim_r is the OUTER radius; the contact land's centreline
    # sits half a land-wall inboard of it
    rim_r = float(g.skirt_rim_r) - 0.5 * P["skirt_wall_land"]
    return np.asarray(g.nozzle_T, float), float(g.skirt_rim_x), rim_r


NOZZLE_T, SKIRT_RIM_X, SKIRT_RIM_R = _nozzle_frame()
# the Ø19 rim centre in the ASSEMBLY frame -- what the seating search lands on the
# canal aperture.  At cant = 0 this is just (SKIRT_RIM_X, 0, 0).
RIM_CENTRE = NOZZLE_T[:3, :3] @ np.array([SKIRT_RIM_X, 0.0, 0.0]) + NOZZLE_T[:3, 3]
NOZZLE_AXIS = NOZZLE_T[:3, 0] / np.linalg.norm(NOZZLE_T[:3, 0])
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

    def __init__(self, patch, n=250_000, seed=0):
        # seeded: the report quotes per-ear millimetres, and an unseeded surface
        # sampling moved them by ~0.1 mm between runs
        pts, fid = trimesh.sample.sample_surface(patch, n, seed=seed)
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

def _rim_circle(n=72, nT=None, rim_x=None, rim_r=None):
    """The sealing rim, built in the nozzle-local frame and canted into the
    assembly frame -- so it follows `nozzle_cant_deg` automatically."""
    nT = NOZZLE_T if nT is None else nT
    rim_x = SKIRT_RIM_X if rim_x is None else rim_x
    rim_r = SKIRT_RIM_R if rim_r is None else rim_r
    th = np.linspace(0, 2 * np.pi, n, endpoint=False)
    loc = np.stack([np.full(n, rim_x), rim_r * np.cos(th), rim_r * np.sin(th)],
                   axis=1)
    return loc @ nT[:3, :3].T + nT[:3, 3]


def iem_points(stl_dir=None, seed=0, cant=None):
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
    # this point set's own nozzle frame, so a cant-0 and a cant-45 build can be
    # scored in the same process without either picking up the other's transform
    nT, rim_x, rim_r = ((NOZZLE_T, SKIRT_RIM_X, SKIRT_RIM_R) if cant is None
                        else _nozzle_frame(cant))

    def load(n):
        return trimesh.load(os.path.join(sdir, n), force="mesh")

    core = load("core.stl")
    face = load("faceplate.stl")
    jw = load("jacket_wing.stl")
    # carrier.stl is written in the NOZZLE-LOCAL frame (see _nozzle_frame); the
    # core, faceplate and jacket/wing are already in the assembly frame.
    car = load("carrier.stl")
    car.apply_transform(nT)

    out["rim"] = _rim_circle(nT=nT, rim_x=rim_x, rim_r=rim_r)
    out["_nozzle_T"] = nT
    out["_rim_centre"] = nT[:3, :3] @ np.array([rim_x, 0.0, 0.0]) + nT[:3, 3]

    jv = jw.vertices
    tipidx = np.argsort(jv[:, 1])[-40:]
    out["wing_tip"] = jv[tipidx]
    ymid = 0.5 * (jv[:, 1].max() + 8.0)
    band = np.where(np.abs(jv[:, 1] - ymid) < 0.6)[0]
    out["wing_mid"] = jv[rng.choice(band, min(40, len(band)), replace=False)]

    # jacket ear face: the -Z skin, excluding the wing span
    body = jv[(jv[:, 1] < 8.0) & (jv[:, 2] < -3.0)]
    out["jacket"] = body[rng.choice(len(body), min(140, len(body)), replace=False)]

    s_core = trimesh.sample.sample_surface(core, 260, seed=seed)[0]
    s_face = trimesh.sample.sample_surface(face, 140, seed=seed)[0]
    s_jw = trimesh.sample.sample_surface(jw, 700, seed=seed)[0]
    out["rigid"] = np.vstack([s_core, s_face, s_jw])
    # `shell` is `rigid` minus the wing.  The wing is *meant* to press into the
    # antihelix, so counting it as interference would make every seated pose
    # look like a crash; the shell is not meant to press anywhere.
    out["shell"] = np.vstack([s_core, s_face, s_jw[s_jw[:, 1] < 8.0]])

    cv = car.vertices
    # radius about the NOZZLE axis, measured in the nozzle-local frame -- hypot on
    # world y,z would only be the skirt radius when the cant is zero
    loc = np.einsum("ij,jk->ik", cv - nT[:3, 3], nT[:3, :3])
    r = np.hypot(loc[:, 1], loc[:, 2])
    out["soft"] = cv[rng.choice(len(cv), 400, replace=False)]
    out["skirt"] = cv[r > 6.0]

    fv = face.vertices
    fp = fv[fv[:, 2] > fv[:, 2].max() - 1.5]
    # decimated: protrusion is a max over the outer band, and 300 well-spread
    # points land within ~0.05 mm of the max over all ~22 k.  The full set is
    # transformed thousands of times inside the seating search, where it was
    # the single most expensive term.
    if len(fp) > 300:
        fp = fp[rng.choice(len(fp), 300, replace=False)]
    out["faceplate"] = fp

    # --- separate sets for the contact contract (tryon.py) ---------------- #
    # the contract judges parts individually, so it needs core and faceplate
    # apart from `shell`, and the nozzle insert, which nothing else samples.
    out["core_s"] = s_core
    out["face_s"] = s_face
    try:
        ins = trimesh.load(os.path.join(sdir, "nozzle_insert_short.stl"),
                           force="mesh")
        ins.apply_transform(nT)          # nozzle-local, like the carrier
        out["nozzle"] = trimesh.sample.sample_surface(ins, 400, seed=seed)[0]
    except Exception:                                            # noqa: BLE001
        out["nozzle"] = np.zeros((0, 3))

    # volume-weighted centroid of the worn assembly, for inertial moments in
    # stability.py.  Ti and silicone differ in density, but the silicone carrier
    # is a small fraction of the mass and this is a 3 g load check, not a modal
    # analysis -- volume weighting is close enough and is stated as such.
    try:
        vols = np.array([abs(m.volume) for m in (core, face, jw, car)], float)
        cens = np.array([m.centroid for m in (core, face, jw, car)], float)
        out["_com"] = (vols[:, None] * cens).sum(axis=0) / max(vols.sum(), 1e-9)
    except Exception:                                            # noqa: BLE001
        out["_com"] = np.vstack([core.vertices, face.vertices]).mean(axis=0)

    # cable-exit proxy: centre of the 2-pin socket pocket.  There is no cable,
    # boot or strain relief in the build, so a tug has to be applied somewhere
    # defensible -- this is the connector, and every use of it is flagged.
    try:
        import generate
        g = generate.G(generate.PARAMS if cant is None
                       else dict(generate.PARAMS, nozzle_cant_deg=float(cant)))
        loc = np.array([0.5 * (g.socket_x0 + g.socket_x1), 0.0,
                        generate.PARAMS["socket_z"]])
        out["_cable_exit"] = nT[:3, :3] @ loc + nT[:3, 3]
    except Exception:                                            # noqa: BLE001
        out["_cable_exit"] = out["_com"]

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
