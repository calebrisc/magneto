#!/usr/bin/env python3
"""
Magneto IEM -- parametric STL generator (v0).

Everything is built as a signed-distance field, sampled on a regular grid and
polygonised with marching cubes.  That is slower than a B-rep kernel but it is
watertight by construction and it lets the gyroid lattice and the solid bodies
live in the same expression tree.

Coordinate convention (see README):
    origin  = nozzle base centre (the point where the nozzle stub leaves the core)
    +X      = nozzle axis, pointing INTO the ear
    +Z      = toward the faceplate (outward from the head)
    +Y      = up, toward the antihelix / cymba concha
    units   = mm, right-handed, RIGHT ear is the master; left is mirrored in X.

Usage:
    python generate.py --all
    python generate.py --part core --part jacket_wing
    python generate.py --list
    python generate.py --all --magnet-preset n52_long
    python generate.py --part jacket_wing --voxel 0.07     # fine, slow
"""

from __future__ import annotations

import argparse
import copy
import math
import os
import sys
import time

import numpy as np

try:
    from skimage import measure
except ImportError:  # pragma: no cover
    sys.exit("scikit-image missing -- activate cad/iem/.venv")

try:
    import trimesh
except ImportError:  # pragma: no cover
    sys.exit("trimesh missing -- activate cad/iem/.venv")


# --------------------------------------------------------------------------
# PARAMETERS
# --------------------------------------------------------------------------

MAGNET_PRESETS = {
    # name             od    id    t     rest_gap  material          F@rest  ratio
    "bonded_compact": dict(od=7.0, id=3.0, t=1.5, gap=2.25,
                           material="bonded NdFeB (Br 0.65 T)",
                           f_rest=0.254, f_lo=0.433, f_hi=0.165, ratio=2.62,
                           mass_g=0.283, stray_mT_at_10mm=3.30),
    "n52_long":       dict(od=7.0, id=3.0, t=1.5, gap=6.40,
                           material="N52 sintered NdFeB (Br 1.45 T)",
                           f_rest=0.200, f_lo=0.263, f_hi=0.155, ratio=1.69,
                           mass_g=0.353, stray_mT_at_10mm=None),
    "n52_clean_bore": dict(od=8.0, id=4.0, t=1.5, gap=6.80,
                           material="N52 sintered NdFeB (Br 1.45 T)",
                           f_rest=0.203, f_lo=0.259, f_hi=0.161, ratio=1.61,
                           mass_g=0.424, stray_mT_at_10mm=8.26),
    "n52_small_pkg":  dict(od=6.0, id=3.0, t=1.0, gap=3.30,
                           material="N52 sintered NdFeB (Br 1.45 T)",
                           f_rest=0.226, f_lo=0.359, f_hi=0.153, ratio=2.35,
                           mass_g=0.159, stray_mT_at_10mm=None),
}


PARAMS = dict(
    # ---- process / manufacturing limits -------------------------------
    min_wall=0.20,             # mm, thinnest printable lattice wall (LPBF Ti / resin)
    min_cell=1.00,             # mm, smallest usable gyroid unit cell
    clearance=0.15,            # mm, jacket-to-core sliding clearance
    press_clearance=0.15,      # mm, press/slip fit clearance on cylindrical joints

    # ---- magnets ------------------------------------------------------
    magnet_preset="bonded_compact",
    magnet_pocket_clear=0.05,  # mm added to magnet OD/thickness for the pocket

    # ---- core shell ---------------------------------------------------
    core_rx=8.5,               # mm half-extent along X (nozzle axis)
    core_ry=7.0,               # mm half-extent along Y (up)
    core_rz=5.0,               # mm half-extent along Z (outward)
    core_wall=1.20,            # mm shell wall
    faceplate_z=1.00,          # mm, Z of the core/faceplate parting plane
    cavity_cap_z=-0.60,        # mm, Z where the acoustic cavity is capped off
    nose_cone_r0=5.0,          # mm nose-cone radius at its rear
    nose_cone_r1=3.2,          # mm nose-cone radius at the core face (x=0)
    nose_cone_x0=-6.0,         # mm nose-cone rear station

    # ---- driver -------------------------------------------------------
    driver_dia=10.0,           # mm dynamic driver
    driver_carrier_od=12.0,    # mm press-fit carrier ring OD
    driver_carrier_id=10.0,    # mm carrier ring ID (driver slip fit)
    driver_carrier_h=3.0,      # mm carrier ring height
    driver_pocket_clear=0.20,  # mm on the carrier OD -> core pocket dia
    driver_pocket_depth=3.0,   # mm, from the faceplate parting plane, -Z

    # ---- nozzle -------------------------------------------------------
    nozzle_bore=4.00,          # mm acoustic bore
    stub_od=5.00,              # mm core nozzle stub OD
    stub_len=3.00,             # mm core nozzle stub length (+X from origin)
    lug_h=0.60,                # mm bayonet lug radial height
    lug_w=1.50,                # mm bayonet lug axial width
    socket_od=8.00,            # mm nozzle-insert socket (slips over the stub) OD
    insert_od=5.00,            # mm nozzle-insert tube OD
    insert_tube_lengths=dict(short=6.0, med=8.0, long=10.0),  # mm beyond the magnet flange
    insert_lug_x0=9.00,        # mm, carrier bayonet lug station (start)
    damper_dia=4.00,           # mm damper disc
    damper_recess=0.30,        # mm damper disc recess depth

    # ---- mag-float carrier --------------------------------------------
    carrier_bore=5.20,         # mm sliding bore on the insert tube
    carrier_od=8.00,           # mm carrier body OD
    carrier_len=8.00,          # mm carrier body length
    carrier_travel=1.50,       # mm allowed axial float
    skirt_flare_deg=35.0,      # deg half-angle of the sealing skirt
    skirt_wall=0.35,           # mm skirt wall
    skirt_max_dia=19.0,        # mm skirt rim diameter (outer)

    # ---- jacket + wing -------------------------------------------------
    jacket_thick=1.60,         # mm total jacket thickness (lattice + skin)
    jacket_x_clip=-2.00,       # mm, jacket stops here so it never fouls the nozzle
    gyroid_cell=1.20,          # mm gyroid unit cell
    wall_face=0.20,            # mm gyroid wall at the outer / ear face
    wall_root=0.40,            # mm gyroid wall at the root
    grade_len=6.00,            # mm over which the wall grades face->root
    skin_t=0.60,               # mm solid skin membrane on the ear-facing surface
    perf_dia=0.40,             # mm sweat perforation diameter
    perf_pitch=1.50,           # mm perforation grid pitch
    solid_root=1.00,           # mm of solid Ti before the lattice starts
    wing_len=14.0,             # mm wing length along its centreline
    wing_thick=4.00,           # mm wing thickness (in the XY plane, across the blade)
    wing_width=7.00,           # mm wing width (in Z, into the concha)
    wing_z_top=-0.20,          # mm top of the wing (just under the parting plane)
    wing_taper_deg=40.0,       # deg overhang of the wing's deep-edge taper (<45 = self-supporting)
    wing_taper=2.60,           # mm of tapered depth on the wing's deep edge
    wing_rise=10.0,            # mm the tip lands above the core rim
    wing_back_deg=30.0,        # deg the tip is angled toward -X
    wing_root_dx=1.00,         # mm, wing root offset from the core centre in X

    # ---- jacket/core interface -----------------------------------------
    gasket_w=0.60,             # mm gasket groove width
    gasket_d=0.40,             # mm gasket groove depth
    jmag_dia=2.00,             # mm jacket magnet diameter
    jmag_depth=1.00,           # mm jacket magnet pocket depth
    pin_dia=1.00,              # mm locating pin diameter
    pin_depth=1.50,            # mm locating pin hole depth

    # ---- electronics ----------------------------------------------------
    socket_w=5.60,             # mm 2-pin socket pocket width  (Y)
    socket_h=2.80,             # mm 2-pin socket pocket height (Z)
    socket_d=6.00,             # mm 2-pin socket pocket depth  (X)
    socket_z=1.20,             # mm socket pocket centre Z (toward the top)
    bone_w=4.00,               # mm bone-conduction sensor pocket (X)
    bone_h=3.00,               # mm bone sensor pocket (Z)
    bone_d=1.50,               # mm bone sensor pocket depth (into -Y)
    wire_dia=1.00,             # mm wire channel diameter
    vent_dia=0.80,             # mm front/rear vent diameter

    # ---- mould ----------------------------------------------------------
    mold_wall=4.5,             # mm mould block wall around the cavity
    mold_pin_dia=3.0,          # mm mould alignment pin diameter
    mold_pin_len=3.0,          # mm mould alignment pin length
    mold_spout_dia=3.0,        # mm pour spout diameter
    mold_vent_dia=1.2,         # mm vent diameter

    # ---- meshing --------------------------------------------------------
    voxel=None,                # mm; None -> derived from the voxel budget
    budget=3.0e6,              # voxels for solid parts
    budget_lattice=6.0e6,      # voxels for the lattice part
)


# --------------------------------------------------------------------------
# SDF PRIMITIVES  (arrays broadcast: X is (nx,1,1), Y is (1,ny,1), Z is (1,1,nz))
# --------------------------------------------------------------------------

def U(*ds):
    """Union."""
    out = ds[0]
    for d in ds[1:]:
        out = np.minimum(out, d)
    return out


def I(*ds):
    """Intersection."""
    out = ds[0]
    for d in ds[1:]:
        out = np.maximum(out, d)
    return out


def S(a, b):
    """Subtract b from a."""
    return np.maximum(a, -b)


def smin(a, b, k):
    """Polynomial smooth union (adds a fillet of roughly k)."""
    h = np.clip(0.5 + 0.5 * (b - a) / k, 0.0, 1.0)
    return b * (1 - h) + a * h - k * h * (1 - h)


def ellipsoid(X, Y, Z, c, r):
    px, py, pz = X - c[0], Y - c[1], Z - c[2]
    k0 = np.sqrt((px / r[0]) ** 2 + (py / r[1]) ** 2 + (pz / r[2]) ** 2)
    k1 = np.sqrt((px / r[0] ** 2) ** 2 + (py / r[1] ** 2) ** 2 + (pz / r[2] ** 2) ** 2)
    return k0 * (k0 - 1.0) / np.maximum(k1, 1e-9)


def box(X, Y, Z, c, h):
    qx = np.abs(X - c[0]) - h[0]
    qy = np.abs(Y - c[1]) - h[1]
    qz = np.abs(Z - c[2]) - h[2]
    outside = np.sqrt(np.maximum(qx, 0) ** 2 + np.maximum(qy, 0) ** 2
                      + np.maximum(qz, 0) ** 2)
    inside = np.minimum(np.maximum(np.maximum(qx, qy), qz), 0.0)
    return outside + inside


def rbox(X, Y, Z, c, h, r):
    """Rounded box."""
    return box(X, Y, Z, c, (h[0] - r, h[1] - r, h[2] - r)) - r


def slab(A, lo, hi):
    """1-D slab lo <= A <= hi."""
    return np.maximum(lo - A, A - hi)


def cyl_x(X, Y, Z, cy, cz, r, x0, x1):
    rr = np.sqrt((Y - cy) ** 2 + (Z - cz) ** 2)
    return np.maximum(rr - r, slab(X, x0, x1))


def cyl_z(X, Y, Z, cx, cy, r, z0, z1):
    rr = np.sqrt((X - cx) ** 2 + (Y - cy) ** 2)
    return np.maximum(rr - r, slab(Z, z0, z1))


def cyl_y(X, Y, Z, cx, cz, r, y0, y1):
    rr = np.sqrt((X - cx) ** 2 + (Z - cz) ** 2)
    return np.maximum(rr - r, slab(Y, y0, y1))


def tube_x(X, Y, Z, cy, cz, ri, ro, x0, x1):
    rr = np.sqrt((Y - cy) ** 2 + (Z - cz) ** 2)
    return np.maximum(np.maximum(ri - rr, rr - ro), slab(X, x0, x1))


def cone_x(X, Y, Z, x0, r0, x1, r1):
    """Truncated cone about the X axis, radius r0 at x0 -> r1 at x1."""
    t = np.clip((X - x0) / (x1 - x0), 0.0, 1.0)
    r = r0 + (r1 - r0) * t
    rr = np.sqrt(Y ** 2 + Z ** 2)
    return np.maximum(rr - r, slab(X, min(x0, x1), max(x0, x1)))


def capsule(X, Y, Z, a, b, r):
    ax, ay, az = a
    bx, by, bz = b
    dx, dy, dz = bx - ax, by - ay, bz - az
    dd = dx * dx + dy * dy + dz * dz
    px, py, pz = X - ax, Y - ay, Z - az
    h = np.clip((px * dx + py * dy + pz * dz) / dd, 0.0, 1.0)
    return np.sqrt((px - dx * h) ** 2 + (py - dy * h) ** 2 + (pz - dz * h) ** 2) - r


def _wrap(a):
    return (a + np.pi) % (2 * np.pi) - np.pi


def arc_slot_x(X, Y, Z, ri, ro, x0, x1, th0, th1):
    """Annular sector about the X axis, angles in radians measured from +Y."""
    rr = np.sqrt(Y ** 2 + Z ** 2)
    th = np.arctan2(Z, Y)
    thc = 0.5 * (th0 + th1)
    half = 0.5 * (th1 - th0)
    dth = (np.abs(_wrap(th - thc)) - half) * np.maximum(rr, 0.2)
    return np.maximum(np.maximum(np.maximum(ri - rr, rr - ro), slab(X, x0, x1)), dth)


def revolve_segment(X, Y, Z, a, b, r):
    """Capsule in the (x, rho) meridian half-plane, revolved about the X axis.

    a and b are (x, rho) pairs; the result is a rounded-ended conical shell of
    radius r.  Used for the sealing skirt.
    """
    rho = np.sqrt(Y ** 2 + Z ** 2)
    ax, ar = a
    bx, br = b
    dx, dr = bx - ax, br - ar
    dd = dx * dx + dr * dr
    px, pr = X - ax, rho - ar
    h = np.clip((px * dx + pr * dr) / dd, 0.0, 1.0)
    return np.sqrt((px - dx * h) ** 2 + (pr - dr * h) ** 2) - r


def gyroid(X, Y, Z, cell, wall):
    """Sheet gyroid with an (approximately) metric distance field.

    g = sin x cos y + sin y cos z + sin z cos x with k = 2 pi / cell.
    RMS |grad g| = k * sqrt(3/2), so the |g| < t band has width 2t/(1.2247 k);
    inverting gives t = wall * 2 pi * 1.2247 / (2 * cell) = wall * 3.848 / cell.
    """
    k = 2.0 * np.pi / cell
    sx, cx_ = np.sin(k * X), np.cos(k * X)
    sy, cy_ = np.sin(k * Y), np.cos(k * Y)
    sz, cz_ = np.sin(k * Z), np.cos(k * Z)
    g = sx * cy_ + sy * cz_ + sz * cx_
    t = wall * 3.848 / cell
    return (np.abs(g) - t) * cell / 7.695


# --------------------------------------------------------------------------
# DERIVED GEOMETRY
# --------------------------------------------------------------------------

class G:
    """Derived geometry -- everything downstream reads from here."""

    def __init__(self, P):
        self.P = P
        m = copy.deepcopy(MAGNET_PRESETS[P["magnet_preset"]])
        self.mag = m
        mc = P["magnet_pocket_clear"]

        self.core_cx = -P["core_rx"]
        self.core_r = (P["core_rx"], P["core_ry"], P["core_rz"])
        self.core_c = (self.core_cx, 0.0, 0.0)
        self.inner_r = tuple(v - P["core_wall"] for v in self.core_r)
        self.z_cut = P["faceplate_z"]

        self.pocket_r = 0.5 * (P["driver_carrier_od"] + P["driver_pocket_clear"])
        self.pocket_z1 = self.z_cut
        self.pocket_z0 = self.z_cut - P["driver_pocket_depth"]
        self.front_wall_x = self.core_cx + self.pocket_r      # +X wall of the driver pocket

        # ---- nozzle / insert stack (all magnets live here, none in the core)
        self.stub_x1 = P["stub_len"]
        self.socket_x0 = -0.20
        self.socket_x1 = P["stub_len"] + 0.20                 # 3.20
        self.mag_pocket_r = 0.5 * m["od"] + mc
        self.flange_x0 = self.socket_x1
        self.flange_x1 = self.flange_x0 + m["t"] + 0.20       # 4.90
        self.fixed_mag_x0 = self.flange_x1 - m["t"]           # 3.40
        self.fixed_mag_face = self.flange_x1                  # +X face of the fixed ring

        # carrier: its magnet floor sits rest_gap in front of the fixed ring face
        self.carrier_mag_face = self.fixed_mag_face + m["gap"]
        self.cbore_depth = m["gap"] + 0.25                    # counterbore over the flange
        self.carrier_x0 = self.carrier_mag_face - self.cbore_depth
        self.carrier_x1 = self.carrier_x0 + P["carrier_len"]
        self.cbore_r = 0.5 * P["socket_od"] + 0.15
        self.collar_r = self.cbore_r + 0.65
        # the carrier ring must have ID >= carrier_bore -- it slides on the tube
        self.carrier_mag_ri = 0.5 * P["carrier_bore"] + 0.10
        self.carrier_mag_ro = 0.5 * m["od"] + mc
        self.carrier_mag_id_eff = 2.0 * self.carrier_mag_ri

        self.insert_lug_x0 = P["insert_lug_x0"]
        self.insert_lug_x1 = self.insert_lug_x0 + P["lug_w"]

        # skirt: rim at the carrier tip so nothing protrudes past the seal plane
        self.skirt_rim_r = 0.5 * P["skirt_max_dia"] - 0.5 * P["skirt_wall"]
        self.skirt_root_r = 0.5 * P["carrier_od"]
        dr = self.skirt_rim_r - self.skirt_root_r
        self.skirt_dx = dr / math.tan(math.radians(P["skirt_flare_deg"]))
        self.skirt_rim_x = self.carrier_x1
        self.skirt_root_x = self.skirt_rim_x - self.skirt_dx

        self.tip_protrusion = self.carrier_x1                 # from the core face
        self.tip_protrusion_max = self.carrier_x1 + P["carrier_travel"]

        # ---- faceplate / jacket magnet stations
        outer_half_x = P["core_rx"] * math.sqrt(max(1e-9, 1 - (self.z_cut / P["core_rz"]) ** 2))
        self.fp_mag_off = 0.5 * (self.pocket_r + 0.5 * P["jmag_dia"] + 0.35
                                 + outer_half_x - 0.5 * P["jmag_dia"] - 0.35)
        self.fp_mag_off = max(self.fp_mag_off, self.pocket_r + 0.5 * P["jmag_dia"] + 0.3)

        # jacket magnets (3) + locating pins (2) on the -Z hemisphere, as (x, y)
        cx = self.core_cx
        self.jacket_mags = [(cx + 5.6, 0.0), (cx - 5.6, 0.0), (cx, 4.4)]
        self.jacket_pins = [(cx + 2.6, -4.0), (cx - 2.6, -4.0)]

        # ---- wing centreline (quadratic Bezier in XY)
        y_root = P["core_ry"] + P["clearance"]
        back = math.radians(P["wing_back_deg"])
        self.wing_p0 = (cx + P["wing_root_dx"], y_root - 1.6)
        self.wing_p1 = (cx + P["wing_root_dx"], y_root + 0.62 * P["wing_rise"])
        self.wing_p2 = (cx + P["wing_root_dx"] - P["wing_rise"] * math.tan(back),
                        y_root + P["wing_rise"])
        self.y_root = y_root

        # a Ti part gets no magnet in the core for the float; jacket magnets only
        self.warnings = []
        if P["wall_face"] < P["min_wall"] - 1e-9:
            self.warnings.append(
                f"wall_face {P['wall_face']} mm < min_wall {P['min_wall']} mm")
        if P["gyroid_cell"] < P["min_cell"] - 1e-9:
            self.warnings.append(
                f"gyroid_cell {P['gyroid_cell']} mm < min_cell {P['min_cell']} mm")
        if m["id"] < P["nozzle_bore"] - 1e-9:
            self.warnings.append(
                f"magnet ID {m['id']} mm < nozzle bore {P['nozzle_bore']} mm -- the "
                f"fixed ring necks the acoustic bore to {m['id']} mm over {m['t']} mm")
        if self.carrier_mag_id_eff > m["id"] + 1e-9:
            self.warnings.append(
                f"carrier ring ID opened {m['id']}->{self.carrier_mag_id_eff:.1f} mm to "
                f"clear the {P['carrier_bore']} mm bore; the pair is asymmetric, so the "
                f"real rest gap will be SHORTER than the modelled {m['gap']} mm")
        if self.front_wall_x > self.fixed_mag_x0:
            pass  # nothing in the core any more; kept for symmetry

    def core_outer(self, X, Y, Z, C):
        """The bare outer surface of the core (no pockets) -- the jacket offsets from this."""
        P = self.P
        e = ellipsoid(X, Y, Z, self.core_c, self.core_r)
        nose = cone_x(X, Y, Z, P["nose_cone_x0"], P["nose_cone_r0"], 0.0, P["nose_cone_r1"])
        return smin(e, nose, 1.2)


# --------------------------------------------------------------------------
# GRID EVALUATION
# --------------------------------------------------------------------------

class Ctx:
    """Per-evaluation context: 1-D axis samples plus a scratch cache."""

    def __init__(self, x, y, z):
        self.x, self.y, self.z = x, y, z
        self.cache = {}


def evaluate(fn, bounds, spacing, slab_mb=48.0):
    """Sample fn on a padded regular grid, in Z slabs, returning (field, origin)."""
    (x0, x1), (y0, y1), (z0, z1) = bounds
    nx = int(math.ceil((x1 - x0) / spacing)) + 1
    ny = int(math.ceil((y1 - y0) / spacing)) + 1
    nz = int(math.ceil((z1 - z0) / spacing)) + 1
    x = x0 + spacing * np.arange(nx)
    y = y0 + spacing * np.arange(ny)
    z = z0 + spacing * np.arange(nz)

    field = np.empty((nx, ny, nz), dtype=np.float32)
    per_slab = max(1, int(slab_mb * 1e6 / (nx * ny * 8)))
    C = Ctx(x, y, z)
    Xb = x.reshape(nx, 1, 1)
    Yb = y.reshape(1, ny, 1)
    for k0 in range(0, nz, per_slab):
        k1 = min(nz, k0 + per_slab)
        Zb = z[k0:k1].reshape(1, 1, k1 - k0)
        field[:, :, k0:k1] = fn(Xb, Yb, Zb, C).astype(np.float32)

    pad = float(max(4.0 * spacing, 1.0))
    field = np.pad(field, 1, mode="constant", constant_values=pad)
    origin = (x0 - spacing, y0 - spacing, z0 - spacing)
    return field, origin, spacing


def polygonise(field, origin, spacing):
    verts, faces, _, _ = measure.marching_cubes(field, level=0.0, spacing=(spacing,) * 3)
    verts = verts + np.asarray(origin)
    # NOTE: marching_cubes already emits a manifold, index-shared surface.  Do NOT
    # run merge_vertices()/nondegenerate_faces() on it -- welding across a 0.2 mm
    # lattice wall creates non-manifold edges and loses watertightness.
    mesh = trimesh.Trimesh(vertices=verts, faces=faces, process=False)
    if mesh.volume < 0:
        mesh.invert()
    return mesh


def spacing_for(bounds, budget, override=None):
    if override:
        return float(override)
    vol = 1.0
    for lo, hi in bounds:
        vol *= max(hi - lo, 1e-6)
    return (vol / budget) ** (1.0 / 3.0)


# --------------------------------------------------------------------------
# PART: core
# --------------------------------------------------------------------------

def part_core(g):
    P, cx = g.P, g.core_cx
    w = P["core_wall"]
    bore_r = 0.5 * P["nozzle_bore"]

    def fn(X, Y, Z, C):
        outer = g.core_outer(X, Y, Z, C)
        d = np.maximum(outer, Z - g.z_cut)                       # cut at the faceplate plane

        # nozzle stub with two external bayonet lugs
        stub = cyl_x(X, Y, Z, 0, 0, 0.5 * P["stub_od"], -1.0, g.stub_x1)
        lug_x0 = 1.0
        lug_x1 = lug_x0 + P["lug_w"]
        for th in (0.0, np.pi):
            stub = U(stub, arc_slot_x(X, Y, Z, 0.0,
                                      0.5 * P["stub_od"] + P["lug_h"],
                                      lug_x0, lug_x1, th - 0.42, th + 0.42))
        d = U(d, np.maximum(stub, Z - g.z_cut))

        # 2-pin socket sleeve (unioned so the pocket never opens into the cavity)
        sx0 = cx - P["core_rx"] - 0.6
        sx1 = sx0 + P["socket_d"] + 0.6
        sl_h = (0.5 * (sx1 - sx0), 0.5 * P["socket_w"] + 0.7, 0.5 * P["socket_h"] + 0.7)
        sl_c = (0.5 * (sx0 + sx1), 0.0, P["socket_z"])
        sleeve = I(rbox(X, Y, Z, sl_c, sl_h, 0.4), outer, Z - g.z_cut)
        d = U(d, sleeve)

        # bone-sensor boss on the -Y (tragus) flank
        bone_y = -(P["core_ry"] + 0.55)
        boss = rbox(X, Y, Z, (cx, bone_y, -1.0),
                    (0.5 * P["bone_w"] + 1.1, 1.5, 0.5 * P["bone_h"] + 1.1), 0.7)
        d = smin(d, boss, 0.8)

        # ---- interior --------------------------------------------------
        cavity = np.maximum(
            ellipsoid(X, Y, Z, g.core_c, g.inner_r), Z - P["cavity_cap_z"])
        pocket = cyl_z(X, Y, Z, cx, 0.0, g.pocket_r, g.pocket_z0, g.pocket_z1 + 1.0)
        bore = cyl_x(X, Y, Z, 0, 0, bore_r, cx - 1.0, g.stub_x1 + 1.0)
        void = U(cavity, pocket, bore)

        # vents
        void = U(void, capsule(X, Y, Z,
                               (g.front_wall_x + 0.6, 0.0, 0.0),
                               (g.front_wall_x + 0.6, -P["core_ry"] - 3.0,
                                -P["core_rz"] - 3.0),
                               0.5 * P["vent_dia"]))
        void = U(void, capsule(X, Y, Z,
                               (cx - 3.0, 0.0, -1.0),
                               (cx - 5.5, -P["core_ry"] - 3.0, -P["core_rz"] - 3.0),
                               0.5 * P["vent_dia"]))

        # 2-pin socket pocket
        void = U(void, box(X, Y, Z,
                           (sx0 + 0.5 * (P["socket_d"] + 1.2), 0.0, P["socket_z"]),
                           (0.5 * (P["socket_d"] + 1.2), 0.5 * P["socket_w"],
                            0.5 * P["socket_h"])))

        # bone-sensor pocket + wire channel to the socket pocket
        pk_y = bone_y - 1.5 + 0.5 * P["bone_d"]
        void = U(void, box(X, Y, Z, (cx, pk_y, -1.0),
                           (0.5 * P["bone_w"], 0.5 * P["bone_d"], 0.5 * P["bone_h"])))
        void = U(void, capsule(X, Y, Z, (cx, pk_y, -1.0),
                               (sx0 + 1.5, -1.4, P["socket_z"]),
                               0.5 * P["wire_dia"]))

        # jacket magnet pockets + locating pin holes on the -Z hemisphere
        for (mx, my) in g.jacket_mags:
            zs = _lower_z(g, mx, my)
            void = U(void, cyl_z(X, Y, Z, mx, my, 0.5 * P["jmag_dia"] + 0.05,
                                 zs - 0.4, zs + P["jmag_depth"]))
        for (px, py) in g.jacket_pins:
            zs = _lower_z(g, px, py)
            void = U(void, cyl_z(X, Y, Z, px, py, 0.5 * P["pin_dia"] + 0.06,
                                 zs - 0.4, zs + P["pin_depth"]))

        # faceplate magnet pockets in the +Z rim
        for sgn in (+1, -1):
            void = U(void, cyl_z(X, Y, Z, cx + sgn * g.fp_mag_off, 0.0,
                                 0.5 * P["jmag_dia"] + 0.05,
                                 g.z_cut - P["jmag_depth"], g.z_cut + 1.0))

        # fixed-magnet counterbore lives in the nozzle INSERT, not here (v0.2 change)

        d = S(d, void)

        # gasket groove around the parting line at z = 0
        groove = I(np.abs(Z) - 0.5 * P["gasket_w"],
                   -(outer + P["gasket_d"]),
                   X - (cx + P["core_rx"] - 2.0))
        d = S(d, groove)
        return d

    b = ((cx - P["core_rx"] - 1.6, g.stub_x1 + 1.2),
         (-P["core_ry"] - 3.2, P["core_ry"] + 1.6),
         (-P["core_rz"] - 1.6, g.z_cut + 1.2))
    return fn, b


def _lower_z(g, x, y):
    """Z of the core's lower (-Z) surface at (x, y), on the bare ellipsoid."""
    P = g.P
    t = 1.0 - ((x - g.core_cx) / P["core_rx"]) ** 2 - (y / P["core_ry"]) ** 2
    return -P["core_rz"] * math.sqrt(max(t, 1e-4))


# --------------------------------------------------------------------------
# PART: faceplate
# --------------------------------------------------------------------------

def part_faceplate(g):
    P, cx = g.P, g.core_cx

    def fn(X, Y, Z, C):
        d = I(g.core_outer(X, Y, Z, C), g.z_cut - Z)
        inner = tuple(v - 1.0 for v in g.core_r)
        d = S(d, I(ellipsoid(X, Y, Z, g.core_c, inner), g.z_cut + 0.9 - Z))
        for sgn in (+1, -1):
            d = S(d, cyl_z(X, Y, Z, cx + sgn * g.fp_mag_off, 0.0,
                           0.5 * P["jmag_dia"] + 0.05,
                           g.z_cut - 1.0, g.z_cut + P["jmag_depth"]))
        return d

    b = ((cx - P["core_rx"] - 1.0, 1.0),
         (-P["core_ry"] - 1.0, P["core_ry"] + 1.0),
         (g.z_cut - 1.0, P["core_rz"] + 1.0))
    return fn, b


# --------------------------------------------------------------------------
# PART: jacket_wing
# --------------------------------------------------------------------------

def _bezier_pts(p0, p1, p2, n=48):
    t = np.linspace(0.0, 1.0, n).reshape(-1, 1)
    p0 = np.asarray(p0); p1 = np.asarray(p1); p2 = np.asarray(p2)
    return (1 - t) ** 2 * p0 + 2 * (1 - t) * t * p1 + t ** 2 * p2


def _bezier_dist2d(C, pts):
    """Min distance from every (x, y) grid node to the polyline, cached."""
    key = ("bez", id(pts))
    if key in C.cache:
        return C.cache[key]
    x = C.x.reshape(-1, 1, 1)
    y = C.y.reshape(1, -1, 1)
    best = None
    for i in range(len(pts) - 1):
        ax, ay = pts[i]
        bx, by = pts[i + 1]
        dx, dy = bx - ax, by - ay
        dd = dx * dx + dy * dy
        px, py = x - ax, y - ay
        h = np.clip((px * dx + py * dy) / dd, 0.0, 1.0)
        d = np.sqrt((px - dx * h) ** 2 + (py - dy * h) ** 2)
        best = d if best is None else np.minimum(best, d)
    C.cache[key] = best
    return best


def part_jacket_wing(g):
    P, cx = g.P, g.core_cx
    pts = _bezier_pts(g.wing_p0, g.wing_p1, g.wing_p2)
    clear = P["clearance"]
    thick = P["jacket_thick"]

    def fn(X, Y, Z, C):
        env = g.core_outer(X, Y, Z, C)

        # --- jacket envelope: offset shell over the -Z hemisphere
        shell = np.maximum(clear - env, env - clear - thick)
        shell = I(shell, Z, X - P["jacket_x_clip"])

        # --- wing envelope: swept blade with a 45 deg taper on its deep edge
        d2 = _bezier_dist2d(C, pts) - 0.5 * P["wing_thick"]
        rate = math.tan(math.radians(P["wing_taper_deg"]))
        shrink = rate * np.maximum(0.0, (-Z) - (P["wing_width"] - P["wing_taper"]))
        wing = np.maximum(d2 + shrink,
                          np.maximum(Z - P["wing_z_top"], -Z - P["wing_width"]))
        wing = S(wing, env - clear)              # keep clear of the core

        envelope = U(shell, wing)

        # --- graded gyroid
        dz = np.maximum(0.0, -Z)
        dy = np.maximum(0.0, Y - g.y_root)
        root_dist = np.sqrt(dz * dz + dy * dy)
        t = np.clip(root_dist / P["grade_len"], 0.0, 1.0)
        wall = P["wall_root"] + (P["wall_face"] - P["wall_root"]) * t
        lat = gyroid(X, Y, Z, P["gyroid_cell"], wall)

        # --- solid root collar
        solid = root_dist - P["solid_root"]

        # --- perforated skin membrane on the ear-facing surface
        skin = np.maximum(shell, (clear + thick - P["skin_t"]) - env)
        mx = np.abs((X + 0.5 * P["perf_pitch"]) % P["perf_pitch"] - 0.5 * P["perf_pitch"])
        my = np.abs((Y + 0.5 * P["perf_pitch"]) % P["perf_pitch"] - 0.5 * P["perf_pitch"])
        perf = np.sqrt(mx * mx + my * my) - 0.5 * P["perf_dia"]
        skin = S(skin, perf)

        d = I(U(lat, solid), envelope)
        d = U(d, skin)

        # --- matching magnet pockets + locating pins
        for (px_, py_) in g.jacket_mags:
            zs = _lower_z(g, px_, py_)
            d = S(d, cyl_z(X, Y, Z, px_, py_, 0.5 * P["jmag_dia"] + 0.05,
                           zs - clear - P["jmag_depth"] - 0.4, zs - clear + 0.4))
        for (px_, py_) in g.jacket_pins:
            zs = _lower_z(g, px_, py_)
            d = U(d, cyl_z(X, Y, Z, px_, py_, 0.5 * P["pin_dia"] - 0.06,
                           zs - clear - P["pin_depth"] - 1.4, zs + P["pin_depth"] - 0.05))
        d = S(d, env - clear)
        return d

    ymax = max(p[1] for p in pts) + 3.0
    xmin = min(min(p[0] for p in pts), cx - P["core_rx"]) - 3.5
    b = ((xmin, P["jacket_x_clip"] + 1.0),
         (-P["core_ry"] - 3.0, ymax),
         (-P["core_rz"] - thick - 2.0, 1.0))
    return fn, b


# --------------------------------------------------------------------------
# PART: nozzle inserts
# --------------------------------------------------------------------------

def _carrier_lslot(g, X, Y, Z, th, radial_lo, radial_hi):
    """L-slot solid (to be subtracted from the carrier bore wall)."""
    P = g.P
    w = 0.5 * (P["lug_w"] + 0.30)
    entry = arc_slot_x(X, Y, Z, radial_lo, radial_hi,
                       g.insert_lug_x0 - 0.35, g.carrier_x1 + 1.0,
                       th - w / 2.6, th + w / 2.6)
    arc = arc_slot_x(X, Y, Z, radial_lo, radial_hi,
                     g.insert_lug_x0 - 0.35, g.insert_lug_x1 + 0.35,
                     th, th + math.radians(85.0))
    th2 = th + math.radians(85.0)
    pocket = arc_slot_x(X, Y, Z, radial_lo, radial_hi,
                        g.insert_lug_x0 - 0.35,
                        g.insert_lug_x1 + 0.35 + P["carrier_travel"],
                        th2 - w / 2.6, th2 + w / 2.6)
    return U(entry, arc, pocket)


def part_nozzle_insert(g, name):
    P = g.P
    Ltube = P["insert_tube_lengths"][name]
    x_end = g.flange_x1 + Ltube
    ro = 0.5 * P["insert_od"]
    bore_r = 0.5 * P["nozzle_bore"]

    def fn(X, Y, Z, C):
        socket = cyl_x(X, Y, Z, 0, 0, 0.5 * P["socket_od"], g.socket_x0, g.socket_x1)
        flange = cyl_x(X, Y, Z, 0, 0, 0.5 * P["socket_od"], g.flange_x0, g.flange_x1)
        tube = cyl_x(X, Y, Z, 0, 0, ro, g.flange_x1 - 0.4, x_end)
        d = U(socket, flange, tube)

        # carrier bayonet lugs on the tube
        for th in (0.0, np.pi):
            d = U(d, arc_slot_x(X, Y, Z, 0.0, ro + P["lug_h"],
                                g.insert_lug_x0, g.insert_lug_x1,
                                th - 0.30, th + 0.30))

        # bore
        d = S(d, cyl_x(X, Y, Z, 0, 0, bore_r, g.flange_x0 - 1.0, x_end + 1.0))
        # socket bore over the core stub + L-slots for the core lugs
        sb = 0.5 * (P["stub_od"] + P["press_clearance"])
        d = S(d, cyl_x(X, Y, Z, 0, 0, sb, g.socket_x0 - 1.0, g.socket_x1))
        for th in (0.0, np.pi):
            d = S(d, U(
                arc_slot_x(X, Y, Z, sb - 0.05, sb + P["lug_h"] + 0.15,
                           g.socket_x0 - 1.0, 2.30, th - 0.34, th + 0.34),
                arc_slot_x(X, Y, Z, sb - 0.05, sb + P["lug_h"] + 0.15,
                           1.55, 2.30, th, th + math.radians(80.0)),
                arc_slot_x(X, Y, Z, sb - 0.05, sb + P["lug_h"] + 0.15,
                           1.55, 2.30,
                           th + math.radians(80.0) - 0.34,
                           th + math.radians(80.0) + 0.34)))

        # fixed ring-magnet counterbore, opening +X, floor at fixed_mag_x0
        d = S(d, cyl_x(X, Y, Z, 0, 0, g.mag_pocket_r,
                       g.fixed_mag_x0, g.flange_x1 + 1.0))
        # damper-disc recess at the ear end
        d = S(d, cyl_x(X, Y, Z, 0, 0, 0.5 * P["damper_dia"] + 0.15,
                       x_end - P["damper_recess"], x_end + 1.0))
        return d

    b = ((g.socket_x0 - 1.2, x_end + 1.2), (-5.5, 5.5), (-5.5, 5.5))
    return fn, b


# --------------------------------------------------------------------------
# PART: mag-float carrier (+ moulds)
# --------------------------------------------------------------------------

def carrier_field(g, X, Y, Z, C):
    P = g.P
    body = cyl_x(X, Y, Z, 0, 0, 0.5 * P["carrier_od"], g.carrier_x0, g.carrier_x1)
    collar = cyl_x(X, Y, Z, 0, 0, g.collar_r, g.carrier_x0, g.carrier_mag_face)
    skirt = revolve_segment(X, Y, Z,
                            (g.skirt_root_x, g.skirt_root_r),
                            (g.skirt_rim_x, g.skirt_rim_r),
                            0.5 * P["skirt_wall"])
    d = U(body, collar, skirt)

    # counterbore that swallows the insert's magnet flange (this is the air gap)
    d = S(d, cyl_x(X, Y, Z, 0, 0, g.cbore_r,
                   g.carrier_x0 - 1.0, g.carrier_mag_face))
    # annular ring-magnet pocket
    d = S(d, tube_x(X, Y, Z, 0, 0, g.carrier_mag_ri, g.carrier_mag_ro,
                    g.carrier_mag_face, g.carrier_mag_face + g.mag["t"] + 0.05))
    # sliding bore
    d = S(d, cyl_x(X, Y, Z, 0, 0, 0.5 * P["carrier_bore"],
                   g.carrier_mag_face, g.carrier_x1 + 1.0))
    # L-slots
    rl = 0.5 * P["carrier_bore"] - 0.05
    rh = 0.5 * P["carrier_bore"] + P["lug_h"] + 0.15
    for th in (0.0, np.pi):
        d = S(d, _carrier_lslot(g, X, Y, Z, th, rl, rh))
    return d


def part_carrier(g):
    P = g.P
    r = g.skirt_rim_r + P["skirt_wall"] + 1.0

    def fn(X, Y, Z, C):
        return carrier_field(g, X, Y, Z, C)

    b = ((g.carrier_x0 - 1.2, g.carrier_x1 + 1.2), (-r, r), (-r, r))
    return fn, b


def mold_geom(g):
    P = g.P
    x0 = g.carrier_x0 - 3.5
    x1 = g.carrier_x1 + 3.5
    r = g.skirt_rim_r + P["skirt_wall"] + P["mold_wall"]
    return x0, x1, r


def mold_core_field(g, X, Y, Z, C):
    """The removable core rod: defines the bore, L-slots and the magnet seat."""
    P = g.P
    x0, x1, r = mold_geom(g)
    rod = cyl_x(X, Y, Z, 0, 0, 0.5 * P["carrier_bore"], x0 - 2.0, x1 + 2.0)
    seat = cyl_x(X, Y, Z, 0, 0, g.carrier_mag_ro,
                 g.carrier_mag_face, g.carrier_mag_face + g.mag["t"] + 0.05)
    plug = cyl_x(X, Y, Z, 0, 0, g.cbore_r, x0 - 2.0, g.carrier_mag_face)
    d = U(rod, seat, plug)
    for th in (0.0, np.pi):
        d = U(d, _carrier_lslot(g, X, Y, Z, th,
                                0.0, 0.5 * P["carrier_bore"] + P["lug_h"] + 0.15))
    grip = cyl_x(X, Y, Z, 0, 0, 4.0, x1 + 2.0, x1 + 8.0)
    return U(d, grip)


def part_mold_core(g):
    x0, x1, r = mold_geom(g)

    def fn(X, Y, Z, C):
        return mold_core_field(g, X, Y, Z, C)

    return fn, ((x0 - 3.0, x1 + 9.5), (-5.0, 5.0), (-5.0, 5.0))


def _mold_half(g, upper):
    P = g.P
    x0, x1, r = mold_geom(g)
    sgn = 1.0 if upper else -1.0
    pinx = (x1 - 2.5, x0 + 2.5)
    pinz = g.skirt_rim_r + 0.5 * P["mold_wall"]

    def fn(X, Y, Z, C):
        block = box(X, Y, Z, (0.5 * (x0 + x1), 0.0, 0.0),
                    (0.5 * (x1 - x0), r, r))
        d = I(block, -sgn * Y)

        void = carrier_field(g, X, Y, Z, C)
        void = U(void, mold_core_field(g, X, Y, Z, C))
        # pour spout into the ring end, from +Z
        void = U(void, cyl_z(X, Y, Z, g.carrier_x1 - 1.6, 0.0,
                             0.5 * P["mold_spout_dia"],
                             0.5 * P["carrier_od"] - 0.4, r + 2.0))
        # vent from the skirt rim, from +Z
        void = U(void, cyl_z(X, Y, Z, g.skirt_rim_x - 0.4, 0.0,
                             0.5 * P["mold_vent_dia"], g.skirt_rim_r - 0.6, r + 2.0))
        d = S(d, void)

        for px in pinx:
            for pz in (pinz, -pinz):
                pin = cyl_y(X, Y, Z, px, pz, 0.5 * P["mold_pin_dia"],
                            -P["mold_pin_len"], P["mold_pin_len"])
                if upper:
                    d = U(d, I(pin, Y))                      # pins on the +Y half
                    d = S(d, I(pin, -Y))
                else:
                    hole = cyl_y(X, Y, Z, px, pz, 0.5 * P["mold_pin_dia"] + 0.12,
                                 -P["mold_pin_len"] - 0.2, 0.2)
                    d = S(d, hole)
        return d

    b = ((x0 - 1.5, x1 + 1.5), (-r - 4.0, r + 4.0), (-r - 1.5, r + 1.5))
    return fn, b


def part_mold_a(g):
    return _mold_half(g, True)


def part_mold_b(g):
    return _mold_half(g, False)


# --------------------------------------------------------------------------
# PART: driver carrier
# --------------------------------------------------------------------------

def part_driver_carrier(g):
    P, cx = g.P, g.core_cx
    z0 = g.pocket_z0 + 0.1
    z1 = z0 + P["driver_carrier_h"]

    def fn(X, Y, Z, C):
        ring = cyl_z(X, Y, Z, cx, 0.0, 0.5 * P["driver_carrier_od"], z0, z1)
        ring = S(ring, cyl_z(X, Y, Z, cx, 0.0, 0.5 * P["driver_carrier_id"],
                             z0 - 1.0, z1 + 1.0))
        # rear vent notch, on the -X flank, cut from the +Z face
        ring = S(ring, box(X, Y, Z,
                           (cx - 0.5 * P["driver_carrier_od"] + 0.2, 0.0, z1 - 0.4),
                           (1.2, 1.0, 0.9)))
        return ring

    r = 0.5 * P["driver_carrier_od"] + 1.2
    return fn, ((cx - r, cx + r), (-r, r), (z0 - 1.2, z1 + 1.2))


# --------------------------------------------------------------------------
# PART: damper jig
# --------------------------------------------------------------------------

def part_damper_jig(g):
    P = g.P
    dd = P["damper_dia"]

    def fn(X, Y, Z, C):
        die = rbox(X, Y, Z, (0, 0, 3.0), (13.0, 8.0, 3.0), 1.2)
        die = S(die, cyl_z(X, Y, Z, -6.0, 0.0, 0.5 * dd + 0.03, -1.0, 7.0))   # punch die
        die = S(die, cyl_z(X, Y, Z, 0.0, 0.0, 0.5 * dd + 0.10, 2.6, 7.0))     # magazine
        die = S(die, cyl_z(X, Y, Z, 6.0, 0.0, 0.5 * dd + 0.35, 4.0, 7.0))     # drop nest
        die = S(die, box(X, Y, Z, (6.0, 0.0, 6.4), (3.4, 1.1, 1.2)))          # tweezer relief
        rod = cyl_z(X, Y, Z, 0.0, 14.0, 0.5 * dd - 0.10, 0.0, 18.0)           # plunger
        rod = U(rod, cyl_z(X, Y, Z, 0.0, 14.0, 4.0, 18.0, 20.0))              # knob
        return U(die, rod)

    return fn, ((-15.0, 15.0), (-10.0, 19.0), (-1.0, 21.5))


# --------------------------------------------------------------------------
# OVERHANG ANALYSIS
# --------------------------------------------------------------------------

def overhang_stats(mesh, growth):
    """growth = unit vector the printer builds along.

    Returns (worst_deg, area_fraction_over_45).  A face whose outward normal
    points along `growth` is a down-facing (unsupported) face; 0 deg is a
    vertical wall, 90 deg a flat ceiling.
    """
    g = np.asarray(growth, dtype=float)
    g /= np.linalg.norm(g)
    n = np.nan_to_num(np.asarray(mesh.face_normals, dtype=float), nan=0.0,
                      posinf=0.0, neginf=0.0)
    a = np.asarray(mesh.area_faces, dtype=float)
    c = np.clip(n @ g, -1.0, 1.0)
    ang = 90.0 - np.degrees(np.arccos(c))
    m = c > 1e-6
    if not m.any():
        return 0.0, 0.0
    worst = float(ang[m].max())
    frac = float(a[m & (ang > 45.0)].sum() / a.sum())
    # area-weighted 99th percentile: robust to marching-cubes staircase facets
    order = np.argsort(ang[m])
    aw = np.cumsum(a[m][order])
    p99 = float(ang[m][order][np.searchsorted(aw, 0.99 * aw[-1])])
    return worst, frac, p99


def best_build_dir(mesh, n=42):
    """Fibonacci-sample build directions; return the one with the least >45 deg area."""
    i = np.arange(n) + 0.5
    phi = np.arccos(1 - 2 * i / n)
    gold = np.pi * (1 + 5 ** 0.5)
    theta = gold * i
    dirs = np.stack([np.cos(theta) * np.sin(phi),
                     np.sin(theta) * np.sin(phi),
                     np.cos(phi)], axis=1)
    best = None
    for d in dirs:
        w, f, p = overhang_stats(mesh, d)
        if best is None or f < best[2]:
            best = (d, w, f, p)
    return best


def wing_envelope_mesh(g, spacing=0.22):
    """The wing's macro envelope (no lattice) -- the meaningful overhang metric."""
    P = g.P
    pts = _bezier_pts(g.wing_p0, g.wing_p1, g.wing_p2)

    def fn(X, Y, Z, C):
        d2 = _bezier_dist2d(C, pts) - 0.5 * P["wing_thick"]
        rate = math.tan(math.radians(P["wing_taper_deg"]))
        shrink = rate * np.maximum(0.0, (-Z) - (P["wing_width"] - P["wing_taper"]))
        wing = np.maximum(d2 + shrink,
                          np.maximum(Z - P["wing_z_top"], -Z - P["wing_width"]))
        # only the free-standing span beyond the jacket rim: inboard of that the
        # wing is fused to (and supported by) the jacket shell.
        return I(wing, g.y_root - Y)

    ymax = max(p[1] for p in pts) + 3.0
    xmin = min(p[0] for p in pts) - 4.0
    b = ((xmin, g.core_cx + P["core_rx"] + 2.0),
         (g.y_root - 1.0, ymax),
         (-P["wing_width"] - 1.5, 1.5))
    f, o, s = evaluate(fn, b, spacing)
    return polygonise(f, o, s)


# --------------------------------------------------------------------------
# DRIVER
# --------------------------------------------------------------------------

PARTS = {
    "core": (part_core, "solid"),
    "faceplate": (part_faceplate, "solid"),
    "jacket_wing": (part_jacket_wing, "lattice"),
    "nozzle_insert_short": (lambda g: part_nozzle_insert(g, "short"), "solid"),
    "nozzle_insert_med": (lambda g: part_nozzle_insert(g, "med"), "solid"),
    "nozzle_insert_long": (lambda g: part_nozzle_insert(g, "long"), "solid"),
    "carrier": (part_carrier, "solid"),
    "carrier_mold_a": (part_mold_a, "solid"),
    "carrier_mold_b": (part_mold_b, "solid"),
    "carrier_mold_core": (part_mold_core, "solid"),
    "driver_carrier": (part_driver_carrier, "solid"),
    "damper_jig": (part_damper_jig, "solid"),
}

# parts that live in the assembly frame (moulds and the jig have their own frames)
ASSEMBLY_PARTS = ["core", "faceplate", "jacket_wing",
                  "nozzle_insert_short", "carrier", "driver_carrier"]


def build(name, g, voxel=None):
    fn, bounds = PARTS[name][0](g)
    budget = g.P["budget_lattice"] if PARTS[name][1] == "lattice" else g.P["budget"]
    sp = spacing_for(bounds, budget, voxel)
    t0 = time.time()
    field, origin, sp = evaluate(fn, bounds, sp)
    mesh = polygonise(field, origin, sp)
    return mesh, sp, time.time() - t0


def mirror(mesh):
    m = mesh.copy()
    m.apply_transform(np.diag([-1.0, 1.0, 1.0, 1.0]))
    m.invert()
    return m


def report_row(name, mesh, sp, dt):
    b = mesh.bounds
    ext = b[1] - b[0]
    return (f"{name:24s} {'YES' if mesh.is_watertight else 'NO ':>4s}  "
            f"vol {mesh.volume:9.1f} mm3  bbox "
            f"{ext[0]:6.2f} x {ext[1]:6.2f} x {ext[2]:6.2f}  "
            f"tris {len(mesh.faces):7d}  voxel {sp:.3f}  {dt:5.1f}s")


def main():
    ap = argparse.ArgumentParser(description="Magneto IEM parametric STL generator")
    ap.add_argument("--all", action="store_true", help="generate every part + assembly")
    ap.add_argument("--part", action="append", default=[], help="generate one part")
    ap.add_argument("--list", action="store_true")
    ap.add_argument("--out", default=None, help="output dir (default ./stl)")
    ap.add_argument("--ear", choices=["right", "left", "both"], default="both")
    ap.add_argument("--voxel", type=float, default=None, help="override voxel size (mm)")
    ap.add_argument("--magnet-preset", default=PARAMS["magnet_preset"],
                    choices=sorted(MAGNET_PRESETS))
    ap.add_argument("--no-assembly", action="store_true")
    args = ap.parse_args()

    if args.list:
        for k in PARTS:
            print(k)
        return 0

    P = copy.deepcopy(PARAMS)
    P["magnet_preset"] = args.magnet_preset
    if args.voxel:
        P["voxel"] = args.voxel
    g = G(P)

    here = os.path.dirname(os.path.abspath(__file__))
    out = args.out or os.path.join(here, "stl")
    os.makedirs(os.path.join(out, "right"), exist_ok=True)
    if args.ear in ("left", "both"):
        os.makedirs(os.path.join(out, "left"), exist_ok=True)

    names = list(PARTS) if args.all else args.part
    if not names:
        ap.error("nothing to do -- pass --all or --part NAME")

    # resolution advisory for the lattice part
    _, jb = part_jacket_wing(g)
    jsp = spacing_for(jb, P["budget_lattice"], P["voxel"])
    if jsp > P["wall_face"] / 2.5:
        print(f"  NOTE: jacket voxel {jsp:.3f} mm vs {P['wall_face']} mm lattice wall "
              f"({P['wall_face']/jsp:.1f} samples/wall).  Good enough for form and fit; "
              f"re-run with --voxel {P['wall_face']/3.0:.2f} before sending to print.")

    m = g.mag
    print("=" * 100)
    print(f"Magneto IEM generator -- magnet preset '{P['magnet_preset']}': "
          f"{m['od']}x{m['id']}x{m['t']} mm, {m['material']}, rest gap {m['gap']} mm, "
          f"F(rest)={m['f_rest']} N, ratio {m['ratio']}")
    print(f"core  {2*P['core_rx']:.1f} x {2*P['core_ry']:.1f} x {2*P['core_rz']:.1f} mm    "
          f"carrier x {g.carrier_x0:.2f}..{g.carrier_x1:.2f} mm    "
          f"tip protrusion {g.tip_protrusion:.2f} mm "
          f"(max {g.tip_protrusion_max:.2f} at full float)")
    for w in g.warnings:
        print(f"  WARNING: {w}")
    print("=" * 100)

    meshes = {}
    rows = []
    for n in names:
        if n not in PARTS:
            print(f"  ?? unknown part {n}")
            continue
        mesh, sp, dt = build(n, g, P["voxel"])
        meshes[n] = mesh
        rows.append(report_row(n, mesh, sp, dt))
        print(rows[-1])
        if args.ear in ("right", "both"):
            mesh.export(os.path.join(out, "right", f"{n}.stl"))
        if args.ear in ("left", "both"):
            mirror(mesh).export(os.path.join(out, "left", f"{n}.stl"))

    # ---- assembly ------------------------------------------------------
    if args.all and not args.no_assembly:
        parts = [meshes[k] for k in ASSEMBLY_PARTS if k in meshes]
        if parts:
            asm = trimesh.util.concatenate(parts)
            rows.append(report_row("assembly", asm, 0.0, 0.0))
            print(rows[-1])
            asm.export(os.path.join(out, "right", "assembly.stl"))
            if args.ear in ("left", "both"):
                mirror(asm).export(os.path.join(out, "left", "assembly.stl"))

    # ---- overhang analysis --------------------------------------------
    if args.all or "jacket_wing" in meshes:
        we = wing_envelope_mesh(g)
        w_rim, f_rim, p_rim = overhang_stats(we, (0, 0, -1))
        d, w_best, f_best, p_best = best_build_dir(we)
        print("-" * 100)
        print(f"wing macro envelope, printed rim-down (build dir 0,0,-1): "
              f"worst overhang {w_rim:.1f} deg, p99 {p_rim:.1f} deg, "
              f"{100*f_rim:.1f}% of area over 45 deg")
        print(f"best sampled build direction ({d[0]:+.2f},{d[1]:+.2f},{d[2]:+.2f}): "
              f"worst {w_best:.1f} deg, p99 {p_best:.1f} deg, "
              f"{100*f_best:.1f}% of area over 45 deg")
        if "jacket_wing" in meshes:
            w2, f2, p2 = overhang_stats(meshes["jacket_wing"], (0, 0, -1))
            print(f"full jacket+wing incl. lattice, rim-down: worst {w2:.1f} deg, "
                  f"p99 {p2:.1f} deg, {100*f2:.1f}% of area over 45 deg "
                  f"(gyroid micro-facets dominate this number; the gyroid is "
                  f"self-supporting in practice)")

    bad = [n for n, mm in meshes.items() if not mm.is_watertight]
    print("-" * 100)
    print(f"{len(meshes)} parts written to {out}/  "
          + ("ALL WATERTIGHT" if not bad else f"NOT WATERTIGHT: {bad}"))
    return 1 if bad else 0


if __name__ == "__main__":
    sys.exit(main())
