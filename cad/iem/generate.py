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
    # Asymmetric pairs -- the moving ring must slide over the 5 mm nozzle tube, so
    # its ID is pinned at 5.4 mm and the pair can never be identical.
    # See docs/MAGFLOAT_MAGNETS.md, "Asymmetric pair (as-built)".
    "asym_as_built": dict(
        fixed=(7.0, 3.0, 1.5), moving=(9.0, 5.4, 2.0), gap=2.25,
        material="bonded NdFeB (Br 0.65 T)",
        f_lo=0.307, f_rest=0.214, f_hi=0.155, ratio=1.98,
        mass_fixed_g=0.283, mass_moving_g=0.489, stray_mT_at_10mm=3.30),
    "asym_light": dict(
        fixed=(7.0, 3.0, 1.5), moving=(9.0, 5.4, 1.5), gap=1.25,
        material="bonded NdFeB (Br 0.65 T)",
        f_lo=0.405, f_rest=0.288, f_hi=0.199, ratio=2.04,
        mass_fixed_g=0.283, mass_moving_g=0.366, stray_mT_at_10mm=3.30),
    "asym_8mm": dict(
        fixed=(7.0, 3.0, 1.5), moving=(8.0, 5.4, 2.0), gap=1.89,
        material="bonded NdFeB (Br 0.65 T)",
        f_lo=0.382, f_rest=0.230, f_hi=0.151, ratio=2.53,
        mass_fixed_g=0.283, mass_moving_g=0.328, stray_mT_at_10mm=3.30),
    # Symmetric presets from the first study -- kept for reference.  They assume a
    # moving ring with ID 3 mm, which physically cannot slide on a 5 mm tube; the
    # generator opens it and warns.
    "n52_long": dict(
        fixed=(7.0, 3.0, 1.5), moving=(7.0, 3.0, 1.5), gap=6.40,
        material="N52 sintered NdFeB (Br 1.45 T)",
        f_lo=0.263, f_rest=0.200, f_hi=0.155, ratio=1.69,
        mass_fixed_g=0.353, mass_moving_g=0.353, stray_mT_at_10mm=None),
    "n52_clean_bore": dict(
        fixed=(8.0, 4.0, 1.5), moving=(8.0, 4.0, 1.5), gap=6.80,
        material="N52 sintered NdFeB (Br 1.45 T)",
        f_lo=0.259, f_rest=0.203, f_hi=0.161, ratio=1.61,
        mass_fixed_g=0.424, mass_moving_g=0.424, stray_mT_at_10mm=8.26),
}


PARAMS = dict(
    # ---- process / manufacturing limits -------------------------------
    min_wall=0.20,             # mm, thinnest printable lattice wall (LPBF Ti / resin)
    min_cell=1.00,             # mm, smallest usable gyroid unit cell
    clearance=0.15,            # mm, jacket-to-core sliding clearance
    press_clearance=0.15,      # mm, press/slip fit clearance on cylindrical joints

    # ---- magnets ------------------------------------------------------
    magnet_preset="asym_as_built",
    magnet_pocket_clear=0.05,  # mm added to magnet OD/thickness for the pocket
    magnet_encap=0.50,         # mm of silicone over the moving ring's axial faces

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
    nozzle_cant_deg=45.0,      # deg: nozzle axis rotated about +Y relative to the
                               #      core/faceplate/jacket, so the body lies in the
                               #      concha plane instead of down the canal axis
    bore_arc_r=9.00,           # mm centreline radius of the curved sound bore
    bore_arc_end=1.20,         # mm behind the nozzle base where the arc ends
    bore_run_in=2.50,          # mm of straight bore into the front volume
    nozzle_bore=4.00,          # mm acoustic bore
    stub_od=5.00,              # mm core nozzle stub OD
    stub_len=3.00,             # mm core nozzle stub length (+X from origin)
    lug_h=0.60,                # mm bayonet lug radial height
    lug_w=1.50,                # mm bayonet lug axial width
    socket_od=8.00,            # mm nozzle-insert socket (slips over the stub) OD
    insert_od=5.00,            # mm nozzle-insert tube OD
    insert_tube_lengths=dict(short=7.0, med=9.0, long=11.0),  # mm beyond the magnet flange
    damper_dia=4.00,           # mm damper disc
    damper_recess=0.30,        # mm damper disc recess depth

    # ---- mag-float carrier --------------------------------------------
    carrier_bore=5.20,         # mm sliding bore on the insert tube
    carrier_od=10.50,          # mm carrier body OD (moving ring OD + 2 x 0.75 wall)
    carrier_len=9.00,          # mm carrier body length
    carrier_travel=1.50,       # mm allowed axial float
    skirt_flare_deg=35.0,      # deg half-angle of the sealing skirt
    skirt_wall=0.35,           # mm skirt wall
    skirt_max_dia=19.0,        # mm skirt rim diameter (outer)
    skirt_land_w=4.50,         # mm slant width of the conical contact land (FEA: >=4.0)
    skirt_wall_neck=0.40,      # mm wall at the skirt root (structural)
    skirt_wall_hinge=0.20,     # mm wall in the compliance groove behind the land
    skirt_wall_land=0.25,      # mm wall through the contact land
    skirt_hinge_w=0.60,        # mm slant width of the compliance groove

    # ---- jacket skin (fine gyroid, structural-rigid by design) ----------
    jacket_thick=1.60,         # mm total jacket thickness (lattice + skin)
    jacket_x_clip=-2.00,       # mm, jacket stops here so it never fouls the nozzle
    gyroid_cell=1.20,          # mm gyroid unit cell, jacket skin
    wall_face=0.20,            # mm gyroid wall at the outer / ear face
    wall_root=0.40,            # mm gyroid wall at the root
    grade_len=6.00,            # mm over which the wall grades face->root
    skin_t=0.60,               # mm solid skin membrane on the ear-facing surface
    perf_dia=0.40,             # mm sweat perforation diameter
    perf_pitch=1.50,           # mm perforation grid pitch
    solid_root=1.00,           # mm of solid Ti before the lattice starts

    # ---- wing: MACRO-scale gyroid shell (a compliant doubly-curved sheet) --
    gyroid_cell_wing=12.00,    # mm wing unit cell -- 1-2 cells across the envelope
    wing_wall_root=0.22,       # mm sheet wall at the root
    wing_wall_tip=0.20,        # mm sheet wall at the tip
    wing_edge_wall=0.40,       # mm rolled/thickened rim on exposed sheet edges
    wing_edge_band=0.70,       # mm over which the wall ramps up to the rolled edge
    wing_root_solid=1.20,      # mm of solid Ti transition into the jacket rim
    wing_len=14.0,             # mm nominal wing length along its centreline
    wing_thick=7.00,           # mm envelope across the press direction (in XY)
    wing_width=5.00,           # mm envelope depth in Z, into the concha
    wing_anchor_w=2.40,        # mm Z-width at the anchor (necked foot; softens the wing)
    wing_anchor_len=7.00,      # mm over which the Z-width opens anchor_w -> wing_width
    wing_z_top=-0.20,          # mm top of the wing, just under the parting plane
    wing_taper_deg=40.0,       # deg overhang of the wing's deep-edge taper
    wing_edge_round=0.85,      # fraction of the half-section used as a corner radius
    wing_taper=1.60,           # mm of tapered depth on the wing's deep edge
    wing_rise=11.0,            # mm the tip lands above the core rim
    wing_back_deg=30.0,        # deg the tip is angled toward -X
    wing_root_dx=1.00,         # mm, wing root offset from the core centre in X
    shell_chi=0.40,            # sheet-orientation factor <cos^2 th> for the k estimate

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
    dd = max(dx * dx + dy * dy + dz * dz, 1e-12)
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
        self.fix_od, self.fix_id, self.fix_t = m["fixed"]
        self.mov_od, self.mov_id, self.mov_t = m["moving"]
        mc = P["magnet_pocket_clear"]

        self.core_cx = -P["core_rx"]
        self.core_r = (P["core_rx"], P["core_ry"], P["core_rz"])
        self.core_c = (self.core_cx, 0.0, 0.0)
        self.inner_r = tuple(v - P["core_wall"] for v in self.core_r)
        self.z_cut = P["faceplate_z"]

        self.pocket_r = 0.5 * (P["driver_carrier_od"] + P["driver_pocket_clear"])
        self.pocket_z1 = self.z_cut
        self.pocket_z0 = self.z_cut - P["driver_pocket_depth"]
        self.front_wall_x = self.core_cx + self.pocket_r

        # ---- canted nozzle frame ------------------------------------------
        # The nozzle axis is rotated about +Y by nozzle_cant_deg, so +X -> -Z.
        # The core, faceplate, jacket and wing do NOT move: only the nozzle,
        # insert, carrier and skirt follow this axis.  At cant = 0 the frame is
        # the identity and the geometry is exactly the old collinear design.
        self.cant = math.radians(P["nozzle_cant_deg"])
        ca, sa = math.cos(self.cant), math.sin(self.cant)
        self.n_ax = np.array([ca, 0.0, -sa])       # nozzle axis, into the ear
        self.n_ay = np.array([0.0, 1.0, 0.0])
        self.n_az = np.array([sa, 0.0, ca])
        # the axis passes through the core centre, so the nozzle root is buried
        # the full ellipsoid support distance -- that is what carries the side load
        self.nozzle_t_exit = 1.0 / math.sqrt((ca / P["core_rx"]) ** 2
                                             + (sa / P["core_rz"]) ** 2)
        self.nozzle_base = np.array(self.core_c) + self.nozzle_t_exit * self.n_ax
        T = np.eye(4)
        T[:3, 0], T[:3, 1], T[:3, 2] = self.n_ax, self.n_ay, self.n_az
        T[:3, 3] = self.nozzle_base
        self.nozzle_T = T                          # nozzle-local -> world

        # ---- nozzle / insert stack: BOTH magnets live here, none in the core
        self.stub_x1 = P["stub_len"]
        self.socket_x0 = -0.20
        self.socket_x1 = P["stub_len"] + 0.20                 # 3.20
        self.mag_pocket_r = 0.5 * self.fix_od + mc
        self.flange_x0 = self.socket_x1
        self.flange_x1 = self.flange_x0 + self.fix_t + 0.20   # 4.90
        self.fixed_mag_x0 = self.flange_x1 - self.fix_t
        self.fixed_mag_face = self.flange_x1                  # +X face of the fixed ring

        # carrier: the moving ring's -X face sits rest_gap in front of the fixed one
        self.carrier_x0 = self.flange_x1 - 0.25               # collar overlaps the flange
        self.carrier_mag_face = self.fixed_mag_face + m["gap"]
        self.cbore_floor = self.carrier_mag_face - P["magnet_encap"]
        self.cbore_depth = self.cbore_floor - self.carrier_x0
        self.carrier_mag_x1 = self.carrier_mag_face + self.mov_t + mc
        self.carrier_x1 = self.carrier_x0 + P["carrier_len"]
        self.cbore_r = 0.5 * P["socket_od"] + 0.15
        self.collar_r = max(self.cbore_r + 0.65, 0.5 * P["carrier_od"])
        # the moving ring's ID face is flush with the sliding bore: ID 5.4 over a
        # 5.2 bore leaves 0.1 mm, which will not cast.  The mould core's post
        # locates the ring instead.
        self.carrier_mag_ri = 0.5 * P["carrier_bore"]
        self.carrier_mag_ro = 0.5 * self.mov_od + mc
        self.carrier_wall_od = 0.5 * P["carrier_od"] - self.carrier_mag_ro

        # bayonet station: just clear of the encapsulated ring
        self.insert_lug_x0 = self.carrier_mag_x1 + 0.35
        self.insert_lug_x1 = self.insert_lug_x0 + P["lug_w"]
        self.lslot_x1 = self.insert_lug_x1 + 0.35 + P["carrier_travel"]

        # skirt.  The OUTER face is one exact 35 deg cone from the carrier OD to the
        # Would-be rim diameter, so the contact land is a true conical band and the
        # rim diameter is exactly skirt_max_dia.  All compliance is cut from the
        # inside, by varying the wall along the slant.
        fl = math.radians(P["skirt_flare_deg"])
        self.skirt_root_r = 0.5 * P["carrier_od"]
        self.skirt_rim_r = 0.5 * P["skirt_max_dia"]
        dr = self.skirt_rim_r - self.skirt_root_r
        self.skirt_dx = dr / math.tan(fl)
        self.skirt_slant = dr / math.sin(fl)
        self.skirt_rim_x = self.carrier_x1
        self.skirt_root_x = self.skirt_rim_x - self.skirt_dx
        self.skirt_u = (self.skirt_dx / self.skirt_slant, dr / self.skirt_slant)
        # wall profile along the slant: neck | compliance groove | contact land
        land = min(P["skirt_land_w"], self.skirt_slant - P["skirt_hinge_w"] - 1.0)
        s_l0 = self.skirt_slant - land
        s_h1 = s_l0 - 0.20
        s_h0 = s_h1 - P["skirt_hinge_w"]
        s_n1 = max(0.05, s_h0 - 0.20)
        self.skirt_land_w = land
        self.skirt_wall_xp = [0.0, s_n1, s_h0, s_h1, s_l0, self.skirt_slant]
        self.skirt_wall_fp = [P["skirt_wall_neck"], P["skirt_wall_neck"],
                              P["skirt_wall_hinge"], P["skirt_wall_hinge"],
                              P["skirt_wall_land"], P["skirt_wall_land"]]
        self.skirt_land_x0 = self.skirt_root_x + s_l0 * self.skirt_u[0]
        self.skirt_land_d0 = 2.0 * (self.skirt_root_r + s_l0 * self.skirt_u[1])

        self.tip_protrusion = self.carrier_x1                 # from the core face
        self.tip_protrusion_max = self.carrier_x1 + P["carrier_travel"]

        # ---- faceplate / jacket magnet stations
        outer_half_x = P["core_rx"] * math.sqrt(max(1e-9, 1 - (self.z_cut / P["core_rz"]) ** 2))
        self.fp_mag_off = 0.5 * (self.pocket_r + 0.5 * P["jmag_dia"] + 0.35
                                 + outer_half_x - 0.5 * P["jmag_dia"] - 0.35)
        self.fp_mag_off = max(self.fp_mag_off, self.pocket_r + 0.5 * P["jmag_dia"] + 0.3)

        # jacket magnets (3) + locating pins (2) on the -Z hemisphere, as (x, y)
        cx = self.core_cx
        self.jacket_mags = [(cx - 6.0, 0.0), (cx + 0.5, 4.8), (cx + 0.5, -4.8)]
        self.jacket_pins = [(cx - 3.0, 4.0), (cx - 3.0, -4.0)]

        # ---- wing centreline (quadratic Bezier in XY)
        y_root = P["core_ry"] + P["clearance"]
        back = math.radians(P["wing_back_deg"])
        self.wing_p0 = (cx + P["wing_root_dx"], y_root - 1.6)
        self.wing_p1 = (cx + P["wing_root_dx"], y_root + 0.62 * P["wing_rise"])
        self.wing_p2 = (cx + P["wing_root_dx"] - P["wing_rise"] * math.tan(back),
                        y_root + P["wing_rise"])
        self.y_root = y_root

        # ---- curved sound bore ---------------------------------------------
        # Straight along the nozzle axis through the stub, then one circular arc
        # of radius bore_arc_r turning through the cant angle back into the body
        # plane, then straight into the front volume under the driver.  Constant
        # Ø nozzle_bore: the solid is a swept sphere along this polyline.
        ax = self.n_ax
        pts = [self.nozzle_base + ax * (P["stub_len"] + 1.5)]
        arc_end = self.nozzle_base - ax * P["bore_arc_end"]
        pts.append(arc_end)
        R, N = P["bore_arc_r"], 24
        for i in range(N + 1):
            phi = self.cant * (1.0 - i / N)
            pts.append(arc_end - R * np.array([math.sin(self.cant) - math.sin(phi),
                                               0.0,
                                               math.cos(self.cant) - math.cos(phi)]))
        pts.append(pts[-1] + np.array([-P["bore_run_in"], 0.0, 0.0]))
        # drop repeated points -- at cant = 0 the whole arc collapses to one point
        keep = [pts[0]]
        for q in pts[1:]:
            if np.linalg.norm(q - keep[-1]) > 1e-6:
                keep.append(q)
        self.bore_path = np.array(keep)

        # a Ti part gets no magnet in the core for the float; jacket magnets only
        self.warnings = []
        if P["wall_face"] < P["min_wall"] - 1e-9:
            self.warnings.append(
                f"wall_face {P['wall_face']} mm < min_wall {P['min_wall']} mm")
        if P["gyroid_cell"] < P["min_cell"] - 1e-9:
            self.warnings.append(
                f"gyroid_cell {P['gyroid_cell']} mm < min_cell {P['min_cell']} mm")
        if self.fix_id < P["nozzle_bore"] - 1e-9:
            self.warnings.append(
                f"fixed ring ID {self.fix_id} mm < nozzle bore {P['nozzle_bore']} mm -- it "
                f"necks the acoustic bore to {self.fix_id} mm over {self.fix_t} mm")
        if self.mov_id < P["carrier_bore"] - 1e-9:
            self.warnings.append(
                f"moving ring ID {self.mov_id} mm < carrier bore {P['carrier_bore']} mm -- "
                f"it cannot slide on the {P['insert_od']} mm tube; pick a preset whose "
                f"moving ring has ID >= {P['carrier_bore']} mm")
        if self.carrier_wall_od < 0.5 - 1e-9:
            self.warnings.append(
                f"only {self.carrier_wall_od:.2f} mm of silicone outside the moving ring; "
                f"raise carrier_od to >= {self.mov_od + 1.1:.1f} mm")
        if self.cbore_depth < 0.3:
            self.warnings.append(
                f"carrier counterbore is only {self.cbore_depth:.2f} mm deep -- the rest "
                f"gap {m['gap']} mm is too small for {P['magnet_encap']} mm encapsulation")
        need = self.insert_lug_x1 + 0.35 + P["carrier_travel"]
        if need > self.carrier_x1 - 0.2:
            self.warnings.append(
                f"L-slot needs the carrier to reach x={need + 0.2:.2f} mm but it ends at "
                f"{self.carrier_x1:.2f}; raise carrier_len to "
                f"{need + 0.2 - self.carrier_x0:.2f} mm")

    def nz(self, X, Y, Z):
        """World -> nozzle-local coordinates (the canted frame)."""
        ca, sa = math.cos(self.cant), math.sin(self.cant)
        bx, by, bz = self.nozzle_base
        px, py, pz = X - bx, Y - by, Z - bz
        return px * ca - pz * sa, py, px * sa + pz * ca

    def core_body(self, X, Y, Z):
        """The core ellipsoid alone -- what the jacket offsets from."""
        return ellipsoid(X, Y, Z, self.core_c, self.core_r)

    def nose_cone(self, X, Y, Z, grow=0.0):
        P = self.P
        Xn, Yn, Zn = self.nz(X, Y, Z)
        return cone_x(Xn, Yn, Zn, P["nose_cone_x0"], P["nose_cone_r0"] + grow,
                      0.0, P["nose_cone_r1"] + grow)

    def core_outer(self, X, Y, Z, C):
        """The bare outer surface of the core (no pockets) -- the jacket offsets from this."""
        P = self.P
        return smin(self.core_body(X, Y, Z), self.nose_cone(X, Y, Z), 1.2)


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


def polygonise(field, origin, spacing, tag=None):
    verts, faces, _, _ = measure.marching_cubes(field, level=0.0, spacing=(spacing,) * 3)
    verts = verts + np.asarray(origin)
    # NOTE: marching_cubes already emits a manifold, index-shared surface.  Do NOT
    # run merge_vertices()/nondegenerate_faces() on it -- welding across a 0.2 mm
    # lattice wall creates non-manifold edges and loses watertightness.
    mesh = trimesh.Trimesh(vertices=verts, faces=faces, process=False)
    if mesh.volume < 0:
        mesh.invert()
    return drop_specks(mesh, tag=tag)


DROPPED = {}


def drop_specks(mesh, min_body=0.02, min_void=1.0, tag=None):
    """Clean up two marching-cubes artefacts, without breaking watertightness.

    * positive components under min_body mm^3 -- sub-voxel islands where a lattice
      wall grazes the sampling grid; they would be loose particles in the print.
    * negative components (enclosed voids) under min_void mm^3 -- trapped-powder
      pockets at gyroid/skin junctions.  Dropping the shell fills the pocket.

    Dropping whole closed components leaves the surface closed, so the result is
    still watertight.
    """
    comps = mesh.split(only_watertight=False)
    if len(comps) <= 1:
        return mesh
    keep = [c for c in comps
            if (c.volume >= min_body) or (c.volume <= -min_void)]
    if not keep or len(keep) == len(comps):
        return mesh
    if tag:
        DROPPED[tag] = (len(comps) - len(keep),
                        sum(abs(c.volume) for c in comps if c not in keep))
    out = trimesh.util.concatenate(keep)
    return out if out.is_watertight else mesh


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

        # nozzle stub with two external bayonet lugs, on the CANTED axis
        Xn, Yn, Zn = g.nz(X, Y, Z)
        stub = cyl_x(Xn, Yn, Zn, 0, 0, 0.5 * P["stub_od"], -3.0, g.stub_x1)
        lug_x0 = 1.0
        lug_x1 = lug_x0 + P["lug_w"]
        for th in (0.0, np.pi):
            stub = U(stub, arc_slot_x(Xn, Yn, Zn, 0.0,
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
        bore = None
        for i in range(len(g.bore_path) - 1):
            c_ = capsule(X, Y, Z, g.bore_path[i], g.bore_path[i + 1], bore_r)
            bore = c_ if bore is None else U(bore, c_)
        void = U(cavity, pocket, bore)

        # vents
        fv = g.bore_path[len(g.bore_path) // 2]
        void = U(void, capsule(X, Y, Z, fv,
                               (fv[0], -P["core_ry"] - 3.0, -P["core_rz"] - 3.0),
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

    tip = g.nozzle_base + g.n_ax * (g.stub_x1 + 0.6)
    b = ((cx - P["core_rx"] - 1.6, max(0.0, float(tip[0])) + 4.2),
         (-P["core_ry"] - 3.2, P["core_ry"] + 1.6),
         (min(-P["core_rz"], float(tip[2]) - 4.2) - 1.6, g.z_cut + 1.2))
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


def _bezier_dist_and_s(C, pts):
    """(min distance, arc-length station of the closest point) for every (x, y)
    grid node, against the wing centreline polyline.  Cached per evaluation."""
    key = ("bez", id(pts))
    if key in C.cache:
        return C.cache[key]
    x = C.x.reshape(-1, 1, 1)
    y = C.y.reshape(1, -1, 1)
    seg = np.linalg.norm(np.diff(pts, axis=0), axis=1)
    cum = np.concatenate([[0.0], np.cumsum(seg)])
    best_d = best_s = None
    for i in range(len(pts) - 1):
        ax, ay = pts[i]
        bx, by = pts[i + 1]
        dx, dy = bx - ax, by - ay
        dd = dx * dx + dy * dy
        px, py = x - ax, y - ay
        h = np.clip((px * dx + py * dy) / dd, 0.0, 1.0)
        d = np.sqrt((px - dx * h) ** 2 + (py - dy * h) ** 2)
        st = cum[i] + h * seg[i]
        if best_d is None:
            best_d, best_s = d, st
        else:
            m = d < best_d
            best_s = np.where(m, st, best_s)
            best_d = np.minimum(d, best_d)
    C.cache[key] = (best_d, best_s, float(cum[-1]))
    return C.cache[key]


def _bezier_dist2d(C, pts):
    return _bezier_dist_and_s(C, pts)[0]


def wing_envelope(g, X, Y, Z, D2):
    """Swept wing envelope: constant `wing_thick` across the press direction, a
    Z-depth that necks down to `wing_anchor_w` at the foot (the softening lever),
    and a self-supporting taper on the deep edge."""
    P = g.P
    wy = np.maximum(0.0, Y - g.y_root)                 # distance beyond the jacket rim
    half = 0.5 * (P["wing_anchor_w"] + (P["wing_width"] - P["wing_anchor_w"])
                  * np.clip(wy / P["wing_anchor_len"], 0.0, 1.0))
    z_c = P["wing_z_top"] - 0.5 * P["wing_width"]
    z_bot = z_c - half
    rate = math.tan(math.radians(P["wing_taper_deg"]))
    shrink = rate * np.maximum(0.0, (-Z) + z_bot + P["wing_taper"])
    # rounded (stadium) cross-section: no flat ceiling on the deep edge
    rb = P["wing_edge_round"] * np.minimum(half, 0.5 * P["wing_thick"])
    qa = D2 + shrink - (0.5 * P["wing_thick"] - rb)
    qb = np.abs(Z - z_c) - (half - rb)
    return (np.sqrt(np.maximum(qa, 0) ** 2 + np.maximum(qb, 0) ** 2)
            + np.minimum(np.maximum(qa, qb), 0.0) - rb)


def part_jacket_wing(g):
    P, cx = g.P, g.core_cx
    pts = _bezier_pts(g.wing_p0, g.wing_p1, g.wing_p2)
    clear = P["clearance"]
    thick = P["jacket_thick"]

    def fn(X, Y, Z, C):
        env = g.core_body(X, Y, Z)                 # ellipsoid only

        # ---- jacket skin: fine gyroid in an offset shell over the -Z hemisphere
        shell = np.maximum(clear - env, env - clear - thick)
        shell = I(shell, Z, X - P["jacket_x_clip"])
        # clear the canted nozzle root where it passes through the shell
        Xn, Yn, Zn = g.nz(X, Y, Z)
        nozzle_clear = U(g.nose_cone(X, Y, Z, clear),
                         cyl_x(Xn, Yn, Zn, 0, 0,
                               0.5 * P["socket_od"] + clear + 0.30, -0.5, 60.0))
        shell = S(shell, nozzle_clear)

        dz = np.maximum(0.0, -Z)
        dy = np.maximum(0.0, Y - g.y_root)
        root_dist = np.sqrt(dz * dz + dy * dy)
        t = np.clip(root_dist / P["grade_len"], 0.0, 1.0)
        wall = P["wall_root"] + (P["wall_face"] - P["wall_root"]) * t
        lat = gyroid(X, Y, Z, P["gyroid_cell"], wall)
        solid = root_dist - P["solid_root"]

        skin = np.maximum(shell, (clear + thick - P["skin_t"]) - env)
        mx = np.abs((X + 0.5 * P["perf_pitch"]) % P["perf_pitch"] - 0.5 * P["perf_pitch"])
        my = np.abs((Y + 0.5 * P["perf_pitch"]) % P["perf_pitch"] - 0.5 * P["perf_pitch"])
        perf = np.sqrt(mx * mx + my * my) - 0.5 * P["perf_dia"]
        skin = S(skin, perf)

        jacket = U(I(U(lat, solid), shell), skin)

        # ---- wing: MACRO gyroid, 1-2 cells, i.e. a doubly-curved 0.2 mm Ti sheet
        D2, S_, L_tot = _bezier_dist_and_s(C, pts)
        env_w = wing_envelope(g, X, Y, Z, D2)
        env_w = S(env_w, env - clear)

        wall_w = (P["wing_wall_root"]
                  + (P["wing_wall_tip"] - P["wing_wall_root"])
                  * np.clip(S_ / max(L_tot, 1e-6), 0.0, 1.0))
        # roll/thicken every exposed sheet edge so the wing has no knife edges
        prox = np.clip(1.0 + env_w / P["wing_edge_band"], 0.0, 1.0)
        wall_w = wall_w + (P["wing_edge_wall"] - wall_w) * prox
        sheet = gyroid(X, Y, Z, P["gyroid_cell_wing"], wall_w)

        # solid transition into the jacket rim
        wy = Y - g.y_root
        root_plug = np.maximum(wy - P["wing_root_solid"], -wy - 0.6)

        wing = I(U(sheet, root_plug), env_w)
        d = smin(jacket, wing, 0.35)

        # ---- matching magnet pockets + locating pins
        for (px_, py_) in g.jacket_mags:
            zs = _lower_z(g, px_, py_)
            d = S(d, cyl_z(X, Y, Z, px_, py_, 0.5 * P["jmag_dia"] + 0.05,
                           zs - clear - P["jmag_depth"] - 0.4, zs - clear + 0.4))
        for (px_, py_) in g.jacket_pins:
            zs = _lower_z(g, px_, py_)
            d = U(d, cyl_z(X, Y, Z, px_, py_, 0.5 * P["pin_dia"] - 0.06,
                           zs - clear - P["pin_depth"] - 1.4, zs + P["pin_depth"] - 0.05))
        d = S(d, env - clear)
        d = S(d, nozzle_clear)
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
                        g.lslot_x1,
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

def skirt_field(g, X, Y, Z):
    """Sealing skirt: an exact 35 deg outer cone with a variable wall cut from the
    inside -- structural neck, compliance groove, then the >=4 mm contact land."""
    P = g.P
    rho = np.sqrt(Y ** 2 + Z ** 2)
    ax, ar = g.skirt_root_x, g.skirt_root_r
    ux, ur = g.skirt_u
    nx, nr = -ur, ux                       # outward cone normal in the meridian plane
    px, pr = X - ax, rho - ar
    sl = px * ux + pr * ur                 # station along the slant
    dn = px * nx + pr * nr                 # +ve outside the cone face
    w = np.interp(sl, g.skirt_wall_xp, g.skirt_wall_fp)
    band = np.maximum(dn, -dn - w)         # between the cone and its inward offset
    d = np.maximum(band, np.maximum(-sl, sl - g.skirt_slant))
    # rolled rim lip: a torus tangent to the cone face, so the rim stays exactly
    # skirt_max_dia and there is no knife edge on the sealing lip
    lw = 0.5 * P["skirt_wall_land"]
    lipx = ax + ux * g.skirt_slant + nx * (-lw)
    lipr = ar + ur * g.skirt_slant + nr * (-lw)
    lip = np.sqrt((X - lipx) ** 2 + (rho - lipr) ** 2) - lw
    return U(d, lip)


def carrier_field(g, X, Y, Z, C):
    P = g.P
    body = cyl_x(X, Y, Z, 0, 0, 0.5 * P["carrier_od"], g.carrier_x0, g.carrier_x1)
    collar = cyl_x(X, Y, Z, 0, 0, g.collar_r, g.carrier_x0, g.carrier_mag_face)
    d = smin(U(body, collar), skirt_field(g, X, Y, Z), 0.35)

    # counterbore that swallows the insert's magnet flange (this is the air gap)
    d = S(d, cyl_x(X, Y, Z, 0, 0, g.cbore_r,
                   g.carrier_x0 - 1.0, g.cbore_floor))
    # annular ring-magnet pocket
    d = S(d, tube_x(X, Y, Z, 0, 0, g.carrier_mag_ri, g.carrier_mag_ro,
                    g.carrier_mag_face, g.carrier_mag_x1))
    # sliding bore -- starts at the counterbore floor, so the 0.5 mm encapsulation
    # web in front of the moving ring survives
    d = S(d, cyl_x(X, Y, Z, 0, 0, 0.5 * P["carrier_bore"],
                   g.cbore_floor, g.carrier_x1 + 1.0))
    # L-slots
    rl = 0.5 * P["carrier_bore"] - 0.05
    rh = 0.5 * P["carrier_bore"] + P["lug_h"] + 0.15
    for th in (0.0, np.pi):
        d = S(d, _carrier_lslot(g, X, Y, Z, th, rl, rh))
    return d


def part_carrier(g):
    P = g.P
    r = g.skirt_rim_r + 1.2

    def fn(X, Y, Z, C):
        return carrier_field(g, X, Y, Z, C)

    b = ((g.carrier_x0 - 1.2, g.carrier_x1 + 1.2), (-r, r), (-r, r))
    return fn, b


def mold_geom(g):
    P = g.P
    x0 = g.carrier_x0 - 3.5
    x1 = g.carrier_x1 + 3.5
    r = g.skirt_rim_r + P["mold_wall"]
    return x0, x1, r


def mold_core_field(g, X, Y, Z, C):
    """The removable core rod: defines the bore, L-slots and the magnet seat."""
    P = g.P
    x0, x1, r = mold_geom(g)
    rod = cyl_x(X, Y, Z, 0, 0, 0.5 * P["carrier_bore"], x0 - 2.0, x1 + 2.0)
    seat = cyl_x(X, Y, Z, 0, 0, g.carrier_mag_ro,
                 g.carrier_mag_face, g.carrier_mag_x1)
    plug = cyl_x(X, Y, Z, 0, 0, g.cbore_r, x0 - 2.0, g.cbore_floor)
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
                             0.5 * P["mold_vent_dia"], g.skirt_rim_r - 0.8, r + 2.0))
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

def overhang_stats(mesh, growth, sel=None):
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
    if sel is not None:
        n, a = n[sel], a[sel]
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


def wing_envelope_mesh(g, spacing=0.18):
    """The wing's macro envelope (no sheet) -- the meaningful overhang metric."""
    P = g.P
    pts = _bezier_pts(g.wing_p0, g.wing_p1, g.wing_p2)

    def fn(X, Y, Z, C):
        D2 = _bezier_dist_and_s(C, pts)[0]
        # only the free-standing span beyond the jacket rim: inboard of that the
        # wing is fused to (and supported by) the jacket shell.
        return I(wing_envelope(g, X, Y, Z, D2), g.y_root - Y)

    ymax = max(p[1] for p in pts) + 3.0
    xmin = min(p[0] for p in pts) - 4.0
    b = ((xmin, g.core_cx + P["core_rx"] + 2.0),
         (g.y_root - 1.0, ymax),
         (-P["wing_width"] - 1.5, 1.5))
    f, o, sp = evaluate(fn, b, spacing)
    return polygonise(f, o, sp)


def wing_report(g, measure=True):
    """Relative density (measured from the SDF) and a shell-bending estimate of
    tip stiffness.  See README section 6 for the model and its assumptions."""
    P = g.P
    E, NU = 110000.0, 0.31                     # Ti-6Al-4V, N/mm^2
    pts = _bezier_pts(g.wing_p0, g.wing_p1, g.wing_p2)

    # free span: arc length of the centreline beyond the jacket rim
    seg = np.linalg.norm(np.diff(pts, axis=0), axis=1)
    cum = np.concatenate([[0.0], np.cumsum(seg)])
    L_tot = float(cum[-1])
    yv = pts[:, 1]
    s_root = float(np.interp(g.y_root, yv, cum)) if yv[0] < g.y_root < yv[-1] else 0.0
    L_free = L_tot - s_root

    # ---- shell-bending stiffness, integrated per station
    n = 800
    sv = np.linspace(0.0, L_free, n)
    tw = (P["wing_wall_root"] + (P["wing_wall_tip"] - P["wing_wall_root"])
          * np.clip((s_root + sv) / L_tot, 0.0, 1.0))
    half = 0.5 * (P["wing_anchor_w"] + (P["wing_width"] - P["wing_anchor_w"])
                  * np.clip(sv / P["wing_anchor_len"], 0.0, 1.0))
    A_cut = P["wing_thick"] * 2.0 * half                    # mm^2 of cut plane
    L_A = (np.pi / 4.0) * (3.09 / P["gyroid_cell_wing"])    # sheet chord per unit area
    chord = L_A * A_cut                                     # mm of sheet in the cut
    D = E * tw ** 3 / (12.0 * (1.0 - NU ** 2))              # plate rigidity, N.mm
    EI = D * chord * P["shell_chi"]
    k = 1.0 / float(np.trapezoid((L_free - sv) ** 2 / EI, sv))

    rho_nom = 3.09 * float(tw.mean()) / P["gyroid_cell_wing"]

    out = dict(L_free=L_free, L_tot=L_tot, k=k, rho_nom=rho_nom,
               chord_root=float(chord[0]), chord_tip=float(chord[-1]),
               stations=(sv, tw, 2 * half, chord, EI))

    if measure:
        # measured as-built density over the free span, sheet region only
        fn, _ = PARTS["jacket_wing"][0](g)

        def wfn(X, Y, Z, C):
            D2 = _bezier_dist_and_s(C, pts)[0]
            return I(wing_envelope(g, X, Y, Z, D2),
                     (g.y_root + P["wing_root_solid"]) - Y)

        pts_x = [p[0] for p in pts]
        b = ((min(pts_x) - 5.0, max(pts_x) + 5.0),
             (g.y_root + P["wing_root_solid"], max(p[1] for p in pts) + 3.0),
             (-P["wing_width"] - 1.0, 0.5))
        sp = 0.055
        fenv, _o, _s = evaluate(wfn, b, sp)
        fpart, _o2, _s2 = evaluate(fn, b, sp)
        inside = fenv < 0
        out["vol_env"] = float(inside.sum()) * sp ** 3
        out["rho_meas"] = (float((inside & (fpart < 0)).sum()) / max(inside.sum(), 1))
    return out


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

# these are modelled about the nozzle axis, so their STLs are axis-aligned in the
# nozzle-local frame and get canted into place only for the assembly
NOZZLE_FRAME_PARTS = {"nozzle_insert_short", "nozzle_insert_med",
                      "nozzle_insert_long", "carrier"}


def build(name, g, voxel=None):
    fn, bounds = PARTS[name][0](g)
    budget = g.P["budget_lattice"] if PARTS[name][1] == "lattice" else g.P["budget"]
    sp = spacing_for(bounds, budget, voxel)
    t0 = time.time()
    field, origin, sp = evaluate(fn, bounds, sp)
    mesh = polygonise(field, origin, sp, tag=name)
    return mesh, sp, time.time() - t0


def mirror(mesh):
    m = mesh.copy()
    m.apply_transform(np.diag([-1.0, 1.0, 1.0, 1.0]))
    m.invert()
    return m


def mesh_holes(mesh):
    """Boundary edges (edges with a single face).  0 == the surface is closed.

    trimesh's is_watertight also demands that no edge has MORE than two faces,
    so a self-touching lattice sheet can be a closed surface and still report
    False.  Holes are what actually break a slicer; pinches do not.
    """
    e = np.sort(mesh.edges, axis=1)
    _, c = np.unique(e, axis=0, return_counts=True)
    return int((c == 1).sum()), int((c > 2).sum())


def report_row(name, mesh, sp, dt):
    b = mesh.bounds
    ext = b[1] - b[0]
    tail = ""
    if not mesh.is_watertight:
        h, pn = mesh_holes(mesh)
        tail = f"  [holes {h}, pinch edges {pn}]"
    return (f"{name:24s} {'YES' if mesh.is_watertight else 'NO ':>4s}  "
            f"vol {mesh.volume:9.1f} mm3  bbox "
            f"{ext[0]:6.2f} x {ext[1]:6.2f} x {ext[2]:6.2f}  "
            f"tris {len(mesh.faces):7d}  voxel {sp:.3f}  {dt:5.1f}s" + tail)


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
    ap.add_argument("--cant", type=float, default=None,
                    help="nozzle cant in degrees about +Y (default 45)")
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
    if args.cant is not None:
        P["nozzle_cant_deg"] = args.cant
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
    print(f"Magneto IEM generator -- magnet preset '{P['magnet_preset']}': fixed "
          f"{g.fix_od}x{g.fix_id}x{g.fix_t} / moving {g.mov_od}x{g.mov_id}x{g.mov_t} mm, "
          f"{m['material']}, rest gap {m['gap']} mm, F {m['f_lo']}->{m['f_rest']}->"
          f"{m['f_hi']} N over {P['carrier_travel']} mm, ratio {m['ratio']}")
    print(f"core  {2*P['core_rx']:.1f} x {2*P['core_ry']:.1f} x {2*P['core_rz']:.1f} mm    "
          f"carrier OD {P['carrier_od']} x {g.carrier_x0:.2f}..{g.carrier_x1:.2f} mm    "
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
        if n in DROPPED:
            cnt, vol = DROPPED[n]
            print(f"{'':26s}cleaned {cnt} marching-cubes artefacts "
                  f"({vol:.3f} mm3: sub-voxel islands + trapped-powder voids)")
        if args.ear in ("right", "both"):
            mesh.export(os.path.join(out, "right", f"{n}.stl"))
        if args.ear in ("left", "both"):
            mirror(mesh).export(os.path.join(out, "left", f"{n}.stl"))

    # ---- assembly ------------------------------------------------------
    if args.all and not args.no_assembly:
        parts = []
        for k in ASSEMBLY_PARTS:
            if k not in meshes:
                continue
            mk = meshes[k]
            if k in NOZZLE_FRAME_PARTS:            # built in the nozzle-local frame
                mk = mk.copy()
                mk.apply_transform(g.nozzle_T)
            parts.append(mk)
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
        print(f"wing SOLID envelope (bounding shape, not the part), rim-down: "
              f"worst overhang {w_rim:.1f} deg, p99 {p_rim:.1f} deg, "
              f"{100*f_rim:.1f}% of area over 45 deg")
        print(f"best sampled build direction ({d[0]:+.2f},{d[1]:+.2f},{d[2]:+.2f}): "
              f"worst {w_best:.1f} deg, p99 {p_best:.1f} deg, "
              f"{100*f_best:.1f}% of area over 45 deg")
        wr = wing_report(g, measure=("jacket_wing" in meshes))
        print(f"wing macro gyroid: cell {P['gyroid_cell_wing']} mm, wall "
              f"{P['wing_wall_root']}->{P['wing_wall_tip']} mm, envelope "
              f"{wr['L_free']:.1f} long x {P['wing_thick']} (press) x {P['wing_width']} "
              f"(deep) mm")
        print(f"  relative density  nominal 3.09t/a = {100*wr['rho_nom']:.1f}%"
              + (f"   as-built (incl. rolled edges) = {100*wr['rho_meas']:.1f}%"
                 if "rho_meas" in wr else ""))
        print(f"  shell-bending tip stiffness k = {wr['k']:.3f} N/mm  ->  "
              f"F(1.0 mm) = {wr['k']:.3f} N, F(1.5 mm) = {1.5*wr['k']:.3f} N   "
              f"(target k 0.15-0.35 N/mm; chi={P['shell_chi']})")
        if "jacket_wing" in meshes:
            jw = meshes["jacket_wing"]
            w2, f2, p2 = overhang_stats(jw, (0, 0, -1))
            sel = jw.triangles_center[:, 1] > g.y_root
            w3, f3, p3 = overhang_stats(jw, (0, 0, -1), sel=sel)
            print(f"as-built wing sheet (y > rim), rim-down: worst {w3:.1f} deg, "
                  f"p99 {p3:.1f} deg, {100*f3:.1f}% of area over 45 deg")
            print(f"whole jacket+wing, rim-down: worst {w2:.1f} deg, "
                  f"p99 {p2:.1f} deg, {100*f2:.1f}% of area over 45 deg")

    # ---- protrusion stack along the nozzle axis (docs/TRYON_REPORT.md rec 1)
    rigid = [meshes[k] for k in ("core", "faceplate", "jacket_wing") if k in meshes]
    if rigid:
        lo = min(float(((mm.vertices - g.nozzle_base) @ g.n_ax).min()) for mm in rigid)
        print("-" * 100)
        print(f"nozzle cant {P['nozzle_cant_deg']:.0f} deg about +Y   "
              f"axis {tuple(round(float(v), 3) for v in g.n_ax)}   "
              f"nozzle base {tuple(round(float(v), 2) for v in g.nozzle_base)}")
        print(f"  rigid body along the nozzle axis: {lo:.2f} .. {g.carrier_x0:.2f} mm "
              f"= {g.carrier_x0 - lo:.2f} mm stack   "
              f"(uncanted baseline 23.40 mm, TRYON_REPORT.md rec 1)")
        print(f"  seal plane at {g.carrier_x1:.2f} mm; everything behind it spans "
              f"{g.carrier_x1 - lo:.2f} mm along the axis")

    # ---- skirt contact land + pressure budget (docs/MECH_VALIDATION.md JOB 2)
    if args.all or "carrier" in meshes:
        fl = math.radians(P["skirt_flare_deg"])
        w = g.skirt_land_w
        print("-" * 100)
        print(f"skirt contact land: {w:.2f} mm slant width at {P['skirt_flare_deg']}deg, "
              f"Ø{g.skirt_land_d0:.1f} -> Ø{P['skirt_max_dia']:.1f} mm, wall "
              f"{P['skirt_wall_land']} mm, compliance groove {P['skirt_wall_hinge']} mm "
              f"x {P['skirt_hinge_w']} mm behind it   (FEA minimum 4.0 mm)")
        f_max = g.mag["f_lo"]
        cells = []
        for dia in (10.0, 13.0, 16.0, 19.0):
            area = math.pi * dia * w * math.sin(fl)      # mm^2, cone-normal
            kpa = 1000.0 * f_max / area
            tag = ("comfortable" if kpa <= 2.15 else
                   "borderline" if kpa <= 4.27 else "TOO MUCH")
            cells.append(f"Ø{dia:.0f}: {kpa:4.2f} kPa {tag}")
        print(f"  cone-normal pressure at F_max {f_max} N   " + " | ".join(cells))

    # verify what a slicer actually sees: re-read the exported binary STL
    print("-" * 100)
    holed = []
    for n in meshes:
        rm = trimesh.load(os.path.join(out, "right", f"{n}.stl"), force="mesh")
        if rm.is_watertight:
            continue
        h, pn = mesh_holes(rm)
        print(f"  re-read {n}: holes {h}, pinch edges {pn}"
              + ("  (closed surface; pinches only)" if h == 0 else "  *** OPEN ***"))
        if h:
            holed.append(n)
    bad = [n for n, mm in meshes.items() if not mm.is_watertight] + holed
    print(f"{len(meshes)} parts written to {out}/  "
          + ("ALL WATERTIGHT" if not bad else f"NOT CLOSED: {sorted(set(bad))}"))
    return 1 if bad else 0


if __name__ == "__main__":
    sys.exit(main())
