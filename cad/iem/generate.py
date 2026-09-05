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

    # ---- posterior-inferior corner roll (docs/TRYON_REPORT.md v2) ----------
    corner_chamfer=4.00,       # mm cut back along corner_dir at the faceplate corner
    corner_dir=(-0.41, -0.65, 0.49),   # protrusion-sensitivity diagonal, (x, y, z)
    corner_roll=1.50,          # mm blend radius -- a roll, not a knife edge
    corner_min_wall=0.90,      # mm of Ti kept over every internal pocket
    corner_ramp=1.20,          # mm over which the cut opens up above the magnet band

    # ---- body trim (docs/TRYON_REPORT.md, recalibrated 13e75b3) ------------
    body_trim_mm=0.00,         # mm of PROTRUSION to remove by shrinking the shell.
                               #   0 = ACCEPT: protrusion settled 2026-08-31, no
                               #   driver change and no further body shrink. The
                               #   solver below is kept for the record; see README 2c.
    body_trim_w=(0.41, 0.65, 0.49),   # mm of protrusion removed per mm cut on X/Y/Z
    body_trim_min_wall=0.90,   # mm of Ti kept around the driver pocket while trimming
    body_trim_min_front_cc=0.50,   # cc of acoustic void that must survive the trim
    body_trim_keep_magnets=True,   # refuse an X trim that worsens the faceplate rim
    body_trim_force=None,      # (tx, ty, tz) override, set by solve_body_trim()
    notch_measured_reach=None, # filled in by calibrate_notch()
    notch_calib_hist=None,

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
    # ---- intertragic-notch sector: COMPLIANCE, not reach ------------------
    # v3 added radial flare; v4 (2d2c491) showed that is net negative -- a radial
    # extension commits the lip toward the notch opening, where there is no flesh
    # to seal against.  The rim is plain O19 all round again; the sector now just
    # gets a thinner land on a longer hinge so the lip can drape onto the
    # cartilage flanking the notch.
    notch_compliance=True,     # enable the compliant inferior sector
    notch_sector_ext=0.00,     # mm radial flare -- REVERTED to 0, do not re-enable
    notch_sector_deg=90.0,     # deg of skirt perimeter treated as the notch sector
    notch_sector_center_deg=180.0,  # deg from +Y_local; 180 = inferior
    notch_sector_trans_deg=20.0,    # deg of smooth azimuthal blend at each edge
    notch_hinge_wall=0.15,     # mm hinge wall in the sector (vs skirt_wall_hinge)
    notch_hinge_w=1.60,        # mm hinge slant width in the sector -- the free
                               #   length the lip rotates over
    notch_sector_wall=0.22,    # mm land wall inside the sector (vs skirt_wall_land)

    # ---- bell tip (2026-09-01 decision, artifact "Aperture Tip") -------------
    # The mag-float carrier + drape skirt above is CANCELLED (8/31 reset).  It is
    # kept behind tip_style="carrier" so the record stays reproducible; the
    # shipped tip is the bell: one nose cone for every size, a hollow rolled lip
    # on an oval footprint in S/M/L, and a thin web between them.
    tip_style="bell",          # "bell" (default) | "carrier" (legacy mag-float)
    bell_asm_size="M",         # which lip size the assembly shows
    bell_nozzle_od=4.00,       # mm Ti nozzle tube the tip seats on (insert 'bell')
    bell_bore=2.60,            # mm sound bore through insert and tip
    bell_nose_tip_d=5.00,      # mm nose cone tip
    bell_nose_base_d=12.00,    # mm nose cone base -- the rim stop
    bell_nose_len=2.40,        # mm axial cone length -> 55.6 deg half-angle
    bell_base_land=0.40,       # mm cylindrical O12 land behind the cone; the web fuses here
    bell_sleeve_wall=0.50,     # mm silicone over the Ti tube behind the cone
    bell_noz_recess=1.00,      # mm the Ti tube stops short of the tip face
    bell_groove_w=1.00,        # mm retaining groove in the Ti tube (ridge in the tip)
    bell_groove_d=0.25,        # mm groove depth
    bell_ridge_interf=0.05,    # mm the silicone ridge is undersize on the groove floor
    bell_lip_tube_d=3.00,      # mm rolled-lip tube diameter
    bell_lip_wall=0.40,        # mm rolled-lip wall (Shore 00-30)
    bell_lip_hollow=True,      # True: C-section rolled edge; False: solid 00-30 bead
    bell_lip_slit_deg=90.0,    # deg of the tube opened toward the axis so the core ring pulls
    bell_sizes=dict(XS=(11.5, 9.5), S=(14.0, 12.0), M=(17.0, 14.0), L=(21.0, 16.0)),  # outer H x W, mm
    # ladder re-spaced 9/4 against 102 scanned apertures (was S 14x11 / M 17x13 / L 20x15):
    # XS 50 / S 18 / M 21 / L 9 / too-big 4.  S and M keep their heights, +1 mm width
    # buys the rocking room the O3 tube needs over the O5 sleeve.
    # XS added 2026-09-04: 63/102 scanned apertures took S under a 1 mm clearance rule
    # (median aperture 7.3 x 4.1 mm).  A O3 lip tube cannot get under W=11 without
    # the lip inner edge hitting the O5 sleeve, so XS runs a O2 tube.
    bell_tube_by_size=dict(XS=2.0),   # mm per-size lip tube override (default bell_lip_tube_d)
    bell_lip_x=2.70,           # mm from the tip's proximal face to the lip centre plane
    bell_web_t=0.50,           # mm web between nose and lip
    bell_web_ax=1.20,          # mm axial run of the web's S-curve
    bell_web_stub=0.30,        # mm straight stub where the web leaves the lip tube
    bell_fillet=0.40,          # mm web-to-nose fillet
    bell_ant_center_deg=90.0,  # deg from +Y toward +Z: +Z_local is anterior (tragus)
    bell_ant_deg=90.0,         # deg of lip treated as the anterior sector
    bell_ant_wall=0.30,        # mm lip wall in the anterior sector
    bell_ant_free=1.00,        # mm extra web free length (lip moved proximally) there
    bell_inf_center_deg=180.0, # deg: -Y_local is inferior (intertragic notch)
    bell_inf_deg=90.0,
    bell_inf_ext=1.50,         # mm radial lip extension over the notch
    bell_sector_trans_deg=20.0,
    bell_vent_dia=0.80,        # mm vent channel at the tip/nozzle interface
    bell_vent_az_deg=180.0,    # deg: the vent runs down the inferior side, exits under the lip
    bell_disc_t=2.00,          # mm mould-core base disc behind the tip's proximal face
    budget_fine=8.0e6,         # voxels for the bell tip and its moulds (0.3-0.4 mm walls)

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

    # ---- wing mechanism: THREE RADIAL MAG-PLUNGERS -------------------------
    # v5 (docs/MECH_VALIDATION.md 8eb6ac1) replaced the Ti spring wing with a
    # mag-plunger; the plan then moved from one rail to three independent radial
    # plungers, one per contact site, each with its own aim, cam preset and stops.
    # "gyroid" keeps the old compliant sheet as a legacy option.
    wing_style="plungers",     # "plungers" | "gyroid"
    # site name -> aim unit vector (x, y, z) in design coords; +Y superior,
    # -Y inferior, -Z ear-facing.  Normalised at load.
    plunger_aims=(("cymba",              (0.30,  0.94, -0.15)),
                  ("antihelix_undercut", (-0.45, 0.85, -0.28)),
                  ("tragus_inner",       (0.82, -0.17, -0.54))),
    # per-leg enable.  tragus_inner is defined but OFF: its aim clashes with the
    # nozzle stack (see README 6), so the shipped build is the two-leg + interlock
    # variant for stability testing.
    plunger_enable={"cymba": True, "antihelix_undercut": True,
                    "tragus_inner": False},
    cymba_lip_bias=7.0,        # deg the cymba aim rotates toward +Y, so the pad
                               #   lands UNDER the cymba's overhanging lip
    cymba_pad_extra=1.50,      # mm of extra pad diameter on the cymba leg only
    plunger_pad_roll=0.40,     # mm rolled edge on the pad shoulder
    plunger_cam_ext_sites=("tragus_inner",),   # legs that get the extended cam
    plunger_cam_ext_steps=6,   # detents on the extended cam
    plunger_cam_ext_range=4.5, # mm of coarse engagement on the extended cam
    plunger_mag_od=5.00,       # mm plunger ring OD   (N35, 5 x 2.5 x 1)
    plunger_mag_id=2.50,       # mm plunger ring ID   -- the guide pin runs through it
    plunger_mag_t=1.00,        # mm plunger ring thickness
    plunger_gap=2.75,          # mm rest gap between the plunger faces
    plunger_travel=0.75,       # mm dynamic travel each way, hard-limited by stops
    plunger_pin_od=2.00,       # mm guide pin OD
    plunger_pin_sleeve=0.30,   # mm polymer sleeve wall in the jacket bore
    plunger_foot_od=9.00,      # mm piston / pad diameter
    plunger_plate=0.80,        # mm Ti plate over the moving magnet
    plunger_pad_t=1.00,        # mm soft silicone contact pad
    plunger_rocker=0.35,       # mm crown on the pad -- the slight rocker
    plunger_boss_od=9.60,      # mm jacket boss OD
    plunger_boss_h=5.50,       # mm boss height: cam + fixed magnet + pin sleeve
    plunger_cam_h=3.20,        # mm cam preset ring height (0.2 base + 3.0 of range)
    plunger_cam_steps=4,       # detents on the cam preset ring
    plunger_cam_range=3.0,     # mm of coarse engagement the cam covers
    plunger_cam_od=7.40,       # mm cam preset ring OD

    # ---- cable exit boot ----------------------------------------------------
    cable_exit="up_back",      # "up_back" | "back" | "none"
    cable_boot_len=6.00,       # mm strain-relief stub
    cable_boot_od0=5.50,       # mm boot OD at the shell
    cable_boot_od1=3.50,       # mm boot OD at the tip
    cable_bore=1.80,           # mm cable bore through the boot
    cable_boot_angle=35.0,     # deg the boot rakes up (+Y) from straight back

    # ---- wing: MACRO-scale gyroid shell (a compliant doubly-curved sheet) --
    gyroid_cell_wing=12.00,    # mm wing unit cell -- 1-2 cells across the envelope
    wing_wall_root=0.20,       # mm sheet wall at the root
    wing_wall_tip=0.20,        # mm sheet wall at the tip
    wing_edge_wall=0.40,       # mm rolled/thickened rim on exposed sheet edges
    wing_edge_band=0.70,       # mm over which the wall ramps up to the rolled edge
    wing_root_solid=1.20,      # mm of solid Ti transition into the jacket rim
    wing_len=14.0,             # mm nominal wing length along its centreline
    wing_thick=7.00,           # mm envelope across the press direction (in XY)
    wing_width=5.00,           # mm envelope depth in Z, into the concha
    wing_anchor_w=1.40,        # mm Z-width at the anchor (necked foot; softens the
                               #   wing).  Narrowed at v4 to hold k in band after the
                               #   2.75 mm shortening stiffened it (k ~ 1/L^3).
    wing_anchor_len=7.00,      # mm over which the Z-width opens anchor_w -> wing_width
    wing_z_top=-0.20,          # mm top of the wing, just under the parting plane
    wing_taper_deg=40.0,       # deg overhang of the wing's deep-edge taper
    wing_edge_round=0.85,      # fraction of the half-section used as a corner radius
    wing_taper=1.60,           # mm of tapered depth on the wing's deep edge
    wing_shorten=2.75,         # mm taken off the free span (v4: 44% overpressed)
    wing_splay_deg=-5.0,       # deg the whole wing rotates about its root in XY,
                               #   which splays the press direction by the same angle
    wing_rise=11.0,            # mm the tip lands above the core rim (pre-shorten)
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
    socket_z=-0.60,            # mm socket pocket centre Z; keeps the full 2.8 mm
                               #     height under the faceplate parting plane
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


def ssub(a, b, k):
    """Smooth subtraction: remove b from a with a blend radius of roughly k."""
    if k <= 1e-9:
        return np.maximum(a, -b)
    h = np.clip(0.5 - 0.5 * (a + b) / k, 0.0, 1.0)
    return a * (1.0 - h) + (-b) * h + k * h * (1.0 - h)


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

        # ---- body trim ------------------------------------------------------
        # Distribute the requested protrusion reduction across X/Y/Z in proportion
        # to each axis's efficiency, then clamp each axis to what the driver
        # pocket and the faceplate/driver stack actually allow.
        w = np.array(P["body_trim_w"], dtype=float)
        pr0 = 0.5 * (P["driver_carrier_od"] + P["driver_pocket_clear"])
        alpha = P["body_trim_mm"] / float(np.dot(w, w))
        req = alpha * w                                   # full-extent cuts, mm
        mw = P["body_trim_min_wall"]
        zfloor = max(P["driver_pocket_depth"] - P["faceplate_z"] + mw,
                     P["faceplate_z"] + 1.5)
        cap = np.array([2.0 * (P["core_rx"] - (pr0 + mw)),
                        2.0 * (P["core_ry"] - (pr0 + mw)),
                        2.0 * (P["core_rz"] - zfloor)])
        cap = np.maximum(cap, 0.0)
        # the faceplate magnets live on the rim between the driver pocket and the
        # shell; if that rim is already tight, an X trim cannot make it worse
        st = math.hypot(6.5, 2.9)
        if P["body_trim_keep_magnets"] and \
                st - pr0 - 0.5 * P["jmag_dia"] < 0.30:
            cap[0] = 0.0
            self.trim_block_x = (f"faceplate magnet rim already "
                                 f"{st - pr0 - 0.5 * P['jmag_dia']:.2f} mm")
        else:
            self.trim_block_x = None
        got = np.minimum(req, cap)
        if P["body_trim_force"] is not None:
            got = np.array(P["body_trim_force"], dtype=float)
        self.trim_req, self.trim_cap, self.trim_got = req, cap, got
        self.trim_protrusion_req = float(np.dot(w, req))
        self.trim_protrusion_got = float(np.dot(w, got))
        rx = P["core_rx"] - 0.5 * got[0]
        ry = P["core_ry"] - 0.5 * got[1]
        rz = P["core_rz"] - 0.5 * got[2]

        self.core_cx = -rx
        self.core_r = (rx, ry, rz)
        self.core_c = (self.core_cx, 0.0, 0.0)
        self.inner_r = tuple(v - P["core_wall"] for v in self.core_r)
        self.core_rx, self.core_ry, self.core_rz = rx, ry, rz
        self.z_cut = P["faceplate_z"]

        self.pocket_r = 0.5 * (P["driver_carrier_od"] + P["driver_pocket_clear"])
        self.pocket_z1 = self.z_cut
        self.pocket_z0 = self.z_cut - P["driver_pocket_depth"]
        self.front_wall_x = self.core_cx + self.pocket_r

        # ---- cable exit boot ------------------------------------------------
        if P["cable_exit"] == "none":
            self.boot = None
        else:
            ang = math.radians(P["cable_boot_angle"]
                               if P["cable_exit"] == "up_back" else 0.0)
            bd = np.array([-math.cos(ang), math.sin(ang), 0.0])
            a0 = np.array([-2.0 * P["core_rx"] + 1.2, 0.0, P["socket_z"]])
            self.boot = (a0, a0 + bd * P["cable_boot_len"],
                         0.5 * P["cable_boot_od0"], 0.5 * P["cable_boot_od1"])

        # ---- three radial mag-plungers -------------------------------------
        self.plungers = []
        self.plungers_all = []
        rr_ = np.array(self.core_r)
        for name, aim in P["plunger_aims"]:
            a = np.array(aim, dtype=float)
            a /= np.linalg.norm(a)
            if name == "cymba" and P["cymba_lip_bias"]:
                # rotate the aim toward +Y in the plane it shares with +Y, so the
                # pad lands under the cymba's overhanging lip and can interlock
                yv = np.array([0.0, 1.0, 0.0])
                perp = yv - a * np.dot(yv, a)
                if np.linalg.norm(perp) > 1e-6:
                    perp /= np.linalg.norm(perp)
                    t = math.radians(P["cymba_lip_bias"])
                    a = a * math.cos(t) + perp * math.sin(t)
                    a /= np.linalg.norm(a)
            ext = name in P["plunger_cam_ext_sites"]
            cam_h = 0.2 + (P["plunger_cam_ext_range"] if ext
                           else P["plunger_cam_range"])
            boss_h = max(P["plunger_boss_h"],
                         cam_h + P["plunger_mag_t"] + 1.20)
            # point on the core where the outward normal is a, then out to the
            # jacket's outer surface, then up the boss to the mount face
            surf = np.array(self.core_c) + (rr_ ** 2 * a) / np.linalg.norm(rr_ * a)
            base = surf + a * (P["clearance"] + P["jacket_thick"])
            mount = base + a * boss_h
            ref = np.array([0.0, 0.0, 1.0])
            if abs(np.dot(ref, a)) > 0.9:
                ref = np.array([1.0, 0.0, 0.0])
            u = np.cross(a, ref); u /= np.linalg.norm(u)
            v = np.cross(a, u)
            self.plungers.append(dict(
                name=name, aim=a, base=base, mount=mount, u=u, v=v,
                ext=ext, cam_h=cam_h, boss_h=boss_h,
                cam_steps=(P["plunger_cam_ext_steps"] if ext
                           else P["plunger_cam_steps"]),
                cam_range=(P["plunger_cam_ext_range"] if ext
                           else P["plunger_cam_range"]),
                pad_extra=(P["cymba_pad_extra"] if name == "cymba" else 0.0),
                enabled=bool(P["plunger_enable"].get(name, True))))
        self.plungers_all = list(self.plungers)
        self.plungers = [q for q in self.plungers_all if q["enabled"]]
        # canonical axial stations, measured from the mount face along the aim
        self.pl_mag_fix = -P["plunger_mag_t"]                      # fixed ring back
        self.pl_mag_mov = P["plunger_gap"]                         # moving ring face
        self.pl_mag_mov1 = self.pl_mag_mov + P["plunger_mag_t"]
        self.pl_plate1 = self.pl_mag_mov1 + P["plunger_plate"]
        self.pl_pad1 = self.pl_plate1 + P["plunger_pad_t"]
        self.pl_depth_stack = self.pl_mag_mov1 - self.pl_mag_fix   # = 4.75
        self.pl_stop_in = self.pl_mag_mov - P["plunger_travel"]
        self.pl_stop_out = self.pl_mag_mov + P["plunger_travel"]

        # ---- posterior-inferior corner roll --------------------------------
        # The v2 try-on found the worst-protruding point is a corner, not a face:
        # median (-11.5, -4.1, +3.7), i.e. the -X/-Y/+Z octant of the faceplate,
        # with protrusion sensitivity 0.65 / 0.49 / 0.41 mm per mm on +Y / +Z / +X.
        # Cut perpendicular to that diagonal, which is the most material-efficient
        # direction available.
        n = np.array(P["corner_dir"], dtype=float)
        self.corner_n = n / np.linalg.norm(n)
        cn = self.corner_n
        rr = np.array(self.core_r)
        self.corner_h = float(np.dot(self.core_c, cn) + np.linalg.norm(rr * cn))

        # How deep can the cut go before it breaks into something?  Only CLOSED
        # internals bound it: the driver pocket and the acoustic cavity.  The
        # 2-pin socket pocket is an intentional opening in the shell, so "keep
        # 0.9 mm of Ti over it" is not a meaningful constraint.
        sup = []
        pr, pz0, pz1 = self.pocket_r, self.pocket_z0, self.pocket_z1
        sup.append(self.core_cx * cn[0] + pr * math.hypot(cn[0], cn[1])
                   + max(pz0 * cn[2], pz1 * cn[2]))                    # driver pocket
        ir = np.array(self.inner_r)
        sup.append(float(np.dot(self.core_c, cn) + np.linalg.norm(ir * cn)))   # cavity
        self.corner_sup = max(sup)
        self.corner_c_core = float(np.clip(
            self.corner_h - self.corner_sup - P["corner_min_wall"],
            0.0, P["corner_chamfer"]))
        # the faceplate's own magnet pockets sit in z = [z_cut, z_cut + jmag_depth];
        # the full-depth cut only opens above them
        self.corner_z_lo = self.z_cut

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
        self.nozzle_t_exit = 1.0 / math.sqrt((ca / self.core_rx) ** 2
                                             + (sa / self.core_rz) ** 2)
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
        # sector profile: same land start, but a thinner hinge over a longer run,
        # so the lip has more free length to rotate through
        n_h1 = s_h1
        n_h0 = max(0.05, n_h1 - P["notch_hinge_w"])
        n_n1 = max(0.02, n_h0 - 0.20)
        self.notch_wall_xp = [0.0, n_n1, n_h0, n_h1, s_l0, self.skirt_slant]
        self.notch_wall_fp = [P["skirt_wall_neck"], P["skirt_wall_neck"],
                              P["notch_hinge_wall"], P["notch_hinge_wall"],
                              P["notch_sector_wall"], P["notch_sector_wall"]]
        self.notch_hinge_span = n_h1 - n_h0
        self.skirt_land_x0 = self.skirt_root_x + s_l0 * self.skirt_u[0]
        self.skirt_land_d0 = 2.0 * (self.skirt_root_r + s_l0 * self.skirt_u[1])

        # ---- bell tip: axial stations in the nozzle-local frame ----------------
        Rt = 0.5 * P["bell_lip_tube_d"]
        self.bell_Rt = Rt                       # shared tube (nose/web stations)
        self.bell_Rt_by = {k: 0.5 * P.get("bell_tube_by_size", {}).get(k, P["bell_lip_tube_d"])
                           for k in P["bell_sizes"]}   # per-size lip tube
        self.bell_x0 = self.socket_x1              # proximal face seats on the socket face
        self.bell_seat_r = 0.5 * P["bell_nozzle_od"]
        self.bell_sleeve_r = self.bell_seat_r + P["bell_sleeve_wall"]
        # the vent groove (O0.8, centred 0.15 mm outside the seat) needs a full
        # wall over it, carried by a local rib on the sleeve at the vent azimuth
        self.bell_rib_r = (self.bell_seat_r + 0.15 + 0.5 * P["bell_vent_dia"]
                           + P["bell_sleeve_wall"])
        self.bell_lip_x = self.bell_x0 + P["bell_lip_x"]
        self.bell_wl_x = self.bell_lip_x + Rt + P["bell_web_stub"]
        self.bell_cb_x = self.bell_wl_x + P["bell_web_ax"]          # cone base = rim stop
        self.bell_land_x0 = self.bell_cb_x - P["bell_base_land"]
        self.bell_web_x1 = self.bell_cb_x - 0.5 * P["bell_base_land"]
        self.bell_base_r = 0.5 * P["bell_nose_base_d"]
        self.bell_nose_r = 0.5 * P["bell_nose_tip_d"]
        self.bell_web_r1 = self.bell_base_r - 0.35   # web centreline enters the land here
        self.bell_tip_x = self.bell_cb_x + P["bell_nose_len"]
        self.bell_noz_end = self.bell_tip_x - P["bell_noz_recess"]
        self.bell_groove_x0 = self.bell_cb_x - 0.2
        self.bell_groove_x1 = self.bell_groove_x0 + P["bell_groove_w"]
        self.bell_ridge_r = self.bell_seat_r - P["bell_groove_d"] - P["bell_ridge_interf"]
        self.bell_half_angle = math.degrees(math.atan2(
            self.bell_base_r - self.bell_nose_r, P["bell_nose_len"]))
        # lip centreline semi-axes (a along Y = superior-inferior, b along Z)
        self.bell_lip_axes = {k: (0.5 * h - self.bell_Rt_by[k], 0.5 * w - self.bell_Rt_by[k])
                              for k, (h, w) in P["bell_sizes"].items()}
        self.bell_lip_rmax = {k: a + self.bell_Rt_by[k] + P["bell_inf_ext"]
                              for k, (a, b) in self.bell_lip_axes.items()}
        self.bell_disc_r = {k: a + P["bell_inf_ext"] + 0.5
                            for k, (a, b) in self.bell_lip_axes.items()}

        if P["tip_style"] == "bell":
            self.tip_x0, self.tip_x1 = self.bell_x0, self.bell_tip_x
            self.seal_x = self.bell_lip_x + Rt           # lip's distal face
            self.tip_protrusion = self.bell_tip_x
            self.tip_protrusion_max = self.bell_tip_x
        else:
            self.tip_x0, self.tip_x1 = self.carrier_x0, self.carrier_x1
            self.seal_x = self.carrier_x1
            self.tip_protrusion = self.carrier_x1             # from the core face
            self.tip_protrusion_max = self.carrier_x1 + P["carrier_travel"]

        # ---- faceplate / jacket magnet stations
        outer_half_x = self.core_rx * math.sqrt(max(1e-9, 1 - (self.z_cut / self.core_rz) ** 2))
        self.fp_mag_off = 0.5 * (self.pocket_r + 0.5 * P["jmag_dia"] + 0.35
                                 + outer_half_x - 0.5 * P["jmag_dia"] - 0.35)
        self.fp_mag_off = max(self.fp_mag_off, self.pocket_r + 0.5 * P["jmag_dia"] + 0.3)
        # A diagonal pair, not an axial one: an axial -X station sits inside the
        # posterior-inferior corner roll.  Both are >= pocket_r + jmag_dia/2 from
        # the driver pocket centre and on the wide part of the faceplate rim.
        self.fp_mags = [(self.core_cx - 6.5, 2.9), (self.core_cx + 6.5, -2.9)]

        # jacket magnets (3) + locating pins (2) on the -Z hemisphere, as (x, y)
        cx = self.core_cx
        self.jacket_mags = [(cx - 6.0, 0.0), (cx + 0.5, 4.8), (cx + 0.5, -4.8)]
        self.jacket_pins = [(cx - 3.0, 4.0), (cx - 3.0, -4.0)]

        # ---- wing centreline (quadratic Bezier in XY) ----------------------
        # v4 (2d2c491): 44% of ears are OVERPRESSED by the wing (median -2.66 mm)
        # and only one is short, so the wing is shortened and the whole blade is
        # rotated about its root, which splays the press direction by the same
        # angle.
        y_root = self.core_ry + P["clearance"]
        back = math.radians(P["wing_back_deg"])
        splay = math.radians(P["wing_splay_deg"])
        dx0 = self.core_cx + P["wing_root_dx"]

        def _pts(rise):
            p0 = np.array([dx0, y_root - 1.6])
            p1 = np.array([dx0, y_root + 0.62 * rise])
            p2 = np.array([dx0 - rise * math.tan(back), y_root + rise])
            c, sn = math.cos(splay), math.sin(splay)
            R = np.array([[c, -sn], [sn, c]])
            return p0, p0 + R @ (p1 - p0), p0 + R @ (p2 - p0)

        def _free(rise):
            q = _bezier_pts(*_pts(rise))
            seg = np.linalg.norm(np.diff(q, axis=0), axis=1)
            cum = np.concatenate([[0.0], np.cumsum(seg)])
            yv = q[:, 1]
            s_r = float(np.interp(y_root, yv, cum)) if yv[0] < y_root < yv[-1] else 0.0
            return float(cum[-1]) - s_r

        target = _free(P["wing_rise"]) - P["wing_shorten"]
        lo, hi = 1.0, P["wing_rise"]
        for _ in range(40):
            mid = 0.5 * (lo + hi)
            if _free(mid) < target:
                lo = mid
            else:
                hi = mid
        self.wing_rise = 0.5 * (lo + hi)
        self.wing_p0, self.wing_p1, self.wing_p2 = (tuple(v) for v in _pts(self.wing_rise))
        self.wing_L_free = _free(self.wing_rise)

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
        if P["tip_style"] == "carrier" and need > self.carrier_x1 - 0.2:
            self.warnings.append(
                f"L-slot needs the carrier to reach x={need + 0.2:.2f} mm but it ends at "
                f"{self.carrier_x1:.2f}; raise carrier_len to "
                f"{need + 0.2 - self.carrier_x0:.2f} mm")
        if P["tip_style"] == "bell":
            if self.bell_groove_x1 > self.bell_noz_end - 0.4:
                self.warnings.append(
                    f"bell retaining groove ends at {self.bell_groove_x1:.2f} mm, only "
                    f"{self.bell_noz_end - self.bell_groove_x1:.2f} mm of Ti tube past it")
            if P["bell_bore"] < P["nozzle_bore"] - 1e-9:
                self.warnings.append(
                    f"bell bore O{P['bell_bore']} necks the core's O{P['nozzle_bore']} bore "
                    f"at the insert socket (spec: O4 Ti nozzle, O2.6 bore)")
            for k, (a, b) in self.bell_lip_axes.items():
                Rt = self.bell_Rt_by[k]
                gap = b - Rt - self.bell_sleeve_r
                if gap < 0.3:
                    self.warnings.append(
                        f"bell {k}: lip inner edge is {gap:+.2f} mm from the sleeve on the "
                        f"minor axis (W {P['bell_sizes'][k][1]:.0f} = 2 x {Rt:.1f} tube + "
                        f"O{2 * self.bell_sleeve_r:.0f} sleeve) -- no rocking room there")

    def corner_cut(self, X, Y, Z):
        """Half-space to remove at the posterior-inferior corner.

        Depth ramps with Z: below the faceplate magnet band it is limited to what
        the driver pocket / cavity / socket sleeve allow; above it, the full
        corner_chamfer.
        """
        P = self.P
        cn = self.corner_n
        t = np.clip((Z - self.corner_z_lo) / P["corner_ramp"], 0.0, 1.0)
        c = self.corner_c_core + (P["corner_chamfer"] - self.corner_c_core) \
            * (t * t * (3.0 - 2.0 * t))
        # negative INSIDE the region to remove, which is what ssub() expects
        return (self.corner_h - c) - (X * cn[0] + Y * cn[1] + Z * cn[2])

    def pl_coords(self, pl, X, Y, Z):
        """(axial station from the mount face, radial distance from the axis)."""
        a = pl["aim"]
        m = pl["mount"]
        dx, dy, dz = X - m[0], Y - m[1], Z - m[2]
        sA = dx * a[0] + dy * a[1] + dz * a[2]
        r2 = (dx * dx + dy * dy + dz * dz) - sA * sA
        u, v = pl["u"], pl["v"]
        th = np.arctan2(dx * v[0] + dy * v[1] + dz * v[2],
                        dx * u[0] + dy * u[1] + dz * u[2])
        return sA, np.sqrt(np.maximum(r2, 0.0)), th

    def nz(self, X, Y, Z):
        """World -> nozzle-local coordinates (the canted frame)."""
        ca, sa = math.cos(self.cant), math.sin(self.cant)
        bx, by, bz = self.nozzle_base
        px, py, pz = X - bx, Y - by, Z - bz
        return px * ca - pz * sa, py, px * sa + pz * ca

    def bone_boss(self, X, Y, Z, grow=0.0):
        P = self.P
        return rbox(X, Y, Z, (self.core_cx, -(self.core_ry + 0.55), -1.0),
                    (0.5 * P["bone_w"] + 1.1 + grow, 1.5 + grow,
                     0.5 * P["bone_h"] + 1.1 + grow), 0.7)

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
        sx0 = cx - g.core_rx - 0.6
        sx1 = sx0 + P["socket_d"] + 0.6
        sl_h = (0.5 * (sx1 - sx0), 0.5 * P["socket_w"] + 0.7, 0.5 * P["socket_h"] + 0.7)
        sl_c = (0.5 * (sx0 + sx1), 0.0, P["socket_z"])
        sleeve = I(rbox(X, Y, Z, sl_c, sl_h, 0.4), outer, Z - g.z_cut)
        d = U(d, sleeve)

        # bone-sensor boss on the -Y (tragus) flank
        bone_y = -(g.core_ry + 0.55)
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
                               (fv[0], -g.core_ry - 3.0, -g.core_rz - 3.0),
                               0.5 * P["vent_dia"]))
        void = U(void, capsule(X, Y, Z,
                               (cx - 3.0, 0.0, -1.0),
                               (cx - 5.5, -g.core_ry - 3.0, -g.core_rz - 3.0),
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
        for (mx, my) in g.fp_mags:
            void = U(void, cyl_z(X, Y, Z, mx, my, 0.5 * P["jmag_dia"] + 0.05,
                                 g.z_cut - P["jmag_depth"], g.z_cut + 1.0))

        # fixed-magnet counterbore lives in the nozzle INSERT, not here (v0.2 change)

        d = S(d, void)

        # gasket groove around the parting line at z = 0
        groove = I(np.abs(Z) - 0.5 * P["gasket_w"],
                   -(outer + P["gasket_d"]),
                   X - (cx + g.core_rx - 2.0))
        d = S(d, groove)
        d = ssub(d, g.corner_cut(X, Y, Z), P["corner_roll"])

        # ---- cable exit boot: strain relief so the try-on can score clearance
        if g.boot is not None:
            a0, a1, r0, r1 = g.boot
            n = 20
            boot = None
            for i in range(n + 1):
                t = i / n
                c = a0 + (a1 - a0) * t
                rr = r0 + (r1 - r0) * t
                sp = np.sqrt((X - c[0]) ** 2 + (Y - c[1]) ** 2 + (Z - c[2]) ** 2) - rr
                boot = sp if boot is None else U(boot, sp)
            d = smin(d, boot, 0.8)
            d = S(d, capsule(X, Y, Z, a0 - (a1 - a0) * 0.6, a1 + (a1 - a0) * 0.15,
                             0.5 * P["cable_bore"]))
        return d

    tip = g.nozzle_base + g.n_ax * (g.stub_x1 + 0.6)
    bx = cx - g.core_rx - 1.6
    by = g.core_ry + 1.6
    if g.boot is not None:
        bx = min(bx, float(g.boot[1][0]) - 3.4)
        by = max(by, float(g.boot[1][1]) + 3.4)
    b = ((bx, max(0.0, float(tip[0])) + 4.2),
         (-g.core_ry - 3.2, by),
         (min(-g.core_rz, float(tip[2]) - 4.2) - 1.6, g.z_cut + 1.2))
    return fn, b


def _lower_z(g, x, y):
    """Z of the core's lower (-Z) surface at (x, y), on the bare ellipsoid."""
    P = g.P
    t = 1.0 - ((x - g.core_cx) / g.core_rx) ** 2 - (y / g.core_ry) ** 2
    return -g.core_rz * math.sqrt(max(t, 1e-4))


# --------------------------------------------------------------------------
# PART: faceplate
# --------------------------------------------------------------------------

def part_faceplate(g):
    P, cx = g.P, g.core_cx

    def fn(X, Y, Z, C):
        outer = ssub(g.core_outer(X, Y, Z, C), g.corner_cut(X, Y, Z), P["corner_roll"])
        d = I(outer, g.z_cut - Z)
        # shell by eroding the CHAMFERED outer, so the roll keeps a uniform wall
        d = S(d, np.maximum(outer + 1.0, g.z_cut + 0.9 - Z))
        for (mx, my) in g.fp_mags:
            d = S(d, cyl_z(X, Y, Z, mx, my, 0.5 * P["jmag_dia"] + 0.05,
                           g.z_cut - 1.0, g.z_cut + P["jmag_depth"]))
        return d

    b = ((cx - g.core_rx - 1.0, 1.0),
         (-g.core_ry - 1.0, g.core_ry + 1.0),
         (g.z_cut - 1.0, g.core_rz + 1.0))
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

        if P["wing_style"] == "gyroid":
            # ---- legacy: MACRO gyroid sheet wing (see README 6)
            D2, S_, L_tot = _bezier_dist_and_s(C, pts)
            env_w = wing_envelope(g, X, Y, Z, D2)
            env_w = S(env_w, env - clear)
            wall_w = (P["wing_wall_root"]
                      + (P["wing_wall_tip"] - P["wing_wall_root"])
                      * np.clip(S_ / max(L_tot, 1e-6), 0.0, 1.0))
            prox = np.clip(1.0 + env_w / P["wing_edge_band"], 0.0, 1.0)
            wall_w = wall_w + (P["wing_edge_wall"] - wall_w) * prox
            sheet = gyroid(X, Y, Z, P["gyroid_cell_wing"], wall_w)
            wy = Y - g.y_root
            root_plug = np.maximum(wy - P["wing_root_solid"], -wy - 0.6)
            d = smin(jacket, I(U(sheet, root_plug), env_w), 0.35)
        else:
            # ---- three radial mag-plunger bosses
            d = jacket
            for pl in g.plungers:
                nstep = int(pl["cam_steps"])
                A, R, TH = g.pl_coords(pl, X, Y, Z)
                boss = np.maximum(R - 0.5 * P["plunger_boss_od"],
                                  slab(A, -pl["boss_h"] - 2.5, 0.0))
                d = smin(d, boss, 0.8)
                # cam + fixed-ring bore
                d = S(d, np.maximum(R - (0.5 * P["plunger_cam_od"] + 0.05),
                                    slab(A, -(pl["cam_h"]
                                              + P["plunger_mag_t"] + 0.15), 0.40)))
                # sleeved guide-pin bore
                d = S(d, np.maximum(
                    R - (0.5 * P["plunger_pin_od"] + P["plunger_pin_sleeve"]),
                    slab(A, -pl["boss_h"] - 3.0, 0.5)))
                # detent notches for the cam bumps
                for i in range(nstep):
                    a_ = -np.pi + (i + 0.5) * (2 * np.pi / nstep)
                    rn = 0.5 * P["plunger_cam_od"] + 0.05
                    dd = np.sqrt((R - rn) ** 2 + (rn * _wrap(TH - a_)) ** 2)
                    d = S(d, np.maximum(dd - 0.42,
                                        slab(A, -(pl["cam_h"] + 1.15) + 0.25,
                                             -(pl["cam_h"] + 1.15) + 1.05)))

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
        d = S(d, g.bone_boss(X, Y, Z, clear))
        return d

    ymax = max(p[1] for p in pts) + 3.0
    xmin = min(min(p[0] for p in pts), cx - g.core_rx) - 3.5
    zmin, zmax = -g.core_rz - thick - 2.0, 1.0
    if P["wing_style"] != "gyroid" and g.plungers:      # no plungers -> plain jacket bounds
        ymax = -1e9
        for pl in g.plungers:
            tip = pl["mount"] + pl["aim"] * (g.pl_pad1 + P["plunger_rocker"])
            pad = 0.5 * (P["plunger_boss_od"] + P["cymba_pad_extra"]) + 1.0
            xmin = min(xmin, float(min(tip[0], pl["base"][0])) - pad)
            ymax = max(ymax, float(max(tip[1], pl["base"][1])) + pad)
            zmin = min(zmin, float(min(tip[2], pl["base"][2])) - pad)
            zmax = max(zmax, float(max(tip[2], pl["base"][2])) + pad)
        ymin = min(-g.core_ry - 3.0,
                   min(float(pl["mount"][1] + pl["aim"][1] *
                             (g.pl_pad1 + P["plunger_rocker"])) for pl in g.plungers)
                   - 0.5 * P["plunger_boss_od"] - 1.0)
        xmax = max(P["jacket_x_clip"] + 1.0,
                   max(float(pl["mount"][0] + pl["aim"][0] *
                             (g.pl_pad1 + P["plunger_rocker"])) for pl in g.plungers)
                   + 0.5 * P["plunger_boss_od"] + 1.0)
    else:
        ymin, xmax = -g.core_ry - 3.0, P["jacket_x_clip"] + 1.0
    b = ((xmin, xmax), (ymin, ymax), (zmin, zmax))
    return fn, b


# --------------------------------------------------------------------------
# PART: mag-plungers (three radial sites)
# --------------------------------------------------------------------------

def _ring(A, R, ri, ro, a0, a1):
    """Annulus about the +A axis, in (axial, radial) coordinates."""
    return np.maximum(np.maximum(ri - R, R - ro), slab(A, a0, a1))


def _seg_seg_dist(p1, q1, p2, q2):
    """Shortest distance between two 3-D segments."""
    d1, d2 = q1 - p1, q2 - p2
    r = p1 - p2
    a, e, f = np.dot(d1, d1), np.dot(d2, d2), np.dot(d2, r)
    c = np.dot(d1, r)
    b = np.dot(d1, d2)
    den = a * e - b * b
    sN = np.clip((b * f - c * e) / den, 0, 1) if den > 1e-12 else 0.0
    tN = np.clip((b * sN + f) / e, 0, 1) if e > 1e-12 else 0.0
    sN = np.clip((b * tN - c) / a, 0, 1) if a > 1e-12 else 0.0
    return float(np.linalg.norm((p1 + d1 * sN) - (p2 + d2 * tN)))


def _seg_pt_dist(p1, q1, pt):
    d = q1 - p1
    t = np.clip(np.dot(pt - p1, d) / max(np.dot(d, d), 1e-12), 0.0, 1.0)
    return float(np.linalg.norm(p1 + d * t - pt))


def nozzle_stack_profile(g):
    """The nozzle/insert/carrier/skirt stack as a swept-sphere profile:
    (point on the canted axis, radius) samples along it."""
    P = g.P
    nb, ax = np.array(g.nozzle_base), np.array(g.n_ax)
    out = []
    if P["tip_style"] == "bell":
        # socket, then the bell tip's lip envelope (conservative: the largest
        # radius of the assembly size), then the nose cone
        rl = g.bell_lip_rmax[P["bell_asm_size"]]
        spans = [(g.socket_x0, g.socket_x1, 0.5 * P["socket_od"]),
                 (g.bell_x0, g.bell_cb_x, rl)]
        for x0, x1, r in spans:
            for t in np.linspace(x0, x1, 12):
                out.append((nb + ax * t, r))
        for t in np.linspace(g.bell_cb_x, g.bell_tip_x, 8):
            f = (t - g.bell_cb_x) / max(g.bell_tip_x - g.bell_cb_x, 1e-9)
            out.append((nb + ax * t, g.bell_base_r + f * (g.bell_nose_r - g.bell_base_r)))
        return out
    spans = [(g.socket_x0, g.flange_x1, 0.5 * P["socket_od"]),
             (g.flange_x1, g.carrier_x0, 0.5 * P["insert_od"]),
             (g.carrier_x0, g.carrier_x1, 0.5 * P["carrier_od"])]
    for x0, x1, r in spans:
        for t in np.linspace(x0, x1, 12):
            out.append((nb + ax * t, r))
    # skirt: a cone, radius growing from the carrier OD to the rim
    for t in np.linspace(g.skirt_root_x, g.skirt_rim_x, 16):
        f = (t - g.skirt_root_x) / max(g.skirt_rim_x - g.skirt_root_x, 1e-9)
        out.append((nb + ax * t, g.skirt_root_r + f * (g.skirt_rim_r - g.skirt_root_r)))
    return out


def clearance_for(g, base, aim, pad_extra=0.0, boss_h=None):
    """Clearance from a hypothetical plunger at (base, aim) to the nozzle stack."""
    P = g.P
    base = np.asarray(base, dtype=float)
    aim = np.asarray(aim, dtype=float) / np.linalg.norm(aim)
    bh = P["plunger_boss_h"] if boss_h is None else boss_h
    mount = base + aim * bh
    tip = mount + aim * (g.pl_pad1 + P["plunger_rocker"] + P["plunger_travel"])
    secs = [(base, mount, 0.5 * P["plunger_boss_od"]),
            (mount, tip, 0.5 * (P["plunger_foot_od"] + pad_extra))]
    prof = nozzle_stack_profile(g)
    return min(_seg_pt_dist(a0, a1, c) - r - ra
               for a0, a1, ra in secs for c, r in prof)


def jacket_point(g, n):
    """Point on the jacket's outer surface whose outward normal is n."""
    P = g.P
    n = np.asarray(n, dtype=float) / np.linalg.norm(n)
    rr = np.array(g.core_r)
    surf = np.array(g.core_c) + (rr ** 2 * n) / np.linalg.norm(rr * n)
    return surf + n * (P["clearance"] + P["jacket_thick"])


def leg3_feasibility(g, nominal=(0.82, -0.17, -0.54), name="tragus_inner",
                     min_clear=0.8, max_aim_dev=35.0, n_dir=4000):
    """Can leg 3 reach the tragus inner wall from ANY base on the jacket?

    The base is a free variable -- only the far end of the aim line is fixed.  So
    sweep bases over the anterior / anterior-superior jacket surface, re-aim each
    at the target the nominal leg reaches, and score clearance to the real stack.
    """
    P = g.P
    nom = np.asarray(nominal, dtype=float)
    nom /= np.linalg.norm(nom)
    pl0 = next((q for q in g.plungers_all if q["name"] == name), None)
    bh = pl0["boss_h"] if pl0 else P["plunger_boss_h"]
    need = bh + g.pl_pad1 + P["plunger_rocker"]
    T = jacket_point(g, nom) + nom * need               # the tragus-wall target
    half = 0.5 * P["plunger_cam_range"]                 # the cam absorbs this much

    i = np.arange(n_dir) + 0.5
    phi = np.arccos(1 - 2 * i / n_dir)
    th = np.pi * (1 + 5 ** 0.5) * i
    dirs = np.stack([np.cos(th) * np.sin(phi), np.sin(th) * np.sin(phi),
                     np.cos(phi)], axis=1)
    best, ok = None, []
    for d in dirs:
        if d[0] < 0.0 or d[1] < -0.7:                   # anterior / ant-superior
            continue
        B = jacket_point(g, d)
        v = T - B
        L = float(np.linalg.norm(v))
        if abs(L - need) > half:                        # cam cannot take up more
            continue
        a = v / L
        dev = math.degrees(math.acos(min(1.0, float(np.dot(a, nom)))))
        if dev > max_aim_dev:
            continue
        cl = clearance_for(g, B, a, boss_h=bh)
        rec = dict(base=B, aim=a, clear=cl, dev=dev, L=L,
                   axis_deg=math.degrees(math.acos(min(1.0, abs(float(
                       np.dot(a, g.n_ax)))))))
        if best is None or cl > best["clear"]:
            best = rec
        if cl >= min_clear:
            ok.append(rec)
    # why: where does the target itself sit relative to the nozzle stack?
    v = T - np.array(g.nozzle_base)
    st = float(np.dot(v, g.n_ax))
    perp = float(np.linalg.norm(v - st * np.array(g.n_ax)))
    prof = nozzle_stack_profile(g)
    stations = [float(np.dot(c - np.array(g.nozzle_base), g.n_ax)) for c, r in prof]
    k = int(np.argmin([abs(x - st) for x in stations]))
    near_st = stations[k]
    # several bodies share a station (carrier body and skirt rim both end at
    # carrier_x1); the widest one is what a leg has to clear
    rad = max(r for (c, r), x in zip(prof, stations) if abs(x - near_st) < 0.5)
    return dict(target=T, need=need, best=best, ok=ok, n_ok=len(ok),
                tgt_station=st, tgt_perp=perp, tgt_stack_r=rad,
                tgt_near_station=near_st)


def plunger_clearance(g, pl):
    """Clearance from a plunger's swept envelope to the nozzle stack.

    The plunger is a capsule from its boss base to the pad tip at full outward
    travel.  Negative means interference.
    """
    P = g.P
    # two sections, because the boss is narrower than the pad: only the
    # protruding envelope counts -- inside the shell the parts cannot meet
    base = np.array(pl["base"])
    mount = np.array(pl["mount"])
    tip = mount + pl["aim"] * (g.pl_pad1 + P["plunger_rocker"] + P["plunger_travel"])
    secs = [(base, mount, 0.5 * P["plunger_boss_od"]),
            (mount, tip, 0.5 * (P["plunger_foot_od"] + pl["pad_extra"]))]
    prof = nozzle_stack_profile(g)
    return min(_seg_pt_dist(a0, a1, c) - r - ra
               for a0, a1, ra in secs for c, r in prof)


def part_plunger_foot(g):
    """Moving piston: magnet pocket, pin bore, and the inward travel stop."""
    P = g.P
    ro = 0.5 * P["plunger_foot_od"]

    def fn(X, Y, Z, C):
        A, R = X, np.sqrt(Y ** 2 + Z ** 2)
        body = np.maximum(R - ro, slab(A, g.pl_mag_mov, g.pl_plate1))
        skirt = _ring(A, R, ro - 0.8, ro, g.pl_stop_in, g.pl_mag_mov)
        d = U(body, skirt)
        d = S(d, np.maximum(R - (0.5 * P["plunger_mag_od"] + 0.05),
                            slab(A, g.pl_mag_mov - 0.05, g.pl_mag_mov1 + 0.05)))
        d = S(d, np.maximum(R - (0.5 * P["plunger_pin_od"] + 0.025),
                            slab(A, g.pl_stop_in - 0.5, g.pl_plate1 + 0.5)))
        return d

    return fn, ((g.pl_stop_in - 0.6, g.pl_plate1 + 0.6), (-ro - 0.6, ro + 0.6),
                (-ro - 0.6, ro + 0.6))


def part_plunger_pad(g, extra=0.0):
    """Soft contact pad: rocker crown plus a rolled shoulder so the edge can tuck
    under an overhanging lip instead of digging into it."""
    P = g.P
    ro = 0.5 * (P["plunger_foot_od"] + extra)
    rl = P["plunger_pad_roll"]

    def fn(X, Y, Z, C):
        A, R = X, np.sqrt(Y ** 2 + Z ** 2)
        crown = P["plunger_rocker"] * np.clip(1.0 - (R / ro) ** 2, 0.0, 1.0)
        # rounded-box profile in (axial, radial): gives a >= rl rolled edge all
        # the way round the shoulder
        qa = np.abs(A - 0.5 * (g.pl_plate1 + g.pl_pad1 + crown)) \
            - (0.5 * (g.pl_pad1 + crown - g.pl_plate1) - rl)
        qr = R - (ro - rl)
        d = (np.sqrt(np.maximum(qa, 0) ** 2 + np.maximum(qr, 0) ** 2)
             + np.minimum(np.maximum(qa, qr), 0.0) - rl)
        d = S(d, np.maximum(R - (0.5 * P["plunger_pin_od"] + 0.10),
                            slab(A, g.pl_plate1 - 0.5, g.pl_pad1 + 1.5)))
        return d

    return fn, ((g.pl_plate1 - 0.8, g.pl_pad1 + P["plunger_rocker"] + 0.8),
                (-ro - 0.8, ro + 0.8), (-ro - 0.8, ro + 0.8))


def part_plunger_pin(g):
    """Guide pin: pressed into the foot, sliding in the jacket's sleeved bore."""
    P = g.P
    x0 = -(P["plunger_boss_h"] + 1.5)

    def fn(X, Y, Z, C):
        A, R = X, np.sqrt(Y ** 2 + Z ** 2)
        shaft = np.maximum(R - 0.5 * P["plunger_pin_od"], slab(A, x0, g.pl_plate1))
        head = np.maximum(R - 1.7, slab(A, x0, x0 + 0.6))
        return U(shaft, head)

    return fn, ((x0 - 0.5, g.pl_plate1 + 0.5), (-2.2, 2.2), (-2.2, 2.2))


def part_plunger_cam(g, steps=None, rng=None):
    """Cam preset ring: a staircase top face gives coarse engagement."""
    P = g.P
    n = int(steps if steps else P["plunger_cam_steps"])
    rng = rng if rng else P["plunger_cam_range"]
    step = rng / max(n - 1, 1)
    ro = 0.5 * P["plunger_cam_od"]
    hmax = 0.2 + rng

    def fn(X, Y, Z, C):
        A, R = X, np.sqrt(Y ** 2 + Z ** 2)
        th = np.arctan2(Z, Y)
        k = np.floor((th + np.pi) / (2 * np.pi / n))
        top = 0.2 + step * np.clip(k, 0, n - 1)
        d = _ring(A, R, 0.5 * P["plunger_mag_id"] + 0.10, ro, 0.0, hmax)
        d = np.maximum(d, A - top)                        # staircase top face
        # seat that locates the fixed ring on whichever step is selected
        d = S(d, np.maximum(R - (0.5 * P["plunger_mag_od"] + 0.05),
                            A - top + 0.55))
        # detent bumps on the OD
        for i in range(n):
            a_ = -np.pi + (i + 0.5) * (2 * np.pi / n)
            cy, cz = ro * math.cos(a_), ro * math.sin(a_)
            d = U(d, np.sqrt((X - 0.6) ** 2 + (Y - cy) ** 2 + (Z - cz) ** 2) - 0.40)
        return d

    return fn, ((-0.5, hmax + 0.5), (-ro - 0.8, ro + 0.8), (-ro - 0.8, ro + 0.8))


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


def _insert_socket_cuts(g, X, Y, Z, d):
    """Socket bore over the core stub + the two L-slots for the core's bayonet
    lugs.  Shared by every nozzle insert."""
    P = g.P
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
    return d


def part_nozzle_insert_bell(g):
    """Fixed nozzle for the bell tip: the same bayonet socket onto the core stub,
    then a plain O4 Ti tube with an O2.6 bore and a retaining groove for the
    tip's internal ridge.  No magnet flange -- the float is gone."""
    P = g.P
    x_end = g.bell_noz_end
    ro = g.bell_seat_r
    bore_r = 0.5 * P["bell_bore"]

    def fn(X, Y, Z, C):
        socket = cyl_x(X, Y, Z, 0, 0, 0.5 * P["socket_od"], g.socket_x0, g.socket_x1)
        tube = cyl_x(X, Y, Z, 0, 0, ro, g.socket_x1 - 0.4, x_end)
        d = U(socket, tube)
        # retaining groove for the tip's ridge
        d = S(d, tube_x(X, Y, Z, 0, 0, ro - P["bell_groove_d"], ro + 1.0,
                        g.bell_groove_x0, g.bell_groove_x1))
        d = S(d, cyl_x(X, Y, Z, 0, 0, bore_r, g.socket_x1 - 1.0, x_end + 1.0))
        d = _insert_socket_cuts(g, X, Y, Z, d)
        # damper-disc recess at the ear end (disc = bore size, not the O4 jig)
        d = S(d, cyl_x(X, Y, Z, 0, 0, bore_r + 0.15,
                       x_end - P["damper_recess"], x_end + 1.0))
        return d

    b = ((g.socket_x0 - 1.2, x_end + 1.2), (-5.5, 5.5), (-5.5, 5.5))
    return fn, b


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
        d = _insert_socket_cuts(g, X, Y, Z, d)

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

    # ---- intertragic-notch sector: compliance only, no radial reach ---------
    if P["notch_compliance"]:
        hw = math.radians(0.5 * P["notch_sector_deg"])
        tr = math.radians(P["notch_sector_trans_deg"])
        th = np.arctan2(Z, Y)
        dth = np.abs(_wrap(th - math.radians(P["notch_sector_center_deg"])))
        tb = np.clip((dth - (hw - tr)) / tr, 0.0, 1.0)
        b = 1.0 - tb * tb * (3.0 - 2.0 * tb)
        wn = np.interp(sl, g.notch_wall_xp, g.notch_wall_fp)
        w = w + (wn - w) * b
        if P["notch_sector_ext"]:                     # kept, but 0 by default
            dn = dn - P["notch_sector_ext"] * b * np.clip(sl / g.skirt_slant, 0, 1) * ux

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
    r = g.skirt_rim_r + P["notch_sector_ext"] + 1.2

    def fn(X, Y, Z, C):
        return carrier_field(g, X, Y, Z, C)

    b = ((g.carrier_x0 - 1.2, g.carrier_x1 + 1.2), (-r, r), (-r, r))
    return fn, b


def mold_geom(g):
    P = g.P
    x0 = g.carrier_x0 - 3.5
    x1 = g.carrier_x1 + 3.5
    r = g.skirt_rim_r + P["notch_sector_ext"] + P["mold_wall"]
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
# PART: bell tip S/M/L (+ moulds) -- replaces the mag-float carrier + skirt
# --------------------------------------------------------------------------
# Decided 2026-09-01 (artifact "Aperture Tip"): one blunt nose cone for every
# size, cast Shore A 10-15, O5 -> O12 over 2.4 mm (55 deg half-angle) so the
# cone's base stops at the aperture rim and insertion is self-limited; a hollow
# rolled lip (tube O3, wall 0.4, Shore 00-30) on an oval footprint in S/M/L; a
# 0.5 mm bell-shaped web between them so the lip can rock; an anterior sector
# thinned + lengthened for jaw motion; an inferior sector extended over the
# intertragic notch; a 0.8 mm vent exiting under the lip.
#
# Frame: nozzle-local, +X along the nozzle into the ear -- the same frame as the
# nozzle inserts and the old carrier.  th = atan2(Z, Y): 0 = +Y = superior,
# 90 deg = +Z = anterior (tragus side), 180 deg = -Y = inferior (notch).
#
# Mould: two halves split on y = 0 pulling +/-Y, exactly like the carrier
# mould, plus a removable core that carries the bore, the seat with its
# retaining ridge, the vent, the plug that fills the hollow under the web, and
# the ring inside the rolled lip.  A closed hollow torus cannot be pulled, so
# the lip is a C-section open toward the axis (bell_lip_slit_deg) and the core's
# ring comes out through the slit -- Shore 00-30 stretches ~1.5x there.

BELL_SIZES = ("XS", "S", "M", "L")
BELL_LEGACY_PARTS = {"carrier", "carrier_mold_a", "carrier_mold_b", "carrier_mold_core",
                     "nozzle_insert_short", "nozzle_insert_med", "nozzle_insert_long"}


def _sector(th, center_deg, width_deg, trans_deg):
    """1 inside a sector of azimuth, 0 outside, smoothstep over trans_deg."""
    hw = math.radians(0.5 * width_deg)
    tr = math.radians(trans_deg)
    dth = np.abs(_wrap(th - math.radians(center_deg)))
    tb = np.clip((dth - (hw - tr)) / tr, 0.0, 1.0)
    return 1.0 - tb * tb * (3.0 - 2.0 * tb)


def _smoothstep(t):
    t = np.clip(t, 0.0, 1.0)
    return t * t * (3.0 - 2.0 * t)


def _polyline_dist(px, pr, pts):
    """Distance from (px, pr) to a polyline of (x, r) points; points may be arrays."""
    best = None
    for (ax, ar), (bx, br) in zip(pts[:-1], pts[1:]):
        dx, dr = bx - ax, br - ar
        dd = np.maximum(dx * dx + dr * dr, 1e-12)
        h = np.clip(((px - ax) * dx + (pr - ar) * dr) / dd, 0.0, 1.0)
        d = np.sqrt((px - ax - dx * h) ** 2 + (pr - ar - dr * h) ** 2)
        best = d if best is None else np.minimum(best, d)
    return best


def bell_lip(g, size, th):
    """Per-azimuth lip: (centreline radius, tube centre x, wall)."""
    P = g.P
    a_c, b_c = g.bell_lip_axes[size]
    r_e = 1.0 / np.sqrt((np.cos(th) / a_c) ** 2 + (np.sin(th) / b_c) ** 2)
    b_inf = _sector(th, P["bell_inf_center_deg"], P["bell_inf_deg"],
                    P["bell_sector_trans_deg"])
    b_ant = _sector(th, P["bell_ant_center_deg"], P["bell_ant_deg"],
                    P["bell_sector_trans_deg"])
    r_lip = r_e + P["bell_inf_ext"] * b_inf
    x_lip = g.bell_lip_x - P["bell_ant_free"] * b_ant
    wall = P["bell_lip_wall"] + (P["bell_ant_wall"] - P["bell_lip_wall"]) * b_ant
    return r_lip, x_lip, wall


def bell_web_r(g, X, r_lip, x_wl):
    """Web centreline radius at station X: an S-curve from the lip radius to the
    land, axial tangent at both ends -- the bell."""
    t = (X - x_wl) / np.maximum(g.bell_web_x1 - x_wl, 1e-6)
    return r_lip + (g.bell_web_r1 - r_lip) * _smoothstep(t)


def bell_web_pts(g, r_lip, x_lip, n=10, Rt=None):
    P = g.P
    Rt = g.bell_Rt if Rt is None else Rt
    x_wl = x_lip + Rt + P["bell_web_stub"]
    pts = [(x_lip + Rt - 0.05, r_lip), (x_wl, r_lip)]     # stub ends in the shell wall
    for t in np.linspace(0.0, 1.0, n + 1)[1:]:
        pts.append((x_wl + (g.bell_web_x1 - x_wl) * t,
                    r_lip + (g.bell_web_r1 - r_lip) * _smoothstep(t)))
    return pts


def bell_vent_field(g, X, Y, Z):
    """The O0.8 vent, as a solid to cut from the tip / add to the mould core.

    canal <- bore -> radial link just past the Ti tube end -> groove along the
    Ti/silicone seat (closed by the tube) -> radial channel across the tip's
    proximal face (closed by the insert's socket face) -> the hollow under the
    web, which is open to the concha under the lip.  Nothing on the core bridges
    the sleeve, so the core pulls -X without tearing it.
    """
    P = g.P
    rv = 0.5 * P["bell_vent_dia"]
    az = math.radians(P["bell_vent_az_deg"])
    uy, uz = math.cos(az), math.sin(az)
    rb = g.bell_seat_r + 0.15                    # a near-full O0.8 groove in the seat
    bead = cyl_x(X, Y, Z, rb * uy, rb * uz, rv, g.bell_x0 - 0.6, g.bell_noz_end + 0.7)
    xl = g.bell_noz_end + 0.3
    link = capsule(X, Y, Z, (xl, 0.0, 0.0), (xl, rb * uy, rb * uz), rv)
    xf = g.bell_x0 + 0.15
    r1 = g.bell_rib_r + 0.3
    face = capsule(X, Y, Z, (xf, 1.5 * uy, 1.5 * uz), (xf, r1 * uy, r1 * uz), rv)
    return U(bead, link, face)


def bell_tip_field(g, size, X, Y, Z, C=None):
    P = g.P
    key = ("bell_tip", size, X.shape, Y.shape, Z.shape,
           float(X.ravel()[0]), float(Y.ravel()[0]), float(Z.ravel()[0]))
    if C is not None and key in C.cache:
        return C.cache[key]
    rho = np.sqrt(Y ** 2 + Z ** 2)
    th = np.arctan2(Z, Y)
    Rt = g.bell_Rt_by[size]
    r_lip, x_lip, wall = bell_lip(g, size, th)

    # ---- nose: O12 land + 55 deg cone + sleeve over the Ti tube + vent rib
    land = cyl_x(X, Y, Z, 0, 0, g.bell_base_r, g.bell_land_x0, g.bell_cb_x)
    cone = cone_x(X, Y, Z, g.bell_cb_x, g.bell_base_r, g.bell_tip_x, g.bell_nose_r)
    sleeve = cyl_x(X, Y, Z, 0, 0, g.bell_sleeve_r, g.bell_x0, g.bell_land_x0 + 0.2)
    az = math.radians(P["bell_vent_az_deg"])
    u = Y * math.cos(az) + Z * math.sin(az)      # radial along the vent azimuth
    v = -Y * math.sin(az) + Z * math.cos(az)
    # a box, not a wedge: its sides are parallel to the mould pull
    rib = box(X, u, v, (0.5 * (g.bell_x0 + g.bell_land_x0),
                        0.5 * (1.5 + g.bell_rib_r), 0.0),
              (0.5 * (g.bell_land_x0 - g.bell_x0) + 0.1,
               0.5 * (g.bell_rib_r - 1.5), 0.9))
    nose = U(land, cone, sleeve, rib)

    # ---- web sheet
    web = _polyline_dist(X, rho, bell_web_pts(g, r_lip, x_lip, Rt=Rt)) - 0.5 * P["bell_web_t"]
    d = smin(nose, web, P["bell_fillet"])

    # ---- rolled lip: tube O3 in every meridian plane, on the oval centreline
    uu, vv = X - x_lip, rho - r_lip
    rr = np.sqrt(uu ** 2 + vv ** 2)
    dt = rr - Rt
    if P["bell_lip_hollow"]:
        shell = np.maximum(dt, -(dt + wall))
        phi = np.arctan2(uu, -vv)                # 0 at the pole facing the axis
        slit = (np.abs(phi) - math.radians(0.5 * P["bell_lip_slit_deg"])) \
            * np.maximum(rr, 0.2)
        lip = S(shell, slit)
    else:
        lip = dt
    d = U(d, lip)

    # ---- cuts: seat on the Ti tube (keeping the ridge), bore, vent
    seat = cyl_x(X, Y, Z, 0, 0, g.bell_seat_r, g.bell_x0 - 1.0, g.bell_noz_end)
    ridge = tube_x(X, Y, Z, 0, 0, g.bell_ridge_r, g.bell_seat_r + 0.2,
                   g.bell_groove_x0, g.bell_groove_x1)
    # the vent groove crosses the ridge: open the ridge there rather than leave
    # a sliver of silicone between the groove and the ridge's inner face
    ridge = S(ridge, arc_slot_x(X, Y, Z, 0.0, g.bell_seat_r + 0.3,
                                g.bell_groove_x0 - 0.2, g.bell_groove_x1 + 0.2,
                                az - 0.32, az + 0.32))
    seat = S(seat, ridge)
    bore = cyl_x(X, Y, Z, 0, 0, 0.5 * P["bell_bore"], g.bell_x0 - 1.0, g.bell_tip_x + 1.0)
    d = S(d, U(seat, bore))
    d = S(d, bell_vent_field(g, X, Y, Z))
    if C is not None:
        C.cache[key] = d
    return d


def bell_mold_geom(g, size):
    P = g.P
    x0 = g.bell_x0 - P["bell_disc_t"] - 3.5
    x1 = g.bell_tip_x + 3.5
    r = g.bell_lip_rmax[size] + P["mold_wall"]
    return x0, x1, r


def bell_core_field(g, size, X, Y, Z, C=None):
    """The removable core = everything the mould halves must not fill that the
    tip does not occupy: rod (seat + O2.6 bore), the ridge groove, the vent bead,
    the plug under the web, the ring inside the lip, and a base disc behind the
    tip's proximal face.  Built as (region - tip) so the cast is exact."""
    P = g.P
    rho = np.sqrt(Y ** 2 + Z ** 2)
    th = np.arctan2(Z, Y)
    Rt = g.bell_Rt_by[size]
    r_lip, x_lip, _ = bell_lip(g, size, th)
    x_wl = x_lip + Rt + P["bell_web_stub"]
    rp = np.where(X < g.bell_x0, g.bell_disc_r[size],
         np.where(X < x_wl, r_lip,
         np.where(X < g.bell_land_x0, bell_web_r(g, X, r_lip, x_wl),
         np.where(X < g.bell_noz_end, g.bell_seat_r, 0.5 * P["bell_bore"]))))
    x0, x1, r = bell_mold_geom(g, size)
    region = np.maximum(rho - rp, slab(X, g.bell_x0 - P["bell_disc_t"], x1 + 2.0))
    region = U(region, np.sqrt((X - x_lip) ** 2 + (rho - r_lip) ** 2) - Rt)
    core = S(region, bell_tip_field(g, size, X, Y, Z, C))
    rod = cyl_x(X, Y, Z, 0, 0, g.bell_seat_r, x0 - 2.0,
                g.bell_x0 - P["bell_disc_t"] + 0.5)
    handle = cyl_x(X, Y, Z, 0, 0, 4.0, x0 - 8.0, x0 - 1.5)
    return U(core, rod, handle, bell_vent_field(g, X, Y, Z))


def bell_void_field(g, size, X, Y, Z, C=None):
    """tip + core: the cavity the two halves close around."""
    return U(bell_tip_field(g, size, X, Y, Z, C), bell_core_field(g, size, X, Y, Z, C))


def part_bell_tip(g, size):
    rr = g.bell_lip_rmax[size] + 1.2

    def fn(X, Y, Z, C):
        return bell_tip_field(g, size, X, Y, Z, C)

    return fn, ((g.bell_x0 - 1.0, g.bell_tip_x + 1.0), (-rr, rr), (-rr, rr))


def part_bell_mold_core(g, size):
    x0, x1, r = bell_mold_geom(g, size)
    rd = g.bell_disc_r[size] + 1.0

    def fn(X, Y, Z, C):
        return bell_core_field(g, size, X, Y, Z, C)

    return fn, ((x0 - 9.0, x1 + 3.0), (-rd, rd), (-rd, rd))


def _bell_mold_half(g, size, upper):
    P = g.P
    x0, x1, r = bell_mold_geom(g, size)
    sgn = 1.0 if upper else -1.0
    pinx = (x1 - 2.5, x0 + 2.5)
    pinz = g.bell_lip_rmax[size] + 0.5 * P["mold_wall"]
    a_c, b_c = g.bell_lip_axes[size]
    Rt = g.bell_Rt_by[size]

    def fn(X, Y, Z, C):
        block = box(X, Y, Z, (0.5 * (x0 + x1), 0.0, 0.0),
                    (0.5 * (x1 - x0), r, r))
        d = I(block, -sgn * Y)
        void = bell_void_field(g, size, X, Y, Z, C)
        # pour spout: from the +X face down the cone's flank (pour nose-up, so
        # the lip fills first -- that is where a metered 00-30 first pour goes)
        void = U(void, cyl_x(X, Y, Z, 0.0, 3.6, 0.5 * P["mold_spout_dia"],
                             g.bell_cb_x + 0.3, x1 + 2.0))
        # vents at the lip's proximal-outer quadrant, both sides, in the parting plane
        xv = g.bell_lip_x - 1.0
        zv = b_c + Rt - 0.6
        void = U(void, cyl_z(X, Y, Z, xv, 0.0, 0.5 * P["mold_vent_dia"], zv, r + 2.0))
        void = U(void, cyl_z(X, Y, Z, xv, 0.0, 0.5 * P["mold_vent_dia"], -r - 2.0, -zv))
        d = S(d, void)

        for px in pinx:
            for pz in (pinz, -pinz):
                pin = cyl_y(X, Y, Z, px, pz, 0.5 * P["mold_pin_dia"],
                            -P["mold_pin_len"], P["mold_pin_len"])
                if upper:
                    d = U(d, I(pin, Y))
                    d = S(d, I(pin, -Y))
                else:
                    hole = cyl_y(X, Y, Z, px, pz, 0.5 * P["mold_pin_dia"] + 0.12,
                                 -P["mold_pin_len"] - 0.2, 0.2)
                    d = S(d, hole)
        return d

    b = ((x0 - 1.5, x1 + 1.5), (-r - 4.0, r + 4.0), (-r - 1.5, r + 1.5))
    return fn, b


def _first_run(inside, step):
    """Length of the first contiguous True run."""
    idx = np.flatnonzero(inside)
    if len(idx) == 0:
        return 0.0
    n = 1
    while n < len(idx) and idx[n] == idx[n - 1] + 1:
        n += 1
    return n * step


def _field_at(fn, pts):
    """Evaluate a broadcasting field function on an (n, 3) list of points."""
    p = np.asarray(pts, dtype=float)
    X = p[:, 0].reshape(-1, 1, 1)
    Y = p[:, 1].reshape(-1, 1, 1)
    Z = p[:, 2].reshape(-1, 1, 1)
    return np.asarray(fn(X, Y, Z)).reshape(-1)


def _wall_along(fn, origins, dirs, tmax, step):
    """March from each origin along its direction and return the length of the
    first negative (inside) run of the field -- the wall thickness there."""
    o = np.asarray(origins, dtype=float)
    d = np.asarray(dirs, dtype=float)
    d = d / np.maximum(np.linalg.norm(d, axis=1, keepdims=True), 1e-12)
    t = np.arange(0.0, tmax + step, step)
    pts = (o[:, None, :] + d[:, None, :] * t[None, :, None]).reshape(-1, 3)
    f = _field_at(fn, pts).reshape(len(o), len(t)) < 0
    out = np.zeros(len(o))
    for i in range(len(o)):
        out[i] = _first_run(f[i], step)
    return out


def bell_measure(g, size, mesh, n=800):
    """Numbers off the built tip: footprint and volume from the mesh, insertion
    depth from the rim stop, the named walls probed through the field, the
    global minimum wall (field marched inward from mesh vertices), and the
    rocking clearance."""
    P = g.P
    Rt = g.bell_Rt_by[size]
    fn = lambda X, Y, Z: bell_tip_field(g, size, X, Y, Z)
    v = np.asarray(mesh.vertices)
    rho = np.hypot(v[:, 1], v[:, 2])
    x_max = float(v[:, 0].max())
    stop = v[rho >= g.bell_base_r - 0.06, 0]
    x_stop = float(stop.max()) if len(stop) else float("nan")
    out = dict(H=float(v[:, 1].max() - v[:, 1].min()),
               W=float(v[:, 2].max() - v[:, 2].min()),
               vol=float(mesh.volume), x_max=x_max, x_stop=x_stop,
               insertion=x_max - x_stop, protrusion=x_max)

    def probe(p0, p1):
        p0, p1 = np.asarray(p0, float), np.asarray(p1, float)
        L = float(np.linalg.norm(p1 - p0))
        return float(_wall_along(fn, [p0], [p1 - p0], L, L / (n - 1))[0])

    def lip_at(th_deg):
        th = np.array([math.radians(th_deg)])
        r_lip, x_lip, wall = bell_lip(g, size, th)
        return float(r_lip[0]), float(x_lip[0]), float(wall[0])

    # lip wall, probed inward from outside the tube's outer pole
    for tag, th_deg in (("lip", 0.0), ("lip_ant", P["bell_ant_center_deg"])):
        rl, xl, wl = lip_at(th_deg)
        c, s = math.cos(math.radians(th_deg)), math.sin(math.radians(th_deg))
        out[tag] = probe((xl, (rl + Rt + 0.3) * c, (rl + Rt + 0.3) * s),
                         (xl, (rl + 0.2) * c, (rl + 0.2) * s))
        out[tag + "_spec"] = wl
    # web, probed along its normal at mid-run (superior azimuth)
    rl, xl, _ = lip_at(0.0)
    x_wl = xl + Rt + P["bell_web_stub"]
    xm = 0.5 * (x_wl + g.bell_web_x1)
    h = 0.05
    r0 = float(bell_web_r(g, np.array([xm - h]), rl, x_wl)[0])
    r1 = float(bell_web_r(g, np.array([xm + h]), rl, x_wl)[0])
    rm = float(bell_web_r(g, np.array([xm]), rl, x_wl)[0])
    ln = math.hypot(2 * h, r1 - r0)
    nx, nr = -(r1 - r0) / ln, 2 * h / ln
    out["web"] = probe((xm + 0.8 * nx, rm + 0.8 * nr, 0.0),
                       (xm - 0.8 * nx, rm - 0.8 * nr, 0.0))
    # sleeve over the Ti tube (superior azimuth, clear of the vent rib)
    xs = 0.5 * (g.bell_x0 + g.bell_land_x0)
    out["sleeve"] = probe((xs, 3.6, 0.0), (xs, 1.6, 0.0))
    # nose: bore wall near the tip face
    xt = g.bell_tip_x - 0.4
    out["nose"] = probe((xt, 3.5, 0.0), (xt, 0.8, 0.0))

    # global minimum wall: march inward from a sample of mesh vertices
    rng = np.random.default_rng(0)
    sel = rng.choice(len(v), size=min(3000, len(v)), replace=False)
    nrm = np.asarray(mesh.vertex_normals)[sel]
    t = _wall_along(fn, v[sel] - nrm * 0.02, -nrm, 2.0, 0.01)
    keep = t > 0.03
    t, pv = t[keep], v[sel][keep]
    out["min_wall"] = float(t.min()) if len(t) else float("nan")
    out["p05_wall"] = float(np.percentile(t, 5)) if len(t) else float("nan")
    if len(t):
        p = pv[int(np.argmin(t))]
        pr, pth = math.hypot(p[1], p[2]), math.degrees(math.atan2(p[2], p[1]))
        dv = abs(_wrap(np.array([math.radians(pth - P["bell_vent_az_deg"])]))[0])
        if pr < g.bell_rib_r + 0.2 and dv < math.radians(30):
            out["min_where"] = "vent"          # cusp where the O0.8 vent meets the bore/seat
        elif abs(pr - _field_lip_r(g, size, pth)) < Rt + 0.2 and \
                abs(p[0] - _field_lip_x(g, size, pth)) < Rt + 0.2:
            out["min_where"] = "lip"           # free edge of the C-section
        else:
            out["min_where"] = "body"
        out["min_at"] = (float(p[0]), float(pr), float(pth))

    # rocking: radial gap from the lip's inner edge to the sleeve, over the azimuth
    th = np.linspace(-np.pi, np.pi, 721)
    r_lip, x_lip, _ = bell_lip(g, size, th)
    gap = r_lip - Rt - g.bell_sleeve_r
    x_wl = x_lip + Rt + P["bell_web_stub"]
    L = np.hypot(g.bell_web_x1 - x_wl, g.bell_web_r1 - r_lip) + P["bell_web_stub"]
    i = int(np.argmin(gap))
    out["rock_gap"] = float(gap[i])
    out["rock_deg"] = float(math.degrees(math.atan2(max(gap[i], 0.0), L[i])))
    out["rock_az"] = float(math.degrees(th[i]))
    out["web_len_min"], out["web_len_max"] = float(L.min()), float(L.max())
    return out


def _field_lip_r(g, size, th_deg):
    return float(bell_lip(g, size, np.array([math.radians(th_deg)]))[0][0])


def _field_lip_x(g, size, th_deg):
    return float(bell_lip(g, size, np.array([math.radians(th_deg)]))[1][0])


def bell_draw_check(g, size, spacing=0.16, n_pts=4000, step=0.08):
    """Can the two halves pull +/-Y off (tip + core)?  From points on the cavity
    surface, march along the pull through the void field: mould material that
    meets cavity again further along is trapped under a silicone overhang.
    Returns (fraction of shadowed points, worst overhang thickness in mm)."""
    rr = g.bell_lip_rmax[size] + 1.0
    b = ((g.bell_x0 - 0.5, g.bell_tip_x + 0.5), (-rr, rr), (-rr, rr))
    fn = lambda X, Y, Z: bell_void_field(g, size, X, Y, Z)
    fld, org, sp = evaluate(lambda X, Y, Z, C: fn(X, Y, Z), b, spacing)
    m = polygonise(fld, org, sp)
    v = np.asarray(m.vertices)
    nrm = np.asarray(m.vertex_normals)
    rng = np.random.default_rng(1)
    shadow, worst, total = 0, 0.0, 0
    for sgn in (1.0, -1.0):
        # cavity-surface points on this half, off the sampling box's cut faces
        sel = np.flatnonzero((sgn * v[:, 1] > 0.3)
                             & (v[:, 0] > b[0][0] + 0.3) & (v[:, 0] < b[0][1] - 0.3))
        if len(sel) == 0:
            continue
        sel = rng.choice(sel, size=min(n_pts, len(sel)), replace=False)
        o = v[sel] + nrm[sel] * 0.04
        ok = _field_at(fn, o) > 0.0            # origin must sit in mould material
        sel, o = sel[ok], o[ok]
        t = np.arange(0.0, 2 * rr, step)
        pts = (o[:, None, :] + np.array([0.0, sgn, 0.0])[None, None, :]
               * t[None, :, None]).reshape(-1, 3)
        f = _field_at(fn, pts).reshape(len(sel), len(t))
        inside = f < -0.02
        total += len(sel)
        for i in range(len(sel)):
            if inside[i].any():
                shadow += 1
                worst = max(worst, _first_run(inside[i], step))
    return shadow / max(total, 1), worst


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
    b = ((xmin, g.core_cx + g.core_rx + 2.0),
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
    "plunger_foot": (part_plunger_foot, "solid"),
    "plunger_pad": (part_plunger_pad, "solid"),
    "plunger_pad_cymba": (lambda g: part_plunger_pad(g, g.P["cymba_pad_extra"]),
                          "solid"),
    "plunger_cam_ext": (lambda g: part_plunger_cam(g, g.P["plunger_cam_ext_steps"],
                                                   g.P["plunger_cam_ext_range"]),
                        "solid"),
    "plunger_pin": (part_plunger_pin, "solid"),
    "plunger_cam": (part_plunger_cam, "solid"),
    "driver_carrier": (part_driver_carrier, "solid"),
    "damper_jig": (part_damper_jig, "solid"),
    "nozzle_insert_bell": (part_nozzle_insert_bell, "solid"),
}
for _sz in BELL_SIZES:
    PARTS[f"bell_tip_{_sz}"] = ((lambda g, s=_sz: part_bell_tip(g, s)), "fine")
    PARTS[f"bell_tip_{_sz}_mold_a"] = ((lambda g, s=_sz: _bell_mold_half(g, s, True)), "fine")
    PARTS[f"bell_tip_{_sz}_mold_b"] = ((lambda g, s=_sz: _bell_mold_half(g, s, False)), "fine")
    PARTS[f"bell_tip_{_sz}_mold_core"] = ((lambda g, s=_sz: part_bell_mold_core(g, s)), "fine")

BELL_PARTS = {n for n in PARTS if n.startswith("bell_tip_")} | {"nozzle_insert_bell"}


def enabled_parts(P):
    """What --all builds: the bell family by default, the mag-float carrier
    family only with tip_style='carrier'."""
    off = BELL_PARTS if P["tip_style"] == "carrier" else BELL_LEGACY_PARTS
    return [n for n in PARTS if n not in off]


def assembly_parts(P):
    """Parts that live in the assembly frame (moulds and the jig have their own)."""
    if P["tip_style"] == "carrier":
        tip = ["nozzle_insert_short", "carrier"]
    else:
        tip = ["nozzle_insert_bell", f"bell_tip_{P['bell_asm_size']}"]
    return ["core", "faceplate", "jacket_wing"] + tip + ["driver_carrier"]


PLUNGER_PARTS = ["plunger_foot", "plunger_pad", "plunger_pin", "plunger_cam"]

# these are modelled about the nozzle axis, so their STLs are axis-aligned in the
# nozzle-local frame and get canted into place only for the assembly
NOZZLE_FRAME_PARTS = {"nozzle_insert_short", "nozzle_insert_med",
                      "nozzle_insert_long", "carrier", "nozzle_insert_bell"} \
    | {f"bell_tip_{s}" for s in BELL_SIZES}


def build(name, g, voxel=None):
    fn, bounds = PARTS[name][0](g)
    budget = {"lattice": g.P["budget_lattice"],
              "fine": g.P["budget_fine"]}.get(PARTS[name][1], g.P["budget"])
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


def acoustic_void(g, spacing=0.20):
    """Volume of the enclosed void inside the core shell, below the parting plane."""
    fn, b = PARTS["core"][0](g)
    f, o, sp = evaluate(fn, b, spacing)
    X = (o[0] + sp * np.arange(f.shape[0])).reshape(-1, 1, 1)
    Y = (o[1] + sp * np.arange(f.shape[1])).reshape(1, -1, 1)
    Z = (o[2] + sp * np.arange(f.shape[2])).reshape(1, 1, -1)
    ins = (g.core_outer(X, Y, Z, None) < 0) & (Z < g.z_cut)
    return float((ins & (f > 0)).sum()) * sp ** 3


def measure_skirt_wall(g, P, theta_deg, n=600):
    """Wall thickness of the sealing land, measured by probing the built field.

    Marches inward along the cone normal at mid-land and returns the length of
    the interval where carrier_field() is negative.
    """
    ux, ur = g.skirt_u
    sl = g.skirt_slant - 0.5 * g.skirt_land_w
    x0 = g.skirt_root_x + ux * sl
    r0 = g.skirt_root_r + ur * sl
    t = np.linspace(-0.15, 1.20, n)
    th = math.radians(theta_deg)
    X = (x0 + ur * t).reshape(-1, 1, 1)
    rr = r0 - ux * t
    Y = (rr * math.cos(th)).reshape(1, -1, 1)
    Z = (rr * math.sin(th)).reshape(1, 1, -1)
    # sample along the ray only (diagonal of the broadcast grid)
    d = np.array([float(carrier_field(g, X[i], Y[0, i], Z[0, 0, i], None))
                  for i in range(n)])
    inside = d < 0
    return float(inside.sum()) * (t[1] - t[0]) if inside.any() else 0.0


def measure_notch_reach(mesh, g, P):
    """Realised radial reach of the notch sector, measured on the built mesh.

    Takes the rim ring of the carrier and compares the rim radius inside the
    sector against the rim radius outside it.  This is the number that matters:
    the cone slant and marching-cubes rounding both eat into the nominal.
    """
    v = np.asarray(mesh.vertices)
    sel = v[:, 0] > g.skirt_rim_x - 0.60
    if sel.sum() < 200:
        sel = np.ones(len(v), dtype=bool)
    vv = v[sel]
    rho = np.hypot(vv[:, 1], vv[:, 2])
    th = np.arctan2(vv[:, 2], vv[:, 1])
    hw = math.radians(0.5 * P["notch_sector_deg"])
    tr = math.radians(P["notch_sector_trans_deg"])
    dth = np.abs(_wrap(th - math.radians(P["notch_sector_center_deg"])))
    inn, out = dth < (hw - tr), dth > (hw + tr)
    if inn.sum() < 20 or out.sum() < 20:
        return float("nan")
    return float(np.percentile(rho[inn], 95) - np.percentile(rho[out], 95))


def calibrate_notch(P, target=3.25, lo=3.0, hi=3.5, iters=5):
    """Drive notch_sector_ext until the MEASURED reach lands in [lo, hi]."""
    P = copy.deepcopy(P)
    hist = []
    for _ in range(iters):
        g = G(P)
        fn, b = PARTS["carrier"][0](g)
        fld, org, spc = evaluate(fn, b, spacing_for(b, P["budget"], P["voxel"]))
        m = polygonise(fld, org, spc)
        got = measure_notch_reach(m, g, P)
        hist.append((P["notch_sector_ext"], got))
        if not math.isfinite(got) or got <= 1e-6:
            break
        if lo <= got <= hi:
            break
        P["notch_sector_ext"] *= target / got
    P["notch_measured_reach"] = hist[-1][1]
    P["notch_calib_hist"] = hist
    return P


def solve_body_trim(P):
    """Largest trim that keeps every internal margin.

    X and Y are bounded analytically (driver pocket wall, faceplate magnet rim).
    Z is bounded by the acoustic void, which has to be measured, so it is found
    by bisection on the real geometry.
    """
    g = G(P)
    tx, ty, tz = g.trim_got
    target = 1000.0 * P["body_trim_min_front_cc"]
    if tz > 0:
        def vol(t):
            Q = copy.deepcopy(P)
            Q["body_trim_force"] = (tx, ty, t)
            return acoustic_void(G(Q))
        if vol(tz) < target:
            lo, hi = 0.0, tz
            for _ in range(6):
                mid = 0.5 * (lo + hi)
                if vol(mid) >= target:
                    lo = mid
                else:
                    hi = mid
            tz = lo
    P = copy.deepcopy(P)
    P["body_trim_force"] = (float(tx), float(ty), float(tz))
    return P


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
    ap.add_argument("--trim", type=float, default=None,
                    help="mm of protrusion to remove by shrinking the shell (default 5)")
    ap.add_argument("--cant", type=float, default=None,
                    help="nozzle cant in degrees about +Y (default 45)")
    ap.add_argument("--no-assembly", action="store_true")
    ap.add_argument("--tip-style", choices=["bell", "carrier"], default=None,
                    help="bell (default) or the cancelled mag-float carrier + skirt")
    ap.add_argument("--bell-size", choices=list(BELL_SIZES), default=None,
                    help="lip size shown in the assembly (default M)")
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
    if args.trim is not None:
        P["body_trim_mm"] = args.trim
    if args.tip_style is not None:
        P["tip_style"] = args.tip_style
    if args.bell_size is not None:
        P["bell_asm_size"] = args.bell_size
    P = solve_body_trim(P)
    if P["notch_sector_ext"] > 0:
        P = calibrate_notch(P)
    g = G(P)

    here = os.path.dirname(os.path.abspath(__file__))
    out = args.out or os.path.join(here, "stl")
    os.makedirs(os.path.join(out, "right"), exist_ok=True)
    if args.ear in ("left", "both"):
        os.makedirs(os.path.join(out, "left"), exist_ok=True)

    names = enabled_parts(P) if args.all else args.part
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
    if P["tip_style"] == "bell":
        print(f"tip = BELL (S/M/L, assembly shows {P['bell_asm_size']}): nose "
              f"Ø{P['bell_nose_tip_d']:.1f}->Ø{P['bell_nose_base_d']:.1f} over "
              f"{P['bell_nose_len']:.2f} mm ({g.bell_half_angle:.1f} deg half-angle), "
              f"lip tube Ø{P['bell_lip_tube_d']:.1f} wall {P['bell_lip_wall']:.2f} "
              f"({'C-section, slit ' + format(P['bell_lip_slit_deg'], '.0f') + ' deg' if P['bell_lip_hollow'] else 'solid bead'}), "
              f"web {P['bell_web_t']:.2f} mm")
        print(f"      seated on the Ø{P['bell_nozzle_od']:.1f} Ti insert tube (bore "
              f"Ø{P['bell_bore']:.1f}), ridge in a {P['bell_groove_w']:.1f} x "
              f"{P['bell_groove_d']:.2f} groove at x = {g.bell_groove_x0:.2f}.."
              f"{g.bell_groove_x1:.2f}; Ti ends {g.bell_noz_end:.2f}, tip face "
              f"{g.bell_tip_x:.2f}, rim stop {g.bell_cb_x:.2f}, lip plane "
              f"{g.bell_lip_x:.2f} (carrier was {g.carrier_x1:.2f} at the seal plane)")
        print(f"      mag-float carrier + skirt: CANCELLED, kept behind --tip-style carrier")
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
        # the bell family also lands flat in stl/tips/ as bell_tip_<size>_<R|L>...
        if n.startswith("bell_tip_"):
            os.makedirs(os.path.join(out, "tips"), exist_ok=True)
            size, _, rest = n[len("bell_tip_"):].partition("_")
            suffix = f"_{rest}" if rest else ""
            if args.ear in ("right", "both"):
                mesh.export(os.path.join(out, "tips", f"bell_tip_{size}{suffix}_R.stl"))
            if args.ear in ("left", "both"):
                mirror(mesh).export(os.path.join(out, "tips",
                                                 f"bell_tip_{size}{suffix}_L.stl"))

    # ---- assembly ------------------------------------------------------
    if args.all and not args.no_assembly:
        parts = []
        for k in assembly_parts(P):
            if k not in meshes:
                continue
            mk = meshes[k]
            if k in NOZZLE_FRAME_PARTS:            # built in the nozzle-local frame
                mk = mk.copy()
                mk.apply_transform(g.nozzle_T)
            parts.append(mk)
        if P["wing_style"] != "gyroid":
            for pl in g.plungers:
                T = np.eye(4)
                T[:3, 0], T[:3, 1], T[:3, 2] = pl["aim"], pl["u"], pl["v"]
                T[:3, 3] = pl["mount"]
                pick = {"plunger_pad": ("plunger_pad_cymba"
                                        if pl["pad_extra"] else "plunger_pad"),
                        "plunger_cam": ("plunger_cam_ext" if pl["ext"]
                                        else "plunger_cam")}
                for k in PLUNGER_PARTS:
                    src = pick.get(k, k)
                    if src not in meshes:
                        continue
                    mk = meshes[src].copy()
                    if k == "plunger_cam":       # cam sits below the fixed ring
                        Tc = T.copy()
                        Tc[:3, 3] = pl["mount"] - pl["aim"] * (
                            pl["cam_h"] + P["plunger_mag_t"] + 0.15)
                        mk.apply_transform(Tc)
                    else:
                        mk.apply_transform(T)
                    parts.append(mk)
        if parts:
            asm = trimesh.util.concatenate(parts)
            rows.append(report_row("assembly", asm, 0.0, 0.0))
            print(rows[-1])
            asm.export(os.path.join(out, "right", "assembly.stl"))
            if args.ear in ("left", "both"):
                mirror(asm).export(os.path.join(out, "left", "assembly.stl"))

    # ---- overhang analysis --------------------------------------------
    if P["wing_style"] != "gyroid" and ("jacket_wing" in meshes or args.all):
        print("-" * 100)
        print(f"wing mechanism = THREE RADIAL MAG-PLUNGERS "
              f"({P['plunger_mag_od']}x{P['plunger_mag_id']}x{P['plunger_mag_t']} mm "
              f"N35 pairs, rest gap {P['plunger_gap']} mm)")
        print(f"  depth stack (fixed ring back -> moving ring face): "
              f"{g.pl_depth_stack:.2f} mm      dynamic travel: "
              f"+/-{P['plunger_travel']:.2f} mm "
              f"(stops at s = {g.pl_stop_in:.2f} and {g.pl_stop_out:.2f} mm, "
              f"{g.pl_stop_out - g.pl_stop_in:.2f} mm total)")
        print(f"  cam preset: {P['plunger_cam_steps']} detents over "
              f"{P['plunger_cam_range']:.1f} mm of coarse engagement; guide pin "
              f"Ø{P['plunger_pin_od']} in a {P['plunger_pin_sleeve']} mm sleeve")
        if "plunger_foot" in meshes:
            fm = meshes["plunger_foot"]
            sk = float(fm.bounds[0][0])
            print(f"  MEASURED on plunger_foot.stl: stop skirt reaches s = "
                  f"{sk:.2f} mm -> {g.pl_mag_mov - sk:.2f} mm of inward travel "
                  f"before it bottoms on the boss")
        clash = []
        for pl in g.plungers_all:
            tip = pl["mount"] + pl["aim"] * (g.pl_pad1 + P["plunger_rocker"])
            cl = plunger_clearance(g, pl)
            ang = math.degrees(math.acos(min(1.0, abs(float(
                np.dot(pl["aim"], g.n_ax))))))
            if cl < 0.3 and pl["enabled"]:
                clash.append((pl["name"], cl, ang))
            print(f"    [{'ON ' if pl['enabled'] else 'OFF'}] "
                  f"{pl['name']:19s} aim ({pl['aim'][0]:+.2f},"
                  f"{pl['aim'][1]:+.2f},{pl['aim'][2]:+.2f})  cam "
                  f"{pl['cam_steps']}x{pl['cam_range']:.1f}mm  boss "
                  f"{pl['boss_h']:.2f} mm  pad Ø"
                  f"{P['plunger_foot_od'] + pl['pad_extra']:.1f}  reach "
                  f"{float(np.linalg.norm(tip - np.array(g.core_c))):5.2f} mm  "
                  f"{ang:4.1f} deg off the nozzle axis  clearance "
                  f"{cl:+6.2f} mm"
                  + ("  *** CLASH ***" if cl < 0.3 else ""))
        if clash:
            print(f"  *** {len(clash)} plunger(s) interfere with the nozzle/carrier/"
                  f"skirt stack: "
                  + ", ".join(f"{n} {c:+.2f} mm ({a:.0f} deg off the nozzle axis)"
                              for n, c, a in clash))
        off = [q["name"] for q in g.plungers_all if not q["enabled"]]
        if off:
            print(f"  BUILD = {len(g.plungers)}-leg variant; disabled via "
                  f"plunger_enable: {', '.join(off)}")
        # ---- leg-3 feasibility: is there ANY base that clears the stack?
        fs = leg3_feasibility(g)
        b = fs["best"]
        print(f"  leg-3 base sweep (base is free, aim must still end at the "
              f"tragus wall): {fs['n_ok']} of the swept anterior bases clear "
              f"0.80 mm")
        print(f"    the target sits {fs['tgt_perp']:.2f} mm off the nozzle axis at "
              f"station {fs['tgt_station']:.2f} mm; the nearest stack station is "
              f"{fs['tgt_near_station']:.2f} mm at {fs['tgt_stack_r']:.2f} mm radius"
              + ("  -- the tragus wall lies inside the skirt's shadow, so no base "
                 "can reach it" if fs["tgt_perp"] < fs["tgt_stack_r"] else ""))
        if b is not None:
            print(f"    best base ({b['base'][0]:+6.2f},{b['base'][1]:+6.2f},"
                  f"{b['base'][2]:+6.2f})  aim ({b['aim'][0]:+.2f},"
                  f"{b['aim'][1]:+.2f},{b['aim'][2]:+.2f})  {b['dev']:4.1f} deg "
                  f"off nominal  {b['axis_deg']:4.1f} deg off the nozzle axis  "
                  f"clearance {b['clear']:+.2f} mm"
                  + ("  -> FEASIBLE" if b["clear"] >= 0.8 else "  -> still clashes"))
        if P["cymba_lip_bias"]:
            print(f"  cymba lip bias: aim rotated {P['cymba_lip_bias']:.0f} deg "
                  f"toward +Y, pad Ø{P['plunger_foot_od']:.1f} -> Ø"
                  f"{P['plunger_foot_od'] + P['cymba_pad_extra']:.1f} with a "
                  f"{P['plunger_pad_roll']:.2f} mm rolled shoulder, so contact "
                  f"happens under the cymba lip")
        if "jacket_wing" in meshes:
            m_ti = meshes["jacket_wing"].volume * 4.43e-3      # g, Ti-6Al-4V
            print(f"  jacket + 3 bosses = {meshes['jacket_wing'].volume:.0f} mm3 "
                  f"= {m_ti:.2f} g Ti.  The bosses dominate it; if mass matters, "
                  f"thin plunger_boss_od or lattice them.")
        if g.boot is not None:
            print(f"  cable exit boot: '{P['cable_exit']}', "
                  f"{P['cable_boot_len']:.1f} mm, Ø{P['cable_boot_od0']}->"
                  f"Ø{P['cable_boot_od1']}, Ø{P['cable_bore']} bore, raked "
                  f"{P['cable_boot_angle'] if P['cable_exit']=='up_back' else 0:.0f} deg up")

    if P["wing_style"] == "gyroid" and (args.all or "jacket_wing" in meshes):
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
        print(f"wing: free span {wr['L_free']:.2f} mm "
              f"(shortened {P['wing_shorten']:.2f} mm, splayed "
              f"{P['wing_splay_deg']:+.0f} deg about the root, rise "
              f"{g.wing_rise:.2f} mm)")
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

    # ---- body trim margins (docs/TRYON_REPORT.md, recalibrated)
    if "core" in meshes:
        print("-" * 100)
        print(f"body trim: requested {P['body_trim_mm']:.1f} mm of protrusion, "
              f"weights {P['body_trim_w']}")
        print(f"  per-axis cut (mm of full extent)  X/Y/Z requested "
              f"{g.trim_req[0]:.2f}/{g.trim_req[1]:.2f}/{g.trim_req[2]:.2f}   "
              f"feasible {g.trim_cap[0]:.2f}/{g.trim_cap[1]:.2f}/{g.trim_cap[2]:.2f}   "
              f"applied {g.trim_got[0]:.2f}/{g.trim_got[1]:.2f}/{g.trim_got[2]:.2f}")
        print(f"  protrusion removed by the trim: {g.trim_protrusion_got:.2f} mm "
              f"of the {g.trim_protrusion_req:.2f} mm asked for"
              + (f"   [X blocked: {g.trim_block_x}]" if g.trim_block_x else ""))
        print(f"  core {2*g.core_rx:.2f} x {2*g.core_ry:.2f} x {2*g.core_rz:.2f} mm")

        # internal margins
        mar = [("driver pocket wall, +/-X", g.core_rx - g.pocket_r),
               ("driver pocket wall, +/-Y", g.core_ry - g.pocket_r),
               ("shell under the driver pocket", g.core_rz + g.pocket_z0),
               ("faceplate cap above the parting plane", g.core_rz - g.z_cut),
               ("socket pocket to shell bottom",
                g.core_rz + (P["socket_z"] - 0.5 * P["socket_h"]))]
        for (mx, my) in g.fp_mags:
            r = math.hypot(mx - g.core_cx, my)
            mar.append((f"faceplate magnet rim at ({mx:.1f},{my:.1f})",
                        r - g.pocket_r - 0.5 * P["jmag_dia"]))
        for (mx, my) in g.jacket_mags:
            zs = _lower_z(g, mx, my)
            mar.append((f"jacket magnet depth at ({mx:.1f},{my:.1f})",
                        -zs - P["jmag_depth"] - P["core_wall"] * 0.0))
        bad = [n for n, v in mar if v < 0.30]
        for n, v in mar:
            print(f"    {'FLAG' if v < 0.30 else '    '} {n:44s} {v:6.2f} mm")

        # acoustic volume actually left inside the shell (same probe the solver used)
        vol = acoustic_void(g)
        print(f"    {'FLAG' if vol < 500 else '    '} "
              f"{'acoustic void inside the shell':44s} {vol:6.1f} mm3 "
              f"(front-volume target 500 mm3 = 0.5 cc)")
        if bad or vol < 500:
            print(f"  *** {len(bad) + (vol < 500)} margin(s) under threshold: "
                  + ", ".join(bad + ([] if vol >= 500 else ["acoustic volume"])))

    # ---- posterior-inferior corner roll (docs/TRYON_REPORT.md v2)
    if ("core" in meshes and "faceplate" in meshes):
        cn = g.corner_n
        after = max(float((mm.vertices @ cn).max())
                    for mm in (meshes["core"], meshes["faceplate"]))
        P0 = copy.deepcopy(P)
        P0["corner_chamfer"] = 0.0
        P0["body_trim_force"] = (0.0, 0.0, 0.0)     # baseline = untrimmed, unrolled
        g0 = G(P0)
        before = 0.0
        for nm in ("core", "faceplate"):
            f0, b0 = PARTS[nm][0](g0)
            fld, org, spc = evaluate(f0, b0, 0.22)
            before = max(before, float((polygonise(fld, org, spc).vertices @ cn).max()))
        print("-" * 100)
        print(f"corner roll {P['corner_chamfer']:.1f} mm along "
              f"({cn[0]:+.2f},{cn[1]:+.2f},{cn[2]:+.2f}) with a "
              f"{P['corner_roll']:.1f} mm blend")
        print(f"  internals allow {g.corner_c_core:.2f} mm below z = "
              f"{g.corner_z_lo:.1f} mm (driver pocket wall {P['corner_min_wall']} mm), "
              f"so the roll lives in the faceplate")
        print(f"  worst point along that diagonal, chamfer + trim COMBINED: "
              f"{before:.2f} -> {after:.2f} mm = {before - after:.2f} mm removed "
              f"(measured on the meshes)")
        print(f"  for comparison, the report's per-axis coefficients predict "
              f"{g.trim_protrusion_got:.2f} mm from the trim alone; they include "
              f"re-seating gains this generator cannot see")
        print(f"  core volume {meshes['core'].volume:.1f} mm3, faceplate "
              f"{meshes['faceplate'].volume:.1f} mm3")

    # ---- protrusion stack along the nozzle axis (docs/TRYON_REPORT.md rec 1)
    rigid = [meshes[k] for k in ("core", "faceplate", "jacket_wing") if k in meshes]
    if rigid:
        lo = min(float(((mm.vertices - g.nozzle_base) @ g.n_ax).min()) for mm in rigid)
        print("-" * 100)
        print(f"nozzle cant {P['nozzle_cant_deg']:.0f} deg about +Y   "
              f"axis {tuple(round(float(v), 3) for v in g.n_ax)}   "
              f"nozzle base {tuple(round(float(v), 2) for v in g.nozzle_base)}")
        print(f"  rigid body along the nozzle axis: {lo:.2f} .. {g.tip_x0:.2f} mm "
              f"= {g.tip_x0 - lo:.2f} mm stack   "
              f"(uncanted baseline 23.40 mm, TRYON_REPORT.md rec 1)")
        print(f"  seal plane at {g.seal_x:.2f} mm; everything behind it spans "
              f"{g.seal_x - lo:.2f} mm along the axis; tip face at {g.tip_x1:.2f} mm")

    # ---- bell tip: sizes, volumes, walls, insertion -- measured on the meshes
    built = [s for s in BELL_SIZES if f"bell_tip_{s}" in meshes]
    if built:
        print("-" * 100)
        print(f"bell tip: nose Ø{P['bell_nose_tip_d']:.1f}->Ø{P['bell_nose_base_d']:.1f} "
              f"x {P['bell_nose_len']:.2f} mm, {g.bell_half_angle:.1f} deg; lip tube "
              f"Ø{P['bell_lip_tube_d']:.1f}, wall {P['bell_lip_wall']:.2f} "
              f"({P['bell_ant_wall']:.2f} anterior, +{P['bell_ant_free']:.1f} mm free "
              f"length); inferior +{P['bell_inf_ext']:.1f} mm; web {P['bell_web_t']:.2f}; "
              f"vent Ø{P['bell_vent_dia']:.1f} at {P['bell_vent_az_deg']:.0f} deg")
        print(f"  {'size':4s} {'spec HxW':>10s} {'built HxW':>12s} {'vol mm3':>8s} "
              f"{'lip':>5s} {'ant':>5s} {'web':>5s} {'slv':>5s} {'nose':>5s} "
              f"{'min':>5s} {'p05':>5s}  {'insert':>6s} {'protr':>6s} "
              f"{'web L':>9s} {'rock':>10s}")
        for s in built:
            m = bell_measure(g, s, meshes[f"bell_tip_{s}"])
            h, w = P["bell_sizes"][s]
            print(f"  {s:4s} {h:5.1f}x{w:4.1f} {m['H']:6.2f}x{m['W']:5.2f} "
                  f"{m['vol']:8.1f} {m['lip']:5.2f} {m['lip_ant']:5.2f} {m['web']:5.2f} "
                  f"{m['sleeve']:5.2f} {m['nose']:5.2f} {m['min_wall']:5.2f} "
                  f"{m['p05_wall']:5.2f}  {m['insertion']:6.2f} {m['protrusion']:6.2f} "
                  f"{m['web_len_min']:4.2f}-{m['web_len_max']:4.2f} "
                  f"{m['rock_deg']:4.1f}deg@{m['rock_az']:4.0f}")
        print(f"  walls are MEASURED through the mesh (mm): lip = tube wall at the "
              f"superior pole (spec {P['bell_lip_wall']:.2f}), ant = anterior sector "
              f"(spec {P['bell_ant_wall']:.2f}), web (spec {P['bell_web_t']:.2f}), "
              f"slv = sleeve over the Ti tube (spec {P['bell_sleeve_wall']:.2f}), "
              f"nose = bore wall at the tip face; min/p05 = ray thickness over 2500 "
              f"surface points.  insert = tip face to the most distal Ø"
              f"{P['bell_nose_base_d']:.0f} station (the rim stop), spec <= "
              f"{P['bell_nose_len']:.1f}.  rock = inward rock before the lip's inner "
              f"edge meets the sleeve, at the worst azimuth (spec +/-15).  H includes "
              f"the inferior extension.")
        for s in built:
            if f"bell_tip_{s}_mold_a" in meshes:
                frac, worst = bell_draw_check(g, s)
                print(f"  {s} mould, split y=0 pull +/-Y: {100 * frac:.1f}% of cavity "
                      f"points shadowed along the pull, worst overhang "
                      f"{worst:.2f} mm of silicone  -> "
                      + ("pulls clean" if frac < 0.01 else
                         "pulls with a flex of the lip (soft part, hard mould)"
                         if worst < 1.0 else "*** UNDERCUT ***"))

    # ---- skirt contact land + pressure budget (docs/MECH_VALIDATION.md JOB 2)
    if "carrier" in meshes:
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
        # intertragic-notch sector + two-part mould draw check
        hw = math.radians(0.5 * P["notch_sector_deg"])
        tr = math.radians(P["notch_sector_trans_deg"])
        th = np.linspace(-np.pi, np.pi, 2001)
        dth = np.abs(_wrap(th - math.radians(P["notch_sector_center_deg"])))
        tb = np.clip((dth - (hw - tr)) / tr, 0.0, 1.0)
        bb = 1.0 - tb * tb * (3.0 - 2.0 * tb)
        rr_ = g.skirt_rim_r + P["notch_sector_ext"] * bb
        yy, zz = rr_ * np.cos(th), np.abs(rr_ * np.sin(th))
        draw = []
        for sgn in (1, -1):
            m = sgn * yy >= 0
            o = np.argsort(sgn * yy[m])
            draw.append(int((np.diff(zz[m][o]) > 1e-6).sum()))
        print(f"  notch sector: {P['notch_sector_deg']:.0f} deg centred "
              f"{P['notch_sector_center_deg']:.0f} deg from +Y (inferior), "
              f"+{P['notch_sector_ext']:.2f} mm radial reach -> Ø"
              f"{2*(g.skirt_rim_r+P['notch_sector_ext']):.1f} mm there, Ø"
              f"{2*g.skirt_rim_r:.1f} mm elsewhere; land wall "
              f"{P['notch_sector_wall']} vs {P['skirt_wall_land']} mm, "
              f"{P['notch_sector_trans_deg']:.0f} deg blends")
        w_in = measure_skirt_wall(g, P, P["notch_sector_center_deg"])
        w_out = measure_skirt_wall(g, P, P["notch_sector_center_deg"] + 180.0)
        mr = measure_notch_reach(meshes["carrier"], g, P) if "carrier" in meshes else 0.0
        print(f"  notch = COMPLIANCE, not reach (v4): rim is plain Ø"
              f"{P['skirt_max_dia']:.1f} mm all round, measured sector reach "
              f"{mr:+.2f} mm")
        print(f"  MEASURED land wall: {w_in:.3f} mm in the sector vs "
              f"{w_out:.3f} mm outside (spec {P['notch_sector_wall']} / "
              f"{P['skirt_wall_land']} mm); hinge {P['notch_hinge_wall']} mm over "
              f"{g.notch_hinge_span:.2f} mm of free length vs "
              f"{P['skirt_wall_hinge']} mm over {P['skirt_hinge_w']:.2f} mm")
        print(f"  two-part mould, split y=0, pull +/-Y: non-monotone rim steps "
              f"{draw[0]} / {draw[1]}  -> "
              + ("demoulds, no undercut" if max(draw) == 0
                 else "*** UNDERCUT, angle the extension ***"))

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
