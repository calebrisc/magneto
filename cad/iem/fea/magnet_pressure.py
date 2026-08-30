"""JOB 2 -- contact-pressure budget for the mag-float skirt and the wing.

No tissue FEA: this is a load-spreading budget.  Take the magnet preload the
mag-float study actually produced, spread it over the skirt's contact band, and
compare the resulting interface pressure against published sustained-skin-load
thresholds.

Geometry comes from ``generate.py`` (skirt rim diameter, flare half-angle, wall);
forces come from ``docs/MAGFLOAT_MAGNETS.md``, preset ``asym_as_built``.

Run:  .venv/bin/python fea/magnet_pressure.py
"""

import json
import math
import os
import sys

import numpy as np

import _common as K

# ---------------------------------------------------------------------------
# thresholds -- see docs/MECH_VALIDATION.md for provenance and caveats
# ---------------------------------------------------------------------------

THRESH = [
    ("shear-derated sustained", 2.15,
     "Landis capillary closing / 2, per Bennett 1979 (shear halves the "
     "occluding pressure).  The number to design to for an all-day device "
     "that micro-slides with jaw motion."),
    ("capillary closing (Landis 1930)", 4.27,
     "32 mmHg arteriolar nailfold capillary pressure.  An ischaemia flag, "
     "not a validated design limit."),
    ("NIV mask injury band (Brill 2018)", 6.34,
     "Low end of measured nasal-bridge mask pressure (47.6 mmHg) that "
     "correlated with discomfort and skin injury in 20 subjects."),
]
P_LANDIS = 4.27      # kPa  (32 mmHg)
P_SHEAR = 2.15       # kPa
P_NIV_LO = 6.34      # kPa  (47.6 mmHg)
P_NIV_HI = 12.25     # kPa  (91.9 mmHg)
MMHG = 7.50062       # mmHg per kPa


def verdict(p_kpa):
    if p_kpa <= P_SHEAR:
        return "comfortable"
    if p_kpa <= P_LANDIS:
        return "borderline"
    if p_kpa <= P_NIV_LO:
        return "too much"
    return "too much (>>)"


def band_pressure(F, D, w, flare_deg=None):
    """Interface pressure over an annular band of diameter D and slant width w.

    Returns (P_nominal, P_normal) in kPa.

    P_nominal  = F / (pi D w)                -- the figure asked for: axial force
                 divided by the wetted band area.
    P_normal   = F / (pi D w sin(flare))     -- the pressure on the surface once
                 the axial force is resolved onto the cone normal.  For a cone
                 whose meridian makes angle `flare` with the axis, only
                 sin(flare) of the surface normal is axial, so the true contact
                 pressure is 1/sin(flare) times the nominal.  This is the
                 conservative number and the physically correct one.
    """
    A = math.pi * D * w                     # mm^2  (slant area of the band)
    p_nom = F / A * 1e3                     # N/mm^2 -> kPa
    if flare_deg is None:
        return p_nom, p_nom
    return p_nom, p_nom / math.sin(math.radians(flare_deg))


def min_band_width(F, D, p_kpa, flare_deg=None):
    """Slant band width (mm) needed to keep the pressure at or under p_kpa."""
    s = 1.0 if flare_deg is None else math.sin(math.radians(flare_deg))
    return F / (math.pi * D * p_kpa * 1e-3 * s)


def main():
    g, P = K.geom()
    m = g.mag
    flare = P["skirt_flare_deg"]
    F = dict(max=m["f_lo"], rest=m["f_rest"], min=m["f_hi"])

    # ---- what the generator actually builds ------------------------------
    slant = (g.skirt_rim_r - g.skirt_root_r) / math.sin(math.radians(flare))
    print("=" * 78)
    print("SKIRT GEOMETRY (from generate.py, preset '%s')" % P["magnet_preset"])
    print("=" * 78)
    print(f"  rim outer dia        {P['skirt_max_dia']:.2f} mm  "
          f"(mid-wall r = {g.skirt_rim_r:.3f} mm)")
    print(f"  root dia             {2*g.skirt_root_r:.2f} mm  (= carrier OD)")
    print(f"  flare half-angle     {flare:.1f} deg from the axis")
    print(f"  axial run            {g.skirt_dx:.2f} mm  (x {g.skirt_root_x:.2f} "
          f"-> {g.skirt_rim_x:.2f})")
    print(f"  available slant      {slant:.2f} mm  -- the widest band the funnel "
          f"can offer")
    print(f"  wall                 {P['skirt_wall']:.2f} mm")
    print(f"  preload              max {F['max']:.3f} N / rest {F['rest']:.3f} N "
          f"/ min {F['min']:.3f} N  ({m['material']}, rest gap {m['gap']} mm)")
    print(f"  axial->normal factor 1/sin({flare:.0f}) = "
          f"{1/math.sin(math.radians(flare)):.3f}")

    # ---- pressure sweep ---------------------------------------------------
    Ds = [10.0, 13.0, 16.0, 19.0]
    ws = [0.5, 1.0, 2.0, 3.0]
    out = {"skirt": [], "wing": [], "min_width": []}

    for case, f in (("max", F["max"]), ("rest", F["rest"]), ("min", F["min"])):
        print("\n" + "-" * 78)
        print(f"CONTACT PRESSURE, F = {f:.3f} N ({case})   "
              f"[nominal kPa | cone-normal kPa | mmHg normal | verdict on normal]")
        print("-" * 78)
        print(f"{'D_contact':>10} " + "".join(f"{'w=%.1f mm' % w:>26}" for w in ws))
        for D in Ds:
            cells = []
            for w in ws:
                pn, pc = band_pressure(f, D, w, flare)
                out["skirt"].append(dict(case=case, F=f, D=D, w=w,
                                         p_nom_kpa=pn, p_norm_kpa=pc,
                                         p_norm_mmhg=pc * MMHG,
                                         verdict=verdict(pc)))
                cells.append(f"{pn:7.2f}|{pc:6.2f}|{verdict(pc)[:9]:>9} ")
            print(f"{D:8.0f} mm " + "".join(f"{c:>26}" for c in cells))

    # ---- minimum band width ----------------------------------------------
    print("\n" + "=" * 78)
    print("MINIMUM SLANT BAND WIDTH at F_max = %.3f N (cone-normal pressure)" % F["max"])
    print("=" * 78)
    print(f"{'D_contact':>10} {'w for 2.15 kPa':>16} {'w for 4.27 kPa':>16} "
          f"{'w for 6.34 kPa':>16}   available")
    for D in Ds:
        row = [min_band_width(F["max"], D, t, flare) for t in (P_SHEAR, P_LANDIS, P_NIV_LO)]
        out["min_width"].append(dict(D=D, w_shear=row[0], w_landis=row[1],
                                     w_niv=row[2], available=slant))
        print(f"{D:8.0f} mm {row[0]:14.2f} mm {row[1]:14.2f} mm {row[2]:14.2f} mm"
              f"   {slant:.2f} mm"
              + ("" if row[0] <= slant else "   <-- shear-derated target NOT reachable"))

    # ---- the wing, for comparison ----------------------------------------
    print("\n" + "=" * 78)
    print("WING, for comparison -- intended 0.3-0.4 N over its contact patch")
    print("=" * 78)
    print("  (patch length along the blade x effective contact width; the blade "
          "is\n   %.1f mm wide in Z, so 7 mm = full-face contact and 1 mm = an "
          "edge line contact)" % P["wing_width"])
    print(f"\n{'F (N)':>6} {'len (mm)':>9} {'width (mm)':>11} {'area (mm2)':>11} "
          f"{'P (kPa)':>9} {'P (mmHg)':>9}  verdict")
    for f in (0.30, 0.40):
        for L in (8.0, 12.0):
            for w in (1.0, 2.0, 4.0, P["wing_width"]):
                A = L * w
                p = f / A * 1e3
                out["wing"].append(dict(F=f, length=L, width=w, area=A,
                                        p_kpa=p, verdict=verdict(p)))
                print(f"{f:6.2f} {L:9.1f} {w:11.1f} {A:11.1f} {p:9.2f} "
                      f"{p*MMHG:9.1f}  {verdict(p)}")

    dest = os.path.join(os.path.dirname(os.path.abspath(__file__)),
                        "pressure_results.json")
    with open(dest, "w") as fh:
        json.dump(out, fh, indent=2)
    print(f"\nwrote {dest}")


if __name__ == "__main__":
    main()
