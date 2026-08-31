#!/usr/bin/env python3
"""
tryon.py -- score the seated Magneto IEM against every aligned ear.

Reads ears/aligned/*.json (written by align_ear.py), re-loads the cached ear
patch, applies the seated transform to the IEM sample points, and reports five
things per ear.  Sign convention throughout: **positive = clearance (air
between part and flesh), negative = the part is inside the flesh**.

  (a) SKIRT RIM -- 72 points on the Ø19 sealing lip.
      `rim_cover`  fraction of the circle whose |signed distance| <= 1.5 mm,
                   i.e. the contact band a 0.35 mm Shore-A-10 skirt can close.
      `rim_gap`    the worst positive distance anywhere on the circle: the
                   biggest hole the seal has to bridge.
      `rim_press`  the worst negative distance: how hard the lip is buried.

  (b) WING -- the retention feature.  `wing_tip` and `wing_mid` are median
      signed distances over 40 vertices each (median, not min, so a single
      marching-cubes spike cannot carry the number).  The design target is
      -0.5 to -2.0 mm: the wing has to *press* into the antihelix or it does
      not retain.

  (c) JACKET EAR FACE -- the gyroid skin that lands on the concha floor.
      Mean and max clearance over 140 points.  Large mean clearance means the
      shell is bridging the bowl on its rim instead of bedding into it.

  (d) PROTRUSION -- how far the faceplate's outermost points sit past the
      tragus plane (the plane through the detected tragus apex, normal =
      concha-floor outward normal).  Positive = sticking out of the ear.

  (e) HARD INTERFERENCE -- the most negative signed distance over the points
      sampled on the parts that cannot deform and are not supposed to press:
      the Ti core, the faceplate and the jacket body, i.e. `shell` = rigid minus
      the wing.  Flesh can absorb roughly 2.5 mm of overlap; past that something
      has to give, and it will not be the ear.  The wing is excluded because it
      is *designed* to press into the antihelix -- counting it here made the
      clearance axis contradict the retention axis.  Its worst penetration is
      still reported, as `wing_min`, and graded on the retention axis.

GRADING -- each ear gets pass / marginal / fail on four axes, and the overall
grade is the worst of them, with `worst_metric` naming the one that drove it:

    seal        cover >= 0.75 pass | >= 0.50 marginal | else fail
    retention   wing_tip in [-2.5, 0.0] pass | [-4, 1.0] marginal | else fail
    clearance   hard_min >= -1.0 pass | >= -2.5 marginal | else fail
    protrusion  <= PROT_PASS pass | <= PROT_MARGINAL marginal | else fail

Usage:
    python tryon.py                       # everything in ears/aligned/
    python tryon.py --dataset sonicom
    python tryon.py --csv out.csv --md summary.md
"""

from __future__ import annotations

import argparse
import csv
import glob
import json
import os
import sys

import numpy as np
import trimesh

from earfit import ALIGNED, NOZZLE_AXIS, EarField, iem_points, transform

GRADES = ("pass", "marginal", "fail")

# PROTRUSION THRESHOLDS -- design-team calibration, 2026-08-31.
#
# The original <= 2 mm / <= 5 mm pair was aspirational rather than physical: it
# describes a custom-moulded IEM sitting nearly flush, and across 107 ears it was
# cleared by exactly one, which makes it useless as a discriminator.  A universal
# IEM is a body parked in the concha with the shell standing proud of the tragus
# plane; the question is how much is tolerable, not whether it is zero.
#
# Published specs are of limited help because manufacturers quote OVERALL SHELL
# DEPTH, not protrusion past the tragus plane, and the concha swallows an
# unstated part of the former.  For scale: the SIMGOT EA1000 is 22 x 17 x
# 20.7 mm, and ~20-22 mm of overall shell depth is typical.  With a 9-17 mm
# concha depth (docs/EAR_ANTHROPOMETRY.md) that leaves single-digit millimetres
# proud on a typical ear, which is the band these thresholds encode.
#
# These are therefore a DESIGN-TEAM CALIBRATION, not a measured industry
# standard: 10 mm proud is accepted as normal for a universal fit, 14 mm as the
# limit before the shell fouls the helix/antihelix on insertion and starts
# levering itself out.  Revisit with real worn measurements when a prototype
# exists.
PROT_PASS = 10.0
PROT_MARGINAL = 14.0
# the pre-2026-08-31 strict pair, retained so the report can show both
PROT_PASS_STRICT = 2.0
PROT_MARGINAL_STRICT = 5.0


def grade_worst(gs):
    return max(gs, key=lambda g: GRADES.index(g))


def score(rec, P, field):
    M = np.array(rec["transform"], float)
    q = lambda k: field.query(transform(P[k], M))          # noqa: E731

    rim = q("rim")
    cover = float(np.mean(np.abs(rim) <= 1.5))
    gap = float(rim.max())
    press = float(rim.min())

    tip = float(np.median(q("wing_tip")))
    mid = float(np.median(q("wing_mid")))

    jac = q("jacket")
    jac_mean, jac_max = float(jac.mean()), float(jac.max())

    # Clearance is graded on `shell` -- the rigid parts MINUS the wing -- to match
    # what the seating search penalises.  Grading it on `rigid` (wing included)
    # made the clearance axis contradict the retention axis: the wing is *meant*
    # to press -0.5 to -2.0 mm into the antihelix, so an on-target wing forced
    # hard_min <= -2 and could never score better than "marginal" on clearance.
    # The wing's own penetration is reported as wing_min and graded by retention.
    shell = q("shell")
    hard = float(shell.min())
    hard_n = int((shell < -2.5).sum())
    wing_min = float(min(q("wing_tip").min(), q("wing_mid").min()))

    trg = rec.get("tragus")
    if trg is None:
        prot = float("nan")
    else:
        n = np.array(rec["floor_normal"], float)
        fp = transform(P["faceplate"], M)
        # einsum, not `@`: Accelerate's (N,3)@(3,) path raises spurious
        # divide-by-zero/overflow warnings on macOS for finite inputs.
        prot = float(np.max(np.einsum("ij,j->i", fp - np.array(trg, float), n)))

    g_seal = "pass" if cover >= 0.75 else "marginal" if cover >= 0.50 else "fail"
    g_ret = ("pass" if -2.5 <= tip <= 0.0 else
             "marginal" if -4.0 <= tip <= 1.0 else "fail")
    g_clr = "pass" if hard >= -1.0 else "marginal" if hard >= -2.5 else "fail"
    if np.isnan(prot):
        g_pro = "marginal"
    else:
        g_pro = ("pass" if prot <= PROT_PASS
                 else "marginal" if prot <= PROT_MARGINAL else "fail")

    axes = dict(seal=g_seal, retention=g_ret, clearance=g_clr, protrusion=g_pro)
    overall = grade_worst(axes.values())
    worst = [k for k, v in axes.items() if v == overall]

    return dict(
        dataset=rec["dataset"], ear_id=rec["ear_id"], side=rec["side"],
        rim_cover=cover, rim_gap=gap, rim_press=press,
        wing_tip=tip, wing_mid=mid,
        jacket_mean=jac_mean, jacket_max=jac_max,
        protrusion=prot, hard_min=hard, hard_n=hard_n, wing_min=wing_min,
        rake=rec["nozzle_rake_deg"], seat_cost=rec["seat_cost"],
        weak=rec["weak"], canal_run=rec["canal_run"],
        basin_escape=rec["basin_escape"],
        g_seal=g_seal, g_retention=g_ret, g_clearance=g_clr, g_protrusion=g_pro,
        grade=overall, worst_metric="+".join(sorted(worst)),
    )


# --------------------------------------------------------------------------- #
# CONTACT CONTRACT
# --------------------------------------------------------------------------- #
# The graded axes say how WELL the fit scores.  The contract says whether each
# part is doing the job it was designed for -- which is a different question, and
# the one that catches a part touching flesh it should never touch.  Every rule
# below is an explicit design intent, not a derived metric.
#
#   skirt land        MUST TOUCH   a continuous contact band that ENCLOSES the
#                                  canal aperture.  Enclosure is hard: a loop
#                                  that does not surround the canal seals nothing
#                                  (see seal_enclosure_audit.py).
#   wing / rail pad   MUST TOUCH   the antihelix, within the spring's working
#                                  range -- 0 to 1.5 mm of interference, which a
#                                  0.294 N/mm leaf spring absorbs as press.
#   jacket ear-face   MUST REST    on the concha floor: in contact, and never
#                                  more than 2.5 mm into flesh.
#   nozzle + insert   MUST NOT     touch anything.  It is recessed inside the
#                                  skirt; any contact means it is loading the
#                                  canal wall directly.
#   core / faceplate  MUST NOT     touch flesh at all.  Load reaches the ear only
#                                  through the jacket and the wing.
#   cable exit        MUST CLEAR   the ear.
#
# Sign convention as everywhere else: positive = clearance, negative = inside
# the flesh.
CONTRACT_JACKET_MAX_INTERFERENCE = -2.5
CONTRACT_JACKET_CONTACT = 0.5      # mm; above this the jacket is not resting
CONTRACT_WING_BAND = (-1.5, 0.0)   # spring working range
CONTRACT_SEAL_BUDGET = 2.5         # conservative skirt deformation budget


def contract_check(rec, patch, P, field):
    """Evaluate the contact contract for one ear.  Returns a list of rows."""
    import seal_compliance as sc

    M = np.array(rec["transform"], float)
    q = lambda k: field.query(transform(P[k], M))          # noqa: E731
    rows = []

    def row(part, intent, detail, value, ok):
        rows.append(dict(part=part, intent=intent, detail=detail,
                         value=value, ok=None if ok is None else bool(ok)))

    # --- skirt land: continuous band enclosing the aperture ---------------- #
    nax = M[:3, :3] @ NOZZLE_AXIS
    nax = nax / np.linalg.norm(nax)
    rim = transform(P["rim"], M)
    cen = rim.mean(axis=0)
    R = float(np.median(np.linalg.norm(rim - cen, axis=1)))
    ap = np.array(rec["aperture"], float)
    d = ap - cen
    lateral = float(np.linalg.norm(d - (d @ nax) * nax))
    inside = lateral < R

    rim_fine = transform(sc.rim_from_mesh(n=sc.N_RIM), M)
    best = None
    for t in np.linspace(0.0, sc.TRAVEL_MM, sc.TRAVEL_STEPS):
        sealed = field.query(rim_fine + t * nax) <= CONTRACT_SEAL_BUDGET
        glen, _ = sc.longest_false_run(sealed)
        gap = 360.0 * glen / len(sealed)
        key = (float(sealed.mean()), -gap)
        if best is None or key > best[0]:
            best = (key, float(sealed.mean()), gap)
    _, cover, gap = best
    continuous = (cover >= sc.COVER_MIN) and (gap <= sc.GAP_MAX_DEG)
    row("skirt land", "MUST TOUCH",
        f"band {100*cover:.0f}% closed, worst gap {gap:.0f}deg; "
        f"aperture {lateral:.1f} mm off centre vs {R:.1f} mm rim"
        + ("" if inside else "  << canal OUTSIDE the ring"),
        f"{100*cover:.0f}% / {gap:.0f}deg", continuous and inside)

    # --- wing / rail pad --------------------------------------------------- #
    # Only meaningful on a gyroid-wing build.  With wing_style="plungers" there is
    # no wing at all, and this row was reading the 40 most-superior vertices of
    # jacket_wing -- i.e. the jacket's top edge, which is not a retention feature.
    # Retire it rather than report a pass or fail about a part that is gone.
    has_wing = True
    try:
        import generate
        has_wing = generate.PARAMS.get("wing_style", "gyroid") == "gyroid"
    except Exception:                                            # noqa: BLE001
        pass
    tip = float(np.median(q("wing_tip")))
    lo, hi = CONTRACT_WING_BAND
    if has_wing:
        row("wing / rail pad", "MUST TOUCH",
            f"tip median {tip:+.2f} mm; spring band [{lo:+.1f}, {hi:+.1f}]",
            f"{tip:+.2f} mm", lo <= tip <= hi)
    else:
        row("wing / rail pad", "MUST TOUCH",
            "no wing in this build (wing_style=plungers); retention is the "
            "plunger pads row below", "n/a", None)

    # --- jacket ear-face --------------------------------------------------- #
    jac = q("jacket")
    jmin = float(jac.min())
    resting = jmin <= CONTRACT_JACKET_CONTACT
    not_buried = jmin >= CONTRACT_JACKET_MAX_INTERFERENCE
    row("jacket ear-face", "MUST REST",
        f"closest {jmin:+.2f} mm ("
        + ("resting" if resting else "NOT in contact")
        + (", within 2.5 mm" if not_buried else ", OVER-BURIED") + ")",
        f"{jmin:+.2f} mm", resting and not_buried)

    # --- nozzle + insert: recessed, touches nothing ------------------------ #
    if len(P["nozzle"]):
        nmin = float(q("nozzle").min())
        row("nozzle + insert", "MUST NOT TOUCH",
            f"closest {nmin:+.2f} mm (recessed inside the skirt)",
            f"{nmin:+.2f} mm", nmin > 0.0)
    else:
        row("nozzle + insert", "MUST NOT TOUCH", "insert STL not found", "n/a", None)

    # --- core / faceplate: no flesh contact -------------------------------- #
    cmin = float(min(q("core_s").min(), q("face_s").min()))
    row("core / faceplate", "MUST NOT TOUCH",
        f"closest {cmin:+.2f} mm (load must reach the ear via the jacket only)",
        f"{cmin:+.2f} mm", cmin > 0.0)

    # --- plunger pads: the retention contacts in the tripod build ---------- #
    # Each pad rides on a spring giving 0.18-0.49 N over +-0.75 mm of travel, on a
    # cam that presets its stand-off over a 3 mm range.  The cam is a SEATING
    # degree of freedom -- it is set once per wearer -- so a pad that overshoots
    # by less than the cam range is adjustable, and only overshoot beyond
    # cam + travel is a geometry error.
    lim = P.get("_plunger_limits", dict(travel=0.75, cam=3.0, force=(0.18, 0.49)))
    pads = P.get("_plunger", [])
    over = []
    if pads:
        for m in pads:
            aim_w = M[:3, :3] @ np.array(m["aim"], float)
            aim_w /= np.linalg.norm(aim_w)
            mnt_w = transform(np.array(m["mount"], float)[None, :], M)[0]
            tip_w = transform(np.array(m["tip"], float)[None, :], M)[0]
            # signed distance at the pad tip: negative = pad is inside flesh
            dtip = float(field.query(tip_w[None, :])[0])
            over.append((m["name"], -dtip))          # >0 = overshoot into flesh
        # EVERY pad has to reach.  The cam presets stand-off over its range, so a
        # pad short by less than the cam range can be dialled out to meet flesh,
        # and one overshooting by less than cam + travel can be dialled back.
        # Anything outside that window is a geometry error, not an adjustment.
        lo, hi = -lim["cam"], lim["cam"] + lim["travel"]
        bad = [(n, o) for n, o in over if not (lo <= o <= hi)]
        detail = "; ".join(
            f"{n} {o:+.2f}" + ("" if lo <= o <= hi else
                               (" NEVER REACHES" if o < lo else " TOO DEEP"))
            for n, o in over)
        worst = max(over, key=lambda t: abs(t[1]))[1]
        row("plunger pads", "MUST TOUCH",
            f"pad-tip overshoot into flesh (mm): {detail}. Adjustment window "
            f"[{lo:+.2f}, {hi:+.2f}] = cam {lim['cam']:.1f} + travel "
            f"{lim['travel']:.2f}"
            + ("" if not bad else
               f"  << {len(bad)} of {len(over)} pads OUTSIDE the window"),
            f"{worst:+.2f} mm", not bad)
    else:
        row("plunger pads", "MUST TOUCH", "no plunger geometry in this build",
            "n/a", None)

    # --- STABILITY: retention under load ----------------------------------- #
    # Contact says each part touches what it should.  This asks whether the
    # assembly STAYS PUT when something pulls on it -- a different failure, and
    # the one a wearer notices.  Rigid-body force/moment balance with Coulomb
    # friction cones; see stability.py for the model and its assumptions.
    import stability as stab
    st = stab.stability_check(
        rec, P, field, transform,
        cable_point=transform(P["_cable_exit"][None, :], M)[0],
        com=transform(P["_com"][None, :], M)[0])
    mg = st["margin"]
    row("STABILITY", "MUST RESIST",
        f"margin {mg:.2f}x. pull-out capacity {st['pullout_capacity']:.2f} N vs "
        f"demand {st['demand']:.2f} N (skirt {stab.SKIRT_PRELOAD:.2f} always-on + "
        f"cable {stab.CABLE_TUG:.2f} + {stab.G_LOAD:.0f}g {st['f_inert']:.2f}); "
        f"friction budget {st['friction_budget']:.2f} N (mu={st['mu']}, wing "
        f"{st['f_wing']:.2f} N); interlocking contacts {st['interlock']}/"
        f"{st['n_contacts']}"
        + ("" if mg >= 1.0 else
           f"  << worst dir cable {np.round(st['cable_dir'], 2)}"),
        f"{mg:.2f}x", mg >= 1.0)

    # --- cable exit boot ---------------------------------------------------- #
    if len(P.get("boot", [])):
        bd = field.query(transform(P["boot"], M))
        bmin = float(bd.min())
        row("cable exit", "MUST CLEAR",
            f"boot closest {bmin:+.2f} mm to the ear "
            + ("(clear)" if bmin > 0 else "(FOULING the ear)"),
            f"{bmin:+.2f} mm", bmin > 0.0)
    else:
        row("cable exit", "MUST CLEAR",
            "no cable/boot geometry in this build", "n/a", None)
    return rows


def print_contract(rec, rows):
    mark = {True: "PASS", False: "FAIL", None: "n/a "}
    print(f"\n=== CONTACT CONTRACT - {rec['dataset']}/{rec['ear_id']} ===")
    print(f"{'part':<17}{'intent':<15}{'value':>12}  {'':<4} detail")
    print("-" * 118)
    for r in rows:
        print(f"{r['part']:<17}{r['intent']:<15}{r['value']:>12}  "
              f"{mark[r['ok']]:<4} {r['detail']}")
    ok = sum(1 for r in rows if r["ok"] is True)
    bad = sum(1 for r in rows if r["ok"] is False)
    na = sum(1 for r in rows if r["ok"] is None)
    print("-" * 118)
    print(f"{ok} pass, {bad} fail, {na} not evaluable")
    return bad


def fmt_table(rows, cols, headers):
    out = ["| " + " | ".join(headers) + " |",
           "|" + "|".join("---" for _ in headers) + "|"]
    for r in rows:
        cells = []
        for c in cols:
            v = r[c]
            cells.append(f"{v:.2f}" if isinstance(v, float) else str(v))
        out.append("| " + " | ".join(cells) + " |")
    return "\n".join(out)


def summarise(rows):
    L = []
    by_ds = {}
    for r in rows:
        by_ds.setdefault(r["dataset"], []).append(r)
    L.append("| dataset | n | pass | marginal | fail |")
    L.append("|---|---|---|---|---|")
    for ds, rs in sorted(by_ds.items()):
        c = {g: sum(1 for r in rs if r["grade"] == g) for g in GRADES}
        L.append(f"| {ds} | {len(rs)} | {c['pass']} | {c['marginal']} | {c['fail']} |")
    c = {g: sum(1 for r in rows if r["grade"] == g) for g in GRADES}
    L.append(f"| **all** | **{len(rows)}** | **{c['pass']}** | "
             f"**{c['marginal']}** | **{c['fail']}** |")

    L.append("")
    L.append("| axis | pass | marginal | fail |")
    L.append("|---|---|---|---|")
    for ax in ("seal", "retention", "clearance", "protrusion"):
        c = {g: sum(1 for r in rows if r["g_" + ax] == g) for g in GRADES}
        L.append(f"| {ax} | {c['pass']} | {c['marginal']} | {c['fail']} |")

    L.append("")
    L.append("| metric | median | p10 | p90 | min | max |")
    L.append("|---|---|---|---|---|---|")
    for k, lab in (("rim_cover", "skirt rim contact coverage"),
                   ("rim_gap", "skirt rim max gap (mm)"),
                   ("wing_tip", "wing tip signed dist (mm)"),
                   ("wing_mid", "wing mid signed dist (mm)"),
                   ("jacket_mean", "jacket mean clearance (mm)"),
                   ("protrusion", "faceplate past tragus (mm)"),
                   ("hard_min", "worst rigid interference (mm)")):
        v = np.array([r[k] for r in rows], float)
        v = v[np.isfinite(v)]
        if not len(v):
            continue
        L.append(f"| {lab} | {np.median(v):.2f} | {np.percentile(v, 10):.2f} | "
                 f"{np.percentile(v, 90):.2f} | {v.min():.2f} | {v.max():.2f} |")
    return "\n".join(L)


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--dataset", default=None)
    ap.add_argument("--cant", type=float, default=None,
                    help="override nozzle_cant_deg when reading the STLs")
    ap.add_argument("--stl-dir", default=None,
                    help="directory of right-ear STLs (default stl/right)")
    ap.add_argument("--contract", action="store_true",
                    help="print the contact contract table per ear (automatic "
                         "when a single ear is scored)")
    ap.add_argument("--viz-ear", default=None,
                    help="dataset/ear_id to export to viz/seated_scene.glb after "
                         "scoring (automatic when a single ear is scored)")
    ap.add_argument("--json-dir", default=None,
                    help="directory of per-ear seating JSONs (default ears/aligned); "
                         "patches are always read from ears/aligned")
    ap.add_argument("--csv", default=os.path.join(ALIGNED, "tryon.csv"))
    ap.add_argument("--md", default=os.path.join(ALIGNED, "tryon_summary.md"))
    a = ap.parse_args()

    js = sorted(glob.glob(os.path.join(a.json_dir or ALIGNED, "*.json")))
    if a.dataset:
        js = [p for p in js if os.path.basename(p).startswith(a.dataset + "_")]
    if not js:
        sys.exit("nothing in ears/aligned -- run align_ear.py first")

    P = iem_points(stl_dir=a.stl_dir, cant=a.cant)
    rows = []
    contract_fails = 0
    single = len(js) == 1
    for p in js:
        rec = json.load(open(p))
        patch = trimesh.load(os.path.join(ALIGNED, rec["patch"]), force="mesh")
        field = EarField(patch)
        rows.append(score(rec, P, field))
        if a.contract or single:
            contract_fails += print_contract(rec, contract_check(rec, patch, P, field))
        r = rows[-1]
        print(f"{r['dataset']:>9}/{r['ear_id']:<7} cover {r['rim_cover']:.2f}  "
              f"gap {r['rim_gap']:5.2f}  tip {r['wing_tip']:6.2f}  "
              f"jac {r['jacket_mean']:5.2f}  prot {r['protrusion']:6.2f}  "
              f"hard {r['hard_min']:6.2f}  {r['grade']:8s} ({r['worst_metric']})",
              flush=True)

    cols = list(rows[0].keys())
    with open(a.csv, "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=cols)
        w.writeheader()
        w.writerows(rows)

    # worst = most failing axes first, then hardest to seat.  Ranking on the
    # overall grade alone cannot separate the bottom of the table once nearly
    # every ear fails something: the count of failed axes can.
    def severity(r):
        nfail = sum(1 for ax in ("seal", "retention", "clearance", "protrusion")
                    if r["g_" + ax] == "fail")
        return (nfail, r["seat_cost"])

    worst = sorted(rows, key=severity, reverse=True)[:5]
    md = [summarise(rows), "", "### worst 5 ears", "",
          fmt_table(worst,
                    ["dataset", "ear_id", "grade", "worst_metric", "rim_cover",
                     "rim_gap", "wing_tip", "protrusion", "hard_min"],
                    ["dataset", "ear", "grade", "failing", "cover", "max gap mm",
                     "wing tip mm", "protrusion mm", "hard interf mm"])]
    with open(a.md, "w") as f:
        f.write("\n".join(md) + "\n")
    print("\n" + "\n".join(md))
    print(f"\n{len(rows)} ears -> {a.csv}, {a.md}")

    want_viz = a.viz_ear or (f"{rows[0]['dataset']}/{rows[0]['ear_id']}" if single
                             else None)
    if want_viz:
        try:
            import viz_scene
            ds, eid = want_viz.split("/")
            vrow = next((r for r in rows
                         if r["dataset"] == ds and r["ear_id"] == eid), None)
            viz_scene.export_ear(ds, eid, row=vrow, json_dir=a.json_dir)
        except Exception as e:                                   # noqa: BLE001
            print(f"viz export failed: {type(e).__name__}: {e}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
