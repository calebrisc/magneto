#!/usr/bin/env python3
"""
align_ear.py -- find the canal aperture on a scanned ear and seat the IEM in it.

Writes, per ear, a JSON in ears/aligned/ holding the detected landmarks, the
concha frame, the seated 4x4 transform from the IEM design frame into the scan
frame, and enough diagnostics to tell a good detection from a bad one.  The
cropped ear patch is cached beside it so tryon.py never re-cuts the head mesh.


LANDMARKING
-----------
1. **Ear window.**  Anchor on the most-lateral vertex of the upper half of the
   head -- on every scan we checked, a point on the helix rim -- and cut an
   85 x 85 x 50 mm box.

2. **Hull-relative depth map.**  Cast medial rays on a 0.5 mm (x, z) grid, once
   against the ear and once against the convex hull of the same patch, and
   subtract.  The skull curving away behind the pinna lies *on* the hull, so it
   flattens to zero; the concha does not.

3. **Basin selection -- the part that actually matters.**  Threshold that map
   at 0.5 x max and label connected components.  The naive "deepest point wins"
   picks the **retroauricular sulcus** (the groove behind the ear) on a large
   fraction of ears -- it is as deep as the concha and can be just as wide.
   The discriminator that does work is **enclosure**: from a point 0.6 mm off
   the surface, cast 400 rays over the outward hemisphere and measure the
   fraction that escape 25 mm without hitting the ear.  A cavum concha is a
   bowl with walls all round it (escape ~0.10-0.25); the sulcus is a groove
   that opens onto open air behind the head (escape ~0.55-0.70).

   Enclosure is used as a **filter, not as the ranking**: basins with escape
   >= 0.42 are dropped.  Depth then selects among the survivors, but measured
   only over each basin's **bowl** -- the cells at least 0.6 x its largest
   inscribed radius from its edge.  That second gate matters as much as the
   first: the cymba concha and the crus-helicis groove are often 1-6 mm
   *deeper* than the cavum floor and score better on enclosure too, and
   without the bowl gate they drag the aperture 8-12 mm superior of the canal
   on about a third of SONICOM ears.  A groove is deep but never wide; the
   cavum is both.  The aperture is finally refined to the most-enclosed of
   twelve well-separated deep cells within 8 mm of the bowl anchor, which
   walks the last few millimetres from the middle of the cavum to the canal
   mouth.

4. **Canal probe (diagnostic only).**  A 60-degree cone of rays fired inward
   from the aperture; `canal_run` is the longest unobstructed run.  Measured
   across SONICOM and HUTUBS it is **0.5-3 mm**, because a structured-light
   scanner cannot see round the first bend -- in these datasets the canal is a
   dimple, not a tube.  We therefore do NOT
   use it to define the nozzle axis.  It stays in the JSON as a data-quality
   flag: `canal_run` under 1.5 mm means the canal was smoothed shut.

   FAILURE MODES, stated plainly: this finds the deepest sufficiently-enclosed
   basin inside the pinna outline.  It fails on (a) ears with a cymba deeper
   than the cavum, which the depth ranking would then pick;
   (b) scans with hair or an earring bridging the concha; (c) badly-cropped
   scans where the window anchor lands on a shoulder rather than the helix.
   Every ear carries `basin_escape`, `basin_area` and `canal_run`; anything
   outside the normal band is flagged and listed in the report.  `--qc-png`
   dumps a per-ear depth map with the pick marked; `--landmarks FILE` overrides
   the aperture and floor normal by hand from a JSON of
   {ear_id: {"aperture": [...], "floor_normal": [...]}}.


PLACING THE IEM
---------------
The concha frame is: **+Z** = outward normal of a plane fitted to the ear
surface within 10 mm of the aperture (the concha floor); **+Y** = world up,
orthogonalised against it; **+X** = Y x Z, which comes out anterior for a right
ear.  That matches the IEM design frame's stated axes (+X toward the canal,
+Y toward the antihelix, +Z toward the faceplate) -- except for one thing the
README already admits: the design treats the nozzle axis and the faceplate
normal as orthogonal, and in a real ear the canal axis and the concha-floor
normal are 30-60 degrees apart, not 90.  So there is no single correct mapping
to hand the optimiser.

We therefore **multi-start**: five nozzle rakes (0-80 degrees of rotation about
the design +Y, sweeping the nozzle from "anterior, lying in the concha" round
to "medial, pointing down the canal") crossed with three rolls about +X, with
the translation seeded each time so the Ø19 skirt rim centre lands on the
aperture.  All fifteen are scored, the best four are refined with a bounded
Powell search (+-10 mm, +-40 degrees), and the winner is kept.  The seating cost
is a physical one, and every term is there because leaving it out produced a
pose a wearer would never adopt:

    c_rim    the rim wants to lie in a +-0.75 mm contact band
    c_pen    the rigid Ti parts want no more than 0.3 mm of penetration,
             graded on the WORST point (a mean lets a corner bury itself)
    c_wing   the wing tip wants about 1 mm of press into the antihelix
    c_jac    the jacket wants to be near the concha floor
    c_soft   the silicone skirt must not be BURIED.  c_rim scores |distance|,
             so a rim 0.5 mm inside the flesh scores as well as one touching
             it; with nothing opposing that, the optimiser bought rim coverage
             by driving the whole Ø19 skirt down the canal -- 13 to 20 mm
             medial of the tragus plane on the synthetic corners, which is
             inside the skull.  A Ø19 skirt cannot enter an aperture that is at
             most 18 x 14 mm; it seals *on* the funnel entrance.
    c_prot   the other half of the same failure.  The assembly is 32.6 mm long
             along its nozzle axis, so burying the rim throws the faceplate
             just as far the other way.  The reported pose has to be one a
             wearer would accept, not one standing 14 mm proud of the tragus.

That local search is what a wearer does when they wiggle the thing in, and it
is what makes the numbers mean something in spite of the frame ambiguity.  The
winning rake is reported per ear as `nozzle_rake_deg`.

STABILITY.  The cost surface is rough -- it is built on a sampled signed-distance
field -- and Powell finds a local minimum, so the per-ear millimetres move a
little between runs.  The sampling is seeded (`--field-seed`) so a given run
reproduces exactly; vary the seed to measure how much any conclusion depends on
it.  `docs/TRYON_REPORT.md` quotes that spread.

Usage:
    python align_ear.py --dataset sonicom
    python align_ear.py --dataset hutubs --limit 20 --qc-png
    python align_ear.py --dataset synthetic --qc-png
"""

from __future__ import annotations

import argparse
import concurrent.futures as cf
import glob
import json
import os
import sys
import time

import numpy as np
import trimesh
from scipy import ndimage, optimize

from earfit import (ALIGNED, DATASETS, EARS, EarField, SKIRT_RIM_X, crop,
                    depth_map, ear_window, iem_points, load_head, rt_matrix,
                    transform)

RAKES = (0.0, 20.0, 40.0, 60.0, 80.0)
ROLLS = (-25.0, 0.0, 25.0)


# --------------------------------------------------------------------------- #
# geometry helpers
# --------------------------------------------------------------------------- #

def cone_dirs(axis, half_deg, n):
    """Roughly-uniform directions inside a cone about `axis` (Fibonacci cap)."""
    axis = np.asarray(axis, float)
    axis = axis / np.linalg.norm(axis)
    a = np.array([0.0, 0.0, 1.0])
    if abs(axis @ a) > 0.9:
        a = np.array([1.0, 0.0, 0.0])
    u = np.cross(axis, a); u /= np.linalg.norm(u)
    v = np.cross(axis, u)
    cmin = np.cos(np.radians(half_deg))
    g = (1 + 5 ** 0.5) / 2
    i = np.arange(n)
    cz = 1 - (1 - cmin) * (i + 0.5) / n
    ph = 2 * np.pi * i / g
    s = np.sqrt(np.clip(1 - cz ** 2, 0, 1))
    return (cz[:, None] * axis + (s * np.cos(ph))[:, None] * u
            + (s * np.sin(ph))[:, None] * v)


def escape_fraction(patch, p, outward, n=400, half=85.0, cap=25.0):
    """Fraction of an outward hemisphere of rays that clears the ear.

    Low = the point sits in a bowl (concha).  High = it sits in a groove that
    opens onto free space (retroauricular sulcus)."""
    start = p + 0.6 * outward
    dirs = cone_dirs(outward, half, n)
    org = np.tile(start, (len(dirs), 1))
    loc, ir, _ = patch.ray.intersects_location(org, dirs, multiple_hits=False)
    if len(loc) == 0:
        return 1.0
    run = np.linalg.norm(loc - start, axis=1)
    hit = np.zeros(len(dirs), bool)
    hit[ir[run < cap]] = True
    return float(1.0 - hit.mean())


def canal_probe(head, p, medial, half=60.0, n=260, cap=30.0):
    start = p - 1.0 * medial
    dirs = cone_dirs(medial, half, n)
    org = np.tile(start, (len(dirs), 1))
    loc, ir, _ = head.ray.intersects_location(org, dirs, multiple_hits=False)
    if len(loc) == 0:
        return 0.0, medial
    run = np.linalg.norm(loc - start, axis=1)
    keep = run < cap
    if not keep.any():
        return 0.0, medial
    run, ir = run[keep], ir[keep]
    j = int(np.argmax(run))
    return float(run[j] - 1.0), dirs[ir[j]]


def surface_point(dm, ix, iz):
    return np.array([dm["xs"][ix], dm["surf"][ix, iz], dm["zs"][iz]])


# --------------------------------------------------------------------------- #
# landmarks
# --------------------------------------------------------------------------- #

def find_aperture(head, patch, dm, win):
    rel = dm["rel"]
    ok = np.isfinite(rel)
    # 99th percentile, not the max: a narrow canal bore can be several mm
    # deeper than the concha floor and would otherwise set the threshold so
    # high that the bowl itself falls below it
    rmax = float(np.nanpercentile(rel, 99.0))
    step = float(dm["xs"][1] - dm["xs"][0])
    lat = win["lat"]
    outward = np.array([0.0, lat, 0.0])
    medial = -outward

    lab, nlab = ndimage.label(ok & (rel > 0.5 * rmax))
    cands = []
    for k in range(1, nlab + 1):
        mask = lab == k
        area = int(mask.sum()) * step * step
        if area < 18.0:
            continue
        # The BOWL gate.  A basin's depth is scored only over the part of it
        # that is wide enough to be a bowl -- cells at least 0.6 x the basin's
        # largest inscribed radius away from its edge.  Without this, the
        # cymba concha and the crus-helicis groove win on raw depth (they can
        # be 1-6 mm deeper than the cavum floor) and drag the aperture 8-12 mm
        # superior of the canal, which happened on about a third of SONICOM
        # ears.  A groove is deep but never wide; the cavum is both.
        dt = ndimage.distance_transform_edt(mask) * step
        wide = mask & (dt >= 0.6 * dt.max())
        sub = np.where(wide, rel, -np.inf)
        ij = np.unravel_index(int(np.argmax(sub)), rel.shape)
        p = surface_point(dm, *ij)
        cands.append(dict(k=k, area=area, depth=float(rel[ij]), p=p, ij=ij,
                          inscribed=float(dt.max()),
                          esc=escape_fraction(patch, p, outward)))
    if not cands:
        raise RuntimeError("no deep basin inside the ear window")

    # Enclosure is a FILTER, bowl-gated depth is the SELECTOR.  Enclosure alone
    # would pick the narrowest groove; depth alone would pick the sulcus.
    pool = [c for c in cands if c["esc"] < 0.42] or cands
    best = max(pool, key=lambda c: c["depth"])

    # refine: most-enclosed of twelve well-separated deep cells within 8 mm of
    # the bowl anchor.  The bowl gate lands us near the centre of the cavum;
    # this walks the last few mm to the canal mouth, which is the most enclosed
    # thing in the bowl.
    ax, az = dm["xs"][best["ij"][0]], dm["zs"][best["ij"][1]]
    near = ((np.abs(dm["xs"] - ax)[:, None] ** 2
             + np.abs(dm["zs"] - az)[None, :] ** 2) < 8.0 ** 2)
    sub = np.where((lab == best["k"]) & near, rel, -np.inf)
    flat = np.argsort(sub, axis=None)[::-1][:200]
    cells = np.array(np.unravel_index(flat, rel.shape)).T
    pick, seen = [], []
    for ix, iz in cells:
        q = np.array([dm["xs"][ix], dm["zs"][iz]])
        if all(np.linalg.norm(q - s) > 2.5 for s in seen):
            seen.append(q); pick.append((ix, iz))
        if len(pick) >= 12:
            break
    scored = [(escape_fraction(patch, surface_point(dm, ix, iz), outward), ix, iz)
              for ix, iz in pick]
    esc, ix, iz = min(scored)
    ap = surface_point(dm, ix, iz)
    run, _ = canal_probe(head, ap, medial)

    return dict(aperture=ap, basin_escape=esc, basin_area=best["area"],
                basin_depth=best["depth"], basin_inscribed=best["inscribed"],
                n_basins=len(cands), canal_run=run)


def floor_normal(patch, ap, outward, radius=10.0):
    """Outward normal of the concha floor: total-least-squares plane through the
    ear surface within `radius` of the aperture -- the annulus the Ø19 skirt rim
    has to seal against."""
    v = patch.vertices
    sel = v[np.linalg.norm(v - ap, axis=1) < radius]
    if len(sel) < 50:
        return outward
    c = sel - sel.mean(0)
    n = np.linalg.svd(c, full_matrices=False)[2][-1]
    if n @ outward < 0:
        n = -n
    return n / np.linalg.norm(n)


def find_tragus(dm, ap, win):
    """Most-lateral ridge point just anterior of the aperture -- the anterior
    wall of the concha, and the plane the faceplate protrusion is measured
    against."""
    lat = win["lat"]
    xs, zs, surf = dm["xs"], dm["zs"], dm["surf"]
    mx = (xs > ap[0] + 2.0) & (xs < ap[0] + 18.0)
    mz = (zs > ap[2] - 9.0) & (zs < ap[2] + 9.0)
    if not (mx.any() and mz.any()):
        return None
    sub = surf[np.ix_(mx, mz)]
    lateral = sub * lat
    lateral[~np.isfinite(lateral)] = -np.inf
    if not np.isfinite(lateral).any():
        return None
    ij = np.unravel_index(int(np.argmax(lateral)), sub.shape)
    return np.array([xs[mx][ij[0]], sub[ij], zs[mz][ij[1]]])


def concha_frame(ap, n_out, lat):
    """4x4 with +X anterior, +Y superior, +Z = concha outward normal, origin
    placed so the skirt rim centre lands on the aperture."""
    ez = n_out / np.linalg.norm(n_out)
    up = np.array([0.0, 0.0, 1.0])
    ey = up - (up @ ez) * ez
    ey /= np.linalg.norm(ey)
    ex = np.cross(ey, ez)
    ex /= np.linalg.norm(ex)
    if ex[0] < 0:                       # +X must point anterior
        ex, ey = -ex, -ey
    M = np.eye(4)
    M[:3, 0], M[:3, 1], M[:3, 2] = ex, ey, ez
    M[:3, 3] = ap - SKIRT_RIM_X * ex
    return M


# --------------------------------------------------------------------------- #
# seating
# --------------------------------------------------------------------------- #

def relu(a):
    return np.maximum(a, 0.0)


def seating_cost(M, P, field, ctx=None):
    if not np.all(np.isfinite(M)):
        return 1e6
    rim = field.query(transform(P["rim"], M))
    rig = field.query(transform(P["shell"], M))
    tip = field.query(transform(P["wing_tip"], M))
    jac = field.query(transform(P["jacket"], M))
    sof = field.query(transform(P["soft"], M))
    c_rim = np.mean(relu(np.abs(rim) - 0.75) ** 2)
    # penetration is graded on the WORST point, not the mean: averaging over
    # 800 samples lets the optimiser bury a corner of the shell 4 mm into the
    # concha wall for the sake of a percent of rim coverage.  Flesh does not
    # work that way -- the deepest single overlap is what stops the fit.
    c_pen = (3.0 * relu(-float(rig.min()) - 0.5) ** 2
             + 2.0 * np.mean(relu(-rig - 0.30) ** 2))
    t = float(np.median(tip))
    c_wing = 0.5 * (relu(t + 0.3) ** 2 + relu(-t - 2.0) ** 2)
    c_jac = 0.15 * np.mean(relu(jac - 2.0) ** 2)

    # SOFT-BODY BURIAL.  c_rim scores |distance|, so a rim sitting 0.5 mm *inside*
    # the flesh scores as well as one touching it.  With nothing else opposing it
    # the optimiser buys rim coverage by driving the whole Ø19 skirt down the ear
    # canal -- on the synthetic corners it put the rim 13-20 mm medial of the
    # tragus plane, which is inside the skull.  A Ø19 mm skirt cannot enter a
    # canal whose aperture is at most 18 x 14 mm; the skirt seals *on* the funnel
    # entrance.  Silicone squashes, so allow ~0.75 mm mean and 1.5 mm worst-point.
    c_soft = (2.0 * relu(-float(sof.min()) - 1.5) ** 2
              + 1.0 * np.mean(relu(-sof - 0.75) ** 2))

    # PROTRUSION.  The other half of the same failure: the assembly is 32.6 mm
    # long along its nozzle axis, so burying the rim throws the faceplate just as
    # far the other way.  Nobody wears an IEM standing 14 mm proud of the tragus,
    # so the pose the optimiser reports has to be one a wearer would accept.
    c_prot = 0.0
    if ctx is not None and ctx.get("tragus") is not None:
        fp = transform(P["faceplate"], M)
        prot = float(np.max(np.einsum("ij,j->i", fp - ctx["tragus"], ctx["normal"])))
        c_prot = 0.08 * relu(prot - 2.0) ** 2

    c = c_rim + c_pen + c_wing + c_jac + c_soft + c_prot
    return float(c) if np.isfinite(c) else 1e6


def start_pose(frame, ap, rake_deg, roll_deg):
    """Rotate the design frame by `rake` about +Y (nozzle anterior -> medial)
    and `roll` about +X, then re-seat so the rim centre stays on the aperture."""
    R = trimesh.transformations.euler_matrix(np.radians(roll_deg),
                                             np.radians(rake_deg), 0.0, "rxyz")
    M = frame.copy()
    M[:3, 3] = 0.0
    M = M @ R
    M[:3, 3] = ap - SKIRT_RIM_X * M[:3, 0]
    return M


def seat(frame, ap, P, field, refine=4, ctx=None):
    starts = []
    for rake in RAKES:
        for roll in ROLLS:
            M = start_pose(frame, ap, rake, roll)
            starts.append((seating_cost(M, P, field, ctx), rake, roll, M))
    starts.sort(key=lambda s: s[0])

    best = None
    for c0, rake, roll, M0 in starts[:refine]:
        def f(z, M0=M0):
            return seating_cost(rt_matrix(z[:3], z[3:], M0), P, field, ctx)
        res = optimize.minimize(f, np.zeros(6), method="Powell",
                                bounds=[(-10, 10)] * 3 + [(-40, 40)] * 3,
                                options=dict(maxfev=1400, xtol=0.05, ftol=1e-3))
        if best is None or res.fun < best[0]:
            best = (float(res.fun), rake, roll, rt_matrix(res.x[:3], res.x[3:], M0),
                    res.x.tolist(), float(c0))
    return best


# --------------------------------------------------------------------------- #
# driver
# --------------------------------------------------------------------------- #

def qc_png(path, dm, rec, tragus):
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    xs, zs = dm["xs"], dm["zs"]
    ap = rec["aperture"]
    fig, ax = plt.subplots(figsize=(5.4, 5.4))
    ax.imshow(dm["rel"].T, origin="lower", cmap="magma",
              extent=[xs[0], xs[-1], zs[0], zs[-1]])
    ax.plot(ap[0], ap[2], "co", ms=10, mfc="none", mew=2, label="aperture")
    if tragus is not None:
        ax.plot(tragus[0], tragus[2], "g^", ms=8, label="tragus")
    ax.set_xlabel("x anterior (mm)"); ax.set_ylabel("z superior (mm)")
    ax.set_title(f"{rec['dataset']} {rec['ear_id']}  esc={rec['basin_escape']:.2f}"
                 f"  rake={rec['nozzle_rake_deg']:.0f}deg  cost={rec['seat_cost']:.2f}",
                 fontsize=9)
    ax.legend(loc="lower left", fontsize=7)
    fig.tight_layout(); fig.savefig(path, dpi=100); plt.close(fig)


def process(path, ds, side, P, want_png, manual):
    cfg = DATASETS[ds]
    eid = cfg["id_from"](path)
    head = load_head(path, cfg["scale"])
    win = ear_window(head, side)
    patch = crop(head, win)
    dm = depth_map(patch, win)
    lm = find_aperture(head, patch, dm, win)
    outward = np.array([0.0, win["lat"], 0.0])
    n_out = floor_normal(patch, lm["aperture"], outward)
    if manual and eid in manual:
        lm["aperture"] = np.array(manual[eid]["aperture"], float)
        n_out = np.array(manual[eid].get("floor_normal", n_out), float)
        lm["manual"] = True

    tragus = find_tragus(dm, lm["aperture"], win)
    frame = concha_frame(lm["aperture"], n_out, win["lat"])
    field = EarField(patch)
    ctx = dict(tragus=None if tragus is None else np.asarray(tragus, float),
               normal=n_out)
    cost, rake, roll, M, z, cost0 = seat(frame, lm["aperture"], P, field, ctx=ctx)

    os.makedirs(ALIGNED, exist_ok=True)
    stem = os.path.join(ALIGNED, f"{ds}_{eid}_{side}")
    patch.export(stem + "_patch.ply")
    rec = dict(
        dataset=ds, ear_id=eid, side=side,
        source=os.path.relpath(path, EARS),
        aperture=np.asarray(lm["aperture"]).tolist(),
        floor_normal=n_out.tolist(),
        outward=outward.tolist(),
        tragus=None if tragus is None else tragus.tolist(),
        basin_escape=lm["basin_escape"], basin_area=lm["basin_area"],
        basin_depth=lm["basin_depth"], basin_inscribed=lm["basin_inscribed"],
        n_basins=lm["n_basins"],
        canal_run=lm["canal_run"], manual=bool(lm.get("manual")),
        weak=bool(lm["basin_escape"] > 0.42 or lm["canal_run"] < 1.5),
        concha_frame=frame.tolist(), transform=M.tolist(),
        nozzle_rake_deg=rake, roll_deg=roll,
        seat_delta=z, seat_cost=cost, seat_cost_start=cost0,
        patch=os.path.basename(stem + "_patch.ply"),
    )
    with open(stem + ".json", "w") as f:
        json.dump(rec, f, indent=1)
    if want_png:
        qc_png(stem + "_qc.png", dm, rec, tragus)
    return rec


def reseat(json_path, P, field_seed=0):
    """Redo only the seating search for an ear that is already landmarked.

    Landmarking (crop, depth map, escape fractions) is the expensive half and
    does not change when the seating cost changes, so tuning the cost does not
    have to re-cut every head mesh.  Reads the cached patch beside the JSON and
    rewrites the transform and seating diagnostics in place.
    """
    rec = json.load(open(json_path))
    patch = trimesh.load(os.path.join(ALIGNED, rec["patch"]), force="mesh")
    ap = np.array(rec["aperture"], float)
    frame = np.array(rec["concha_frame"], float)
    n_out = np.array(rec["floor_normal"], float)
    trg = rec.get("tragus")
    ctx = dict(tragus=None if trg is None else np.array(trg, float), normal=n_out)
    cost, rake, roll, M, z, cost0 = seat(frame, ap, P,
                                         EarField(patch, seed=field_seed), ctx=ctx)
    rec.update(transform=M.tolist(), nozzle_rake_deg=rake, roll_deg=roll,
               seat_delta=z, seat_cost=cost, seat_cost_start=cost0)
    with open(json_path, "w") as f:
        json.dump(rec, f, indent=1)
    return rec


_P = None


def _reseat_worker(json_path, field_seed=0):
    global _P
    if _P is None:
        _P = iem_points()
    return reseat(json_path, _P, field_seed)


def _worker(path, ds, side, want_png, manual):
    """Process pool entry point; the IEM point sets are built once per worker."""
    global _P
    if _P is None:
        _P = iem_points()
    return process(path, ds, side, _P, want_png, manual)


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--dataset", choices=list(DATASETS), default=None)
    ap.add_argument("--field-seed", type=int, default=0,
                    help="seed for the ear-surface point sampling; vary it to "
                         "measure how much a result depends on the sampling")
    ap.add_argument("--reseat", action="store_true",
                    help="re-run only the seating search over ears already in "
                         "ears/aligned/, reusing their cached patch and "
                         "landmarks (for tuning the seating cost)")
    ap.add_argument("--side", default="right", choices=("right", "left"))
    ap.add_argument("--limit", type=int, default=0)
    ap.add_argument("--qc-png", action="store_true")
    ap.add_argument("--landmarks", default=None,
                    help='manual overrides: {ear_id: {"aperture": [x,y,z], '
                         '"floor_normal": [x,y,z]}}')
    ap.add_argument("--jobs", type=int, default=max(1, (os.cpu_count() or 2) - 1))
    ap.add_argument("--skip-existing", action="store_true",
                    help="leave ears that already have a JSON in ears/aligned/ "
                         "alone; makes the run resumable after an interrupt")
    a = ap.parse_args()

    def report(r):
        print(f"{r['dataset']:>9}/{r['ear_id']:<7} esc {r['basin_escape']:.2f}  "
              f"area {r['basin_area']:6.1f}  run {r['canal_run']:4.1f}  "
              f"rake {r['nozzle_rake_deg']:3.0f}  cost {r['seat_cost_start']:6.2f}"
              f" -> {r['seat_cost']:6.3f}{'  WEAK' if r['weak'] else ''}", flush=True)

    if a.reseat:
        pat = f"{a.dataset}_*.json" if a.dataset else "*.json"
        js = sorted(glob.glob(os.path.join(ALIGNED, pat)))
        if not js:
            sys.exit(f"nothing to reseat in {ALIGNED}")
        t0, ok, bad = time.time(), 0, 0
        with cf.ProcessPoolExecutor(max_workers=a.jobs) as ex:
            futs = {ex.submit(_reseat_worker, p, a.field_seed): p for p in js}
            for fut in cf.as_completed(futs):
                try:
                    report(fut.result()); ok += 1
                except Exception as e:                           # noqa: BLE001
                    bad += 1
                    print(f"{os.path.basename(futs[fut])}  FAILED: {e}", flush=True)
        print(f"\n{ok} reseated, {bad} failed, {time.time()-t0:.0f}s")
        return 0

    if not a.dataset:
        ap.error("--dataset is required unless --reseat is given")
    root = os.path.join(EARS, a.dataset)
    files = sorted(glob.glob(os.path.join(root, DATASETS[a.dataset]["glob"])))
    if a.skip_existing:
        idf = DATASETS[a.dataset]["id_from"]
        files = [p for p in files if not os.path.exists(
            os.path.join(ALIGNED, f"{a.dataset}_{idf(p)}_{a.side}.json"))]
    if a.limit:
        files = files[:a.limit]
    if not files:
        if a.skip_existing and glob.glob(os.path.join(root, DATASETS[a.dataset]["glob"])):
            print(f"{a.dataset}: every mesh already aligned in {ALIGNED}, nothing to do")
            return 0
        sys.exit(f"no meshes under {root} -- run fetch_ears.py / make_synthetic_ear.py")

    manual = json.load(open(a.landmarks)) if a.landmarks else None
    ok = bad = 0
    t0 = time.time()

    if a.jobs > 1:
        with cf.ProcessPoolExecutor(max_workers=a.jobs) as ex:
            futs = {ex.submit(_worker, p, a.dataset, a.side, a.qc_png, manual): p
                    for p in files}
            for fut in cf.as_completed(futs):
                try:
                    report(fut.result()); ok += 1
                except Exception as e:                           # noqa: BLE001
                    bad += 1
                    print(f"{a.dataset}/{os.path.basename(futs[fut])}  FAILED: {e}",
                          flush=True)
    else:
        P = iem_points()
        for p in files:
            try:
                report(process(p, a.dataset, a.side, P, a.qc_png, manual)); ok += 1
            except Exception as e:                               # noqa: BLE001
                bad += 1
                print(f"{a.dataset}/{os.path.basename(p)}  FAILED: {e}", flush=True)

    print(f"\n{ok} aligned, {bad} failed, {time.time()-t0:.0f}s -> {ALIGNED}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
