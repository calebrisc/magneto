#!/usr/bin/env python3
"""
size_bands.py -- S/M/L clamp-pad extensions instead of one universal stroke.

Question (Cale, 2026-09-04): if the cymba clamp pad comes in three lengths,
one body, does the lock close on the ears it could not close on universally?

MODEL  (same geometry as hook_config.py, cymba site only -- the antihelix
stays out of v1 per CONCHA_VARIANCE_102)
  reach   along-aim distance from the jacket outer face to the centroid of the
          undercut patch under the cymba lip / crus overhang.  This is the
          length the pad extension has to span -- the sizing variable.
  wrap    extent of that undercut along pull-out; how much lip there is to hook.
  seated  L_c <= reach <= L_c + S  for the size's (compacted, stroke).

SIZING RULE
  universal   one window, hook_config defaults (4.5 + 6.0)
  banded      sort ears by reach, split into 3 contiguous bands minimising the
              widest band; each band gets L_c = min - MARGIN, S = span + 2 MARGIN,
              and the stack must physically fit: (n_gap+1) magnets + pad + plate
              <= L_c with 1.5 mm usable stroke per repelling gap.

STABILITY  stability.stability_check, all pull directions incl. a straight
  outward yank, pad patch on the undercut normals, mu = 0.6 (dry skin-silicone,
  Zhang & Mak 1999), cable tug 0.5 N and 0.2 N (over-ear).  Run for every ear
  under the universal window and under its own band.

SIZE PICKING  reach is what the part must span; a wearer cannot measure it.
  Report how well ear scale (basin inscribed radius) and crus distance predict
  the band, i.e. whether "pick S/M/L by ear size" actually lands people in
  the right band.

Usage:
    python size_bands.py --reach            # phase 1, ~40 s, caches
    python size_bands.py --stability        # phase 2, slow, uses the cache
    python size_bands.py --report           # tables from the caches
"""
from __future__ import annotations

import argparse, glob, json, os, sys, time
import numpy as np
import trimesh

import stability as stab
from earfit import ALIGNED, HERE, NOZZLE_AXIS, EarField, iem_points, transform
from hook_config import hook_at, gen_bits

AIM = np.array([0.20, 0.98, -0.10]); AIM /= np.linalg.norm(AIM)   # cymba, post-bias
EXCLUDE = {"P0023"}                     # scanner hole-fill concha (9/1 audit)
MARGIN = 1.0                            # mm each end of a band
GAP_STROKE, MAG_T, PAD_T, PLATE_T = 1.5, 1.0, 1.0, 0.8
UNIVERSAL = (4.5, 6.0)                  # hook_config defaults (compacted, stroke)
MU = 0.6
TUGS = (0.5, 0.2)
CACHE_R = os.path.join(ALIGNED, "size_bands_reach.json")
CACHE_S = os.path.join(ALIGNED, "size_bands_stability.json")
OUT_MD = os.path.join(HERE, "..", "..", "docs", "CLAMP_SIZE_BANDS.md")


def ears():
    js = sorted(glob.glob(os.path.join(ALIGNED, "*.json")))
    for p in js:
        try:
            rec = json.load(open(p))
        except ValueError:
            continue                      # a cache mid-write, not an ear
        if not isinstance(rec, dict) or "patch" not in rec:
            continue
        if rec.get("dataset", "").startswith("synth") or rec["ear_id"].startswith("synthetic"):
            continue
        if rec["ear_id"] in EXCLUDE:
            continue
        yield rec


def phase_reach():
    cc, cr, off, _ = gen_bits()
    gav = {e["id"]: e for e in json.load(open(os.path.join(ALIGNED, "gyro_arm_variance.json")))}
    out = []
    for rec in ears():
        patch = trimesh.load(os.path.join(ALIGNED, rec["patch"]), force="mesh")
        field = EarField(patch, seed=0)
        h = hook_at(rec, patch, field, AIM, cc, cr, off, 0.0, 1e3)
        g = gav.get(rec["ear_id"], {})
        crus = g.get("crus")
        out.append(dict(ear=rec["ear_id"], dataset=rec.get("dataset"),
                        reach=(None if not np.isfinite(h["reach"]) else h["reach"]),
                        wrap=h["wrap"], n_pts=h["n_pts"],
                        patch_design=(None if h["patch"] is None else h["patch_design"].tolist()),
                        scale=g.get("scale"),
                        crus=(None if crus is None else float(np.linalg.norm(crus)))))
        print(f"{rec['ear_id']:<10} reach {out[-1]['reach'] if out[-1]['reach'] is None else round(out[-1]['reach'],2)}  wrap {h['wrap']:.2f}")
    json.dump(out, open(CACHE_R, "w"), indent=1)
    print(f"\n{len(out)} ears -> {CACHE_R}")


def bands_from(reach):
    """3 contiguous bands over sorted reach minimising the widest span."""
    r = np.sort(reach); n = len(r); best = None
    for i in range(1, n - 1):
        for j in range(i + 1, n):
            spans = (r[i-1]-r[0], r[j-1]-r[i], r[-1]-r[j])
            key = (max(spans), sum(spans))
            if best is None or key < best[0]:
                best = (key, (r[0], r[i-1]), (r[i], r[j-1]), (r[j], r[-1]))
    return best[1:]


def sizing(lo, hi):
    L_c = lo - MARGIN; S = (hi - lo) + 2 * MARGIN
    n_gap = int(np.ceil(S / GAP_STROKE))
    compact = (n_gap + 1) * MAG_T + PAD_T + PLATE_T
    return dict(compacted=L_c, stroke=S, n_gap=n_gap, stack_len=compact, fits=compact <= L_c)


def configs(R):
    reach = np.array([e["reach"] for e in R if e["reach"] is not None])
    b = bands_from(reach)
    cfg = {"universal": dict(window=UNIVERSAL, members=[e["ear"] for e in R])}
    for name, (lo, hi) in zip("SML", b):
        sz = sizing(lo, hi)
        mem = [e["ear"] for e in R if e["reach"] is not None and lo - 1e-9 <= e["reach"] <= hi + 1e-9]
        cfg[name] = dict(window=(sz["compacted"], sz["stroke"]), band=(lo, hi), sizing=sz, members=mem)
    return cfg


def seated(reach, window):
    return reach is not None and window[0] <= reach <= window[0] + window[1]


def phase_stability():
    R = json.load(open(CACHE_R)); cfg = configs(R)
    P = iem_points()
    done = json.load(open(CACHE_S)) if os.path.exists(CACHE_S) else {}
    byid = {e["ear"]: e for e in R}
    recs = {rec["ear_id"]: rec for rec in ears()}
    t0 = time.time()
    for e in R:
        eid = e["ear"]; rec = recs[eid]
        patch = trimesh.load(os.path.join(ALIGNED, rec["patch"]), force="mesh")
        field = EarField(patch, seed=0)
        M = np.array(rec["transform"], float)
        cp = transform(P["_cable_exit"][None, :], M)[0]
        com = transform(P["_com"][None, :], M)[0]
        band = next((k for k in "SML" if eid in cfg[k]["members"]), None)
        for label in ("universal", band):
            if label is None: continue
            key = f"{eid}|{label}"
            if key in done: continue
            w = cfg[label]["window"]
            st = seated(e["reach"], w)
            Pl = dict(P)
            if st and e["patch_design"] is not None:
                Pl["plunger"] = np.array(e["patch_design"]); Pl["_plunger"] = [dict(name="cymba")]
            else:
                Pl["plunger"] = np.zeros((0, 3)); Pl["_plunger"] = []
            row = dict(ear=eid, config=label, seated=bool(st))
            for tug in TUGS:
                x = stab.stability_check(rec, Pl, field, transform, cable_point=cp, com=com,
                                         mu=MU, cable_tug=tug, cable_mode="sphere")
                row[f"tug{tug}"] = round(float(x["margin"]), 3)
                row[f"pullout{tug}"] = round(float(x.get("pullout_capacity", 0)), 3)
                x = stab.stability_check(rec, Pl, field, transform, cable_point=cp, com=com,
                                         mu=MU, cable_tug=tug, cable_mode="cone")
                row[f"cone{tug}"] = round(float(x["margin"]), 3)
            done[key] = row
            json.dump(done, open(CACHE_S, "w"), indent=1)
            print(f"{eid:<10} {label:<9} seated {str(st):<5} "
                  f"{row['tug0.5']:.2f}x {row['tug0.2']:.2f}x   [{time.time()-t0:.0f}s]", flush=True)



VERDICT = """## Verdict (2026-09-04)

1. **Sizing fixes reach, not retention.** Three pad lengths seat the cymba pad on 96/102 ears
   vs 51/102 with one universal 6 mm stroke. But the along-aim reach spans 25 mm and has **no
   correlation with ear size** (r = -0.04 vs basin radius, -0.19 vs crus distance), so a wearer
   cannot pick S/M/L by how big their ear is. Sizes would be picked by trying the three pads.
2. **A single pad does not hold, at any size, in either load case.** Cone (cable down/back) and
   sphere (any yank) both fail on ~97% of ears at mu 0.6. Pad area (60-200 mm2) and skirt
   push-out (0.31-0.15 N) barely move the margin. The binding constraint is **moment balance**: one
   interlock point in the cymba is a lever, and a pull on the nozzle rotates the body about it.
   The two-pad build (cymba + antihelix) reaches 0.57x on the same ear where one pad gives 0.00x.
3. **The 9/1 clamp is not what this harness tested.** The clamp's fixed jaw is the body edge tucked
   under the crus; that is the second reaction that balances the lever. In the cached seated poses
   (plunger-era optimiser) the jacket touches the concha floor at 0-4 points and is not under the
   crus, so the model cannot credit the jaw. Before anything else: re-seat by construction (body
   edge under the crus, 8/31 item), then rerun `size_bands.py --stability` and `crus_bands.py`.
4. **Decision this does support:** one body, pad extension in three lengths chosen by fit, stroke
   ~7-8 mm per size (5-6 gaps), sunk pocket if the S size is to fit. Crus-travel bands from
   `crus_bands.py`: S 4-9, M 10-16, L 17-22 mm (M and L stacks fit; S needs a recess).
5. **Open question the FEA cannot answer:** whether the crus jaw + one pad gives >= 1x. That is the
   first thing the Form 3 fit shell should tell us on Cale's ear.
"""


def report():
    R = json.load(open(CACHE_R)); cfg = configs(R)
    S = json.load(open(CACHE_S)) if os.path.exists(CACHE_S) else {}
    byid = {e["ear"]: e for e in R}
    reach = np.array([e["reach"] for e in R if e["reach"] is not None])
    L = []
    L.append("# Clamp pad S/M/L bands vs one universal stroke (2026-09-04)\n")
    L.append(VERDICT)
    L.append(f"{len(R)} real ears (P0023 excluded), cymba site only, aim {np.round(AIM,2).tolist()}. "
             "`reach` = along-aim distance from the jacket face to the undercut patch under the cymba lip; "
             "it is the length the pad extension must span. Script: `cad/iem/size_bands.py`.\n")
    L.append(f"reach over {len(reach)} ears with an undercut: min {reach.min():.2f}  p10 {np.percentile(reach,10):.2f}  "
             f"median {np.median(reach):.2f}  p90 {np.percentile(reach,90):.2f}  max {reach.max():.2f} mm  "
             f"(span {reach.max()-reach.min():.2f}); no undercut found on {len(R)-len(reach)} ears.\n")
    L.append("## Windows\n")
    L.append("| config | band (reach, mm) | compacted | stroke | gaps | stack length | fits in compacted? | ears | seated |")
    L.append("|---|---|---|---|---|---|---|---|---|")
    for k in ("universal", "S", "M", "L"):
        c = cfg[k]; w = c["window"]
        mem = [byid[m] for m in c["members"]]
        n_seat = sum(seated(m["reach"], w) for m in mem)
        if k == "universal":
            sz = sizing(w[0] + MARGIN, w[0] + w[1] - MARGIN)
            L.append(f"| universal | all | {w[0]:.2f} | {w[1]:.2f} | {sz['n_gap']} | {sz['stack_len']:.2f} | "
                     f"{'yes' if sz['fits'] else 'NO'} | {len(mem)} | {n_seat} ({100*n_seat/len(mem):.0f}%) |")
        else:
            sz = c["sizing"]; lo, hi = c["band"]
            L.append(f"| {k} | {lo:.2f}–{hi:.2f} | {sz['compacted']:.2f} | {sz['stroke']:.2f} | {sz['n_gap']} | "
                     f"{sz['stack_len']:.2f} | {'yes' if sz['fits'] else 'NO'} | {len(mem)} | {n_seat} ({100*n_seat/max(len(mem),1):.0f}%) |")
    tot = sum(sum(seated(byid[m]["reach"], cfg[k]["window"]) for m in cfg[k]["members"]) for k in "SML")
    L.append(f"\nBanded: {tot}/{len(R)} ears seated with the right size vs "
             f"{sum(seated(e['reach'], UNIVERSAL) for e in R)}/{len(R)} universal.\n")
    # stability
    if S:
        L.append("## Stability, mu = 0.6\n")
        L.append("Two load cases. **cone** = cable hanging down and back (the 8/31 matrix convention). "
                 "**sphere** = every pull direction including a straight outward yank along the nozzle axis, "
                 "the case a hook has to survive. Margin >= 1x passes.\n")
        L.append("| config | ears | seated | cone pass @0.5 N | cone pass @0.2 N | cone median 0.5 N | sphere pass @0.5 N | sphere pass @0.2 N | median straight pull-out capacity |")
        L.append("|---|---|---|---|---|---|---|---|---|")
        for k in ("universal", "S", "M", "L", "banded"):
            rows = [v for v in S.values() if (v["config"] == k if k != "banded" else v["config"] in "SML")]
            if not rows: continue
            a = np.array([r["tug0.5"] for r in rows]); b = np.array([r["tug0.2"] for r in rows])
            c = np.array([r.get("cone0.5", np.nan) for r in rows]); d = np.array([r.get("cone0.2", np.nan) for r in rows])
            po = np.array([r.get("pullout0.5", np.nan) for r in rows]); ns = sum(r["seated"] for r in rows)
            L.append(f"| {k} | {len(rows)} | {ns} | {(c>=1).sum()} ({100*(c>=1).mean():.0f}%) | {(d>=1).sum()} ({100*(d>=1).mean():.0f}%) | "
                     f"{np.nanmedian(c):.2f}x | {(a>=1).sum()} ({100*(a>=1).mean():.0f}%) | {(b>=1).sum()} ({100*(b>=1).mean():.0f}%) | {np.nanmedian(po):.2f} N |")
        L.append("")
    # size picking
    sc = np.array([[e["reach"], e["scale"] if e["scale"] is not None else np.nan,
                    e["crus"] if e["crus"] is not None else np.nan] for e in R if e["reach"] is not None])
    L.append("## Can a wearer pick the size by ear size?\n")
    for j, nm in ((1, "ear scale (basin inscribed radius)"), (2, "crus helicis distance")):
        ok = np.isfinite(sc[:, j])
        if ok.sum() < 5: continue
        r = np.corrcoef(sc[ok, 0], sc[ok, j])[0, 1]
        # assign by terciles of predictor, compare to true band
        pred = sc[ok, j]; true = sc[ok, 0]
        q = np.percentile(pred, [100/3, 200/3])
        pb = np.digitize(pred, q)
        edges = [cfg[k]["band"] for k in "SML"]
        tb = np.array([next(i for i, (lo, hi) in enumerate(edges) if lo - 1e-9 <= t <= hi + 1e-9) for t in true])
        hit = (pb == tb).mean()
        L.append(f"- **{nm}**: r = {r:.2f} with reach (n={ok.sum()}); picking S/M/L by terciles of it lands "
                 f"{100*hit:.0f}% of ears in their correct reach band.")
    L.append("\nReach is a body-to-lip distance in the seated pose, so it depends on how the body sits, not only on "
             "how big the ear is. If the correlation is weak, size selection has to be by fit (try the three pads), "
             "not by a tape measure.\n")
    md = "\n".join(L)
    open(OUT_MD, "w").write(md); print(md); print(f"\n-> {OUT_MD}")


if __name__ == "__main__":
    ap = argparse.ArgumentParser()
    ap.add_argument("--sweep", action="store_true"); ap.add_argument("--sweep-report", action="store_true")
    ap.add_argument("--reach", action="store_true")
    ap.add_argument("--stability", action="store_true")
    ap.add_argument("--report", action="store_true")
    a = ap.parse_args()
    if a.reach: phase_reach()
    if a.stability: phase_stability()
    if a.report or not (a.reach or a.stability or a.sweep or a.sweep_report): report()


# ---------------------------------------------------------------------------- #
# Sweep: what the pad has to be for the seated cases to hold.  The stability
# model caps the pad's normal force at PLUNGER_FORCE[0] per pad (0.18 N, the
# repelling-pair minimum) -- a DEPLOYMENT number.  A hooked pad is an interlock:
# the cartilage reacts whatever is asked, bounded by pain pressure (4 kPa).
# So cap = 4 kPa x pad area.  And the skirt push-out (0.31 N, mag-float) is the
# sustained demand; the bell lip's preload is whatever the clamp supplies.
# ---------------------------------------------------------------------------- #
SWEEP = [  # (label, pad cap N, skirt preload N)
    ("pad60mm2_skirt0.31", 0.24, 0.31),
    ("pad120mm2_skirt0.31", 0.48, 0.31),
    ("pad120mm2_skirt0.15", 0.48, 0.15),
    ("pad200mm2_skirt0.15", 0.80, 0.15),
]
CACHE_W = os.path.join(ALIGNED, "size_bands_sweep.json")


def phase_sweep():
    R = json.load(open(CACHE_R)); cfg = configs(R); P = iem_points()
    recs = {rec["ear_id"]: rec for rec in ears()}
    done = json.load(open(CACHE_W)) if os.path.exists(CACHE_W) else {}
    t0 = time.time()
    for e in R:
        eid = e["ear"]; band = next((k for k in "SML" if eid in cfg[k]["members"]), None)
        if band is None or e["patch_design"] is None: continue
        rec = recs[eid]
        patch = trimesh.load(os.path.join(ALIGNED, rec["patch"]), force="mesh"); field = EarField(patch, seed=0)
        M = np.array(rec["transform"], float)
        cp = transform(P["_cable_exit"][None, :], M)[0]; com = transform(P["_com"][None, :], M)[0]
        Pl = dict(P); Pl["plunger"] = np.array(e["patch_design"]); Pl["_plunger"] = [dict(name="cymba")]
        for label, cap, skirt in SWEEP:
            key = f"{eid}|{label}"
            if key in done: continue
            stab.PLUNGER_FORCE = (cap, cap); stab.SKIRT_PRELOAD = skirt
            row = dict(ear=eid, band=band, sweep=label)
            for mode in ("cone", "sphere"):
                x = stab.stability_check(rec, Pl, field, transform, cable_point=cp, com=com,
                                         mu=MU, cable_tug=0.5, cable_mode=mode)
                row[mode] = round(float(x["margin"]), 3)
            done[key] = row
        json.dump(done, open(CACHE_W, "w"), indent=1)
        print(f"{eid:<10} " + "  ".join(f"{l[:12]} {done[f'{eid}|{l}']['cone']:.2f}/{done[f'{eid}|{l}']['sphere']:.2f}" for l, _, _ in SWEEP) + f"  [{time.time()-t0:.0f}s]", flush=True)


def report_sweep():
    W = json.load(open(CACHE_W))
    L = ["\n\n## What the pad has to be (sweep over pad cap and skirt push-out)\n",
         "Pad placed on the cymba overhang, banded sizes (every ear seated), mu 0.6, 0.5 N cable tug. "
         "Pad normal-force cap = 4 kPa pain-onset pressure x pad area; skirt = sustained push-out of the seal.\n",
         "| pad area | skirt push-out | ears | cone pass | sphere pass | cone median | sphere median |", "|---|---|---|---|---|---|---|"]
    for label, cap, skirt in SWEEP:
        rows = [v for v in W.values() if v["sweep"] == label]
        if not rows: continue
        c = np.array([r["cone"] for r in rows]); s = np.array([r["sphere"] for r in rows])
        area = label.split("mm2")[0].replace("pad", "")
        L.append(f"| {area} mm² (cap {cap:.2f} N) | {skirt:.2f} N | {len(rows)} | {(c>=1).sum()} ({100*(c>=1).mean():.0f}%) | "
                 f"{(s>=1).sum()} ({100*(s>=1).mean():.0f}%) | {np.median(c):.2f}x | {np.median(s):.2f}x |")
    md = "\n".join(L); open(OUT_MD, "a").write(md); print(md)


if __name__ == "__main__" and "--sweep" in sys.argv:
    phase_sweep()
if __name__ == "__main__" and "--sweep-report" in sys.argv:
    report_sweep()
