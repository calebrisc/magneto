#!/usr/bin/env python3
"""
ear_sizes.py -- one row per ear: which S/M/L it takes on each sized part.

Three sized parts, three independent size picks per ear:
  bell lip     by canal-aperture footprint (funnel traced 0.8 mm off the
               surface, PCA extents in the plane normal to the canal axis).
               Rule: smallest lip whose contact centreline (outer - 3 mm tube)
               clears the aperture by LIP_MARGIN on both axes.
  cymba pad    by along-aim reach to the cymba overhang (size_bands.py bands)
  crus pad     by travel to the crus overhang (crus_bands.py bands)

Writes docs/EAR_SIZE_MATRIX.md and ears/aligned/ear_sizes.csv.
"""
from __future__ import annotations
import csv, json, os
import numpy as np, trimesh
from earfit import ALIGNED, HERE
from size_bands import ears, configs, bands_from, CACHE_R
from viz_marked import aperture_ring
import generate

LIP_MARGIN = 1.0
BELL = generate.PARAMS["bell_sizes"]           # outer H x W
SIZES = list(generate.BELL_SIZES)              # smallest first
TUBE = {k: generate.PARAMS.get("bell_tube_by_size", {}).get(k, generate.PARAMS["bell_lip_tube_d"]) for k in BELL}
CRUS_CACHE = os.path.join(ALIGNED, "crus_bands.json")
OUT_MD = os.path.join(HERE, "..", "..", "docs", "EAR_SIZE_MATRIX.md")
OUT_CSV = os.path.join(ALIGNED, "ear_sizes.csv")


def aperture_dims(rec, patch):
    ap = np.array(rec["aperture"], float); outward = np.array(rec["outward"], float)
    axis = -np.array(rec["floor_normal"], float); axis /= np.linalg.norm(axis)
    ring = aperture_ring(patch, ap, axis, outward)
    if len(ring) < 12:
        return None, None, 0
    X = ring - ring.mean(axis=0); X -= np.outer(X @ axis, axis)
    _, _, Vt = np.linalg.svd(X, full_matrices=False)
    a = X @ Vt[0]; b = X @ Vt[1]
    return float(a.max() - a.min()), float(b.max() - b.min()), len(ring)


def bell_size(H, W):
    if H is None: return "?"
    for k in SIZES:
        h, w = BELL[k]; hc, wc = h - TUBE[k], w - TUBE[k]
        if hc >= H + 2 * LIP_MARGIN and wc >= W + 2 * LIP_MARGIN:
            return k
    return "L+"


def main():
    R = json.load(open(CACHE_R)); cfg = configs(R); reach = {e["ear"]: e for e in R}
    C = json.load(open(CRUS_CACHE)); crus = {r["ear"]: r for r in C}
    ct = np.array([r["travel"] for r in C if r["travel"] is not None]); cb = bands_from(ct)
    def crus_band(t):
        if t is None: return "–"
        for k, (lo, hi) in zip("SML", cb):
            if lo - 1e-9 <= t <= hi + 1e-9: return k
        return "–"
    def cymba_band(eid):
        return next((k for k in "SML" if eid in cfg[k]["members"]), "–")
    rows = []
    for rec in ears():
        patch = trimesh.load(os.path.join(ALIGNED, rec["patch"]), force="mesh")
        H, W, n = aperture_dims(rec, patch)
        e = reach[rec["ear_id"]]; c = crus.get(rec["ear_id"], {})
        rows.append(dict(ear=rec["ear_id"], dataset=rec["dataset"], scale=round(e["scale"], 2) if e.get("scale") else None,
                         ap_H=None if H is None else round(H, 1), ap_W=None if W is None else round(W, 1),
                         bell=bell_size(H, W), reach=None if e["reach"] is None else round(e["reach"], 1),
                         cymba=cymba_band(rec["ear_id"]), travel=None if c.get("travel") is None else round(c["travel"], 1),
                         crus=crus_band(c.get("travel"))))
        print(rows[-1])
    with open(OUT_CSV, "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=list(rows[0].keys())); w.writeheader(); w.writerows(rows)
    # summary
    def count(key): 
        d = {}; [d.__setitem__(r[key], d.get(r[key], 0) + 1) for r in rows]; return d
    sized = [r for r in rows if r["bell"] in SIZES and r["cymba"] in "SML"]
    same = sum(r["bell"] == r["cymba"] for r in sized)
    three = [r for r in rows if r["bell"] in SIZES and r["cymba"] in "SML" and r["crus"] in "SML"]
    same3 = sum(r["bell"] == r["cymba"] == r["crus"] for r in three)
    L = ["# Ear size matrix — which S/M/L each ear takes on each sized part (2026-09-04)\n",
         f"{len(rows)} real ears. Sizes are picked per PART, independently; an ear can be S on the lip and L on the pad. "
         "Script: `cad/iem/ear_sizes.py`; CSV at `cad/iem/ears/aligned/ear_sizes.csv`.\n",
         f"- bell lip: {count('bell')}  (rule: smallest lip whose centreline (outer minus tube) clears the aperture by {LIP_MARGIN} mm; "
         f"sizes {', '.join(f'{k} {BELL[k][0]:.0f}x{BELL[k][1]:.0f} O{TUBE[k]:.0f} tube' for k in SIZES)}; L+ = aperture too big for the L lip)",
         f"- cymba pad (reach band): {count('cymba')}",
         f"- crus pad (travel band, ears with a detected crus): {count('crus')}",
         f"- bell size == cymba pad size on {same}/{len(sized)} ears ({100*same/max(len(sized),1):.0f}%); all three agree on {same3}/{len(three)} "
         f"({100*same3/max(len(three),1):.0f}%). **Sizes do not travel together** — each part is fitted on its own.\n",
         "| ear | set | scale | aperture H×W | bell | reach | cymba pad | crus travel | crus pad |", "|---|---|---|---|---|---|---|---|---|"]
    for r in sorted(rows, key=lambda r: (r["bell"], r["cymba"], r["ear"])):
        L.append(f"| {r['ear']} | {r['dataset']} | {r['scale']} | {r['ap_H']}×{r['ap_W']} | **{r['bell']}** | {r['reach']} | **{r['cymba']}** | {r['travel']} | **{r['crus']}** |")
    open(OUT_MD, "w").write("\n".join(L)); print("\n".join(L[:8])); print("->", OUT_MD)


if __name__ == "__main__":
    main()
