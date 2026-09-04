#!/usr/bin/env python3
"""viz_test.py -- picture of what the stability test is: ear, seated body, contacts, loads."""
import json, os, sys, numpy as np, trimesh
import matplotlib; matplotlib.use("Agg")
import matplotlib.pyplot as plt
from earfit import ALIGNED, HERE, NOZZLE_AXIS, EarField, iem_points, transform
import seal_compliance as sc
import stability as stab

OUT = os.path.join(HERE, "viz")
R = {e["ear"]: e for e in json.load(open(os.path.join(ALIGNED, "size_bands_reach.json")))}
files = {e: f for e, f in ((r["ear_id"], os.path.basename(p)) for p, r in
         ((p, json.load(open(p))) for p in __import__("glob").glob(os.path.join(ALIGNED, "*_right.json"))))}
P = iem_points()

def draw(eid, label, ax, elev, azim):
    rec = json.load(open(os.path.join(ALIGNED, files[eid])))
    patch = trimesh.load(os.path.join(ALIGNED, rec["patch"]), force="mesh"); field = EarField(patch, seed=0)
    M = np.array(rec["transform"], float); Rm = M[:3, :3]
    ap = np.array(rec["aperture"]); nax = Rm @ NOZZLE_AXIS; nax /= np.linalg.norm(nax)
    com = transform(P["_com"][None, :], M)[0]; cp = transform(P["_cable_exit"][None, :], M)[0]
    V = np.asarray(patch.vertices); near = np.linalg.norm(V - ap, axis=1) < 22
    V = V[near][::max(1, near.sum() // 6000)]
    ax.scatter(*V.T, s=1.5, c="#b9a58a", alpha=0.35, linewidths=0)
    body = transform(np.vstack([P["shell"], P["faceplate"]]), M)[::3]
    ax.scatter(*body.T, s=2, c="#7d8590", alpha=0.6, linewidths=0)
    ring = transform(sc.rim_from_mesh(n=120), M)
    ax.plot(*np.vstack([ring, ring[:1]]).T, c="#2f6fd6", lw=1.6)
    e = R[eid]
    if e["patch_design"] is not None:
        pad = transform(np.array(e["patch_design"]), M); _, idx = field.tree.query(pad); n = field.nrm[idx]
        ax.scatter(*pad.T, s=40, c="#d6336c", edgecolors="k", linewidths=0.4, zorder=5)
        ax.quiver(*pad.T, *(n * 2.5).T, color="#d6336c", lw=1.2, arrow_length_ratio=0.35)
    # loads: skirt push-out at rim centre, cable tug cone at cable exit, straight yank at COM
    rc = ring.mean(axis=0)
    ax.quiver(*rc, *(-nax * 6), color="#2f6fd6", lw=2.2, arrow_length_ratio=0.25)
    ax.text(*(rc - nax * 7), "skirt push-out 0.31 N", fontsize=7, color="#2f6fd6")
    for d in stab._cone_dirs(np.array([-1.0, 0.0, -1.0]), 45.0, 5):
        ax.quiver(*cp, *(d * 5), color="#e8a33d", lw=1.1, arrow_length_ratio=0.3, alpha=0.9)
    ax.text(*(cp + np.array([0, 0, -7])), "cable tug 0.5 N (cone)", fontsize=7, color="#b97a12")
    ax.quiver(*com, *(-nax * 8), color="#111", lw=2.0, arrow_length_ratio=0.22)
    ax.text(*(com - nax * 9.5), "straight yank", fontsize=7, color="#111")
    ax.scatter(*com, s=30, c="k", marker="x")
    c = ap; r = 16
    ax.set_xlim(c[0]-r, c[0]+r); ax.set_ylim(c[1]-r, c[1]+r); ax.set_zlim(c[2]-r, c[2]+r)
    ax.view_init(elev=elev, azim=azim); ax.set_axis_off()
    ax.set_title(f"{label}: {eid}   reach {e['reach']:.1f} mm, wrap {e['wrap']:.1f} mm", fontsize=9, loc="left")

def main():
    from size_bands import configs
    cfg = configs(list(R.values()))
    picks = []
    for k in "SML":
        mem = [R[m] for m in cfg[k]["members"] if R[m]["patch_design"] is not None]
        med = np.median([m["reach"] for m in mem])
        picks.append((f"{k} band", min(mem, key=lambda m: abs(m["reach"] - med))["ear"]))
    print(picks)
    fig = plt.figure(figsize=(13, 8.6), dpi=130)
    for i, (lab, eid) in enumerate(picks):
        for j, (el, azm) in enumerate(((20, -60), (75, -90))):
            ax = fig.add_subplot(2, 3, 1 + i + 3 * j, projection="3d"); draw(eid, lab if j == 0 else "top view", ax, el, azm)
    fig.suptitle("What the retention test is: rigid body seated in a scanned ear. Tan = ear, grey = body, blue = sealing ring, "
                 "pink = pad on the cymba overhang (arrows = cartilage normals). Loads: skirt push-out, cable cone, straight yank.", fontsize=8.5)
    fig.tight_layout(); out = os.path.join(OUT, "stability_test_SML.png"); fig.savefig(out); print(out)

if __name__ == "__main__":
    main()
