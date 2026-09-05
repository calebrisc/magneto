#!/usr/bin/env python3
"""
seat_construct.py -- seat the IEM in a scanned ear by anchors, not by optimiser.

Four anchor pairs fix the pose, one contact event fixes the slide:
  lip centreline centre  <->  aperture centroid       translation
  nozzle axis            <->  canal axis (datum ez)    two rotations
  jaw edge (+Y of body)  <->  crus helicis point       roll   (tragus if no crus)
  faceplate normal (+Z)  <->  concha-floor normal      VERIFY only
  slide along the canal axis until the lip or the nose first touches the ear.

Reports per ear what touched first (lip vs nose), lip contact fraction, shell
interference, jaw-to-crus distance and faceplate angle.  Writes a copy of the
ear record with the new transform to ears/seated_construct/ so nothing else
changes under the old caches.  Rerun any analysis with --json-dir there.

Usage:  python seat_construct.py --ears P0007,P0046 [--size XS] [--render]
"""
from __future__ import annotations
import argparse, glob, json, os, shutil, sys
import numpy as np, trimesh
from earfit import ALIGNED, HERE, EarField, iem_points, transform
from gyro_arm_variance import datum
import generate

OUT = os.path.join(HERE, "ears", "seated_construct")
CONTACT = 0.20          # mm: first-touch threshold
STEP = 0.05
JAW_DIR = np.array([0.0, 1.0, 0.0])          # design +Y = superior edge of the body


def rot_a_to_b(a, b):
    a = a / np.linalg.norm(a); b = b / np.linalg.norm(b)
    v = np.cross(a, b); c = float(a @ b)
    if np.linalg.norm(v) < 1e-9:
        return np.eye(3) if c > 0 else -np.eye(3)
    vx = np.array([[0, -v[2], v[1]], [v[2], 0, -v[0]], [-v[1], v[0], 0]])
    return np.eye(3) + vx + vx @ vx * (1.0 / (1.0 + c))


def rot_axis(n, ang):
    n = n / np.linalg.norm(n); K = np.array([[0, -n[2], n[1]], [n[2], 0, -n[0]], [-n[1], n[0], 0]])
    return np.eye(3) + np.sin(ang) * K + (1 - np.cos(ang)) * K @ K


def tip_samples(g, size):
    """Lip centreline ring (distal face) and nose-cone surface, design frame."""
    n = g.n_ax / np.linalg.norm(g.n_ax); e1 = np.array([0.0, 1.0, 0.0]); e2 = np.cross(n, e1)
    a, b = g.bell_lip_axes[size]; Rt = g.bell_Rt_by[size]
    th = np.linspace(0, 2 * np.pi, 72, endpoint=False)
    lip_c = g.nozzle_base + n * g.seal_x
    lip = lip_c + np.outer(a * np.cos(th), e1) + np.outer(b * np.sin(th), e2)
    nose = []
    for x, r in zip(np.linspace(g.bell_cb_x, g.bell_tip_x, 10), np.linspace(g.bell_base_r, g.bell_nose_r, 10)):
        c = g.nozzle_base + n * x
        nose.append(c + np.outer(r * np.cos(th[::3]), e1) + np.outer(r * np.sin(th[::3]), e2))
    return lip_c, lip, np.vstack(nose)


def seat(rec, patch, field, g, size, crus_d=None):
    """Body-first seat.

    The scans have no canal cavity, so nothing can be seated by the nozzle.
    Seat the BODY instead, which is what the ear actually holds:
      base       concha_frame: design +X anterior, +Y superior, +Z out of head
      roll       body +Y toward the crus (tragus if no crus), in the floor plane
      aim        translate so the nozzle axis line passes through the aperture
      slide      along the nozzle axis until the jacket ear-face touches the
                 concha floor OR the lip touches the funnel, whichever first
    What touched first is the design diagnostic: jacket first with the lip
    short of the funnel = nozzle too short / body rides high; lip first with
    the jacket floating = nozzle too long.
    """
    ap, B = datum(rec)
    A = np.array(rec["aperture"], float)
    F = np.array(rec["concha_frame"], float); R0 = F[:3, :3]
    floor_n = np.array(rec["floor_normal"], float); floor_n /= np.linalg.norm(floor_n)
    crus_w = None if crus_d is None else ap + B @ np.asarray(crus_d, float)
    trag_w = np.array(rec["tragus"], float) if rec.get("tragus") is not None else None
    P = iem_points()
    lip_c, lip, nose = tip_samples(g, size)
    ez = R0 @ np.array([0.0, 0.0, 1.0])                 # out of head
    # roll about ez so body +Y points at the crus in the floor plane
    target = crus_w if crus_w is not None else trag_w
    roll_ref = "crus" if crus_w is not None else "tragus"
    def perp(v): v = v - (v @ ez) * ez; return v / max(np.linalg.norm(v), 1e-9)
    core_w0 = F[:3, 3] + R0 @ np.array(g.core_c)
    d = perp(R0 @ JAW_DIR); v = perp(target - core_w0)
    ang = np.arctan2(np.cross(d, v) @ ez, d @ v)
    R = rot_axis(ez, ang) @ R0
    # insertion: nozzle in first (lip centre ON the aperture), then rotate the
    # body about the lip centre, pitching it down into the bowl until the jacket
    # ear-face meets the concha floor without the shell entering the walls.
    n_w = R @ (g.n_ax / np.linalg.norm(g.n_ax))
    pitch_axis = np.cross(n_w, ez); pitch_axis /= max(np.linalg.norm(pitch_axis), 1e-9)
    jacket_d = P["jacket"]; shell_d = P["shell"]
    best = None
    yaw_axis = ez
    for yaw in np.radians(np.arange(-30.0, 31.0, 5.0)):
        Ry = rot_axis(yaw_axis, yaw) @ R
        pa = np.cross(Ry @ (g.n_ax / np.linalg.norm(g.n_ax)), ez); pa /= max(np.linalg.norm(pa), 1e-9)
        for th in np.radians(np.arange(-60.0, 61.0, 2.0)):
            Rp = rot_axis(pa, th) @ Ry
            Mp = np.eye(4); Mp[:3, :3] = Rp; Mp[:3, 3] = A - Rp @ lip_c
            ds = field.query(transform(shell_d, Mp))
            pen = float((ds < -0.5).mean())
            if best is not None and best[0][0] is False and pen > 0.02:
                continue
            dj = field.query(transform(jacket_d, Mp))
            jc = float((np.abs(dj) <= 1.0).mean())
            key = (pen > 0.02, -jc if pen <= 0.02 else pen, abs(float(dj.min())))
            if best is None or key < best[0]:
                best = (key, th, yaw, Rp, Mp, dj, ds)
    _, th, yaw, R, M, dj, ds = best
    first = "jacket" if dj.min() <= CONTACT else "none"
    s = float(np.degrees(th)); yaw_deg = float(np.degrees(yaw))
    dl = field.query(transform(lip, M))
    dn = field.query(transform(nose, M))
    lip_w = transform(lip, M); lip_cw = lip_w.mean(axis=0)
    jaw_d = np.array(g.core_c) + (np.array(g.core_r) ** 2 * JAW_DIR) / np.linalg.norm(np.array(g.core_r) * JAW_DIR)
    jaw_w = transform(jaw_d[None, :], M)[0]
    shell = transform(P["shell"], M)
    face_n = R @ np.array([0.0, 0.0, 1.0])
    rep = dict(ear=rec["ear_id"], size=size, roll_ref=roll_ref, jacket_touch=first, pitch_deg=round(s, 1), yaw_deg=round(yaw_deg, 1),
               lip_centre_to_aperture=round(float(np.linalg.norm(lip_cw - A)), 2),
               lip_contact_frac=round(float((np.abs(dl) <= 1.0).mean()), 2), lip_min=round(float(dl.min()), 2),
               jacket_contact_frac=round(float((np.abs(dj) <= 1.0).mean()), 2), jacket_min=round(float(dj.min()), 2),
               nose_min=round(float(dn.min()), 2),
               shell_pen_frac=round(float((ds < -0.5).mean()), 3), shell_min=round(float(ds.min()), 2),
               jaw_to_crus=None if crus_w is None else round(float(np.linalg.norm(jaw_w - crus_w)), 2),
               faceplate_vs_floor_deg=round(float(np.degrees(np.arccos(np.clip(face_n @ floor_n, -1, 1)))), 1),
               transform=M.tolist())
    return rep, M, dict(A=A, n=n_w, crus=crus_w, tragus=trag_w, lip=lip_w, nose=transform(nose, M), jaw=jaw_w, shell=shell,
                        jacket=transform(jacket_d, M))


def render(rec, patch, geo, rep, out):
    import matplotlib; matplotlib.use("Agg"); import matplotlib.pyplot as plt
    V = np.asarray(patch.vertices); A = geo["A"]; near = np.linalg.norm(V - A, axis=1) < 22
    V = V[near][::max(1, near.sum() // 7000)]
    fig = plt.figure(figsize=(13, 6.5), dpi=130)
    for j, (el, az) in enumerate(((20, -60), (75, -90), (0, 0))):
        ax = fig.add_subplot(1, 3, 1 + j, projection="3d")
        ax.scatter(*V.T, s=1.5, c="#b9a58a", alpha=0.35, linewidths=0)
        ax.scatter(*geo["shell"][::2].T, s=2.5, c="#7d8590", alpha=0.7, linewidths=0)
        ax.scatter(*geo["jacket"].T, s=6, c="#3f6fd0", alpha=0.9, linewidths=0)
        ax.plot(*np.vstack([geo["lip"], geo["lip"][:1]]).T, c="#2f6fd6", lw=1.8)
        ax.scatter(*geo["nose"][::2].T, s=3, c="#1b1f24", alpha=0.6, linewidths=0)
        ax.scatter(*A, s=60, c="#d6336c", marker="o", zorder=6); ax.text(*A, " aperture", fontsize=7, color="#d6336c")
        ax.quiver(*A, *(geo["n"] * 6), color="#d6336c", lw=1.5, arrow_length_ratio=0.3)
        if geo["crus"] is not None:
            ax.scatter(*geo["crus"], s=60, c="#3b8f5c", marker="^", zorder=6); ax.text(*geo["crus"], " crus", fontsize=7, color="#3b8f5c")
        if geo["tragus"] is not None:
            ax.scatter(*geo["tragus"], s=40, c="#b97a12", marker="s", zorder=6); ax.text(*geo["tragus"], " tragus", fontsize=7, color="#b97a12")
        ax.scatter(*geo["jaw"], s=50, c="#6a4fb3", marker="D", zorder=6); ax.text(*geo["jaw"], " jaw", fontsize=7, color="#6a4fb3")
        r = 15; ax.set_xlim(A[0]-r, A[0]+r); ax.set_ylim(A[1]-r, A[1]+r); ax.set_zlim(A[2]-r, A[2]+r)
        ax.view_init(elev=el, azim=az); ax.set_axis_off()
    fig.suptitle(f"{rec['dataset']}/{rec['ear_id']}  seated by construction, {rep['size']} lip  |  pitch {rep['pitch_deg']}°, jacket touch: {rep['jacket_touch']}  "
                 f"lip→aperture {rep['lip_centre_to_aperture']} mm  jacket contact {int(100*rep['jacket_contact_frac'])}%  shell min {rep['shell_min']} mm (pen {rep['shell_pen_frac']})  "
                 f"jaw→crus {rep['jaw_to_crus']} mm  faceplate vs floor {rep['faceplate_vs_floor_deg']}°", fontsize=8.5)
    fig.tight_layout(); fig.savefig(out); plt.close(fig)


def main():
    ap = argparse.ArgumentParser(); ap.add_argument("--ears", default=None); ap.add_argument("--size", default="XS")
    ap.add_argument("--all", action="store_true")
    ap.add_argument("--render", action="store_true"); a = ap.parse_args()
    os.makedirs(OUT, exist_ok=True)
    g = generate.G(generate.PARAMS)
    gav = {e["id"]: e for e in json.load(open(os.path.join(ALIGNED, "gyro_arm_variance.json")))}
    if a.all:
        from size_bands import ears as _ears
        ids = [r["ear_id"] for r in _ears()]
    else:
        ids = a.ears.split(",")
    summary = []
    for eid in ids:
        f = glob.glob(os.path.join(ALIGNED, f"*_{eid}_right.json"))[0]; rec = json.load(open(f))
        patch = trimesh.load(os.path.join(ALIGNED, rec["patch"]), force="mesh"); field = EarField(patch, seed=0)
        crus = gav.get(eid, {}).get("crus")
        rep, M, geo = seat(rec, patch, field, g, a.size, crus)
        rec2 = dict(rec); rec2["transform"] = M.tolist(); rec2["seated_by"] = "construct"; rec2["seat_report"] = {k: v for k, v in rep.items() if k != "transform"}
        json.dump(rec2, open(os.path.join(OUT, os.path.basename(f)), "w"), indent=1)
        src = os.path.join(ALIGNED, rec["patch"]); dst = os.path.join(OUT, rec["patch"])
        if not os.path.exists(dst): os.symlink(os.path.abspath(src), dst)
        summary.append({k: v for k, v in rep.items() if k != "transform"})
        print(json.dumps(summary[-1]), flush=True)
        if a.render:
            out = os.path.join(HERE, "viz", f"seat_{eid}.png"); render(rec, patch, geo, rep, out); print("  ->", out)
    if a.all:
        import csv
        with open(os.path.join(OUT, "seat_summary.csv"), "w", newline="") as fh:
            w = csv.DictWriter(fh, fieldnames=list(summary[0].keys())); w.writeheader(); w.writerows(summary)
        ok = [r for r in summary if r["shell_pen_frac"] <= 0.02 and r["jacket_contact_frac"] >= 0.1]
        print(f"\nbody fits (shell penetration <= 2%, jacket touching) on {len(ok)}/{len(summary)} ears")
        for r in sorted(ok, key=lambda r: -r["jacket_contact_frac"])[:8]:
            print("  ", r["ear"], "jacket", r["jacket_contact_frac"], "lip", r["lip_contact_frac"], "pitch", r["pitch_deg"], "yaw", r["yaw_deg"], "jaw->crus", r["jaw_to_crus"])


if __name__ == "__main__":
    main()
