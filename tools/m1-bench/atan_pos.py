#!/usr/bin/env python3
"""Position from sin/cos: atan2 + phase unwrap = distance travelled along the
strip. The 1D tracking proof — run on a sled-slide capture.

Usage: atan_pos.py runs/<file>.csv [--lambda 4.0]
Prints travel, peak speed, and writes <file>.pos.csv (t_ms, pos_mm).
"""
import math, os, sys

lam = 4.0
if "--lambda" in sys.argv:
    lam = float(sys.argv[sys.argv.index("--lambda") + 1])
path = sys.argv[1]

rows = []
for line in open(path):
    p = line.strip().split(",")
    if len(p) >= 5:
        try:
            rows.append((float(p[0]), float(p[2]), float(p[4])))
        except ValueError:
            continue
if len(rows) < 10:
    sys.exit("not enough rows")

# centre both channels (mid-rail offset) using median
s0m = sorted(r[1] for r in rows)[len(rows) // 2]
s1m = sorted(r[2] for r in rows)[len(rows) // 2]

pos, prev_ph, unwrapped = [], None, 0.0
for t, a, b in rows:
    ph = math.atan2(a - s0m, b - s1m)
    if prev_ph is not None:
        d = ph - prev_ph
        if d > math.pi:
            d -= 2 * math.pi
        elif d < -math.pi:
            d += 2 * math.pi
        unwrapped += d
    prev_ph = ph
    pos.append((t, unwrapped / (2 * math.pi) * lam))

out = path.rsplit(".csv", 1)[0] + ".pos.csv"
with open(out, "w") as f:
    f.write("t_ms,pos_mm\n")
    for t, x in pos:
        f.write(f"{t:.1f},{x:.4f}\n")

travel = max(x for _, x in pos) - min(x for _, x in pos)
# peak speed over 20 ms windows
peak = 0.0
j = 0
for i in range(len(pos)):
    while pos[i][0] - pos[j][0] > 20.0:
        j += 1
    if i > j:
        dt = (pos[i][0] - pos[j][0]) / 1000.0
        if dt > 0:
            peak = max(peak, abs(pos[i][1] - pos[j][1]) / 1000.0 / dt)
print(f"rows {len(rows)} | travel {travel:.1f} mm | "
      f"peak speed {peak:.2f} m/s | wrote {out}")
print("sanity: hand-fling target >= 1.5 m/s tracked without glitches "
      "(no position jumps > lambda/2 between samples)")
jumps = sum(1 for i in range(1, len(pos))
            if abs(pos[i][1] - pos[i - 1][1]) > lam / 2)
print(f"glitch check: {jumps} inter-sample jumps > {lam/2:.1f} mm"
      + ("  <-- LOST LOCK, slow down or raise sample rate" if jumps else "  (clean)"))
