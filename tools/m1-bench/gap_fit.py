#!/usr/bin/env python3
"""Gap characterization fit: peak-to-peak signal vs shim height against the
exponential decay law  Vpp(z) = V0 * exp(-2*pi*z/lambda).

Reads runs/h<height>.csv files (from capture.py, e.g. h0.5.csv h1.0.csv ...),
takes robust p-p (p2..p98) per channel, fits ln(Vpp) vs z by least squares,
reports the measured decay length vs theory, and where the curve crosses the
sensor's saturation floor — that crossing IS the gap budget.

Usage: gap_fit.py [--lambda 4.0] [--floor-mv 40]
  --lambda   pole PERIOD in mm (2 mm poles = 4.0)
  --floor-mv bridge output at the AAT saturation floor, in mV p-p
             (from the datasheet number you wrote down)
"""
import glob, math, os, re, sys

BASE = os.path.dirname(os.path.abspath(__file__))
lam = 4.0
floor_mv = None
args = sys.argv[1:]
for i, a in enumerate(args):
    if a == "--lambda":
        lam = float(args[i + 1])
    if a == "--floor-mv":
        floor_mv = float(args[i + 1])

pts = []
for path in sorted(glob.glob(os.path.join(BASE, "runs", "h*.csv"))):
    m = re.match(r"h([\d.]+)", os.path.basename(path))
    if not m:
        continue
    z = float(m.group(1))
    s0, s1 = [], []
    for line in open(path):
        parts = line.strip().split(",")
        if len(parts) >= 5:
            try:
                s0.append(float(parts[2]))
                s1.append(float(parts[4]))
            except ValueError:
                pass
    if len(s0) < 50:
        continue

    def pp(v):
        v = sorted(v)
        return v[int(0.98 * len(v))] - v[int(0.02 * len(v))]
    pts.append((z, pp(s0), pp(s1)))

if len(pts) < 2:
    sys.exit("need >=2 runs named h<height>.csv in runs/")

print(f"{'z mm':>6} {'sin p-p mV':>11} {'cos p-p mV':>11}")
for z, a, b in pts:
    print(f"{z:>6.2f} {a:>11.1f} {b:>11.1f}")

# fit ln(best-channel p-p) vs z
xs = [p[0] for p in pts]
ys = [math.log(max(p[1], p[2])) for p in pts]
n = len(xs)
sx, sy = sum(xs), sum(ys)
sxx = sum(x * x for x in xs)
sxy = sum(x * y for x, y in zip(xs, ys))
slope = (n * sxy - sx * sy) / (n * sxx - sx * sx)
intercept = (sy - slope * sx) / n
decay_len = -1.0 / slope if slope < 0 else float("inf")
theory = lam / (2 * math.pi)
print(f"\nmeasured decay length: {decay_len:.2f} mm per 1/e "
      f"(theory lambda/2pi = {theory:.2f} mm at lambda={lam} mm)")
print(f"signal halves every {decay_len * math.log(2):.2f} mm "
      f"(theory {theory * math.log(2):.2f} mm)")
v0 = math.exp(intercept)
print(f"extrapolated surface p-p: {v0:.0f} mV")
if floor_mv:
    z_cross = (math.log(v0) - math.log(floor_mv)) / (1.0 / decay_len)
    print(f"crossing vs floor ({floor_mv} mV): z = {z_cross:.2f} mm  "
          f"<-- THE GAP BUDGET (design target: >= 1.2 mm with margin)")
else:
    print("pass --floor-mv <datasheet saturation p-p> to get the gap budget")
