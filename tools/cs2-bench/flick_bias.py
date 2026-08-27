#!/usr/bin/env python3
"""Flick correction-bias analysis: does the post-flick correction go the SAME
direction (undershoot — flick landed short) or REVERSE (overshoot)?

A sens slower than trained muscle memory predicts same-direction bias.
Usage: flick_bias.py <cap.csv> [<cap.csv> ...]
"""
import csv, sys
from report_gen import cap_time_col

BIN_S = 0.010          # 10 ms bins of net dx
FLICK_MIN = 250        # counts — only real flicks (several degrees)
FLICK_PEAK = 80        # counts in one bin — must be ballistic, not a drag
FLICK_MAX_S = 0.14     # ballistic phase duration cap
CORR_S = 0.15          # correction window after the flick
CORR_FRAC = 0.05       # correction must be >=5% of flick to count as biased


def analyze(path):
    rd = csv.DictReader(open(path))
    tcol, scale = cap_time_col(rd.fieldnames)
    bins = {}
    t0 = None
    for r in rd:
        if r["kind"] != "m":
            continue
        try:
            t = int(r[tcol]) * scale
            dx = int(r["dx"])
        except Exception:
            continue
        if t0 is None:
            t0 = t
        bins[int((t - t0) / BIN_S)] = bins.get(int((t - t0) / BIN_S), 0) + dx
    if not bins:
        return None
    idxs = sorted(bins)
    under = over = clean = 0
    i = 0
    while i < len(idxs):
        b = idxs[i]
        v = bins[b]
        if abs(v) >= FLICK_PEAK:
            sign = 1 if v > 0 else -1
            total = 0
            j = b
            while j in bins and bins[j] * sign > 0 and (j - b) * BIN_S <= FLICK_MAX_S:
                total += bins[j]
                j += 1
            if abs(total) >= FLICK_MIN:
                corr = sum(bins.get(k, 0)
                           for k in range(j, j + int(CORR_S / BIN_S)))
                ratio = corr / total  # >0 same direction, <0 reversal
                if ratio >= CORR_FRAC:
                    under += 1
                elif ratio <= -CORR_FRAC:
                    over += 1
                else:
                    clean += 1
                while i < len(idxs) and idxs[i] < j + int(CORR_S / BIN_S):
                    i += 1
                continue
        i += 1
    n = under + over + clean
    return {"flicks": n, "undershoot": under, "overshoot": over, "clean": clean,
            "under_pct": round(100 * under / n) if n else None,
            "over_pct": round(100 * over / n) if n else None}


if __name__ == "__main__":
    for p in sys.argv[1:]:
        r = analyze(p)
        name = p.rsplit("\\", 1)[-1].rsplit("/", 1)[-1]
        if r:
            print(f"{name:<44} flicks{r['flicks']:>5}  "
                  f"under {r['under_pct']:>3}%  over {r['over_pct']:>3}%  "
                  f"clean {100 - r['under_pct'] - r['over_pct']:>3}%")
        else:
            print(f"{name:<44} no data")
