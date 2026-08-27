#!/usr/bin/env python3
"""Warmup-vs-match gap: same input metrics, split by session type.

Classifies every Val capture in sessions/ by mode (from its .val.jsonl:
queueId / range) and prints the side-by-side that answers "does my warmup
form survive matches?" Usage: warmup_gap.py [days_back]
"""
import glob, json, os, sys, time
from round_watch import input_stats
from flick_bias import analyze as flick_analyze

BASE = os.path.dirname(os.path.abspath(__file__))
SESS = os.path.join(BASE, "sessions")
DAYS = float(sys.argv[1]) if len(sys.argv) > 1 else 7.0
cutoff = time.time() - DAYS * 86400

CLASSES = {"range": "practice", "custom": "practice", "unknown": "practice",
           "deathmatch": "dm", "hurm": "dm", "competitive": "match",
           "unrated": "match", "swiftplay": "match"}

buckets = {}
for vpath in glob.glob(os.path.join(SESS, "*val_*.val.jsonl")):
    base = vpath[:-len(".val.jsonl")]
    try:
        anchor = float(open(base + ".anchor").read().strip())
    except Exception:
        continue
    if anchor < cutoff:
        continue
    mode = "unknown"
    try:
        for line in open(vpath):
            p = json.loads(line)["private"]
            q = p.get("queueId") or ""
            mm = (p.get("matchMap") or "").lower()
            if "range" in mm:
                mode = "range"
            elif q:
                mode = q
            break
    except Exception:
        pass
    cls = CLASSES.get(mode, "match" if mode not in ("range",) else "practice")
    if "range" in os.path.basename(base):
        cls = "practice"
    buckets.setdefault(cls, []).append(base)

print(f"{'':<10} {'sessions':>8} {'shots':>6} {'cstrafe%':>9} {'held%':>6} "
      f"{'moving%':>8} {'early%':>7} {'flicks':>7} {'clean%':>7} {'over%':>6}")
for cls in ("practice", "dm", "match"):
    bases = buckets.get(cls, [])
    if not bases:
        continue
    tot = {"shots": 0, "meas": 0, "cs": 0, "held": 0, "mov": 0, "early": 0}
    fl = {"n": 0, "clean": 0, "over": 0}
    for b in bases:
        s = input_stats(b + ".cap.csv", float(open(b + ".anchor").read()),
                        0, 9e12)
        tot["shots"] += s.get("shots", 0)
        tot["meas"] += s.get("cs_measured", 0)
        tot["cs"] += s.get("cstrafe_shots", 0)
        tot["held"] += s.get("cstrafe_held_at_shot", 0)
        tot["mov"] += s.get("moving_shots", 0)
        tot["early"] += s.get("cs_early", 0)
        f = flick_analyze(b + ".cap.csv")
        if f:
            fl["n"] += f["flicks"]
            fl["clean"] += f["clean"]
            fl["over"] += f["overshoot"]
    def pct(a, b):
        return f"{round(100 * a / b)}%" if b else "-"
    print(f"{cls:<10} {len(bases):>8} {tot['shots']:>6} "
          f"{pct(tot['cs'], tot['meas']):>9} {pct(tot['held'], max(1, tot['cs'])):>6} "
          f"{pct(tot['mov'], tot['shots']):>8} {pct(tot['early'], tot['meas']):>7} "
          f"{fl['n']:>7} {pct(fl['clean'], fl['n']):>7} {pct(fl['over'], fl['n']):>6}")
print("\nNote: on util-heavy agents, ability clicks inflate match shot counts"
      " until cast-click exclusion lands.")
