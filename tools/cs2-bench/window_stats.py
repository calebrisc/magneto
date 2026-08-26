#!/usr/bin/env python3
"""Warmup window stats: aim metrics for the last N minutes of the newest capture.

Usage: window_stats.py <sessions_dir> <minutes>
"""
import glob, json, os, sys, time
from round_watch import input_stats

SESS, MINS = sys.argv[1], float(sys.argv[2])
files = sorted(glob.glob(os.path.join(SESS, "*.gsi.jsonl")), key=os.path.getmtime)
if not files:
    print(json.dumps({"error": "no capture"})); sys.exit(0)
base = files[-1][:-len(".gsi.jsonl")]
anchor = float(open(base + ".anchor").read().strip())
now = time.time()
stats = input_stats(base + ".cap.csv", anchor, now - MINS * 60, now)
stats["window_min"] = MINS
stats["capture"] = os.path.basename(base)
print(json.dumps(stats))
