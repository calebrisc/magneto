#!/usr/bin/env python3
"""Network health probe: 1 Hz pings to gateway / public DNS / Riot edge.

Writes sessions/net_<YYYYMMDD>.csv rows: t_wall,target,rtt_ms (-1 = loss).
Joins to any capture by wall time — flags loss bursts and latency spikes
during rounds. Windows-only (uses ping.exe). Usage: net_probe.py [hours]
"""
import os, re, subprocess, sys, time

BASE = os.path.dirname(os.path.abspath(__file__))
SESS = os.path.join(BASE, "sessions"); os.makedirs(SESS, exist_ok=True)
HOURS = float(sys.argv[1]) if len(sys.argv) > 1 else 8.0
TARGETS = []


def gateway():
    try:
        out = subprocess.run(["route", "print", "0.0.0.0"], capture_output=True,
                             text=True, timeout=5).stdout
        m = re.search(r"0\.0\.0\.0\s+0\.0\.0\.0\s+(\d+\.\d+\.\d+\.\d+)", out)
        return m.group(1) if m else None
    except Exception:
        return None


gw = gateway()
if gw:
    TARGETS.append(("gateway", gw))
TARGETS += [("dns", "8.8.8.8"), ("riot", "glz-na-1.na.a.pvp.net")]


def ping_ms(host):
    try:
        out = subprocess.run(["ping", "-n", "1", "-w", "1000", host],
                             capture_output=True, text=True, timeout=3).stdout
        m = re.search(r"[<=](\d+)ms", out)
        return int(m.group(1)) if m else -1
    except Exception:
        return -1


print(f"net probe up: {[t[0] for t in TARGETS]}", flush=True)
end = time.time() + HOURS * 3600
path = None
f = None
while time.time() < end:
    day = time.strftime("%Y%m%d")
    p = os.path.join(SESS, f"net_{day}.csv")
    if p != path:
        if f: f.close()
        new = not os.path.exists(p)
        f = open(p, "a", buffering=1)
        if new:
            f.write("t_wall,target,rtt_ms\n")
        path = p
    t0 = time.time()
    for name, host in TARGETS:
        f.write(f"{time.time():.3f},{name},{ping_ms(host)}\n")
    time.sleep(max(0.0, 1.0 - (time.time() - t0)))
print("net probe down", flush=True)
