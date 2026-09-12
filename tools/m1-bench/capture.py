#!/usr/bin/env python3
"""M1 bench capture: auto-finds the DK's VCOM port, logs the 200 Hz CSV
stream, and prints live per-second peak-to-peak so you can see signal at the
bench without opening anything else.

Usage: capture.py <label>            e.g. capture.py h1.0   (shim height)
Stop with Ctrl+C. Output: runs/<label>.csv (+ .anchor wall-clock file).
"""
import os, sys, time

try:
    import serial
    from serial.tools import list_ports
except ImportError:
    sys.exit("pip install pyserial")

BASE = os.path.dirname(os.path.abspath(__file__))
RUNS = os.path.join(BASE, "runs")
os.makedirs(RUNS, exist_ok=True)
label = sys.argv[1] if len(sys.argv) > 1 else time.strftime("%H%M%S")

port = None
for p in list_ports.comports():
    blob = f"{p.description} {p.manufacturer or ''}".lower()
    if "jlink" in blob or "j-link" in blob or "segger" in blob:
        port = p.device
        break
if port is None:
    ports = list(list_ports.comports())
    if len(ports) == 1:
        port = ports[0].device
    else:
        sys.exit("DK VCOM port not found. Ports: "
                 + ", ".join(f"{p.device} ({p.description})" for p in ports))

out = os.path.join(RUNS, f"{label}.csv")
open(os.path.join(RUNS, f"{label}.anchor"), "w").write(str(time.time()))
print(f"capturing {port} -> {out}   (Ctrl+C to stop)")
n = 0
window = []          # (t, mv0, mv1) for live peak-to-peak
last_print = time.time()
with serial.Serial(port, 115200, timeout=2) as ser, open(out, "w") as f:
    try:
        while True:
            line = ser.readline().decode(errors="ignore").strip()
            if not line:
                continue
            f.write(line + "\n")
            n += 1
            parts = line.split(",")
            if len(parts) >= 5:
                try:
                    window.append((time.time(), float(parts[2]), float(parts[4])))
                except ValueError:
                    pass
            now = time.time()
            if now - last_print >= 1.0:
                window = [w for w in window if now - w[0] <= 1.0]
                if window:
                    s0 = [w[1] for w in window]
                    s1 = [w[2] for w in window]
                    print(f"  {n:>7} rows | 1s p-p: sin {max(s0)-min(s0):7.1f} mV"
                          f"   cos {max(s1)-min(s1):7.1f} mV", flush=True)
                last_print = now
    except KeyboardInterrupt:
        pass
print(f"done: {n} rows -> {out}")
