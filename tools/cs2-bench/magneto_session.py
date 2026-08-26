#!/usr/bin/env python3
"""Magneto bench session supervisor — the only process that runs.

- Embeds the GSI HTTP listener (port 3202).
- Arms the input tap ONLY while a match is live (menus = tap off, zero work).
- Disables pointer acceleration during capture, restores it after (and on exit).
- Detects match end (map phase 'gameover', or map exit after live play),
  then writes a full HTML report to reports/ and prints REPORT:<path>.

Files per match under sessions/: <stamp>_<map>.gsi.jsonl, .cap.csv, .anchor
Usage: magneto_session.py [hours_to_run]
"""
import json, time, os, sys, csv, threading, subprocess, atexit
from http.server import HTTPServer, BaseHTTPRequestHandler
IS_WIN = os.name == "nt"
if IS_WIN:
    import win_input_tap
else:
    import Quartz
import report_gen

BASE = os.path.dirname(os.path.abspath(__file__))
SESS = os.path.join(BASE, "sessions"); os.makedirs(SESS, exist_ok=True)
REPORTS = os.path.join(BASE, "reports"); os.makedirs(REPORTS, exist_ok=True)
RUN_HOURS = float(sys.argv[1]) if len(sys.argv) > 1 else 6.0

state = {
    "phase": None, "map": None, "mode": None,
    "gsi_file": None, "gsi_path": None,
    "cap_file": None, "cap_path": None, "anchor_path": None,
    "tap_on": False, "live_seen": False, "stamp": None,
    "accel_saved": None,
}
lock = threading.Lock()

# ---------- pointer acceleration ----------
# macOS: toggled off during capture, restored after. Windows: raw input deltas
# are pre-acceleration so capture is clean either way, but the game feel isn't —
# warn once if "Enhance pointer precision" is on instead of silently rewriting
# a system setting.
def accel_off():
    if IS_WIN:
        try:
            import ctypes
            params = (ctypes.c_int * 3)()
            ctypes.windll.user32.SystemParametersInfoW(0x0003, 0, params, 0)  # SPI_GETMOUSE
            if params[2]:
                print("WARN: 'Enhance pointer precision' is ON "
                      "(Settings > Bluetooth & devices > Mouse > Additional mouse settings "
                      "> Pointer Options) — turn it off for consistent aim.", flush=True)
        except Exception: pass
        return
    try:
        out = subprocess.run(["defaults","read","-g","com.apple.mouse.scaling"],
                             capture_output=True, text=True).stdout.strip()
        state["accel_saved"] = out or "0.875"
        subprocess.run(["defaults","write","-g","com.apple.mouse.scaling","-1"])
    except Exception: pass

def accel_restore():
    if IS_WIN: return
    if state["accel_saved"] is not None:
        subprocess.run(["defaults","write","-g","com.apple.mouse.scaling",
                        "-float", state["accel_saved"]])
        state["accel_saved"] = None
atexit.register(accel_restore)

# ---------- platform tap control ----------
def tap_set(on):
    if IS_WIN:
        win_input_tap.set_enabled(on)
    elif tap_holder[0] is not None:
        Quartz.CGEventTapEnable(tap_holder[0], on)

# ---------- input tap (runs in its own thread/runloop) ----------
WASD = {0:"a",1:"s",2:"d",13:"w",49:"j"}  # j = jump (space)
MODS = {56:"walkmod",59:"crouchmod"}  # shift, ctrl via flags
mod_state = {"walkmod":False,"crouchmod":False}
tap_holder = [None]

def tap_cb(proxy, etype, event, refcon):
    if etype in (Quartz.kCGEventTapDisabledByTimeout, Quartz.kCGEventTapDisabledByUserInput):
        if tap_holder[0] is not None:
            Quartz.CGEventTapEnable(tap_holder[0], True)
        return event
    f = state["cap_file"]
    if f is None: return event
    t = Quartz.CGEventGetTimestamp(event)
    try:
        if etype in (Quartz.kCGEventMouseMoved, Quartz.kCGEventLeftMouseDragged,
                     Quartz.kCGEventRightMouseDragged):
            dx = Quartz.CGEventGetIntegerValueField(event, Quartz.kCGMouseEventDeltaX)
            dy = Quartz.CGEventGetIntegerValueField(event, Quartz.kCGMouseEventDeltaY)
            f.write(f"{t},m,{dx},{dy}\n")
        elif etype == Quartz.kCGEventLeftMouseDown: f.write(f"{t},L,0,0\n")
        elif etype == Quartz.kCGEventLeftMouseUp:   f.write(f"{t},Lu,0,0\n")
        elif etype == Quartz.kCGEventRightMouseDown: f.write(f"{t},R,0,0\n")
        elif etype in (Quartz.kCGEventKeyDown, Quartz.kCGEventKeyUp):
            kc = Quartz.CGEventGetIntegerValueField(event, Quartz.kCGKeyboardEventKeycode)
            k = WASD.get(kc)
            if k:
                f.write(f"{t},{k if etype==Quartz.kCGEventKeyDown else k.upper()+'u'},0,0\n")
        elif etype == Quartz.kCGEventFlagsChanged:
            kc = Quartz.CGEventGetIntegerValueField(event, Quartz.kCGKeyboardEventKeycode)
            name = MODS.get(kc)
            if name:
                flags = Quartz.CGEventGetFlags(event)
                bit = Quartz.kCGEventFlagMaskShift if name=="walkmod" else Quartz.kCGEventFlagMaskControl
                down = bool(flags & bit)
                if down != mod_state[name]:
                    mod_state[name] = down
                    tag = ("W2" if name=="walkmod" else "C2") + ("" if down else "u")
                    f.write(f"{t},{tag},0,0\n")
    except Exception:
        pass
    return event

def tap_thread():
    mask = 0
    for et in (Quartz.kCGEventMouseMoved, Quartz.kCGEventLeftMouseDragged,
               Quartz.kCGEventRightMouseDragged, Quartz.kCGEventLeftMouseDown,
               Quartz.kCGEventLeftMouseUp, Quartz.kCGEventRightMouseDown,
               Quartz.kCGEventKeyDown, Quartz.kCGEventKeyUp,
               Quartz.kCGEventFlagsChanged):
        mask |= Quartz.CGEventMaskBit(et)
    tap = Quartz.CGEventTapCreate(Quartz.kCGHIDEventTap, Quartz.kCGHeadInsertEventTap,
                                  Quartz.kCGEventTapOptionListenOnly, mask, tap_cb, None)
    if tap is None:
        print("TAP_FAILED", flush=True); return
    tap_holder[0] = tap
    Quartz.CGEventTapEnable(tap, False)          # off until a match goes live
    src = Quartz.CFMachPortCreateRunLoopSource(None, tap, 0)
    Quartz.CFRunLoopAddSource(Quartz.CFRunLoopGetCurrent(), src, Quartz.kCFRunLoopCommonModes)
    while True:
        Quartz.CFRunLoopRunInMode(Quartz.kCFRunLoopDefaultMode, 2.0, False)
        with lock:
            f = state["cap_file"]
            if f is not None:
                f.flush()
                if not Quartz.CGEventTapIsEnabled(tap):
                    Quartz.CGEventTapEnable(tap, True)

def capture_start(map_name):
    stamp = time.strftime("%Y%m%d_%H%M")
    state["stamp"] = stamp
    base = os.path.join(SESS, f"{stamp}_{map_name}")
    state["cap_path"] = base + ".cap.csv"
    state["anchor_path"] = base + ".anchor"
    f = open(state["cap_path"], "w", buffering=1)
    f.write(("t_ns" if IS_WIN else "t_ticks") + ",kind,dx,dy\n")
    open(state["anchor_path"], "w").write(str(time.time()))
    with lock:
        state["cap_file"] = f
    tap_set(True)
    accel_off()
    print(f"CAPTURE_START:{map_name}", flush=True)

def capture_stop_and_report():
    with lock:
        f = state["cap_file"]; state["cap_file"] = None
    tap_set(False)
    accel_restore()
    if f: f.close()
    gsi_f = state["gsi_file"]
    if gsi_f: gsi_f.flush()
    if not (state["cap_path"] and state["gsi_path"]): return
    try:
        anchor = float(open(state["anchor_path"]).read().strip())
        out = os.path.join(REPORTS, f"{state['stamp']}_{state['map']}.html")
        report_gen.generate(state["gsi_path"], state["cap_path"], anchor, out)
        print(f"REPORT:{out}", flush=True)
    except Exception as e:
        print(f"REPORT_FAILED:{e}", flush=True)

# ---------- GSI listener ----------
class H(BaseHTTPRequestHandler):
    def do_POST(self):
        body = self.rfile.read(int(self.headers.get('Content-Length', 0)))
        try: g = json.loads(body)
        except Exception: g = {}
        self.send_response(200); self.end_headers()
        handle_gsi(g)
    def log_message(self, *a): pass

def handle_gsi(g):
    m = g.get("map") or {}
    name, mphase, mode = m.get("name"), m.get("phase"), m.get("mode")
    in_map = bool(name)
    # new match: map appears (or changes) while no capture running
    if in_map and state["cap_file"] is None and mphase in ("warmup","live","intermission"):
        state["map"], state["mode"] = name, mode
        base = os.path.join(SESS, f"{time.strftime('%Y%m%d_%H%M')}_{name}")
        state["gsi_path"] = base + ".gsi.jsonl"
        state["gsi_file"] = open(state["gsi_path"], "a", buffering=1)
        capture_start(name)
        state["live_seen"] = False
    if state["gsi_file"] is not None:
        state["gsi_file"].write(json.dumps({"t_recv": time.time(), "gsi": g}) + "\n")
    if mphase == "live": state["live_seen"] = True
    # match end: gameover phase, or we leave the map after having seen live play
    ended = (mphase == "gameover") or (not in_map and state["live_seen"] and state["cap_file"] is not None)
    if ended and state["cap_file"] is not None:
        state["live_seen"] = False
        capture_stop_and_report()
        gf = state["gsi_file"]
        if gf: gf.close()
        state["gsi_file"] = None

if IS_WIN:
    win_input_tap.start(lambda: state["cap_file"])
else:
    threading.Thread(target=tap_thread, daemon=True).start()
server = HTTPServer(("127.0.0.1", 3202), H)
threading.Thread(target=server.serve_forever, daemon=True).start()
print("supervisor up — waiting for CS2", flush=True)
end = time.time() + RUN_HOURS*3600
try:
    while time.time() < end:
        time.sleep(5)
finally:
    if state["cap_file"] is not None:
        capture_stop_and_report()
    accel_restore()
    print("supervisor down", flush=True)
