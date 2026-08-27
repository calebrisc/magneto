#!/usr/bin/env python3
"""Valorant session supervisor — presence-polling counterpart of magneto_session.

Valorant has no GSI, so match boundaries come from the Riot Client's local API:
the lockfile (%LOCALAPPDATA%/Riot Games/Riot Client/Config/lockfile) gives
port+password for a localhost HTTPS server, and our own presence's base64
`private` blob carries sessionLoopState (MENUS/PREGAME/INGAME), map, queue and
live score. Read-only localhost polling — no injection, no overlay, no game
memory. Unofficial API: treat every field as best-effort and fail soft.

Arms the same Raw Input tap as the CS2 supervisor, only while INGAME.
Files per match under sessions/:
  <stamp>_val_<map>.cap.csv   input events (t_ns header, same vocabulary)
  <stamp>_val_<map>.anchor    wall clock at capture start
  <stamp>_val_<map>.val.jsonl presence snapshots (~2.5 s cadence: score, state)
Report at match end: input-only metrics (no live kill feed exists; kill-timeline
enrichment comes later from match history). Windows-only.

Usage: val_session.py [hours_to_run]
"""
import base64, html, json, os, ssl, sys, time, urllib.request
import win_input_tap
from round_watch import input_stats

BASE = os.path.dirname(os.path.abspath(__file__))
SESS = os.path.join(BASE, "sessions"); os.makedirs(SESS, exist_ok=True)
REPORTS = os.path.join(BASE, "reports"); os.makedirs(REPORTS, exist_ok=True)
try:
    RUN_HOURS = float(sys.argv[1]) if len(sys.argv) > 1 else 6.0
except ValueError:
    RUN_HOURS = 6.0
LOCKFILE = os.path.join(os.environ.get("LOCALAPPDATA", ""),
                        "Riot Games", "Riot Client", "Config", "lockfile")
POLL_S = 2.5

_ssl = ssl.create_default_context()
_ssl.check_hostname = False
_ssl.verify_mode = ssl.CERT_NONE  # local self-signed cert

state = {"cap_file": None, "cap_path": None, "anchor_path": None,
         "val_file": None, "val_path": None, "stamp": None,
         "map": None, "mode": None, "start_w": None, "last_priv": None}


def read_lockfile():
    try:
        name, pid, port, password, proto = open(LOCKFILE).read().strip().split(":")
        return int(port), password
    except Exception:
        return None


def api_get(port, pw, path):
    req = urllib.request.Request(f"https://127.0.0.1:{port}{path}")
    tok = base64.b64encode(f"riot:{pw}".encode()).decode()
    req.add_header("Authorization", f"Basic {tok}")
    with urllib.request.urlopen(req, timeout=3, context=_ssl) as r:
        return json.loads(r.read())


def own_presence(port, pw, cache):
    if "puuid" not in cache:
        cache["puuid"] = api_get(port, pw, "/chat/v1/session")["puuid"]
    for p in api_get(port, pw, "/chat/v4/presences").get("presences", []):
        # 2026 client: no `product` field, puuid only inside `pid` ("<puuid>@na1...")
        pu = p.get("puuid") or (p.get("pid") or "").split("@")[0]
        if pu != cache["puuid"]:
            continue
        try:
            priv = json.loads(base64.b64decode(p.get("private") or ""))
        except Exception:
            continue
        # 2026 client nests match state under matchPresenceData; flatten so the
        # rest of the code sees one schema (old top-level keys win nothing here)
        mpd = priv.get("matchPresenceData")
        if isinstance(mpd, dict):
            priv = {**priv, **mpd}
        if "sessionLoopState" in priv:
            return priv
    return None


def map_name(priv):
    mm = (priv or {}).get("matchMap") or ""
    return (mm.rsplit("/", 1)[-1] or "unknown").lower()


def mode_name(priv):
    q = (priv or {}).get("queueId") or ""
    if q:
        return q
    return "range" if "range" in map_name(priv) else "custom"


def capture_start(priv):
    stamp = time.strftime("%Y%m%d_%H%M")
    state["stamp"] = stamp
    state["map"], state["mode"] = map_name(priv), mode_name(priv)
    base = os.path.join(SESS, f"{stamp}_val_{state['map']}")
    state["cap_path"] = base + ".cap.csv"
    state["anchor_path"] = base + ".anchor"
    state["val_path"] = base + ".val.jsonl"
    f = open(state["cap_path"], "w", buffering=1)
    f.write("t_ns,kind,dx,dy\n")
    state["start_w"] = time.time()
    open(state["anchor_path"], "w").write(str(state["start_w"]))
    state["val_file"] = open(state["val_path"], "a", buffering=1)
    state["score"] = (0, 0)
    state["cap_file"] = f
    win_input_tap.set_enabled(True)
    print(f"CAPTURE_START:{state['map']}:{state['mode']}", flush=True)


def fmt_stats(s):
    if not s or s.get("shots", 0) == 0:
        return "<p class=note>no shots in window</p>"
    rows = [("shots", s["shots"]),
            ("strafe-release shots measured", s.get("cs_measured", 0)),
            ("fast release&rarr;shot (&lt;60 ms)", s.get("cs_early", 0)),
            ("60–130 ms release&rarr;shot", s.get("cs_inwin", 0)),
            ("counter-strafed shots (CS habit)", s.get("cstrafe_shots", 0)),
            ("&hellip;opposite key still held at shot", s.get("cstrafe_held_at_shot", 0)),
            ("moving shots", s.get("moving_shots", 0)),
            ("sprays (&ge;0.4 s)", s.get("sprays", 0)),
            ("peak flick (counts/s)", s.get("peak_flick_cps", 0)),
            ("total |dx| (counts)", s.get("total_dx", 0))]
    tr = "".join(f"<tr><td>{k}</td><td>{v}</td></tr>" for k, v in rows)
    return f"<table>{tr}</table>"


def write_report(cap_path, anchor, a, b, map_, mode, score, stamp):
    whole = input_stats(cap_path, anchor, a, b)
    mid = (a + b) / 2
    h1 = input_stats(cap_path, anchor, a, mid)
    h2 = input_stats(cap_path, anchor, mid, b)
    out = os.path.join(REPORTS, f"{stamp}_val_{map_}.html")
    dur = (b - a) / 60
    doc = f"""<!doctype html><meta charset=utf-8>
<title>val {html.escape(map_)} {stamp}</title>
<style>body{{font:14px system-ui;margin:2em auto;max-width:640px;color:#222}}
table{{border-collapse:collapse;margin:.5em 0}}td{{border:1px solid #ccc;
padding:3px 10px}}h2{{margin-top:1.4em}}.note{{color:#777}}</style>
<h1>Valorant — {html.escape(map_)}</h1>
<p class=note>{html.escape(mode)} · {dur:.1f} min ·
score {score[0]}–{score[1]} (ally–enemy, last seen)</p>
<h2>Whole match</h2>{fmt_stats(whole)}
<h2>First half</h2>{fmt_stats(h1)}
<h2>Second half</h2>{fmt_stats(h2)}
<p class=note>Input-only metrics (Valorant has no live event feed).
Val stops on key release — fast release&rarr;shot is fine here (unlike CS);
the habit to unlearn is counter-strafing: it's neutral at best, and shots
with the opposite key still held are fired while self-inflicted-inaccurate.
Timing buckets kept CS-identical for cross-game comparison.</p>"""
    open(out, "w", encoding="utf-8").write(doc)
    return out


def capture_stop_and_report():
    win_input_tap.set_enabled(False)
    f = state["cap_file"]; state["cap_file"] = None
    if f: f.close()
    vf = state["val_file"]; state["val_file"] = None
    if vf: vf.close()
    score = state.get("score") or (0, 0)
    try:
        out = write_report(state["cap_path"], state["start_w"], state["start_w"],
                           time.time(), state["map"], state["mode"], score,
                           state["stamp"])
        print(f"REPORT:{out}", flush=True)
    except Exception as e:
        print(f"REPORT_FAILED:{e}", flush=True)
    state["last_priv"] = None


def regen(base):
    """Rebuild the report for an existing capture: regen <sessions/...base>
    (path with or without the .cap.csv suffix)."""
    base = base[:-len(".cap.csv")] if base.endswith(".cap.csv") else base
    anchor = float(open(base + ".anchor").read().strip())
    end = os.path.getmtime(base + ".cap.csv")
    name = os.path.basename(base)                      # <stamp>_val_<map>
    stamp = name.split("_val_")[0]
    map_ = name.split("_val_", 1)[1]
    priv, mode, score = {}, "unknown", (0, 0)
    try:
        for line in open(base + ".val.jsonl"):
            priv = json.loads(line)["private"]
            sc = (priv.get("partyOwnerMatchScoreAllyTeam") or 0,
                  priv.get("partyOwnerMatchScoreEnemyTeam") or 0)
            if sum(sc) > sum(score):
                score = sc
        mode = mode_name(priv)
    except Exception:
        pass
    out = write_report(base + ".cap.csv", anchor, anchor, end, map_, mode,
                       score, stamp)
    print(f"REPORT:{out}", flush=True)


def main():
    win_input_tap.start(lambda: state["cap_file"])
    print("val supervisor up — waiting for Riot client / VALORANT", flush=True)
    end = time.time() + RUN_HOURS * 3600
    seen_lock = False
    cache = {}
    try:
        while time.time() < end:
            lock = read_lockfile()
            if lock is None:
                if state["cap_file"] is not None:
                    capture_stop_and_report()  # client closed mid-match
                if seen_lock:
                    print("riot client gone", flush=True); seen_lock = False
                cache.clear()
                time.sleep(5); continue
            if not seen_lock:
                print("riot client detected", flush=True); seen_lock = True
            try:
                priv = own_presence(*lock, cache)
            except Exception:
                priv = None  # client up but chat not ready / VAL not running
            loop = (priv or {}).get("sessionLoopState")
            if loop == "INGAME" and state["cap_file"] is None:
                capture_start(priv)
            if state["cap_file"] is not None:
                if priv is not None:
                    state["last_priv"] = priv
                    # scores zero out on the end screen — keep the max seen
                    sc = (priv.get("partyOwnerMatchScoreAllyTeam") or 0,
                          priv.get("partyOwnerMatchScoreEnemyTeam") or 0)
                    if sum(sc) > sum(state.get("score") or (0, 0)):
                        state["score"] = sc
                    state["val_file"].write(json.dumps(
                        {"t_recv": time.time(), "private": priv}) + "\n")
                if loop is not None and loop != "INGAME":
                    capture_stop_and_report()
            time.sleep(POLL_S)
    finally:
        if state["cap_file"] is not None:
            capture_stop_and_report()
        print("val supervisor down", flush=True)


if __name__ == "__main__":
    if len(sys.argv) > 2 and sys.argv[1] == "regen":
        regen(sys.argv[2])
    else:
        main()
