#!/usr/bin/env python3
"""Join a Val capture with Riot match-details: per-round + per-fight context.

Usage: val_enrich.py <sessions/...base> <own_puuid>
Expects <base>.match.json (fetched via client tokens), .cap.csv, .anchor.
Prints a JSON summary and writes reports/<base>_rounds.html.

Answers the questions the input tap alone can't: was a bad half game-state
(left in 1vX) or discipline? What does the input look like in the 4 s before
kills vs before deaths?
"""
import html, json, os, sys
from round_watch import input_stats

BASE_DIR = os.path.dirname(os.path.abspath(__file__))
REPORTS = os.path.join(BASE_DIR, "reports")

base = sys.argv[1]
base = base[:-len(".cap.csv")] if base.endswith(".cap.csv") else base
ME = sys.argv[2]
m = json.load(open(base + ".match.json"))
anchor = float(open(base + ".anchor").read().strip())
cap = base + ".cap.csv"

info = m["matchInfo"]
g0 = info["gameStartMillis"] / 1000.0
players = {p["subject"]: p for p in m["players"]}
my_team = players[ME]["teamId"]
mates = {s for s, p in players.items() if p["teamId"] == my_team and s != ME}
foes = {s for s, p in players.items() if p["teamId"] != my_team}

rounds_out = []
kill_windows, death_windows = [], []
for rr in m.get("roundResults", []):
    n = rr.get("roundNum")
    all_kills, seen_victims = [], set()
    for ps in rr.get("playerStats", []):
        for k in ps.get("kills", []):
            v = k.get("victim")
            if v in seen_victims:
                continue  # kills can appear duplicated across playerStats;
                          # one death per player per round in non-respawn modes
            seen_victims.add(v)
            all_kills.append(k)
    all_kills.sort(key=lambda k: k.get("roundTime") or 0)
    my_kills = [k for k in all_kills if k.get("killer") == ME]
    my_death = next((k for k in all_kills if k.get("victim") == ME), None)
    died_vs = None
    if my_death is not None:
        dt = my_death.get("roundTime") or 0
        foes_dead = sum(1 for k in all_kills if k.get("victim") in foes
                        and (k.get("roundTime") or 0) < dt)
        mates_dead = sum(1 for k in all_kills if k.get("victim") in mates
                         and (k.get("roundTime") or 0) < dt)
        died_vs = {"enemies_alive": 5 - foes_dead,
                   "mates_alive": 4 - mates_dead,
                   "at_s": round(dt / 1000, 1)}
        w = g0 + (my_death.get("gameTime") or 0) / 1000.0
        death_windows.append(input_stats(cap, anchor, w - 4, w))
    for k in my_kills:
        w = g0 + (k.get("gameTime") or 0) / 1000.0
        kill_windows.append(input_stats(cap, anchor, w - 4, w))
    fb = all_kills[0] if all_kills else None
    rounds_out.append({
        "round": n + 1 if n is not None else None,
        "won": rr.get("winningTeam") == my_team,
        "plant": rr.get("plantSite") or "",
        "kills": len(my_kills),
        "died": my_death is not None,
        "died_vs": died_vs,
        "first_blood": (fb.get("killer") == ME if fb else False),
        "first_death": (fb.get("victim") == ME if fb else False),
    })


def agg(windows):
    keys = ("shots", "cs_measured", "cstrafe_shots", "cstrafe_held_at_shot",
            "moving_shots")
    tot = {k: sum(w.get(k, 0) for w in windows) for k in keys}
    tot["fights"] = len(windows)
    return tot


# util conversion, bind-aware: an ability EQUIP (per sessions/binds.json) that
# gets CONFIRMED by a fire click within the window — without a cancel key or
# another equip stowing it first — counts as a cast. A cast converts if a
# friendly kill lands within 5 s (mine vs a teammate's). Casts that hit nobody
# still count in the denominator; that's the point. Known miss: abilities that
# activate instantly on the equip key with no confirm click.
def util_conversion():
    try:
        binds = json.load(open(os.path.join(BASE_DIR, "sessions", "binds.json")))
    except Exception:
        binds = {"equips": {"c": "ability", "q": "ability", "e": "ability",
                            "x": "ult"}, "cancels": [], "confirm_window_s": 4}
    equips, cancels = binds["equips"], set(binds.get("cancels", []))
    win = binds.get("confirm_window_s", 4)
    kills_wall = []
    for rr in m.get("roundResults", []):
        for ps in rr.get("playerStats", []):
            for k in ps.get("kills", []):
                w = g0 + (k.get("gameTime") or 0) / 1000.0
                kills_wall.append((w, k.get("killer")))
    kills_wall.sort()
    casts, pending, charging = [], None, None
    # pending = (slot, equip_time); charging = (slot, Ldown_time) — hold-to-
    # charge abilities (Sova darts) fire on RELEASE, so cast time = Lu when
    # one follows within 5 s; instant clicks only shift by ~100 ms, harmless
    try:
        import csv as _csv
        from report_gen import cap_time_col
        rd = _csv.DictReader(open(cap))
        tcol, scale = cap_time_col(rd.fieldnames)
        t0 = None
        for row in rd:
            try:
                t = int(row[tcol]) * scale
            except Exception:
                continue
            if t0 is None:
                t0 = t
            w = anchor + (t - t0)
            kind = row["kind"]
            if pending and w - pending[1] > win:
                pending = None
            if charging and w - charging[1] > 5:
                casts.append((charging[0], charging[1]))  # never released: use Ldown
                charging = None
            if kind in equips:
                pending = (equips[kind], w)  # new equip replaces any pending
            elif kind in cancels:
                pending = charging = None
            elif kind == "R" and pending:
                casts.append((pending[0], w))  # no Ru in captures: R fires now
                pending = None
            elif kind == "L" and pending:
                charging = (pending[0], w)
                pending = None
            elif kind == "Lu" and charging:
                casts.append((charging[0], w))
                charging = None
    except Exception:
        pass
    out = {"casts": len(casts), "converted_self": 0, "converted_team": 0,
           "by_slot": {}}
    for slot, ct in casts:
        out["by_slot"][slot] = out["by_slot"].get(slot, 0) + 1
        window = [k for w, k in kills_wall if ct < w <= ct + 5]
        if any(k == ME for k in window):
            out["converted_self"] += 1
        elif any(k in mates for k in window):
            out["converted_team"] += 1
    out["unconverted"] = out["casts"] - out["converted_self"] - out["converted_team"]
    return out


mystats = players[ME].get("stats") or {}
summary = {
    "map": info.get("mapId", "").rsplit("/", 1)[-1],
    "queue": info.get("queueID"),
    "kda": [mystats.get("kills"), mystats.get("deaths"), mystats.get("assists")],
    "rounds": len(rounds_out),
    "won_rounds": sum(1 for r in rounds_out if r["won"]),
    "first_bloods": sum(1 for r in rounds_out if r["first_blood"]),
    "first_deaths": sum(1 for r in rounds_out if r["first_death"]),
    "deaths_outnumbered": sum(1 for r in rounds_out if r["died_vs"]
                              and r["died_vs"]["enemies_alive"]
                              > r["died_vs"]["mates_alive"] + 1),
    "pre_kill_input": agg(kill_windows),
    "pre_death_input": agg(death_windows),
    "util": util_conversion(),
    "per_round": rounds_out,
}
print(json.dumps(summary, indent=1))

rows = ""
for r in rounds_out:
    dv = r["died_vs"]
    ctx = (f"died vs {dv['enemies_alive']} (with {dv['mates_alive']} mates up, "
           f"{dv['at_s']}s)") if dv else ("survived" if not r["died"] else "died")
    tags = ("FB " if r["first_blood"] else "") + ("FD" if r["first_death"] else "")
    rows += (f"<tr><td>{r['round']}</td><td>{'W' if r['won'] else 'L'}</td>"
             f"<td>{r['kills']}</td><td>{html.escape(ctx)}</td>"
             f"<td>{html.escape(r['plant'])}</td><td>{tags}</td></tr>")
name = os.path.basename(base)
out = os.path.join(REPORTS, name + "_rounds.html")
k, d = summary["pre_kill_input"], summary["pre_death_input"]
doc = f"""<!doctype html><meta charset=utf-8><title>{html.escape(name)} rounds</title>
<style>body{{font:14px system-ui;margin:2em auto;max-width:720px;color:#222}}
table{{border-collapse:collapse}}td,th{{border:1px solid #ccc;padding:3px 9px}}
.note{{color:#777}}</style>
<h1>{html.escape(summary['map'])} &mdash; round context</h1>
<p class=note>KDA {summary['kda'][0]}/{summary['kda'][1]}/{summary['kda'][2]}
&middot; FB {summary['first_bloods']} &middot; first-deaths {summary['first_deaths']}
&middot; outnumbered deaths {summary['deaths_outnumbered']}/{sum(1 for r in rounds_out if r['died'])}</p>
<h2>Input in the 4 s before fights</h2>
<table><tr><th></th><th>fights</th><th>shots</th><th>measured</th>
<th>counter-strafed</th><th>opp held at shot</th><th>moving</th></tr>
<tr><td>before kills</td><td>{k['fights']}</td><td>{k['shots']}</td>
<td>{k['cs_measured']}</td><td>{k['cstrafe_shots']}</td>
<td>{k['cstrafe_held_at_shot']}</td><td>{k['moving_shots']}</td></tr>
<tr><td>before deaths</td><td>{d['fights']}</td><td>{d['shots']}</td>
<td>{d['cs_measured']}</td><td>{d['cstrafe_shots']}</td>
<td>{d['cstrafe_held_at_shot']}</td><td>{d['moving_shots']}</td></tr></table>
<h2>Rounds</h2>
<table><tr><th>#</th><th>W/L</th><th>kills</th><th>death context</th>
<th>plant</th><th></th></tr>{rows}</table>"""
open(out, "w", encoding="utf-8").write(doc)
print(f"REPORT:{out}", file=sys.stderr)
