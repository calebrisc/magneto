#!/usr/bin/env python3
"""Cross-match review over downloaded match-details (sessions/valmatches/*.json).

Per death, classifies context the input tap can't see:
  FD        first death of the round (split by attack/defense)
  traded    my killer died within 5 s of my death
  outnum    enemies alive > mates alive + 1 when I died
  anchor    died to a 4-5 man hit while 3+ mates were still alive elsewhere
  clutch    died with 0 mates alive (1vX attempt)
Side (attack/defense) is derived from who planted (planter team = attackers),
falling back to regulation-half inference; OT uses second-half parity.

Usage: val_review.py <puuid>   (prints JSON: per-match rows + aggregates)
"""
import glob, json, os, sys, urllib.request

BASE = os.path.dirname(os.path.abspath(__file__))
MDIR = os.path.join(BASE, "sessions", "valmatches")
ME = sys.argv[1]
TRADE_S = 5.0

# agent uuid -> (name, role) via community mirror; cached next to the matches
AG_CACHE = os.path.join(MDIR, "_agents.json")
try:
    agents = json.load(open(AG_CACHE))
except Exception:
    raw = json.loads(urllib.request.urlopen(
        "https://valorant-api.com/v1/agents?isPlayableCharacter=true",
        timeout=10).read())["data"]
    agents = {a["uuid"].lower(): [a["displayName"], a["role"]["displayName"]]
              for a in raw}
    json.dump(agents, open(AG_CACHE, "w"))

# mapId -> display name (Riot uses internal codenames: Bonsai=Split, Triad=Haven)
MAP_CACHE = os.path.join(MDIR, "_maps.json")
try:
    map_names = json.load(open(MAP_CACHE))
except Exception:
    raw = json.loads(urllib.request.urlopen(
        "https://valorant-api.com/v1/maps", timeout=10).read())["data"]
    map_names = {mm["mapUrl"]: mm["displayName"] for mm in raw if mm.get("mapUrl")}
    json.dump(map_names, open(MAP_CACHE, "w"))


def map_display(map_id):
    return map_names.get(map_id, (map_id or "").rsplit("/", 1)[-1])

matches = []
for path in glob.glob(os.path.join(MDIR, "*.json")):
    if os.path.basename(path).startswith("_"):
        continue
    m = json.load(open(path))
    if not isinstance(m, dict) or "matchInfo" not in m:
        continue
    matches.append(m)
matches.sort(key=lambda m: m["matchInfo"].get("gameStartMillis") or 0)

rows, deaths_all = [], []
for m in matches:
    info = m["matchInfo"]
    players = {p["subject"]: p for p in m["players"]}
    if ME not in players:
        continue
    me = players[ME]
    my_team = me["teamId"]
    mates = {s for s, p in players.items() if p["teamId"] == my_team and s != ME}
    foes = {s for s, p in players.items() if p["teamId"] != my_team}
    team_row = next((t for t in m.get("teams", []) if t["teamId"] == my_team), {})
    ag = agents.get((me.get("characterId") or "").lower(), ["?", "?"])

    # who attacks first: majority vote of planter teams in regulation rounds
    votes = {"mine": 0, "theirs": 0}
    for rr in m.get("roundResults", []):
        n = rr.get("roundNum", 99)
        pl = rr.get("bombPlanter")
        if pl and n < 24:
            first_half = n < 12
            planter_mine = players.get(pl, {}).get("teamId") == my_team
            attacking_first = planter_mine if first_half else not planter_mine
            votes["mine" if attacking_first else "theirs"] += 1
    my_attack_first = votes["mine"] >= votes["theirs"]

    def my_side(n, rr):
        pl = rr.get("bombPlanter")
        if pl:
            return "attack" if players.get(pl, {}).get("teamId") == my_team \
                else "defense"
        if n < 12:
            return "attack" if my_attack_first else "defense"
        if n < 24:
            return "defense" if my_attack_first else "attack"
        flip = (n - 24) % 2 == 1
        second = "defense" if my_attack_first else "attack"
        return ("attack" if second == "defense" else "defense") if flip else second

    st = me.get("stats") or {}
    r = {"map": map_display(info.get("mapId", "")),
         "start": info.get("gameStartMillis"),
         "agent": ag[0], "role": ag[1],
         "won": bool(team_row.get("won")),
         "score": f"{team_row.get('roundsWon')}-"
                  f"{(team_row.get('roundsPlayed') or 0) - (team_row.get('roundsWon') or 0)}",
         "k": st.get("kills"), "d": st.get("deaths"), "a": st.get("assists"),
         "rounds": st.get("roundsPlayed"),
         "dmg": 0, "fb": 0, "fd": 0, "fd_attack": 0, "fd_traded": 0,
         "deaths_traded": 0, "deaths_outnum": 0, "deaths_anchor": 0,
         "deaths_clutch": 0, "kast": 0, "multi": 0}

    for rr in m.get("roundResults", []):
        n = rr.get("roundNum", 0)
        side = my_side(n, rr)
        my_ps = next((p for p in rr.get("playerStats", [])
                      if p["subject"] == ME), None)
        if my_ps:
            r["dmg"] += sum(d.get("damage", 0) for d in my_ps.get("damage", []))
        all_kills, seen = [], set()
        for ps in rr.get("playerStats", []):
            for k in ps.get("kills", []):
                if k.get("victim") in seen:
                    continue
                seen.add(k.get("victim"))
                all_kills.append(k)
        all_kills.sort(key=lambda k: k.get("roundTime") or 0)
        my_kills = [k for k in all_kills if k.get("killer") == ME]
        if len(my_kills) >= 3:
            r["multi"] += 1
        my_death = next((k for k in all_kills if k.get("victim") == ME), None)
        fb = all_kills[0] if all_kills else None
        if fb and fb.get("killer") == ME:
            r["fb"] += 1
        participated = bool(my_kills) or my_death is None
        if my_ps and not participated:
            participated = any(k.get("assistants") and ME in k["assistants"]
                               for k in all_kills)
        if my_death is not None:
            dt = my_death.get("roundTime") or 0
            foes_dead = sum(1 for k in all_kills if k.get("victim") in foes
                            and (k.get("roundTime") or 0) < dt)
            mates_dead = sum(1 for k in all_kills if k.get("victim") in mates
                             and (k.get("roundTime") or 0) < dt)
            ena, maa = 5 - foes_dead, 4 - mates_dead
            traded = any(k.get("victim") == my_death.get("killer")
                         and 0 <= (k.get("roundTime") or 0) - dt <= TRADE_S * 1000
                         for k in all_kills)
            is_fd = fb is not None and fb.get("victim") == ME
            d = {"map": r["map"], "agent": ag[0], "role": ag[1], "round": n + 1,
                 "side": side, "at_s": round(dt / 1000, 1), "fd": is_fd,
                 "traded": traded, "enemies_alive": ena, "mates_alive": maa,
                 "outnum": ena > maa + 1,
                 "anchor": ena >= 4 and maa >= 3,
                 "clutch": maa == 0,
                 "kills_before": len([k for k in my_kills
                                      if (k.get("roundTime") or 0) < dt])}
            deaths_all.append(d)
            if is_fd:
                r["fd"] += 1
                if side == "attack":
                    r["fd_attack"] += 1
                if traded:
                    r["fd_traded"] += 1
            if traded:
                r["deaths_traded"] += 1
                participated = True
            if d["outnum"]:
                r["deaths_outnum"] += 1
            if d["anchor"]:
                r["deaths_anchor"] += 1
            if d["clutch"]:
                r["deaths_clutch"] += 1
        if participated:
            r["kast"] += 1
    r["adr"] = round(r["dmg"] / max(1, r["rounds"] or 0))
    r["kast_pct"] = round(100 * r["kast"] / max(1, r["rounds"] or 0))
    # match verdict: separates "you underperformed" from "team lost around you"
    if r["adr"] <= 70 or r["kast_pct"] <= 55:
        r["verdict"] = "no-show"            # you weren't in this one, W or L
    elif not r["won"] and (r["adr"] >= 120 or r["kast_pct"] >= 75):
        r["verdict"] = "did-your-job-lost"  # performance there, result wasn't
    elif r["won"] and r["adr"] >= 120:
        r["verdict"] = "carried"
    else:
        r["verdict"] = "standard"
    rows.append(r)

# session position: nth comp game that calendar day (fatigue/tilt signal)
import datetime
day_counts = {}
for r in rows:
    day = datetime.datetime.fromtimestamp((r["start"] or 0) / 1000).strftime("%Y%m%d")
    day_counts[day] = day_counts.get(day, 0) + 1
    r["session_game"] = day_counts[day]


def pct(a, b):
    return round(100 * a / b) if b else None


tot_d = len(deaths_all)
by_role = {}
for r in rows:
    b = by_role.setdefault(r["role"], {"m": 0, "w": 0, "k": 0, "d": 0,
                                       "fb": 0, "fd": 0, "adr": []})
    b["m"] += 1; b["w"] += r["won"]; b["k"] += r["k"] or 0
    b["d"] += r["d"] or 0; b["fb"] += r["fb"]; b["fd"] += r["fd"]
    b["adr"].append(r["adr"])
for b in by_role.values():
    b["adr"] = round(sum(b["adr"]) / len(b["adr"]))

agg = {
    "matches": len(rows),
    "won": sum(1 for r in rows if r["won"]),
    "kd": [sum(r["k"] or 0 for r in rows), sum(r["d"] or 0 for r in rows)],
    "total_deaths_classified": tot_d,
    "deaths_traded_pct": pct(sum(1 for d in deaths_all if d["traded"]), tot_d),
    "deaths_outnum_pct": pct(sum(1 for d in deaths_all if d["outnum"]), tot_d),
    "deaths_anchor_pct": pct(sum(1 for d in deaths_all if d["anchor"]), tot_d),
    "deaths_clutch_pct": pct(sum(1 for d in deaths_all if d["clutch"]), tot_d),
    "fd_total": sum(r["fd"] for r in rows),
    "fd_on_defense": sum(1 for d in deaths_all if d["fd"] and d["side"] == "defense"),
    "fd_on_attack": sum(1 for d in deaths_all if d["fd"] and d["side"] == "attack"),
    "fd_traded_pct": pct(sum(1 for d in deaths_all if d["fd"] and d["traded"]),
                         sum(1 for d in deaths_all if d["fd"])),
    "fb_total": sum(r["fb"] for r in rows),
    "deaths_by_side": {
        "attack": sum(1 for d in deaths_all if d["side"] == "attack"),
        "defense": sum(1 for d in deaths_all if d["side"] == "defense")},
    "untraded_even_deaths": sum(1 for d in deaths_all
                                if not d["traded"] and not d["outnum"]
                                and not d["clutch"]),
    "by_role": by_role,
}
print(json.dumps({"matches": rows, "deaths": deaths_all, "agg": agg}, indent=1))
