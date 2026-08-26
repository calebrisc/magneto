#!/usr/bin/env python3
"""Round watcher for Magneto live coaching.

Tails the newest *.gsi.jsonl in the supervisor's sessions/ dir, tracks one
round, and on round end prints a single `ROUND {json}` line and exits so the
harness wakes the coach. State (file offset, round counter, per-round log)
persists in watch_state.json / rounds.jsonl so relaunching resumes cleanly.

Exit markers on stdout:
  MATCH_START {json}   new capture file seen (map/mode) -- then keeps running
  ROUND {json}         round ended -> exit 0
  MATCH_END {json}     map gameover / capture rotated -> exit 0
"""
import csv, json, math, os, sys, time, bisect, glob

SESS = sys.argv[1] if len(sys.argv) > 1 else None
BASE = os.path.dirname(os.path.abspath(__file__))
STATE_PATH = os.path.join(BASE, "watch_state.json")
ROUNDS_PATH = os.path.join(BASE, "rounds.jsonl")

def load_state():
    try:
        return json.load(open(STATE_PATH))
    except Exception:
        return {"gsi_path": None, "offset": 0, "announced": False,
                "prev": {}, "round_open": None, "totals": {}}

def save_state(st):
    json.dump(st, open(STATE_PATH, "w"))

def newest_gsi():
    files = sorted(glob.glob(os.path.join(SESS, "*.gsi.jsonl")), key=os.path.getmtime)
    return files[-1] if files else None

def input_stats(cap_path, anchor, a, b):
    """Aim/movement stats for wall-clock window [a, b]."""
    clicks, ups, strafe_ups, moves, moving_flags = [], [], [], [], []
    t0 = None
    key_state = {}
    try:
        rd = csv.DictReader(open(cap_path))
        from report_gen import cap_time_col
        tcol, scale = cap_time_col(rd.fieldnames)
        for r in rd:
            try:
                t = int(r[tcol]) * scale
            except Exception:
                continue
            if t0 is None: t0 = t
            w = anchor + (t - t0)
            k = r["kind"]
            if k == "m":
                if a <= w <= b: moves.append((w, abs(int(r["dx"]))))
            elif k == "L":
                if a <= w <= b:
                    clicks.append((w, bool(key_state.get("crouch"))))
                    moving_flags.append(any(key_state.get(x) for x in "wasd"))
            elif k == "Lu":
                if a <= w <= b + 2: ups.append(w)
            elif k in ("a", "d", "w", "s"):
                key_state[k] = True
            elif k in ("Au", "Du", "Wu", "Su"):
                key_state[k[0].lower()] = False
                if k in ("Au", "Du") and a - 1 <= w <= b: strafe_ups.append(w)
            elif k == "C2":
                key_state["crouch"] = True
            elif k == "C2u":
                key_state["crouch"] = False
    except Exception:
        pass
    if not clicks:
        return {"shots": 0}
    crouch_shots = sum(1 for _, cr in clicks if cr)
    delays = []
    for c, cr in clicks:
        if cr: continue  # crouch spray-downs pollute duel-timing stats
        i = bisect.bisect_right(strafe_ups, c) - 1
        if i >= 0 and c - strafe_ups[i] < 0.4:
            delays.append((c - strafe_ups[i]) * 1000)
    inwin = sum(1 for d in delays if 60 <= d < 130)
    early = sum(1 for d in delays if d < 60)
    moving = sum(1 for (_, cr), f in zip(clicks, moving_flags) if f and not cr)
    sprays = 0; longest = 0
    for c, _ in clicks:
        u = [x for x in ups if x > c]
        if u:
            dur = u[0] - c
            if dur >= 0.4:
                sprays += 1
                longest = max(longest, round(dur / 0.1) + 1)
    # peak flick speed: max |dx| summed over 20ms windows (counts/s)
    peak = 0; total_dx = 0
    if moves:
        j = 0
        for i in range(len(moves)):
            total_dx += moves[i][1]
        win = 0.02; j = 0; acc = 0
        for i in range(len(moves)):
            acc += moves[i][1]
            while moves[j][0] < moves[i][0] - win:
                acc -= moves[j][1]; j += 1
            peak = max(peak, acc / win)
    return {"shots": len(clicks), "crouch_shots": crouch_shots,
            "cs_measured": len(delays), "cs_inwin": inwin,
            "cs_early": early, "moving_shots": moving, "sprays": sprays,
            "longest_burst": longest, "peak_flick_cps": int(peak),
            "total_dx": total_dx}

def main():
    st = load_state()
    print("watcher up", flush=True)
    while True:
        path = newest_gsi()
        if path is None:
            time.sleep(2); continue
        if path != st["gsi_path"]:
            st = {"gsi_path": path, "offset": 0, "announced": False,
                  "prev": {}, "round_open": None, "totals": {}}
            save_state(st)
        base = path[:-len(".gsi.jsonl")]
        cap_path, anchor_path = base + ".cap.csv", base + ".anchor"
        try:
            anchor = float(open(anchor_path).read().strip())
        except Exception:
            time.sleep(2); continue

        f = open(path)
        f.seek(st["offset"])
        prev = st["prev"]
        ro = st["round_open"]
        while True:
            line = f.readline()
            if not line or not line.endswith("\n"):
                break  # EOF or partial write; re-read next pass
            st["offset"] = f.tell()  # advance per line so a round-end exit
                                     # never skips the lines after it
            try:
                row = json.loads(line); g = row["gsi"]; t = row["t_recv"]
            except Exception:
                continue
            m = g.get("map") or {}
            rnd = g.get("round") or {}
            p = g.get("player") or {}
            phase = rnd.get("phase")
            mphase = m.get("phase")
            # GSI 'player' follows the spectated teammate after death --
            # only trust player data when it is actually the user.
            own_id = (g.get("provider") or {}).get("steamid")
            is_own = (p.get("steamid") == own_id) if own_id else True
            if is_own:
                st["own"] = {"stt": p.get("state") or {},
                             "ms": p.get("match_stats") or {},
                             "team": p.get("team"),
                             "weapons": p.get("weapons") or {}}
            own = st.get("own") or {}
            stt = own.get("stt") or {}
            ms = own.get("ms") or {}

            if not st["announced"] and m.get("name"):
                st["announced"] = True
                save_state(st)
                print("MATCH_START " + json.dumps(
                    {"map": m.get("name"), "mode": m.get("mode"),
                     "team": p.get("team")}), flush=True)

            # round start
            if phase == "live" and prev.get("phase") != "live":
                weap = ""
                for w in (p.get("weapons") or {}).values():
                    if w.get("state") == "active":
                        weap = w.get("name", "").replace("weapon_", "")
                ro = {"n": (m.get("round") or 0) + 1, "start": t,
                      "money0": prev.get("money", stt.get("money")),
                      "equip": stt.get("equip_value"),
                      "k0": ms.get("kills") or 0, "d0": ms.get("deaths") or 0,
                      "a0": ms.get("assists") or 0,
                      "weapon": weap, "death_t": None, "flash_s": 0.0,
                      "last_t": t}

            if ro is not None:
                if is_own:
                    # active-weapon tracking (last non-knife active gun)
                    for w in (p.get("weapons") or {}).values():
                        if w.get("state") == "active":
                            n = w.get("name", "").replace("weapon_", "")
                            if n not in ("knife", "knife_t", "c4"):
                                ro["weapon"] = n
                    if (stt.get("flashed") or 0) > 128:
                        ro["flash_s"] += max(0.0, t - ro["last_t"])
                # death: own health hits 0, or the feed switches to a teammate
                if ro["death_t"] is None and (
                        (is_own and ((p.get("state") or {}).get("health") or 0) == 0)
                        or (not is_own and phase == "live")):
                    ro["death_t"] = t
                    ro["k_at_death"] = stt.get("round_kills") or 0
                    ro["hs_at_death"] = stt.get("round_killhs") or 0
                ro["last_t"] = t

            # round end
            if ro is not None and phase == "over" and prev.get("phase") == "live":
                end_t = ro["death_t"] or t
                stats = input_stats(cap_path, anchor, ro["start"], end_t)
                death_win = None
                if ro["death_t"] is not None:
                    death_win = input_stats(cap_path, anchor,
                                            ro["death_t"] - 4, ro["death_t"])
                won = rnd.get("win_team") and own.get("team") and \
                      rnd.get("win_team") == own.get("team")
                out = {"round": ro["n"], "map": m.get("name"),
                       "mode": m.get("mode"),
                       "score": {"ct": (m.get("team_ct") or {}).get("score"),
                                 "t": (m.get("team_t") or {}).get("score")},
                       "won": bool(won), "win_team": rnd.get("win_team"),
                       "bomb": rnd.get("bomb"),
                       "kills": (ro.get("k_at_death") if ro["death_t"] is not None
                                 else stt.get("round_kills")) or 0,
                       "hs": (ro.get("hs_at_death") if ro["death_t"] is not None
                              else stt.get("round_killhs")) or 0,
                       "died": ro["death_t"] is not None,
                       "survived_s": round((end_t - ro["start"]), 1),
                       "assists": (ms.get("assists") or 0) - ro["a0"],
                       "money_start": ro["money0"], "equip": ro["equip"],
                       "weapon": ro["weapon"], "flashed_s": round(ro["flash_s"], 1),
                       "kd_total": [ms.get("kills"), ms.get("deaths")],
                       "input": stats, "pre_death_4s": death_win}
                with open(ROUNDS_PATH, "a") as rf:
                    rf.write(json.dumps(out) + "\n")
                st["round_open"] = None
                st["prev"] = {"phase": phase, "money": stt.get("money")}
                save_state(st)
                print("ROUND " + json.dumps(out), flush=True)
                return

            # match end
            if mphase == "gameover" and prev.get("mphase") != "gameover":
                out = {"map": m.get("name"),
                       "score": {"ct": (m.get("team_ct") or {}).get("score"),
                                 "t": (m.get("team_t") or {}).get("score")},
                       "team": own.get("team"),
                       "kd": [ms.get("kills"), ms.get("deaths")],
                       "mvps": ms.get("mvps"), "score_pts": ms.get("score")}
                st["round_open"] = None
                st["prev"] = {"phase": phase, "mphase": mphase}
                save_state(st)
                print("MATCH_END " + json.dumps(out), flush=True)
                return

            prev = {"phase": phase, "mphase": mphase,
                    "money": stt.get("money")}

        f.close()
        st["prev"] = prev
        st["round_open"] = ro
        save_state(st)
        time.sleep(1.0)

if __name__ == "__main__":
    main()
