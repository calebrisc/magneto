#!/usr/bin/env python3
"""Live match coach: spoken rundown on each death and at round end.

Usage: match_coach.py capture.csv anchor_wall.txt gsi_log.jsonl report_out.txt

Tails the GSI log; on death or round phase 'over', computes aim stats for
the just-finished segment from the input capture and speaks a short
summary via `say`. Everything is also appended to report_out.txt.
"""
import csv, json, math, sys, time, bisect, subprocess, os

cap, anchor_path, gsi_path, out_path = sys.argv[1:5]
anchor = float(open(anchor_path).read().strip())
report = open(out_path, "a")

VOICE = os.environ.get("COACH_VOICE") == "1"   # silent by default; novelty on request

def speak(text):
    if VOICE:
        subprocess.Popen(["say", "-r", "210", text])
    report.write(f"[{time.strftime('%H:%M:%S')}] {text}\n")
    report.flush()

def load_inputs():
    clicks=[]; ups=[]; keyups=[]; t0=None
    try:
        rd=csv.DictReader(open(cap))
        from report_gen import cap_time_col
        tcol, scale = cap_time_col(rd.fieldnames)
        for r in rd:
            t=int(r[tcol])*scale
            if t0 is None: t0=t
            w=anchor+(t-t0)
            k=r['kind']
            if k=='L': clicks.append(w)
            elif k=='Lu': ups.append(w)
            elif k in ('Au','Du'): keyups.append(w)
    except Exception:
        pass
    return clicks, ups, keyups

def segment_stats(a, b):
    clicks, ups, keyups = load_inputs()
    cs=[c for c in clicks if a<=c<=b]
    if not cs: return None
    delays=[]
    for c in cs:
        i=bisect.bisect_right(keyups,c)-1
        if i>=0 and c-keyups[i]<0.4:
            delays.append((c-keyups[i])*1000)
    inwin=sum(1 for d in delays if 60<=d<130)
    early=sum(1 for d in delays if d<60)
    sprays=0; longest=0
    for c in cs:
        u=[x for x in ups if x>c]
        if u:
            dur=u[0]-c
            if dur>=0.4:
                sprays+=1; longest=max(longest, round(dur/0.1)+1)
    return dict(shots=len(cs), inwin=inwin, early=early,
                nd=len(delays), sprays=sprays, longest=longest)

def phrase(stats, kills, hs, died):
    if stats is None:
        return "No shots that round." if not died else "Died without shooting."
    bits=[]
    if kills: bits.append(f"{kills} kill{'s'*(kills!=1)}" + (f", {hs} headshot{'s'*(hs!=1)}" if hs else ""))
    else: bits.append("no kills")
    if stats['nd']:
        if stats['early']==0: bits.append("timing clean")
        else: bits.append(f"{stats['early']} of {stats['nd']} shots rushed")
    if stats['sprays']: bits.append(f"{stats['sprays']} spray{'s'*(stats['sprays']!=1)}, longest {stats['longest']} bullets")
    lead="Down. " if died else "Round over. "
    return lead + ", ".join(bits) + "."

# ---- tail the GSI log ----
pos=0
prev_deaths=None
prev_phase=None
seg_start=time.time()
kills_seg_start=0
reported_this_round=False
print("coach armed", flush=True)
while True:
    time.sleep(1.0)
    try:
        with open(gsi_path) as f:
            f.seek(pos)
            new=f.read()
            pos=f.tell()
    except FileNotFoundError:
        continue
    for line in new.splitlines():
        try: row=json.loads(line)
        except Exception: continue
        g=row['gsi']; t=row['t_recv']
        p=g.get('player',{})
        ms=p.get('match_stats') or {}
        rnd=g.get('round',{}) or {}
        stt=p.get('state') or {}
        k,d = ms.get('kills'), ms.get('deaths')
        phase=rnd.get('phase')
        rk=stt.get('round_kills') or 0
        rhs=stt.get('round_killhs') or 0

        died = d is not None and prev_deaths is not None and d>prev_deaths
        over = phase=='over' and prev_phase=='live'

        if died and not reported_this_round:
            stats=segment_stats(seg_start, t)
            speak(phrase(stats, rk-kills_seg_start, rhs, died=True))
            reported_this_round=True
            seg_start=t
        elif over and not reported_this_round:
            stats=segment_stats(seg_start, t)
            speak(phrase(stats, rk, rhs, died=False))
            reported_this_round=True

        if phase=='live' and prev_phase!='live':   # new round started
            seg_start=t
            kills_seg_start=0
            reported_this_round=False

        if d is not None: prev_deaths=d
        if phase is not None: prev_phase=phase
