#!/usr/bin/env python3
"""Generate a post-match aim report (HTML) from GSI jsonl + input capture.

Usage: report_gen.py gsi.jsonl capture.csv anchor_wall out.html
Segments the match into rounds via GSI phases, computes per-round and
overall input metrics, peak/worst sections, writes a self-contained page.
"""
import csv, json, math, sys, bisect, statistics as st, html

CPI = 1314.0

def load(gsi_path, cap_path, anchor):
    rows=[json.loads(l) for l in open(gsi_path) if l.strip()]
    moves=[]; clicks=[]; ups=[]; keyups=[]
    t0=None
    for r in csv.DictReader(open(cap_path)):
        t=int(r['t_ticks'])*125/3/1e9
        if t0 is None: t0=t
        w=anchor+(t-t0)
        k=r['kind']
        if k=='m': moves.append((w,int(r['dx']),int(r['dy'])))
        elif k=='L': clicks.append(w)
        elif k=='Lu': ups.append(w)
        elif k in ('Au','Du'): keyups.append(w)
    return rows, moves, clicks, ups, keyups

def segment_rounds(rows):
    rounds=[]; cur=None; prev_phase=None
    meta={'map':'?','mode':'?','final_kills':0,'final_deaths':0,'final_mvps':0,'score':None}
    prev_k=prev_d=None; prev_hs=0
    for r in rows:
        g=r['gsi']; t=r['t_recv']
        m=g.get('map',{}) or {}
        p=g.get('player',{}) or {}
        ms=p.get('match_stats') or {}
        stt=p.get('state') or {}
        phase=(g.get('round',{}) or {}).get('phase')
        if m.get('name'): meta['map']=m['name']
        if m.get('mode'): meta['mode']=m['mode']
        if ms.get('kills') is not None: meta['final_kills']=ms['kills']
        if ms.get('deaths') is not None: meta['final_deaths']=ms['deaths']
        if ms.get('mvps') is not None: meta['final_mvps']=ms['mvps']
        if phase=='live' and prev_phase!='live':
            cur={'start':t,'end':None,'kills':0,'hs':0,'died':False,
                 'num':(m.get('round') or 0)+1, 'kill_times':[]}
        if cur:
            k=ms.get('kills'); d=ms.get('deaths')
            if k is not None and prev_k is not None and k>prev_k:
                cur['kills']+=k-prev_k; cur['kill_times'].append(t)
            if d is not None and prev_d is not None and d>prev_d: cur['died']=True
            rhs=stt.get('round_killhs')
            if rhs is not None: cur['hs']=max(cur['hs'],rhs)
        if phase=='over' and prev_phase=='live' and cur:
            cur['end']=t; rounds.append(cur); cur=None
        if ms.get('kills') is not None: prev_k=ms['kills']
        if ms.get('deaths') is not None: prev_d=ms['deaths']
        if phase is not None: prev_phase=phase
    if cur: cur['end']=rows[-1]['t_recv']; rounds.append(cur)
    return meta, rounds

def input_metrics(a, b, moves, clicks, ups, keyups, times):
    cs=[c for c in clicks if a<=c<=b]
    out={'shots':len(cs),'delays':[],'sprays':[],'flick_peak':0}
    for c in cs:
        i=bisect.bisect_right(keyups,c)-1
        if i>=0 and c-keyups[i]<0.4: out['delays'].append((c-keyups[i])*1000)
        u=[x for x in ups if x>c]
        if u and u[0]-c>=0.4: out['sprays'].append(round((u[0]-c)/0.1)+1)
    i=bisect.bisect_left(times,a); j=bisect.bisect_right(times,b)-1
    w=moves[i:j+1]
    best=0; wi=0
    for wj in range(len(w)):
        while w[wj][0]-w[wi][0]>0.01: wi+=1
        dx=sum(m[1] for m in w[wi:wj+1]); dy=sum(m[2] for m in w[wi:wj+1])
        v=math.hypot(dx,dy)/CPI/0.01
        if v>best: best=v
    out['flick_peak']=best*0.0254  # m/s
    return out

def fmt_ms(x): return f"{x:.0f}&thinsp;ms"

def generate(gsi_path, cap_path, anchor, out_path):
    rows, moves, clicks, ups, keyups = load(gsi_path, cap_path, anchor)
    times=[m[0] for m in moves]
    meta, rounds = segment_rounds(rows)
    per=[]
    for rd in rounds:
        m=input_metrics(rd['start'], rd['end'], moves, clicks, ups, keyups, times)
        stop_kill=None
        for kt in rd['kill_times']:
            i=bisect.bisect_right(keyups, kt)-1
            if i>=0 and kt-keyups[i]<3:
                v=(kt-keyups[i])*1000
                stop_kill=v if stop_kill is None else min(stop_kill,v)
        per.append({**rd, **m, 'stop_kill':stop_kill})

    all_delays=[d for p in per for d in p['delays']]
    all_sprays=[s for p in per for s in p['sprays']]
    tot_kills=sum(p['kills'] for p in per)
    tot_hs=sum(p['hs'] for p in per)
    deaths=sum(1 for p in per if p['died'])
    rushed=sum(1 for d in all_delays if d<60)
    inwin=sum(1 for d in all_delays if 60<=d<130)

    def med(x): return st.median(x) if x else 0
    # peaks / worsts
    kr=[p for p in per if p['kills']]
    peak_round=max(per,key=lambda p:(p['kills'],p['hs']),default=None)
    fastest_kill=min((p['stop_kill'],p['num']) for p in per if p['stop_kill'] is not None) if any(p['stop_kill'] is not None for p in per) else None
    fastest_flick=max((p['flick_peak'],p['num']) for p in per) if per else None
    longest_spray=max(((max(p['sprays']),p['num']) for p in per if p['sprays']), default=None)
    most_rushed=max(per,key=lambda p:sum(1 for d in p['delays'] if d<60),default=None)
    worst_round=max((p for p in per if p['died'] and p['kills']==0),key=lambda p:p['shots'],default=None)
    streak=0; best_streak=0
    for p in per:
        if p['kills'] and not p['died']: streak+=p['kills']; best_streak=max(best_streak,streak)
        elif p['died']: streak=0

    R=[]
    R.append(f"""<title>Aim Report — {html.escape(meta['map'])}</title>
<style>
:root{{--bg:#F4F5F1;--card:#fff;--ink:#15191A;--mut:#77807A;--sig:#0C6E77;--warn:#9A5B08;--rule:#D5D8D0}}
@media (prefers-color-scheme:dark){{:root:not([data-theme=light]){{--bg:#0E1211;--card:#171C1B;--ink:#E7EAE5;--mut:#7C8681;--sig:#4FC3C9;--warn:#E3A94F;--rule:#2A3230}}}}
:root[data-theme=dark]{{--bg:#0E1211;--card:#171C1B;--ink:#E7EAE5;--mut:#7C8681;--sig:#4FC3C9;--warn:#E3A94F;--rule:#2A3230}}
body{{background:var(--bg);color:var(--ink);font:16px/1.55 "Helvetica Neue",system-ui,sans-serif;margin:0;padding:24px 16px 80px}}
.wrap{{max-width:720px;margin:0 auto}} h1{{font-size:30px;letter-spacing:-.02em;margin:.2em 0}}
h2{{font-size:13px;text-transform:uppercase;letter-spacing:.12em;color:var(--mut);margin:34px 0 10px}}
.grid{{display:grid;grid-template-columns:repeat(auto-fit,minmax(140px,1fr));gap:10px}}
.stat{{background:var(--card);border:1px solid var(--rule);border-radius:6px;padding:12px 14px}}
.stat b{{display:block;font-size:24px;font-variant-numeric:tabular-nums}} .stat span{{font-size:11px;text-transform:uppercase;letter-spacing:.08em;color:var(--mut)}}
table{{width:100%;border-collapse:collapse;font-size:14px;font-variant-numeric:tabular-nums}}
th{{text-align:left;font-size:10.5px;text-transform:uppercase;letter-spacing:.1em;color:var(--mut);padding:6px 8px 6px 0;border-bottom:1px solid var(--ink)}}
td{{padding:7px 8px 7px 0;border-bottom:1px solid var(--rule)}}
.good{{color:var(--sig);font-weight:600}} .bad{{color:var(--warn);font-weight:600}}
.tw{{overflow-x:auto}} .note{{color:var(--mut);font-size:13px}}
</style><div class=wrap>
<h1>Aim Report — {html.escape(meta['map'])}</h1>
<p class=note>{html.escape(meta['mode'])} · {len(per)} rounds · {tot_kills}K / {deaths} deaths · {meta['final_mvps']} MVPs</p>
<h2>Session Averages</h2><div class=grid>
<div class=stat><b>{med(all_delays):.0f} ms</b><span>median fire delay</span></div>
<div class=stat><b>{100*inwin/max(len(all_delays),1):.0f}%</b><span>shots in window</span></div>
<div class=stat><b>{100*rushed/max(len(all_delays),1):.0f}%</b><span>shots rushed (&lt;60ms)</span></div>
<div class=stat><b>{100*tot_hs/max(tot_kills,1):.0f}%</b><span>headshot kills</span></div>
<div class=stat><b>{len(all_sprays)}</b><span>sprays ({med(all_sprays):.0f} bullets med)</span></div>
<div class=stat><b>{best_streak}</b><span>best kill run</span></div>
</div>""")
    R.append("<h2>Peak Performance</h2><div class=grid>")
    if fastest_kill: R.append(f"<div class=stat><b class=good>{fastest_kill[0]:.0f} ms</b><span>fastest stop→kill (R{fastest_kill[1]})</span></div>")
    if fastest_flick: R.append(f"<div class=stat><b class=good>{fastest_flick[0]:.2f} m/s</b><span>fastest flick (R{fastest_flick[1]})</span></div>")
    if peak_round: R.append(f"<div class=stat><b class=good>{peak_round['kills']}K</b><span>best round (R{peak_round['num']})</span></div>")
    if longest_spray: R.append(f"<div class=stat><b class=good>{longest_spray[0]}</b><span>longest spray, bullets (R{longest_spray[1]})</span></div>")
    R.append("</div><h2>Needs Work</h2><div class=grid>")
    if most_rushed and any(d<60 for d in most_rushed['delays']):
        n=sum(1 for d in most_rushed['delays'] if d<60)
        R.append(f"<div class=stat><b class=bad>{n} rushed</b><span>worst timing (R{most_rushed['num']})</span></div>")
    if worst_round: R.append(f"<div class=stat><b class=bad>R{worst_round['num']}</b><span>{worst_round['shots']} shots, 0 kills, died</span></div>")
    slow=[d for d in all_delays if d>200]
    R.append(f"<div class=stat><b class=bad>{len(slow)}</b><span>hesitant shots (&gt;200ms)</span></div>")
    R.append("</div>")
    R.append("<h2>Round by Round</h2><div class=tw><table><tr><th>R</th><th>K</th><th>HS</th><th>Died</th><th>Shots</th><th>Fire delay</th><th>Rushed</th><th>Sprays</th><th>Fast kill</th></tr>")
    for p in per:
        dmed=med(p['delays'])
        rn=sum(1 for d in p['delays'] if d<60)
        R.append(f"<tr><td>{p['num']}</td><td>{p['kills'] or '·'}</td><td>{p['hs'] or '·'}</td>"
                 f"<td>{'✕' if p['died'] else '·'}</td><td>{p['shots'] or '·'}</td>"
                 f"<td>{fmt_ms(dmed) if p['delays'] else '·'}</td>"
                 f"<td>{f'<span class=bad>{rn}</span>' if rn else '·'}</td>"
                 f"<td>{len(p['sprays']) or '·'}</td>"
                 f"<td>{fmt_ms(p['stop_kill']) if p['stop_kill'] else '·'}</td></tr>")
    R.append("</table></div>")
    R.append("<p class=note>Window = 60–130 ms after strafe release. Generated locally from your own input journal + Game State Integration; no screen capture, no game memory.</p></div>")
    open(out_path,'w').write('\n'.join(R))
    return out_path, len(per), tot_kills

if __name__=='__main__':
    g,c,a,o=sys.argv[1:5]
    print(generate(g,c,float(a),o))
