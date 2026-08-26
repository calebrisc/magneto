# cs2-bench — $0 aim-analytics rig (Magneto data-mouse groundwork)

Live CS2 aim/movement analytics with no game hooks: an OS-level input tap
(mouse deltas + movement keys) merged with CS2 Game State Integration over
localhost. VAC-clean by construction — nothing reads or touches game memory.

## Pieces

| File | Role |
|---|---|
| `magneto_session.py` | Supervisor. Embeds the GSI listener (port 3202), arms the input tap only while a match is live, disables pointer accel during capture, writes per-match capture files, emits an HTML report at match end. |
| `report_gen.py` | Match report generator (used by the supervisor). |
| `round_watch.py` | Round watcher. Tails the newest capture, emits one `ROUND {json}` line per finished round (aim timing, moving shots, sprays, pre-death snapshot) then exits — built to wake a coaching agent per round. Steamid-gated: ignores spectated-teammate data after death. |
| `window_stats.py` | Ad-hoc stats for the last N minutes (warmup use). |
| `match_coach.py` | Older spoken per-death/round coach (macOS `say`). |
| `gamestate_integration_magneto.cfg` | Drop into `csgo/cfg/` in the CS2 install. |
| `sessions/` | Raw captures: `.cap.csv` (input events, mach ticks), `.gsi.jsonl` (GSI stream), `.anchor` (wall-clock at capture start), `.rounds.jsonl` (per-round metrics). |
| `reports/` | Generated match reports. |

## Data captured 2026-08-23/24

- `20260823_2156_training` — workshop warmup (~5.6 min)
- `20260823_2203_de_inferno` — full competitive match on Inferno, 28+ rounds
  incl. OT (13-15 at capture end; per-round metrics in `.rounds.jsonl`)

## Platform notes

- The GSI listener and all analysis are cross-platform Python.
- The input tap (`Quartz` event tap) and pointer-accel toggle (`defaults
  write com.apple.mouse.scaling`) are **macOS-only**. A Windows port needs a
  Raw Input (or `pynput`) logger writing the same `t_ticks,kind,dx,dy` CSV
  schema and can skip the accel step (set "Enhance pointer precision" off
  manually). Timestamp conversion `t_ticks*125/3/1e9` is mach-specific —
  replace with QPC seconds on Windows and write seconds directly.
- Deps: Python 3.12, `pyobjc-framework-Quartz` (macOS tap only).

## Counter-strafe metric

For each shot, delay since the last A/D key release; 60–130 ms = in-window,
<60 ms = early (fired before the counter-strafe planted). Crouch-held shots
are excluded (spray-down drills pollute duel timing).
