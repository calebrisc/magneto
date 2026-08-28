# Magneto Val Tracker — setup (Windows)

Privacy-first aim analytics: everything runs and stays on YOUR machine. The
tracker reads your own Riot client session (read-only, localhost — no
injection, no overlays, no game memory; the same category of local-API use
Riot has historically tolerated for community tools).

## What it does

- Detects your matches automatically (Riot client presence) and records your
  raw mouse/keyboard input only while in-game (foreground-gated, chat-gated)
- After each match, fetches YOUR match details from Riot with your own login
  and produces two HTML reports in `reports/`: input metrics (movement
  discipline, flick quality, ability-cast conversion) and round context
  (death taxonomy: traded / outnumbered / anchor / clutch, first bloods)
- `val_review.py` aggregates across matches: verdicts per match
  (carried / standard / did-your-job-lost / no-show), HS%, ACS, DDΔ, RR

## Setup (one time, ~3 minutes)

1. Install Python 3.12 (Microsoft Store or winget: `winget install Python.Python.3.12`)
2. Copy the `cs2-bench` folder anywhere (or `git clone https://github.com/calebrisc/magneto` — it's `tools/cs2-bench`)
3. Sign into VALORANT (Riot client must be running), then:

```
python fetch_binds.py
```

This auto-reads your keybinds + sensitivity from your Riot cloud settings —
no manual configuration. Re-run it whenever you rebind.

## Run (each play session)

```
python val_session.py 12
```

Leave it running. Play. Reports appear in `reports\` on their own after every
match (Riot's indexer can lag a couple minutes). That's the whole workflow.

## Aggregate review (any time)

```
python val_review.py <your-puuid>
```

Your puuid is printed as "player" in `sessions/binds.json`.

## Honest notes

- Unofficial Riot local API: read-only and historically tolerated, but Riot
  can change it any time. Nothing here touches the game process.
- Input capture only records while VALORANT is the foreground window and
  drops keys while in-game chat is open.
- If `fetch_binds.py` warns about "uncaptured ability keys", your binds use
  keys the tap doesn't record yet — casts on those keys won't be tracked
  (everything else still works). Tell Cale; it's a small firmware^W script fix.
