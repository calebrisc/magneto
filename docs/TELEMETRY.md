# Telemetry & Aim Analytics — architecture

Decided 2026-08-26 (project chat, Windows PC). Companion to the CS2 bench work in
`tools/cs2-bench/`. Covers how aim data is captured, stored, labeled, and viewed
for the production mice. Privacy-first: raw data never leaves the user's hardware
unless they explicitly export it.

## The core decision: the dongle is the recorder

The nRF54L15 in the mouse has **no USB peripheral** — the mouse's USB-C is
charge-only. But the dongle (nRF52840, custom board) already receives every
motion/click report as a side effect of the mouse working. So the archive lives
on the **dongle**, not the mouse:

- **Mouse:** streams reports as it already does. Keeps only a small RAM/NVM ring
  buffer to cover radio dropouts; the dongle requests backfill for gaps. No added
  weight, no battery cost, possibly zero added parts (internal NVM may suffice —
  sizing is an open item).
- **Dongle:** appends the received report stream (+ analog click channels,
  session markers, anchors) to onboard SPI NAND. Always USB-powered whenever the
  mouse is in use, by definition — the recorder cannot miss a session.
- **Web app (Wooting-style):** reads the archive over WebUSB (bulk) / WebHID
  (config), pulls match history from game APIs, aligns and labels client-side.
  No installed software required. Offload is a non-event — the dongle is already
  plugged in; opening the site drains whatever is new.

Rejected alternatives, for the record:
- *NAND in the mouse, offload when plugged in:* impossible as imagined — the
  mouse port is charge-only (no USB on 54L15). Also worse: weight, battery,
  wear-leveling firmware, multi-minute offloads.
- *Required host software as system of record:* works, but loses the
  "no software, mouse just works" story and any match played without the agent
  running. The dongle recorder strictly dominates it.

## Storage sizing

- ~10–20 MB per **active** hour after delta-encoding + idle suppression
  (motion + analog click channels; nothing logged when still). Biggest knob is
  the analog click logging rate — raw waveform vs. per-click summary (open item).
- Dongle NAND: 1–4 Gbit SPI NAND (WSON-8, ~$1–3) holds days-to-weeks as a ring
  buffer. Show a retention gauge in the web app ("storing your last ~N hours")
  so wraparound is never a surprise.
- NAND program bursts are microseconds at ~30 mA — irrelevant on USB power.

## Radio budget

- 1 kHz ACK protocol (ship mode): dongle sees every report by construction;
  archive = log of received stream + backfill for link gaps.
- 8 kHz streaming mode (no per-packet ACKs, redundant deltas, per spec v0.4):
  recorder logs what arrives; losses recovered via backfill from the mouse's
  dropout buffer during idle airtime. Verify backfill airtime fits at 8K (open
  item, alongside the existing 8K protocol work on the Radio tab).

## Game integration

**CS2** — two sources:
- Live: GSI (existing `tools/cs2-bench` pipeline, being ported Mac → Windows).
  Round phases, kills, rich state. Optional — enriches, not load-bearing.
- Retroactive: match history APIs for segmentation when GSI wasn't running.

**Valorant** — no GSI equivalent; two sources:
- Retroactive (primary): match-history APIs (Riot official or HenrikDev-style
  proxy) fetched by the web app itself — match start/end, per-kill timestamps
  and positions. Enough to segment the archive into matches and fights after
  the fact. No local software needed.
- Live (optional, needs a native agent): Riot Client local API — lockfile at
  `%LocalAppData%\Riot Games\Riot Client\Config\lockfile`, presence endpoint's
  `sessionLoopState` (MENUS/PREGAME/INGAME) + map/mode/score. Read-only
  localhost polling — the category Riot has historically tolerated (no
  injection/overlay/memory reads) — but unofficial and breakable; never
  load-bearing. Docs: valapidocs.techchrism.me.
- Mouse-motion heuristics ("is he gaming?") considered and rejected as a
  trigger: fuzzy boundaries, no context, a classifier to maintain. Kept only as
  a last-ditch fallback concept if all APIs vanish.

Both games share one adapter shape: emit match-state events into the label
stream; CS2-GSI and Val-presence are two adapters on the same schema.

## Clock alignment (three layers)

The dongle has no RTC; timestamps are monotonic-since-power. Alignment to wall
time (and thence to match timelines):

1. **Anchor on host contact** — whenever the web app (or any tool) talks to the
   dongle, record (dongle tick ↔ host wall time). Same pattern as the
   `.anchor` files already produced by `tools/cs2-bench`. Anchors apply
   retroactively across the whole monotonic span (i.e., back to dongle
   power-on), so a weekly site visit anchors everything since the last PC boot.
2. **Session markers** — logged at dongle power-up and mouse connect/disconnect,
   so spans across power cycles stay distinguishable.
3. **Kill-feed correlation** — per-match refinement: logged click bursts
   correlate tightly with the match's kill timestamps; snaps residual drift
   (crystal ~seconds/day) to well under 100 ms.

Edge case to bench-verify: a span with no anchor at all (dongle power-cycled,
site never opened) must be recoverable by layer 3 alone. Expected to work; not
yet proven.

## Charge passthrough (dongle/puck)

Problem: Razer-style habit is one cable — unplug the dongle, plug the mouse in
to charge. On our stack that kills the mouse (no wired mode, see below).

Fix: the dongle (or its desk puck — the Radio tab already wants the receiver
near the pad) carries a **USB-C power passthrough port**. The charge cable
plugs into the dongle, not a second PC port. One-cable habit preserved; dongle
never unplugged; mouse plays wirelessly while charging; recording never stops.
(Contrast: Razer's wired mode silently bypasses any dongle-side recording.)

- Straight 5 V VBUS passthrough — no data, no hub silicon, no PD/QC.
- Rp pull-ups on the output advertise 1.5 A; charges any basic 5 V USB device
  (phones charge slowly — fine).
- **Current limit / polyfuse on the passthrough rail** so a greedy or faulty
  device browns out the passthrough port, never the radio. The mouse must keep
  working no matter what's plugged into the charge port.

Why nobody ships this: mainstream mice have USB in the mouse (wired mode) so
the problem doesn't exist; certified products won't daisy-chain draw past USB
2.0's guaranteed 500 mA; nano-receiver form factor. None of these bind a
6-unit build with a custom puck.

## Accepted limitation: no wired data mode

The nRF54L15 has no USB — with the dongle lost/broken, the mouse is inert.
Consciously accepted 2026-08-26:

- Mitigations: charge passthrough (above) removes the reason to unplug the
  dongle; stock a spare dongle per unit (~$10).
- Would require different silicon (nRF54H20, or a USB bridge chip in the
  mouse) — decide **before** mouse PCB layout if ever revisited.
- Revisit trigger: tournament/LAN rules requiring wired mice.

## Open engineering items

1. Dongle NAND part selection; add to dongle board design (Parts tab).
2. Log format spec: record types, session markers, anchor records, backfill
   protocol framing.
3. Mouse-side dropout buffer sizing (longest radio gap to survive; internal
   NVM vs. small external flash).
4. Analog click logging rate: full waveform vs. per-click feature summary.
5. Valorant adapter in the web app (match-history fetch + alignment), sibling
   to the CS2 GSI adapter; finish the GSI Windows port first.
6. Backfill airtime check in 8K streaming mode.
7. Bench-verify layer-3-only clock recovery (no-anchor edge case).
