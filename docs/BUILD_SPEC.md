# Haptic Mouse — Build Spec

**Canonical spec: `docs/spec.html`** (published artifact, v0.4). This file is the changelog + quick summary; the HTML carries the diagrams and full rationale. The teardown reference lives in `docs/anatomy.html`.

## v0.4 — 2026-08-13

- **PAW3950 selected** (0.7 mm LoD, 750 IPS, native 8K); PAW3395 demoted to sourcing fallback. Cost/unit ≈$90–115.
- **LDC2114 disqualified** — datasheet check closed the open item: 160 SPS max scan rate vs 1.5 ms click budget. Inductive fork now **LDC1101 ×2** (~180 kSPS, SPI, single-channel, raw output). Hall's case strengthened; bench-off still decides.
- **8K polling committed** (was "frontier maybe"). Staged: 1 kHz ship → 4 kHz on HS-USB dongle (same ACK protocol fits 250 µs slot) → 8 kHz streaming protocol (no per-packet ACKs, redundant deltas, faster nRF54 PHY). Mouse hardware already qualifies (54L15 + 3950). 8K = per-profile mode; 1 kHz default (battery). New section on Radio tab.

## v0.3 — 2026-08-13

Added three tabs + a handoff block to spec.html: **Parts** (full mouse parts catalog — pick + alternates + qty-6/100k est. pricing, with Superstrike and Apple Magic Mouse USB-C comparison rows per part), **Keyboard** (reference-only sibling-project catalog vs Wooting 80HE — Lekker V2 switches by undisclosed mfr, 8K wired via HS-USB, gasket PC plate), **Radio** (ESB-vs-BLE rationale, full latency chain — USB service interval dominates at 1 kHz, motion sync, retry/channel strategy, antenna options + VNA tuning with hand+shell on, one-dongle-two-devices via ESB pipes, acceptance targets table). Bottom of page: always-visible **continuity brief** so a fresh no-context session can resume the project from the HTML alone.

## v0.2 — 2026-08-13 (decisions from reference-building sessions)

**Scope:** 6 units, not 1. Cost-optimized where it doesn't hurt feel (~$80–100/unit at qty 6 + $10 dongle; program $800–2,000). Open to selling later — way down the line.

**Locked stack:**
- nRF54L15 (Zephyr / nRF Connect SDK), ESB → nRF52840 dongle @ 1 kHz. 4K/8K ceiling is the dongle's USB Full-Speed limit, not the mouse — frontier polling = later High-Speed-USB receiver.
- PAW3395 **module** (matched lens; Z-height is the classic DIY failure). PAW3950 if sourceable.
- Power: 300–400 mAh LiPo, BQ25100, MAX17048, TPS61099 boost.
- Scroll: **electromagnetic detents** (MagSpeed-style) — toothed steel rim + 2 coils; same-polarity = ratchet, opposed = freespin; toggle = current flip, no clutch. Brass-rim flywheel (5–8 g), magnetic encoder, coil braking = software coast curves. Middle click is Hall + haptic too — zero microswitches anywhere.
- Chassis + skin shell architecture: one rigid chassis ×6, per-owner printed skins (clay grip capture → LiDAR → CadQuery parameters). Undercut grips legal in print.

**Bench forks (P1 rig decides, blind press-off + measured attack/latency):**
- Sensing: Hall (TMAG5273) vs inductive (LDC2114). Tie → Hall (transfers to keyboard project).
- Actuator: **piezo favored** (BOS1901/DRV8662-driven) vs LRA (DRV2605L) vs DIY voice coil. Key physics: fingertip band ~50–500 Hz peaks ~250 Hz; LRA "roundness" = resonant rise time, piezo attacks <2 ms and uses 4–10× less power. Distinct release waveform is cheap and probably beats the reference.

**Layout (see §05 of spec.html):** sensor at grip centroid (frozen first), MCU radio-side to the nose, antenna keep-out in the nose (no copper/metal; plate flexures stay aft), HV piezo corner diagonal-opposite the antenna on its own pour (**the one EMI problem this build owns** — spectrum-analyzer check at P3), battery on midline aft of sensor, actuators bonded to plates + soft-mounted, scroll on centerline, side buttons on FFC daughterboard, rear of shell deliberately empty. Weight: ~58 g projected vs 75 g cap.

**Cost levers:** TI samples via myTI (.edu — engineer is a 2nd/3rd-yr student; also: lab access, capstone budget, PCB sponsorship), LCSC over DigiKey, panelize every board order, DIY feet/coils/jigs/CNC metal.

## v0.1 — 2026-08-13

Initial architecture. Superseded: LRA-first haptics (now piezo-favored bench-off), DRV2605L ×2 as the assumed driver, single-unit scope, mechanical-encoder scroll assumption.
