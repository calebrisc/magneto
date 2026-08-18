# Pending changes — approved but NOT yet on the website

## Queued 2026-08-18 (design audit session — apply to spec.html as v0.5)

- **Two-track structure made explicit:** the six-unit mouse build is the optical (PAW3950) track per v0.4. Magnetic-pad tracking ("Magneto") is a parallel research track with its own go/no-go plan — now on paper in `docs/MAGNETIC_TRACKING.md` (previously existed only in claude.ai planning chats). If its M2 gate passes, it targets mouse v2. spec.html should state this in §02 so "tracking" is unambiguous.
- **Scroll detent revision (approved by Cale 8/18): permanent-magnet-biased detents.** v0.4 design held detents with continuously-driven coils — a constant battery drain in the default mode. Revised: a small permanent magnet provides passive detent flux through the toothed rim (zero idle power, identical feel — same flux-through-teeth physics), coils now only (a) cancel PM flux for freespin (brief, only while spinning), (b) brake for coast curves, (c) trim detent strength per profile. Feel is unchanged; the always-on battery cost moves to the rarely-used mode. Update §04 + anatomy §03.
- **Fork B sourcing note: include a bender/amplified piezo candidate.** A flat piezo patch bonded to the plate produces only microns of displacement — commercial piezo haptics use mechanically amplified configurations (benders/domes). Sourcing only flat discs would rig the bench-off toward the LRA. Add to §10 open items: at least one bender-style piezo in the P1 order.
- **ADC throughput flag (magnetic track):** 3 TMR sites × sin/cos = 6 channels pushes nRF54L15 SAADC near its ceiling at flick-speed sample rates; validation phases dodge this (speed proven on 2ch, geometry on 6ch slow), production answer may be an external ADC. Lives in MAGNETIC_TRACKING.md; spec.html mention only if M2 passes.

## Applied 2026-08-13 (v0.4)
- PAW3950 selected as tracking sensor (fallback 3395 only if module unsourceable)
- Fork A revised: LDC2114 disqualified (160 SPS), LDC1101 ×2 is the inductive contender
- 8K polling committed — staged path (1 kHz ship → 4 kHz on HS dongle, same protocol → 8 kHz streaming); per-profile mode, 1 kHz default; new Radio-tab section
