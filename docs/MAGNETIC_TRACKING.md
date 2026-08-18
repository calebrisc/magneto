# Magnetic Tracking ("Magneto") — Design Record + Go/No-Go Plan

**Status: research track, NOT part of mouse v1.** The six-unit mouse build uses the PAW3950 optical module per spec.html v0.4 — that decision stands regardless of what happens here. This track answers one question: *can a passive magnetically-patterned pad + TMR sensors in the mouse replace the optical sensor in a future version?* Until 2026-08-18 this design lived only in claude.ai planning chats; this file makes it part of the repo so both projects share one source of truth.

## The architecture (as converged in planning)

- **Pad:** passive, magnetically patterned. Inverted stack: foam base → magnetic sheet → thin textile skin, so the compliant layer sits *below* the reference surface and the sensor rides close to the magnetic layer.
- **Mouse:** three TMR (tunnel magnetoresistance) sensor bridge sites in a triangle, measuring X, Y, and yaw simultaneously. Yaw is used internally to de-rotate motion samples before HID reporting — the host sees ordinary corrected dX/dY, no protocol changes.
- **Pole pitch:** 4 mm design center.
- **Signal path:** raw analog TMR bridges → amplification → MCU ADC → firmware atan2 interpolation. Deliberately bypasses integrated encoder ICs (AS5311-class) whose internal pipelines cap speed around ~650 mm/s; raw-bridge + firmware math has 10–20× margin over human flick speeds.
- **Pad manufacturing:** magnetization write jig built on the existing 3D printer gantry (electromagnet write head, coil + steel core, gantry provides positioning).
- **Prior art:** survey found no existing products, patents, or published research on this exact combination (passive patterned pad + in-mouse TMR read). Note: the repo went public 2026-08-18, which starts the US 1-year patent filing clock and forfeits most foreign filing rights.

## Independent physics audit (2026-08-18)

Checked from first principles, not from the planning chats:

- **Air gap is the dominant constraint.** Field from a periodic magnetization pattern decays exponentially with height above the surface: amplitude ∝ e^(−2πz/λ) where λ is the spatial period. At 4 mm pitch, roughly **half the signal survives a 1 mm gap and only ~20% survives 2 mm**. Consequences: sensors must sit at the very bottom of the mouse, feet must be thin, and the pad's textile skin must be thin. The inverted pad stack exists exactly for this reason. VERDICT: works, but every mechanical decision near the mouse floor must protect the gap.
- **Speed / signal frequency: holds.** At 750 IPS (~19 m/s) over a 4 mm pitch the bridge signal is in the low-kHz range (~2.5–5 kHz depending on whether you count poles or pole-pairs). Firmware atan2 at MCU speeds has enormous margin. The claimed ~1 µm theoretical resolution from 12-bit interpolation over a 4 mm period is arithmetic-true; realistic noise-limited resolution of 10–30 µm still equals or beats ~1000–2500 DPI optical. VERDICT: not the risk.
- **ADC throughput: flag for the engineer.** Three sensor sites × sin/cos = 6 analog channels. At full flick speed with meaningful oversampling this pushes the nRF54L15 SAADC near its ceiling. Likely fine for validation (see phase split below — speed is proven on 2 channels, geometry on 6 slow channels); the production answer may be an external ADC. Board-design input, not a blocker.
- **THE open question: 2D tracking.** Stripes solve one axis beautifully. Reading X *and* Y from a single surface requires a 2D magnetization pattern (checkerboard or crossed gratings), and the two axes' fields superpose at the sensor — separating them cleanly is the unproven part. Yaw compounds it. **This is the make-or-break, and nothing short of a bench test answers it.** Everything in the plan below is sequenced to reach this question for the least money and time.

## Go/No-Go validation plan (phases M0–M2)

Decision rule: if M2 demonstrates clean 2D+yaw tracking on a written checkerboard pad, magnetic tracking becomes the planned sensor for mouse v2 and gets a full spec section. If M2 fails after honest effort, the concept is shelved with the findings written up, and the optical path continues unaffected. Either outcome is a success for the *project* — the point is to spend ~$200, not $2,000, learning which world we're in.

### M0 — 1D geometry with off-the-shelf parts (one weekend after parts arrive)

AS5311 eval board + its matched 2 mm-pole magnetic strip, read by the nRF54L15-DK. Proves: the air-gap sensitivity curve in real life, interpolation behavior, and the integrated IC's speed ceiling (which motivates M1). Zero custom anything. **You can do this phase yourself** — it's wiring a dev kit to an eval board and logging numbers.

### M1 — raw TMR + firmware atan2, 1D speed proof (1–2 weekends)

One TMR bridge pair over the same strip plus a cheap multipole magnet sheet (ordinary extruded "fridge magnet" sheeting is factory-striped at a few-mm pitch — verify pitch with viewing film). Bridge outputs amplified (instrumentation amps) into the DK's ADC; firmware computes atan2 position. Validation: hand-flung 3D-printed sled with an optical mouse sensor module riding alongside as ground truth. Success = tracking stays locked at flick speeds the AS5311 couldn't hold. Only 2 ADC channels needed — this phase proves *speed*, deliberately on 1D.

### M2 — 2D + yaw go/no-go (2–4 weekends; the whole reason this track exists)

Build the write jig (electromagnet head on the 3D printer gantry), write a checkerboard/crossed pattern onto blank magnet sheet, mount three bridge sites in a triangle, and attempt real 2D tracking with yaw de-rotation. 6 ADC channels at modest sample rates — 2D *geometry* can be validated at low speed because M1 already proved speed. Quick-and-dirty precursor experiment (near-free): stack two striped sheets orthogonally and see whether the two axes are separable at the sensor at all — messy fields, but a cheap early read on the separation problem before the jig exists.

## Purchase list — everything for M0 through M2

| # | Item | Purpose | Qty | Est. | Source | Notes |
|---|------|---------|-----|------|--------|-------|
| 1 | nRF54L15-DK | MCU/ADC for all phases; shared with mouse P0 | 1 | $55 | Nordic / DigiKey | Buys into the mouse build too — not a pad-only cost |
| 2 | AS5311 adapterboard (AS5311-TS_EK_AB) | M0 integrated-IC baseline | 1 | $20 | DigiKey / Mouser | |
| 3 | Multipole magnetic strip for AS5311 (MS10-class, 2 mm poles) | M0/M1 known-good target | 1–2 | $25 | ams-OSRAM channel / DigiKey | |
| 4 | TMR bridge sensors — NVE AAT003-10E (or AAT001) | The actual sensor candidate | 6 | $60 | DigiKey | **Engineer verifies part + stock first** — need analog sin/cos bridge output; if NVE unavailable, any raw-output TMR/AMR bridge (Sensitec, MDT) substitutes |
| 5 | Instrumentation amps (INA333 or AD8226-class) | Bridge mV output → ADC range | 8 | $30 | TI / DigiKey | **TI part — myTI sampleable, likely $0** |
| 6 | Flexible magnet sheet, plain + adhesive-back roll | M1 cheap target, M2 write-jig blanks | 2 rolls | $20 | Amazon / craft store | Extruded sheet is factory-striped — check with viewing film |
| 7 | Magnetic viewing film | See pole patterns; the essential diagnostic | 1 | $8 | Amazon | Buy first, it makes everything else debuggable |
| 8 | Optical mouse sensor module (PMW3389-class) | Ground-truth comparator on the sled | 1 | $12 | AliExpress / mouse-mod shops | 3389 is fine for truth; don't burn a $30 PAW3950 here |
| 9 | Write-jig parts: magnet wire, steel bolt core, MOSFET, flyback diode, misc | M2 pad writing | 1 set | $20 | Amazon / existing stock | Magnet wire doubles for scroll-coil experiments later |
| 10 | Pad stack materials: EVA foam sheet, thin textile, spray adhesive | M2 inverted-stack prototype | 1 set | $15 | craft store | |
| 11 | Proto supplies: perfboard, headers, hookup wire, PTFE tape for sled feet | all phases | — | $15 | Amazon | |

**Total: ~$280, or ~$225 after the TI in-amp samples — and $55 of that is the DK the mouse build needs anyway.** Longest lead items: TMR sensors (order first) and the AS5311 board.

## Relationship to the mouse build

Nothing here blocks or changes mouse P0–P6. The DK is shared. If M2 succeeds, magnetic tracking targets mouse v2 with a proper spec section, the ADC question gets a real answer, and the pad becomes its own product line. If it fails, the writeup goes in this file and the optical mouse never notices.
