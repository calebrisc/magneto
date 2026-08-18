# Magneto

A ground-up gaming mouse project, built in the open. Two halves:

**Switchless haptic clicks** — no mechanical switches. Analog sensing of finger pressure (LDC1101 inductance-to-digital), click *feel* generated in firmware through a haptic driver (DRV2605L), which means adjustable actuation force, zero debounce, and no switch to wear out. Six hand-fit shells: each owner squeezes a blank in their natural grip, the grip gets scanned, and the scan drives that shell's parameters.

**Magnetic tracking** (in planning) — replacing the optical sensor entirely: a passive magnetically-patterned mousepad read by TMR sensor bridges in the mouse. Three bridge sites in a triangle measure X, Y, and yaw simultaneously; yaw is used to de-rotate motion before HID reporting, so games just see corrected dX/dY. Prior art survey came up empty — as far as I can tell nobody has shipped this.

## Docs

- [`docs/spec.html`](docs/spec.html) — the build spec, 5 tabs, start here
- [`docs/BUILD_SPEC.md`](docs/BUILD_SPEC.md) — same content, markdown
- [`docs/anatomy.html`](docs/anatomy.html) — mouse anatomy reference
- [`docs/PENDING_CHANGES.md`](docs/PENDING_CHANGES.md) — decision queue

Hardware status: spec converged, bench validation phases planned, boards not yet fabbed. This repo is the source of truth as the build progresses.
