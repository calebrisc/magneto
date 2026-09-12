# m1_logger — bench SAADC logger for the M1 sled

Streams one sin/cos TMR bridge pair (2× AAT003-10E-EVB01) as CSV over the
DK's USB virtual COM port. No amps — bridges drive the SAADC directly
(differential, 40 µs acquisition time for the 40 kΩ source impedance).

## Wiring (one EVB01 board — 4 wires total)

The EVB01's 2×3 header exposes SINGLE-ENDED outputs (per NVE's card:
pin 1 = square pad). Pins 1/6 are both Vdd, 3/4 both GND — use either.

| EVB01 pin       | nRF54L15-DK  | SAADC |
|-----------------|--------------|-------|
| Pin 1 (or 6) Vdd| VDD (3.3 V)  | —     |
| Pin 3 (or 4) GND| GND          | —     |
| Pin 5 Sin       | P1.04        | AIN0  |
| Pin 2 Cos       | P1.06        | AIN2  |

Header layout (component side, pin 1 = square): rows [1 2] / [3 4] / [5 6].
AIN↔P1.xx mapping confirmed from Nordic's adc_dt sample overlay
(AIN0–3 = P1.04–07). Outputs are mid-rail-centered; the overlay uses
gain 1/4 so the 0–3.6 V single-ended range covers them.

## Build & flash

```sh
nrfutil sdk-manager toolchain launch --ncs-version v3.4.0 -- \
  west build -b nrf54l15dk/nrf54l15/cpuapp firmware/m1_logger
nrfutil sdk-manager toolchain launch --ncs-version v3.4.0 -- west flash
```

## Capture

115200 baud on the DK VCOM port:

```sh
screen /dev/tty.usbmodem* 115200        # eyeball it
(stty 115200; cat) < /dev/tty.usbmodem* > run1.csv   # log it
```

Output: `t_ms,ch0_raw,ch0_mv,ch1_raw,ch1_mv`, 200 Hz.

## M1 task #1 — gap characterization

For each shim stack height (0.5 → ~4 mm in steps): set the gap with the
sled procedure, slide the sled slowly one full stroke, log ~5 s of CSV.
Peak-to-peak mv per height should follow e^(−2πz/λ) with λ = 4 mm
(2 mm poles). One channel is enough for this; the second is for atan2
position once amplitude looks sane.
