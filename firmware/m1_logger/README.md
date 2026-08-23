# m1_logger — bench SAADC logger for the M1 sled

Streams one sin/cos TMR bridge pair (2× AAT003-10E-EVB01) as CSV over the
DK's USB virtual COM port. No amps — bridges drive the SAADC directly
(differential, 40 µs acquisition time for the 40 kΩ source impedance).

## Wiring (one EVB01 board — the AAT003 outputs sin AND cos)

| EVB01 pin  | nRF54L15-DK        | SAADC   |
|------------|--------------------|---------|
| V+ / VCC   | VDD (3.3 V)        | —       |
| V− / GND   | GND                | —       |
| SIN+       | P1.04              | AIN0    |
| SIN−       | P1.05              | AIN1    |
| COS+       | P1.06              | AIN2    |
| COS−       | P1.07              | AIN3    |

AIN↔P1.xx mapping confirmed from Nordic's adc_dt sample overlay
(AIN0–3 = P1.04–07, AIN4–7 = P1.11–14). Match the EVB01 end by its
silkscreen labels, not position.

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
