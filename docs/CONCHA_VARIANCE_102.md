# Concha variance, 102 real ears (2026-09-03)

gyro_arm_variance.py over every aligned real ear except P0023 (scanner hole-fill concha, excluded). Skirt-datum frame: origin = aperture centroid, +Z = canal axis. Feeds the cymba-fold clamp: body edge under the crus, one magnetic pad expanding into the cymba.

```
102 ears measured; 0 synthetic corners excluded from the statistics (parametric bowls, not anatomy)

n = 102 ears, skirt-datum frame (origin = aperture centroid, +Z = canal axis)
NOTE: this is the 13-ear short list; the same script runs over all ~103 aligned ears unchanged.

measured rim radius: mean 9.53, SD 1.69 mm
size proxy (basin inscribed radius, independent): mean 7.29, SD 0.86, range 5.00–9.01 mm (12 % CV)

--- RAW (mm) ---
| feature | mean | SD | range | n |
|---|---|---|---|---|
| cymba floor, distance from origin | 7.97 | 4.25 | 4.02–15.83 | 100 |
| cymba pocket depth | 8.82 | 3.30 | 2.47–17.76 | 100 |
| cymba lip centroid, distance | 15.00 | 2.92 | 5.82–18.97 | 86 |
| antihelix lip, closest point | 10.16 | 4.46 | 4.06–18.59 | 61 |
| antihelix arc radius | 2.93 | 1.97 | 0.37–7.32 | 61 |
| antihelix arc extent (deg) | 238.20 | 91.69 | 68.05–358.14 | 61 |
| crus helicis, distance from origin | 14.46 | 1.56 | 11.68–18.39 | 41 |

--- SCALE-NORMALISED (x / ear scale) ---
| feature | mean | SD | range | n |
|---|---|---|---|---|
| cymba floor, distance from origin | 1.10 | 0.58 | 0.48–2.28 | 100 |
| cymba pocket depth | 1.22 | 0.47 | 0.38–2.28 | 100 |
| cymba lip centroid, distance | 2.11 | 0.48 | 0.83–3.42 | 86 |
| antihelix lip, closest point | 1.42 | 0.69 | 0.48–3.21 | 61 |
| antihelix arc radius | 0.41 | 0.29 | 0.05–1.24 | 61 |
| antihelix arc extent (deg) | 238.20 | 91.69 | 68.05–358.14 | 61 |
| crus helicis, distance from origin | 1.99 | 0.26 | 1.45–2.63 | 41 |

--- does normalising by size remove the variance? (CV = SD/mean) ---
| feature | CV raw | CV scaled | verdict |
|---|---|---|---|
| cymba floor dist | 0.53 | 0.53 | mixed |
| cymba pocket depth | 0.37 | 0.39 | mixed |
| cymba lip dist | 0.19 | 0.23 | shape |
| antihelix closest | 0.44 | 0.48 | mixed |
| antihelix arc radius | 0.67 | 0.70 | mixed |
| crus dist | 0.11 | 0.13 | shape |

--- concha rim path, point-wise across-ear SD (RAW) ---
  overall: mean 2.32, min 1.74, max 2.96 mm
  azimuths with SD < 1.0: 0/36
    anterior/crus   mean SD 1.99
    cymba           mean SD 2.37
    antihelix       mean SD 2.80
    posterior       mean SD 2.77

--- concha rim path, point-wise across-ear SD (SCALE-NORMALISED) ---
  overall: mean 0.35, min 0.25, max 0.45
  azimuths with SD < 1.0: 36/36  -> -180deg:0.42, -170deg:0.42, -160deg:0.42, -150deg:0.37, -140deg:0.37, -130deg:0.37, -120deg:0.38, -110deg:0.36, -100deg:0.33, -90deg:0.30
    anterior/crus   mean SD 0.29
    cymba           mean SD 0.37
    antihelix       mean SD 0.42
    posterior       mean SD 0.41

```

## Read-out for the clamp
- Crus helicis distance 14.46 ± 1.56 mm, range 11.7–18.4, CV 0.11 (n=41 of 102 detected). Tightest feature in the ear; the fixed jaw references it.
- Cymba lip (helix-root overhang) centroid 15.0 ± 2.9 mm (n=86) — sits just above the crus, so the pad that pinches the crus also lands under the overhang.
- Cymba pocket depth 8.8 ± 3.3, range 2.5–17.8 (n=100): do NOT size the stroke to reach the cymba roof; size it to the crus range (≈7 mm of travel, 11.7→18.4) plus margin.
- Antihelix: CV 0.44–0.67, detected on 61/102 — stays out of v1.
- Rim path SD lowest in the anterior/crus sector (1.99 mm) in this run too.
- **9/4 follow-up:** sized pad extensions and per-size stability in `CLAMP_SIZE_BANDS.md` (`cad/iem/size_bands.py`, `crus_bands.py`). Reach to the overhang does not scale with ear size; one pad alone fails moment balance — the crus jaw must be in the seated pose before the clamp can be scored.
