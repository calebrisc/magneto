# Ear anthropometry for the IEM fit system

Design envelope for a one-size, no-S/M/L fit system. Compiled 2026-08-29 from published population studies; numbers are means ± SD with observed ranges where reported. Use the **ranges**, not the means — universality means the mechanism has to work at both ends.

## Concha (where the wing / jacket lives)

| Dimension | Korean 3D scans (n=200) | Caucasian 3D scans (n=96) | Cross-study range |
|---|---|---|---|
| Concha length (cymba+cavum, vertical) | 17.2 ± 1.3 mm (12.9–21.2) | 14.8 ± 1.3 mm | — |
| Concha width | 16.8 ± 1.8 mm (11.6–21.5) | 16.7 ± 1.8 mm | — |
| Concha depth | — | — | 9–16.8 mm across studies |
| Male vs female (Korean) | L 17.5 vs 16.8 · W 17.2 vs 16.5 | — | ~0.7 mm gender gap |

Source: Lee, Jung, You et al., *Ergonomics* 61(11), 2018 — 3D scans, 25 dimensions; concha depth range from the Hearing Health & Technology Matters review of multiple studies.

## Ear-hole / canal aperture (where the cup tip seals)

| Dimension | Korean (n=193–200) | Caucasian (n=96) | Review of canal studies |
|---|---|---|---|
| Aperture major axis (height) | 12.9 ± 2.0 mm (7.0–17.1) | 12.8 ± 2.1 mm | 8.9–12.5 mm, ~8.5 typical |
| Aperture minor axis (width) | 9.2 ± 2.4 mm (4.7–20.0) | 8.3 ± 2.1 mm | 5.7–9.1 mm, ~6 typical |
| Ear-hole depth (aperture to first bend region) | 8.3 ± 1.2 mm (5.1–11.3) | 8.3 ± 1.1 mm | — |
| Ear-hole length | 14.7 ± 1.3 mm (10.1–18.3) | 12.3 ± 1.3 mm | — |
| Aperture azimuth angle | 24.8° ± 8.8 (4.6–53.2) | — | — |
| Aperture elevation angle | −7.7° ± 23.9 (−70.5 to +49) | — | — |
| Canal length (to eardrum) | — | — | 23–32 mm; male ≈ 2 mm longer |
| Isthmus (narrowest, near 2nd bend) | — | — | ~9 × 5.7 mm |

Notes: aperture is **oval, height > width**, always. Ethnic ordering reported in the canal review: European > Asian > African for canal size. The "earhole minor axis" upper value of 20 mm is likely an outlier/landmarking artifact; treat 4.7–14 mm as the working range.

## What the numbers say for the design

- **Wing / jacket reach**: must be comfortable from a ~12 mm concha to a ~21 mm one (length), 11.6–21.5 mm wide, 9–17 mm deep. That is nearly 2:1 on every axis — confirms that a fixed wing can't do it and a graded-plateau spring + reach adjuster can.
- **Cup tip seal**: must seal an oval anywhere from about **7 × 5 mm to 17 × 13 mm**. A conformable cone with ≥ 2.5:1 diameter range covers it; a flange with a fixed diameter covers a third of people.
- **Nozzle**: 4 mm OD clears the small-aperture tail (5 mm minor axis) with the cup around it. 5.5–6.5 mm nozzles (industry standard) are why small-canal people get hurt.
- **Angles matter**: aperture azimuth varies 5°–53°, elevation −70° to +49°. The nozzle/cup joint needs a few degrees of compliance (soft cup neck), not a fixed angle.
- **Cale's case**: probably at the small-aperture tail — minor axis ≤ 6 mm. Measure it.

## 3D mesh datasets (for shape, not just numbers)

| Dataset | Subjects | Access | Notes |
|---|---|---|---|
| SYMARE-1 | 10 (MRI, full head + ear) | public | best free concha geometry |
| York Ear Model (YEM) | 500 synthesized from 10 scans + 605 landmarked images | academic user agreement — OU student can request | statistical shape model; PCA parameters let you sweep the population |
| Lee et al. Korean/Caucasian scans | 326 | not public (POSTECH) | the numbers above |
| Chinese canal casting study | 700 | paper only | 23 canal variables to the second bend |
| Zhang et al. Chinese cavum concha shape analysis | 1195 scans | paper only | statistical shape of cavum concha + meatus |

Plan: use the ranges here to set the mechanism envelope, use SYMARE-1 meshes to sanity-check the wing/cup against real geometry, have the OU engineer request YEM to sweep the population in software, and validate on our own scans (iPhone photogrammetry or impression putty).

## Sources

- Lee W., Jung H., Bok I., Kim C., Kwon O., Choi T., You H. "Anthropometric analysis of 3D ear scans of Koreans and Caucasians for ear product design." *Ergonomics* 61(11), 2018. https://www.tandfonline.com/doi/full/10.1080/00140139.2018.1493150 — slides: https://www.slideshare.net/WonsupLee1/3d-ear-anthropometry-for-earphone-design
- Staab W. "The Human Ear Canal" series, Hearing Health & Technology Matters, 2023. https://hearinghealthmatters.org/waynesworld/2023/human-ear-canal-v/
- Staab W. "Human Ear Concha Dimensions – Part 4", 2017. https://hearinghealthmatters.org/waynesworld/2017/human-ear-concha-dimensions-part-4/
- "Anthropometric Measurements of the External Auditory Canal for Hearing Protection Earplug" (700 Chinese subjects). https://www.researchgate.net/publication/304702608
- York Ear Model. https://www-users.york.ac.uk/~np7/research/YEM/
- SYMARE. https://pure.york.ac.uk/portal/en/datasets/sydney-york-morphological-and-recording-of-ears-database-symare
