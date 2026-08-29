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
| **SONICOM HRTF dataset** | **300** (200 + 100 extended) | **public download** — https://transfer.ic.ac.uk:9090/#/2022_SONICOM-HRTF-DATASET/ | full-head + ear **STL meshes** (watertight and raw), depth photos, demographics. The primary shape dataset for us. |
| **HUTUBS** (TU Berlin) | **96** | **public** — https://depositonce.tu-berlin.de/items/dc2a3076-a291-417e-97f0-7697e332c960 | head mesh + **pinna scanned at 0.05 mm** with an Artec Space Spider; anthropometric features included. Best concha/aperture detail. |
| 119-subject ear mesh + PRTF set | 119 | public (arXiv 2010.04546) | structured-light ear scans, 18,176 vertices in correspondence across subjects — easy to compute population stats on |
| Notre Dame UND-J2 | 415 (1,800 range images) | license agreement | biometrics-oriented profile range scans; lower detail than the above |
| SYMARE-1 | 10 (MRI, full head + ear) | public | canal included (MRI), tiny n |
| York Ear Model (YEM) | 500 synthesized from 10 scans + 605 landmarked images | academic user agreement — OU student can request | statistical shape model; PCA parameters let you sweep the population |
| Lee et al. Korean/Caucasian scans | 326 | not public (POSTECH) | the numbers above |
| Chinese canal casting study | 700 | paper only | 23 canal variables to the second bend |
| Zhang et al. Chinese cavum concha shape analysis | 1195 scans | paper only | statistical shape of cavum concha + meatus |

Plan: use the ranges here to set the mechanism envelope; download SONICOM + HUTUBS (~400 real ears, no paperwork) and fit the wing/cup against them in software; validate on our own scans (iPhone photogrammetry or impression putty). YEM/UND are optional.

## Cross-population check (added same day)

Caveat first: "concha length" is defined differently between studies. 3D-scan studies (Korean/Caucasian) measure the **cavum** (the bowl) ≈ 15–17 mm. Photo/caliper studies (Indian, Nigerian, Vietnamese, Turkish, Italian) measure the **whole concha** top-of-cymba to bottom ≈ 25–29 mm. Compare only within a method.

| Population | n | Method | Concha length | Concha width | Aperture H × W |
|---|---|---|---|---|---|
| Korean | 200 | 3D scan (cavum) | 17.2 ± 1.3 | 16.8 ± 1.8 | 12.9 × 9.2 (landmark def. generous) |
| Caucasian | 96 | 3D scan (cavum) | 14.8 ± 1.3 | 16.7 ± 1.8 | 12.8 × 8.3 |
| Taiwanese | 38 | CT | — | — | F 9.1 × 6.3 · M 9.6 × 6.8 |
| Multi-study canal review (mostly Western) | — | casts/CT | — | — | 8.9–12.5 × 5.7–9.1 |
| Central Indian | — | caliper (whole) | M 27.5 ± 1.8 · F 25.2 ± 2.3 | M 19.3 ± 2.0 · F 17.9 ± 1.9 | — |
| Nigerian (Hausa/Igbo/Yoruba) | — | photo (whole) | M 29.1 ± 2.0 · F 29.2 ± 2.2 | — | — |
| Vietnamese | 2000 | photo (whole) | F 26.9 ± 2.3 | F 15.9 ± 2.4 | — |
| Maharashtrian Indian | 505 | photo (whole) | (tables in paper; smaller than Central Indian) | — | — |

Ethnic ordering for canal size (Thomas et al., cited in the canal review): European > Asian > African, but "dimensional trends are fairly similar overall."

### The math that matters

Between-population differences in the **means** are ≤ 3–4 mm on every dimension where a same-method comparison exists (Korean vs Caucasian concha width: 0.1 mm; Indian vs Nigerian whole-concha length: ~2–4 mm; Taiwanese vs Western aperture: ~1–2 mm).

Within-population spread is larger: ±2 SD ≈ **±3.6 mm** on concha width, **±4–5 mm** on aperture minor axis, **±2.6 mm** on cavum length. Observed ranges (min–max) are wider still: concha width 11.6–21.5 mm, aperture minor axis 4.7–14 mm.

So: a mechanism that covers one population's observed **range** already contains every other population's **mean**. Shift the envelope outward ~2 mm on each end for the ethnic mean offsets and it covers roughly the 5th–95th percentile of every group measured. That is the design envelope to use:

| Dimension | Design envelope |
|---|---|
| Cavum concha length | 11–23 mm |
| Concha width | 10–24 mm |
| Concha depth | 8–18 mm |
| Aperture height | 7–18 mm |
| Aperture width | 4.5–14 mm |
| Aperture azimuth | 0°–55° |
| Aperture elevation | −70° to +50° (one study; treat as soft) |

### What is *not* verified

- No cross-ethnic data on aperture **angles** or concha **rim shape** (notch width/depth, antihelix undercut). These are what a rigid part would care about; a compliant part cares less.
- African-descent aperture size is only known as "smaller than Asian" from one citation; no numbers found. The 4.5 mm low end already assumes a small tail.
- Public 3D meshes (SONICOM, HUTUBS) are mostly European-recruited university populations; they cover *shape* well but don't fix the ethnic-tail question on their own.

### Is 3D required?

For the **numbers**: no. Ranges above are enough to set spring travel, plateau force, cup cone range, nozzle OD.
For the **shape**: partially. The wing hooks under the antihelix and the cup meets an oval at an angle — those contours aren't in any table. A compliant lattice tolerates shape error that a rigid part wouldn't, so 10–20 scans of *our own* people (varied ears, deliberately including small and large) is the real validation, not a 500-ear dataset. York Ear Model is a nice-to-have for sweeping shapes in software.

Additional sources: Taiwanese CT study (Inter-Noise 2014, https://www.acoustics.asn.au/conference_proceedings/INTERNOISE2014/papers/p121.pdf); Nigerian ear morphometry (https://link.springer.com/article/10.1186/s42269-021-00665-0); Maharashtrian study (https://pmc.ncbi.nlm.nih.gov/articles/PMC6018292/); Vietnamese 2000-subject study (J Craniofac Surg 2022); Central Indian values as cited in the Maharashtrian paper.

## Sources

- Lee W., Jung H., Bok I., Kim C., Kwon O., Choi T., You H. "Anthropometric analysis of 3D ear scans of Koreans and Caucasians for ear product design." *Ergonomics* 61(11), 2018. https://www.tandfonline.com/doi/full/10.1080/00140139.2018.1493150 — slides: https://www.slideshare.net/WonsupLee1/3d-ear-anthropometry-for-earphone-design
- Staab W. "The Human Ear Canal" series, Hearing Health & Technology Matters, 2023. https://hearinghealthmatters.org/waynesworld/2023/human-ear-canal-v/
- Staab W. "Human Ear Concha Dimensions – Part 4", 2017. https://hearinghealthmatters.org/waynesworld/2017/human-ear-concha-dimensions-part-4/
- "Anthropometric Measurements of the External Auditory Canal for Hearing Protection Earplug" (700 Chinese subjects). https://www.researchgate.net/publication/304702608
- York Ear Model. https://www-users.york.ac.uk/~np7/research/YEM/
- SYMARE. https://pure.york.ac.uk/portal/en/datasets/sydney-york-morphological-and-recording-of-ears-database-symare
