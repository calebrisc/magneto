# Magneto IEM — virtual fit validation ("try-on") · v2

Seats the generated IEM assembly in 107 real and synthetic ears and scores the
fit. **v2** scores `cad/iem/generate.py` at commit **`d740ddc`** — the canted
nozzle (`nozzle_cant_deg = 45`) added in response to v1's primary
recommendation — against the collinear baseline generated from the same
checkout with `--cant 0`.

Both sides of the comparison were re-run from scratch under an identical,
corrected seating cost, so the numbers below are directly comparable. They are
**not** comparable to the v1 tables; see [what changed in v2](#what-changed-in-v2).

## Verdict

**The cant is a real improvement and should be kept. It does not fix the
protrusion failure, and body shortening is still required.**

| | cant 0 (baseline) | cant 45 (default) |
|---|---|---|
| pass / marginal / fail | 0 / 1 / 106 | **0 / 0 / 107** |
| protrusion median | 14.19 mm | **12.88 mm** |
| protrusion p90 | 21.08 mm | **17.83 mm** |
| protrusion max | 24.53 mm | **20.26 mm** |

The cant bought **1.31 mm of median protrusion** — it helped 72 ears and hurt
35 — while clearly improving seal and clearance. That is nowhere near the
~12 mm needed. On the overall grade the canted build actually scores *worse*
(107 fails vs 106), because the one baseline ear that squeaked through on
protrusion no longer does: the cant **compresses the distribution** rather than
shifting it, raising the floor from −1.13 mm to +5.61 mm while pulling the
ceiling down from 24.53 mm to 20.26 mm (IQR 7.71 → 5.19 mm).

**Shortening still needed: yes — about 12 mm more.** Details and the sting in
the tail (no single axis is an efficient place to cut) in
[the shortening verdict](#shortening-verdict-yes-12-mm-and-no-single-axis-delivers-it).

---

## Side by side

Same 107 ears, same landmarks, same cost function; only `nozzle_cant_deg`
differs.

### Overall

| dataset | n | cant 0 (P/M/F) | cant 45 (P/M/F) |
|---|---|---|---|
| hutubs | 58 | 0 / 0 / 58 | 0 / 0 / 58 |
| sonicom | 45 | 0 / 1 / 44 | 0 / 0 / 45 |
| synthetic | 4 | 0 / 0 / 4 | 0 / 0 / 4 |
| **all** | **107** | **0 / 1 / 106** | **0 / 0 / 107** |

### Per axis

| axis | cant 0 pass | marg | fail | cant 45 pass | marg | fail | verdict |
|---|---|---|---|---|---|---|---|
| seal | 5 | 38 | **64** | 10 | 53 | **44** | **much better** |
| retention | 57 | 37 | 13 | 48 | 53 | **6** | fails halved, passes down |
| clearance | 76 | 27 | 4 | **101** | 5 | 1 | **much better** |
| protrusion | 1 | 3 | 103 | 0 | 0 | **107** | no real change |

Seal failures drop by 20 and clearance passes rise by 25 — the canted body sits
lower in the bowl and presents the skirt to the aperture at a far better angle.
Retention fails halve (13 → 6) but passes fall (57 → 48), i.e. the wing moves
from "not touching" into "touching a bit too much or not quite enough": median
`wing_tip` goes −0.29 → −1.16 mm, straight through the −0.5 to −2.0 mm target.

### Metric medians

| metric | cant 0 | cant 45 | Δ |
|---|---|---|---|
| skirt rim contact coverage | 0.46 | **0.54** | +0.08 |
| skirt rim max gap (mm) | 4.22 | **2.76** | −1.46 |
| wing tip signed dist (mm) | −0.29 | −1.16 | −0.87 |
| jacket mean clearance (mm) | **3.01** | 3.94 | +0.93 |
| faceplate past tragus (mm) | 14.19 | **12.88** | −1.31 |
| worst rigid interference (mm) | −0.82 | **−0.27** | +0.55 |

The one regression is the jacket, which moves 0.93 mm *further* off the concha
floor — the canted body pivots the gyroid skin away from the bowl. Minor next
to what the cant buys elsewhere, but it is the thing to watch if the cant angle
is tuned further.

### Where the cant actually pays: the tail

| protrusion ≤ | cant 0 | cant 45 |
|---|---|---|
| 2 mm | 1 % | 0 % |
| 5 mm | 4 % | 0 % |
| 8 mm | 9 % | 6 % |
| 10 mm | 17 % | 20 % |
| 12 mm | 35 % | 38 % |
| 15 mm | 55 % | **67 %** |
| 18 mm | 74 % | **93 %** |
| 20 mm | 82 % | **99 %** |

Below ~10 mm the cant is neutral-to-slightly-negative; above 15 mm it is
decisive. **It is a worst-case fix, not a median fix** — which is exactly what
you want from it, and exactly why it cannot substitute for shortening.

### Worst 5 ears

Ranked by number of failing axes, then seating cost.

**cant 0 (baseline)**

| ear | failing axes | cover | max gap | wing tip | protrusion | hard interf |
|---|---|---|---|---|---|---|
| sonicom P0044 | clearance + protrusion + retention + seal | 0.31 | 5.26 | +1.24 | 14.27 | −2.84 |
| sonicom P0036 | clearance + protrusion + seal | 0.31 | 4.08 | −2.86 | **19.05** | −2.61 |
| hutubs pp77 | protrusion + retention + seal | 0.35 | 5.02 | +1.23 | **18.92** | +0.14 |
| hutubs pp69 | protrusion + retention + seal | 0.44 | 6.01 | +1.27 | 17.82 | −0.85 |
| hutubs pp10 | protrusion + retention + seal | 0.46 | 2.84 | +1.55 | 16.86 | −1.40 |

**cant 45 (default)**

| ear | failing axes | cover | max gap | wing tip | protrusion | hard interf |
|---|---|---|---|---|---|---|
| synthetic xl_shallow | protrusion + retention + seal | 0.39 | 0.66 | +3.05 | 11.87 | −1.25 |
| sonicom P0048 | clearance + protrusion + seal | 0.46 | 5.87 | −0.16 | 9.50 | −2.60 |
| sonicom P0009 | protrusion + retention + seal | 0.46 | 4.17 | +1.05 | 13.24 | −1.03 |
| hutubs pp67 | protrusion + retention + seal | 0.36 | 1.47 | +1.24 | 12.87 | −0.42 |
| sonicom P0044 | protrusion + seal | 0.43 | 2.15 | +0.17 | **16.15** | +0.62 |

The tail is visibly less bad: no canted ear fails four axes (the baseline's
P0044 does), the worst canted protrusion in the top 5 is 16.15 mm against the
baseline's 19.05 mm, and P0044 — worst on the baseline with four failing axes —
drops to two. `xl_shallow`, the small-aperture/shallow-concha synthetic corner,
is now the hardest ear in the set: a 4.5 mm-deep bowl simply cannot swallow this
shell at any angle.

---

## Shortening verdict: yes, ~12 mm — and no single axis delivers it

Further reduction needed **on top of** the cant, measured on the canted build:

| along-normal stack removed | protrusion ≤ 2 mm | ≤ 5 mm |
|---|---|---|
| 0 mm (today) | 0.0 % | 0.0 % |
| 4 mm | 2.8 % | 12.1 % |
| 8 mm | 19.6 % | 50.5 % |
| 10 mm | 38.3 % | 67.3 % |
| **12 mm** | **64.5 %** | **84.1 %** |
| **14 mm** | **73.8 %** | **96.3 %** |

**But there is no cheap axis to cut.** Sensitivity of protrusion to 1 mm removed
from the shell along each design axis (1.00 would be perfectly efficient):

| design axis | what it is | mm of protrusion removed per mm cut |
|---|---|---|
| +Y | superior, toward the wing | 0.65 |
| +Z | faceplate normal / stack height | 0.49 |
| +X | body long axis, toward the nozzle | 0.41 |

The reason is that the worst-protruding point is a **corner**, not a face: in
design coordinates it sits at a median of (−11.5, −4.1, +3.7) — the
posterior-inferior corner of the faceplate, which is diagonal to all three axes.
Cutting 12 mm of protrusion out of the +Z stack alone would need 24 mm of Z, and
the entire rigid Z stack is 13.2 mm. Even an aggressive combined trim — 6 mm of
X, 4 mm of Y, 4 mm of Z — buys only ≈ 7 mm.

So the honest recommendation is **not** "shave the shell":

1. **Chamfer the posterior-inferior faceplate corner first.** It is the single
   protruding point on most ears and it is diagonal, so a corner break is worth
   more per mm of material than any face cut. Cheap, and worth re-measuring
   before anything structural.
2. **Then attack the driver/mag-float stack itself.** The core is
   17.0 × 16.0 × 9.8 mm, 394 mm³. Getting 12 mm of protrusion out means the
   shell has to get materially smaller in *volume*, not just thinner in one
   direction — this is a driver-selection and magnet-layout question, not a
   surfacing one.
3. **Re-check the ≤2 mm pass threshold.** It is strict. Even the baseline only
   ever cleared it on 1 ear of 107, and real universal IEMs do stand somewhat
   proud of the tragus plane. If 8–10 mm is in fact acceptable in the industrial
   design, the canted build is already at 6 % / 20 % rather than 0 %, and the
   shortening target drops by more than half. **This threshold should be settled
   before committing to a driver change** — it moves the requirement more than
   any geometry change available.

Tune the cant angle only after the above: at 45° the jacket has already started
lifting off the concha floor (+0.93 mm), so more cant trades seal and bedding
for a protrusion gain that the tail data says is nearly exhausted.

---

## What changed in v2

Three pipeline corrections, all of which affect the numbers. **v1's tables are
superseded and not comparable to these.**

1. **The canted geometry was being read in the wrong frame.** `carrier.stl` and
   the nozzle inserts are written in the **nozzle-local** frame; only
   `assembly.stl` applies `nozzle_T`. Loading the carrier and treating its
   coordinates as assembly-frame coordinates put the entire seal 45° away from
   where it physically is. `earfit.py` now imports the transform from
   `generate.py` rather than hardcoding it, so the cant can never drift out of
   sync, and `iem_points(cant=…)` lets one process score a cant-0 and a cant-45
   build side by side. Verified: every sampled point set now lies within
   0.23 mm of the real `assembly.stl` surface.

2. **The seating search had an unconstrained roll, and was using it.** Nothing
   in the cost said which way was "out", so the optimiser rolled the shell until
   the faceplate pointed sideways or *into the skull* — **23 of 107 baseline
   ears** seated with the faceplate normal more than 90° from the out-of-head
   direction. That is not a fit any wearer could adopt, and it quietly corrupted
   the protrusion metric. Added `c_face`, which leaves 45° of slack free and is
   quadratic beyond it. Flipped poses: 23 → **2** on the baseline, 2 on the
   canted build; median faceplate-vs-out-of-head angle 73° → 60° (baseline) and
   51° (canted).

3. **A concurrent-run hazard was found and closed.** Two `--reseat` jobs writing
   `ears/aligned/*.json` at once silently interleave their poses. The affected
   runs were discarded and re-run serially. Landmark fields survive this (reseat
   only rewrites pose fields), which is why re-landmarking was not needed. If
   you run these by hand, run one at a time.

Carried over from v1: `c_soft` (skirt burial) and `c_prot` (protrusion) in the
seating cost; clearance graded on `shell` rather than `rigid` so it stops
contradicting retention; seeded sampling; Powell refinement over the best 4
starts.

---

## Method

Unchanged from v1 except where noted above. Summarised; the full rationale for
the landmark detector lives in `align_ear.py`'s docstring.

1. **Ear window** — anchor on the most-lateral vertex of the upper half of the
   head, cut an 85 × 85 × 50 mm box.
2. **Hull-relative depth map** — medial rays on a 0.5 mm grid against the ear
   and against its own convex hull, subtracted, so the skull behind the pinna
   flattens to zero and the concha does not.
3. **Basin selection** — an *enclosure* filter (fraction of an outward ray
   hemisphere that escapes; a cavum is a bowl at 0.10–0.25, the retroauricular
   sulcus a groove at 0.55–0.70, cut at 0.42) plus a *bowl* gate (depth scored
   only over cells ≥0.6 × the basin's largest inscribed radius from its edge,
   which stops the cymba and crus-helicis grooves winning on raw depth). Across
   107 ears `basin_escape` ran median 0.16, max 0.33 — the sulcus was never
   selected, and 107/107 ears landmarked with no manual override.
4. **Seating** — 5 nozzle rakes × 3 rolls seeded so the rim centre lands on the
   aperture, best 4 refined by bounded Powell (±10 mm, ±40°). Cost terms:
   `c_rim` (±0.75 mm contact band), `c_pen` (rigid penetration, worst point),
   `c_wing`, `c_jac`, `c_soft`, `c_prot`, `c_face`.

Metrics and grading are unchanged (seal / retention / clearance / protrusion;
overall = worst axis). Sign convention: **positive = clearance, negative = inside
the flesh.** Per-ear rows: `cad/iem/ears/aligned/tryon.csv` (gitignored).

### Dataset

| source | ears | provenance |
|---|---|---|
| **HUTUBS** (TU Berlin) | 58 | CC BY 4.0; all 58 head meshes released of 96 subjects; pinna at 0.05 mm |
| **SONICOM** (Imperial) | 45 | research use; `_preprocessed.stl`, Frankfurt-aligned, canal mouth open |
| **synthetic corners** | 4 | envelope corners from `EAR_ANTHROPOMETRY.md` |
| **total** | **107** | right ears |

### Run-to-run stability

The cost surface is rough and Powell finds a local minimum. Re-seating the
canted build under a different sampling seed (`--field-seed 7`):

| metric | median \|Δ\| | p90 \|Δ\| |
|---|---|---|
| `rim_cover` | 0.10 | 0.32 |
| `rim_gap` | 0.73 mm | 2.38 mm |
| `wing_tip` | 1.02 mm | 3.35 mm |
| `jacket_mean` | 0.87 mm | 2.72 mm |
| `protrusion` | 0.98 mm | 4.18 mm |
| `hard_min` | 0.51 mm | 2.31 mm |

| axis | same grade | fail count, seed 0 → 7 |
|---|---|---|
| protrusion | **107 / 107** | 107 → 107 |
| clearance | 98 / 107 | 1 → 0 |
| seal | 57 / 107 | 44 → 34 |
| retention | 55 / 107 | 6 → 10 |
| **overall** | **107 / 107** | **107 → 107** |

**The protrusion conclusion and the overall verdict are exactly reproducible**
(107/107 on both). The **seal and retention axis counts are noisy** — agreement
only 57/107 and 55/107, because most of the population sits near those
thresholds, and the wing tip rests on the antihelix ridge where a sampled
signed-distance field's sign flips within ~0.3 mm. Read "seal improved
substantially, most ears remain marginal" and "retention fails roughly halved"
as the findings; **do not quote the individual seal/retention counts as
precise.** The side-by-side conclusions above rest on protrusion and clearance,
which are stable.

Tables in this report are `--field-seed 0`.

### Data-quality flags

37 of 107 ears carry `weak`, **all** for `canal_run < 1.5 mm`, none for
enclosure. Structured light cannot see round the first canal bend, so in both
datasets the canal is a dimple, not a tube (median `canal_run` 1.8 mm). Fine for
this design, which seals on the funnel entrance — but this dataset cannot
validate anything that enters the canal, including the three nozzle inserts.

## Known limitations

- **Right ears only** (`--side left` is supported but was not run).
- **Static geometry only** — no skin deformation, no silicone FEA, no insertion
  path, no jaw movement. Soft compliance is fixed millimetre allowances, not
  simulated. Mechanical side: `MECH_VALIDATION.md`.
- **The synthetic XL corners have an unrealistically wide funnel mouth** (built
  at 2× aperture, so 36 × 28 mm), large enough to swallow the Ø19 skirt; their
  seal numbers are optimistic. The 103 real ears carry the conclusions.
- **Both scan populations are European-recruited university cohorts**; the four
  synthetic corners are the only coverage of the envelope tails.
- **The ≤2 mm protrusion pass threshold is unvalidated** against industrial
  design intent — see recommendation 3 above.

## Reproducing

```bash
cd cad/iem && python -m venv .venv && .venv/bin/pip install -r requirements.txt
.venv/bin/python fetch_ears.py --source all --n 45     # ~1.7 GB, gitignored
.venv/bin/python make_synthetic_ear.py
for d in sonicom hutubs synthetic; do
  .venv/bin/python align_ear.py --dataset $d --qc-png --skip-existing
done

# canted default (this report's headline)
.venv/bin/python generate.py --all
.venv/bin/python align_ear.py --reseat --jobs 8 && .venv/bin/python tryon.py

# collinear baseline, same checkout
.venv/bin/python generate.py --all --cant 0 --out stl_cant0 --ear right
.venv/bin/python align_ear.py --reseat --jobs 8 --cant 0 --stl-dir stl_cant0/right
.venv/bin/python tryon.py --cant 0 --stl-dir stl_cant0/right
```

Run reseats **one at a time** — concurrent reseats interleave their poses.
`--field-seed` varies the sampling for a stability check; `--qc-png` writes a
per-ear depth map with the detected aperture and tragus marked; `--landmarks`
overrides a bad pick by hand.
