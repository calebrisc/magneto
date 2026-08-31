# Magneto IEM — virtual fit validation ("try-on") · v2

Seats the generated IEM assembly in 107 real and synthetic ears and scores the
fit. **v2** scores `cad/iem/generate.py` at commit **`d740ddc`** — the canted
nozzle (`nozzle_cant_deg = 45`) added in response to v1's primary
recommendation — against the collinear baseline generated from the same
checkout with `--cant 0`.

Both sides of the comparison were re-run from scratch under an identical,
corrected seating cost, so the numbers below are directly comparable. They are
**not** comparable to the v1 tables; see [what changed in v2](#what-changed-in-v2).

> **Scoring note.** The protrusion thresholds were recalibrated on 2026-08-31
> from ≤2 mm / ≤5 mm to **≤10 mm pass / ≤14 mm marginal**. The old pair
> described a custom-moulded IEM sitting flush and was cleared by 1 ear in 107,
> which made it useless as a discriminator — and, as it turns out, it was
> masking the cant's benefit. Justification:
> [threshold calibration](#protrusion-threshold-calibration). Both scorings are
> shown throughout.

## Verdict

**Keep the cant — under realistic thresholds it is worth far more than it first
appeared. Shortening is still wanted, but ~4–6 mm, not ~12 mm. And the top
failure axis is no longer protrusion: it is the seal.**

| | cant 0 (baseline) | cant 45 (default) |
|---|---|---|
| pass / marginal / fail — **recalibrated** | 0 / 20 / 87 | **0 / 40 / 67** |
| pass / marginal / fail — strict (old) | 0 / 1 / 106 | 0 / 0 / 107 |
| protrusion fails — **recalibrated** | 56 | **38** |
| protrusion fails — strict (old) | 103 | 107 |
| protrusion median | 14.19 mm | **12.88 mm** |
| protrusion p90 / max | 21.08 / 24.53 mm | **17.83 / 20.26 mm** |

Three things follow, and the first two only become visible once the threshold is
realistic:

1. **The cant is worth much more than the strict scoring suggested.** Under
   ≤2 mm it looked like the cant did nothing (103 → 107 fails, nominally
   *worse*). Under ≤10 mm it removes **18 protrusion failures** (56 → 38) and
   **20 overall failures** (87 → 67). The strict threshold was saturated — it
   failed nearly everything either way, so it could not register the
   improvement. This is the clearest argument for the recalibration.
2. **Shortening is now cheap.** Because the distribution's centre sits just above
   the new pass line rather than 12 mm above the old one, **4 mm** of reduction
   takes protrusion to 64 % pass / 93 % within marginal, and **6 mm** to 74 % /
   99 %. The earlier "~12 mm, and no efficient axis to cut it from" problem
   largely dissolves — a corner chamfer plus a modest trim is now a plausible
   route, and **no driver-stack redesign is required**.
3. **The seal is now the binding constraint.** Among the 67 remaining failures,
   **seal drives 44** against protrusion's 38. Worse, seal caps the whole design:
   set protrusion to a perfect score on every ear and the build still only
   reaches **2 pass / 58 marginal / 47 fail**, because just 10 ears pass seal at
   all. Fixing protrusion alone cannot get this design to a pass.

**Priority order is now: seal first, then ~4–6 mm of shortening.** Details in
[the shortening verdict](#shortening-verdict-revised-46-mm-and-seal-is-now-the-binding-constraint).

---

## Protrusion threshold calibration

The original ≤2 mm pass / ≤5 mm marginal pair was aspirational rather than
physical. It describes a **custom-moulded** IEM, which is built to the ear
impression and sits nearly flush; a **universal** IEM is a fixed body parked in
the concha, and standing proud of the tragus plane is its normal condition, not
a defect. Across 107 ears the old pass line was cleared exactly once, so it
carried no information — every ear failed, and a metric that always fails cannot
tell a good change from a bad one. That is precisely what happened to the cant.

**What the published numbers can and cannot tell us.** Manufacturers quote
*overall shell depth*, not protrusion past the tragus plane, and the concha
swallows an unstated portion of the former — so specs cannot be converted into
this metric directly. They do fix the order of magnitude: the SIMGOT EA1000 is
22 × 17 × **20.7 mm**, and ~20–22 mm of overall shell depth is typical for a
universal IEM. Against the 9–17 mm concha depth in
`docs/EAR_ANTHROPOMETRY.md`, that leaves **single-digit millimetres standing
proud on a typical ear**, which is the band these thresholds encode.

**These thresholds are therefore a design-team calibration, not a measured
industry standard**, and the report labels them as such:

| | pass | marginal | fail |
|---|---|---|---|
| **recalibrated** (2026-08-31) | ≤ 10 mm | ≤ 14 mm | > 14 mm |
| strict (original, retained for reference) | ≤ 2 mm | ≤ 5 mm | > 5 mm |

10 mm proud is taken as normal for a universal fit; 14 mm as the limit before
the shell fouls the helix/antihelix on insertion and starts levering itself out.
Both live in `tryon.py` as `PROT_PASS` / `PROT_MARGINAL`, with the strict pair
kept alongside. **Revisit once a prototype exists and can be measured worn** —
this is the single assumption with the most leverage over the programme, and it
is currently a judgement call rather than a measurement.

---

## Side by side

Same 107 ears, same landmarks, same cost function; only `nozzle_cant_deg`
differs.

### Overall

**Recalibrated thresholds (≤10 / ≤14 mm):**

| dataset | n | cant 0 (P/M/F) | cant 45 (P/M/F) |
|---|---|---|---|
| hutubs | 58 | 0 / 6 / 52 | 0 / 22 / 36 |
| sonicom | 45 | 0 / 12 / 33 | 0 / 18 / 27 |
| synthetic | 4 | 0 / 2 / 2 | 0 / 0 / 4 |
| **all** | **107** | **0 / 20 / 87** | **0 / 40 / 67** |

HUTUBS carries most of the gain (6 → 22 marginal). The synthetic corners go the
other way (2 → 0): they are the envelope extremes, and the cant's compression of
the distribution costs the two easy corners more than it wins back on the two
hard ones.

**Strict thresholds (≤2 / ≤5 mm), for reference:**

| dataset | n | cant 0 (P/M/F) | cant 45 (P/M/F) |
|---|---|---|---|
| **all** | **107** | **0 / 1 / 106** | **0 / 0 / 107** |

### Per axis

Only the protrusion row changes with the threshold; seal, retention and
clearance are identical under both scorings.

| axis | cant 0 pass | marg | fail | cant 45 pass | marg | fail | verdict |
|---|---|---|---|---|---|---|---|
| seal | 5 | 38 | **64** | 10 | 53 | **44** | **much better — but now the worst axis** |
| retention | 57 | 37 | 13 | 48 | 53 | **6** | fails halved, passes down |
| clearance | 76 | 27 | 4 | **101** | 5 | 1 | **much better** |
| protrusion — **recalibrated** | 18 | 33 | **56** | 21 | 48 | **38** | **18 fewer fails** |
| protrusion — strict | 1 | 3 | 103 | 0 | 0 | 107 | saturated, no signal |

Seal failures drop by 20 and clearance passes rise by 25 — the canted body sits
lower in the bowl and presents the skirt to the aperture at a far better angle.
Retention fails halve (13 → 6) but passes fall (57 → 48), i.e. the wing moves
from "not touching" into "touching a bit too much or not quite enough": median
`wing_tip` goes −0.29 → −1.16 mm, straight through the −0.5 to −2.0 mm target.

The two protrusion rows are the case for the recalibration in miniature. The
strict row says the cant made things *worse* (103 → 107); the recalibrated row
says it removed a third of the protrusion failures (56 → 38). The underlying
millimetres are identical — only the strict row's saturation differs.

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

Ranked by number of failing axes, then seating cost, under the **recalibrated**
thresholds.

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
| synthetic xl_shallow | retention + seal | 0.39 | 0.66 | +3.05 | 11.87 | −1.25 |
| sonicom P0044 | protrusion + seal | 0.43 | 2.15 | +0.17 | 16.15 | +0.62 |
| hutubs pp63 | protrusion + seal | 0.26 | 0.76 | −1.37 | 17.67 | +3.00 |
| hutubs pp77 | protrusion + retention | 0.64 | 4.45 | +1.49 | **19.32** | +1.65 |
| hutubs pp10 | protrusion + retention | 0.65 | 3.14 | +3.32 | 17.93 | −0.43 |

The tail is markedly less bad. The baseline's worst ear fails **four** axes; the
canted build's worst fails **two**, and no canted ear fails more than two.
P0044 — worst on the baseline with all four axes failing — drops to two.

Two things to read off this table. **Seal appears in three of the five**, which
is the pattern that makes seal the priority. And `xl_shallow`, the
small-aperture/shallow-concha synthetic corner, is now the hardest ear in the
set *without failing protrusion at all* (11.87 mm, inside the marginal band): an
8 mm-deep bowl on a 7 × 4.5 mm aperture defeats the skirt and the wing, not the
shell depth. That is the envelope tail this design is furthest from serving.

---

## Shortening verdict, revised: ~4–6 mm — and seal is now the binding constraint

Further reduction needed **on top of** the cant, measured on the canted build.
Under the recalibrated thresholds the economics change completely:

| stack removed | ≤ 10 mm (pass) | ≤ 14 mm (marginal) | > 14 mm (fail) |
|---|---|---|---|
| 0 mm (today) | 20 % | 64 % | 36 % |
| 2 mm | 38 % | 74 % | 26 % |
| 3 mm | 50 % | 84 % | 16 % |
| **4 mm** | **64 %** | **93 %** | **7 %** |
| **6 mm** | **74 %** | **99 %** | **1 %** |
| 8 mm | 93 % | 100 % | 0 % |

For reference, the same sweep against the strict ≤2 mm line needed **12 mm** for
64 % and **14 mm** for 74 %. The requirement fell by two thirds purely because
the target moved to a defensible place — which is why the threshold, not the
geometry, was the highest-leverage thing on this list.

**4–6 mm is reachable; 12 mm was not.** Sensitivity of protrusion to 1 mm removed
from the shell along each design axis (1.00 would be perfectly efficient):

| design axis | what it is | mm of protrusion removed per mm cut |
|---|---|---|
| +Y | superior, toward the wing | 0.65 |
| +Z | faceplate normal / stack height | 0.49 |
| +X | body long axis, toward the nozzle | 0.41 |

The worst-protruding point is a **corner**, not a face: in design coordinates it
sits at a median of (−11.5, −4.1, +3.7), the posterior-inferior corner of the
faceplate, diagonal to all three axes. That was fatal when 12 mm was needed —
24 mm of +Z out of a 13.2 mm stack. At a 4–6 mm target it is merely a
constraint on *where* to cut: a combined trim of 6 mm X + 4 mm Y + 4 mm Z buys
≈ 7 mm, comfortably past the 6 mm mark, and a corner chamfer alone plausibly
covers the 4 mm case.

### Revised priority

1. **Fix the seal — it is now the ceiling, not protrusion.** Seal drives **44 of
   the 67** remaining failures, against protrusion's 38. And it binds
   independently: give every ear a perfect protrusion score and the build still
   only reaches **2 pass / 58 marginal / 47 fail**, because only 10 ears pass
   seal at all. No amount of shortening rescues this design on its own. The
   v1 seal analysis still stands — median rim gap 2.76 mm against a
   0.25 mm-wall land, i.e. the skirt is asked for travel it does not have — and
   the reworked 4.5 mm contact land moved seal passes only 5 → 10. Treat skirt
   compliance (or a smaller rim that lands inside the cavum floor) as the
   primary open problem.
2. **Chamfer the posterior-inferior faceplate corner.** Highest value per mm of
   material removed, since it is the actual protruding point on most ears, and
   it is cheap. Re-measure after; it may cover the 4 mm case by itself.
3. **Trim 4–6 mm of stack if the chamfer is not enough** — a combined X/Y/Z trim,
   not a single-axis cut. **A driver-stack redesign is no longer indicated**:
   that recommendation was an artefact of the 12 mm target.

Tune the cant angle only after the above: at 45° the jacket has already started
lifting off the concha floor (+0.93 mm), so more cant trades seal and bedding
for a protrusion gain the tail data says is nearly exhausted.

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

Grade agreement under the **recalibrated** thresholds:

| axis | same grade | fail count, seed 0 → 7 |
|---|---|---|
| clearance | 98 / 107 | 1 → 0 |
| protrusion | 81 / 107 | 38 → 36 |
| seal | 57 / 107 | 44 → 34 |
| retention | 55 / 107 | 6 → 10 |
| **overall** | **78 / 107** | **0/40/67 → 1/42/64** |

Under the strict thresholds the protrusion row read 107/107 — but only because
every ear failed under both seeds, which is agreement without information. The
recalibrated row (81/107) is the honest figure, and it is the one to trust.

What survives reseeding, and what does not:

- **Robust — the conclusions this report rests on.** Aggregate counts move by
  only a few ears (protrusion fails 38 → 36, overall fails 67 → 64), so the
  cant's ~18-failure improvement and the ~4–6 mm shortening estimate are far
  larger than the noise. Clearance is near-exact (98/107).
- **Noisy — do not quote per-ear.** Seal (57/107) and retention (55/107) agree
  barely more than half the time, because most of the population sits close to
  those thresholds and the wing tip rests on the antihelix ridge, where a
  sampled signed-distance field's sign flips within ~0.3 mm. Individual seal and
  retention counts are directional only.
- **Seal's primacy holds either way**: it is the largest fail count under both
  seeds (44 and 34), and under both it exceeds protrusion.

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
- **The protrusion thresholds are a design-team calibration, not a measurement.**
  ≤10 / ≤14 mm is a judgement about what a universal fit may stand proud of the
  tragus plane, anchored only on published *overall shell depths* (a different
  quantity) and the concha-depth range. It is the highest-leverage assumption in
  this report — it moved the shortening requirement from ~12 mm to ~4–6 mm — and
  it should be replaced with worn measurements from a prototype.

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
