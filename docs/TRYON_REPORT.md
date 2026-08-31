# Magneto IEM — virtual fit validation ("try-on") · v4

Seats the generated IEM assembly in 107 real and synthetic ears and scores the
fit. **v4** scores `cad/iem/generate.py` at commit **`ff5b52a`** — notch sector
deepened to a measured 3.25 mm, core restored to full size (protrusion accepted,
no trim) — with fresh 6-DOF seating optimisations on all 107 ears at two sampling
seeds, against v3 (`1ae768a`).

Scoring: protrusion **≤10 mm pass / ≤14 mm marginal**
([calibration](#protrusion-threshold-calibration)) and the **compliance-aware
seal** ([method](#compliance-aware-seal-scoring)).

## Verdict

**v4 regressed, and the cause is a recommendation this report got wrong. The
3.25 mm notch flare makes the seal worse, not better — it breaks more ears than
it fixes, on both seeds. Revert it. Separately, the retention diagnosis is the
opposite of what was assumed: the wing is too LONG, not too short.**

| | v3 (`1ae768a`) | v4 (`ff5b52a`) |
|---|---|---|
| overall P / M / F — seed 0 | 21 / 59 / 27 | **7 / 64 / 36** |
| overall P / M / F — seed 7 | 13 / 68 / 26 | 17 / 62 / 28 |
| seal pass @ 2.5 / 4.0 mm | 76 / 103 | **67 / 98** |
| protrusion median | 10.52 mm | 11.26 mm |
| protrusion fails | 17 | 24 |

The protrusion regression is expected and was accepted: the core went back to
full size. The seal regression was not.

### The notch flare is harmful — and that is my error

v3 measured, per ear, the **additional deformation budget** needed to seal:
median 2.0 mm, p90 3.4 mm. It then recommended delivering that as "3.0–3.5 mm of
realised **radial reach**". **Those are not the same quantity, and the difference
is the whole problem.** A deformation budget is direction-free — it asks how far
the silicone must stretch toward *wherever the nearest flesh is*. A radial
extension commits the lip to one direction. At the intertragic notch that
direction points into the notch, and **the notch is an opening in the cartilage**:
there is no flesh there to seal against. The lip now juts into free air.

The evidence is a fully controlled, same-pose comparison — the identical v4
seatings scored with the as-built flared rim versus a plain Ø19 circle:

| budget | seed 0: circle → flared | seed 7: circle → flared |
|---|---|---|
| 1.5 mm | 41 → 39 (broke 7, fixed 5) | 45 → 38 (broke 10, fixed 3) |
| **2.5 mm** | **75 → 67 (broke 12, fixed 4)** | **81 → 79 (broke 8, fixed 6)** |
| 4.0 mm | 102 → 98 (broke 5, fixed 1) | 102 → 101 (broke 3, fixed 2) |

Net negative at every budget on both seeds. And the ears it breaks were not
marginal — many had a **perfect** seal before:

| ear | gap arc, circle → flared | coverage |
|---|---|---|
| hutubs pp31 | 0° → **56°** | 1.00 → 0.84 |
| sonicom P0025 | 0° → 51° | 1.00 → 0.86 |
| hutubs pp95 | 6° → 64° | 0.98 → 0.81 |
| sonicom P0053 | 0° → 40° | 1.00 → 0.89 |
| synthetic xl_shallow | 6° → 37° | 0.98 → 0.88 |

10 of the 12 ears broken at the 2.5 mm budget have their new gap **at the notch
sector itself** — the flare drills its own hole in the seal it was added to
close. Measured across the flared samples, the flare improves the *median* gap
in that sector by 0.70 mm but worsens the p90 (+3.05 → +3.45 mm); since the seal
criterion is continuity, the tail is what decides, and the tail got worse.

**Recommendation: revert the notch sector to a plain Ø19 rim** (or at most the
~1 mm of v3, which was net-neutral). Then pursue notch sealing as **compliance,
not geometry** — thinner wall, longer free length, or lower durometer over that
90° so the lip can *drape toward* the cartilage on either side of the notch
instead of reaching into the void between them. The budget sweep that motivated
this remains valid as a *compliance* target; it was never a geometry target.

### Retention: the wing is too long, not too short

The wing was assumed to need more reach. It does not.

| outcome (spring band, tip ∈ [−1.5, 0] mm) | ears | share |
|---|---|---|
| retained | 31 | 29 % |
| **overpressed** (tip < −1.5 mm) | **47** | **44 %** |
| misdirected (gap, aim ≥ 60° off) | 18 | 17 % |
| blocked (already overlapping ≥1.5 mm elsewhere) | 10 | 9 % |
| short (gap, aiming correctly) | **1** | 1 % |

**Exactly one ear in 107 is short.** The dominant mode is *overpressed* — median
tip overlap −2.66 mm, well past the 1.5 mm the 0.25 N/mm spring can absorb, so
the wing is bottoming and levering the shell out. Full analysis and the wing
change that fixes the most ears: [retention](#retention-analysis).

---

## v3 → v4, side by side

### Overall (compliance seal + recalibrated protrusion)

| | seed 0 P/M/F | seed 7 P/M/F |
|---|---|---|
| v3 (`1ae768a`) | 21 / 59 / 27 | 13 / 68 / 26 |
| **v4 (`ff5b52a`)** | **7 / 64 / 36** | **17 / 62 / 28** |
| v4 with a plain circle rim | 9 / 66 / 32 | — |

The seed spread at v4 is wide (27 → 36 fails between seeds is smaller than the
7 → 17 swing in passes), so **do not over-read the headline pair**. The robust
statements are the paired, same-pose ones: the flare is net-negative on both
seeds, and protrusion worsened as expected with the core restored.

### Per axis (seed 0)

| axis | v3 | v4 | note |
|---|---|---|---|
| protrusion (≤10 / ≤14 mm) | 46 / 44 / 17 | 39 / 44 / 24 | expected — core untrimmed |
| seal (compliance 2.5 / 4.0) | 76 / 27 / 4 | **67 / 31 / 9** | **regression, from the flare** |
| clearance | 105 / 2 / 0 | 105 / 2 / 0 | unchanged, fully clear |
| retention | 55 / 44 / 8 | 50 / 50 / 7 | flat |

### Metric medians (seed 0)

| metric | v3 | v4 |
|---|---|---|
| faceplate past tragus | 10.52 mm | 11.26 mm |
| skirt rim contact coverage (rigid) | 0.61 | 0.56 |
| skirt rim max gap (rigid) | 2.96 mm | 3.03 mm |
| wing tip signed distance | −0.73 mm | −1.03 mm |
| jacket mean clearance | 4.27 mm | 4.63 mm |
| worst rigid interference | +0.26 mm | +0.30 mm |

The jacket keeps drifting off the concha floor (3.01 → 3.94 → 4.27 → 4.63 mm
across v2/v3/v4). It still costs no grades, but it is the surface the wing reacts
against, and it is a plausible contributor to the retention picture below.

### Notch attribution (v4, as built)

| budget | failing ears | gap at the notch | elsewhere | notch share |
|---|---|---|---|---|
| 1.5 mm | 68 | 50 | 18 | 74 % |
| 2.5 mm | 40 | 34 | 6 | **85 %** |
| 4.0 mm | 9 | 7 | 2 | 78 % |

Unchanged in character from v2 and v3: the notch is a quarter of the perimeter
and takes three-quarters to five-sixths of the failures. Deepening the flare did
not move this — it moved failures *into* the sector.

---

## Retention analysis

`retention_analysis.py`, at each ear's frozen seated pose. The wing is a
compliant leaf spring at ~0.25 N/mm, so the physically meaningful target is a
light interference fit — tip signed distance in **[−1.5, 0] mm**, the spring's
working range. That is stricter than `tryon.py`'s grading band, which asks only
for a passing grade.

### Why each graded failure fails

Only 7 ears fail the graded retention axis, and **none of them is a simple reach
problem** except one:

| ear | tip | worst overlap | aim angle | diagnosis |
|---|---|---|---|---|
| synthetic xs_deep | +1.05 | −2.41 | 117° | blocked |
| hutubs pp47 | +1.09 | +0.42 | 59° | short |
| synthetic xs_shallow | +1.24 | −2.47 | 80° | blocked |
| hutubs pp41 | +1.32 | −1.25 | 156° | misdirected |
| hutubs pp29 | +1.33 | +0.62 | 83° | misdirected |
| sonicom P0028 | +1.33 | −2.12 | 164° | blocked |
| sonicom P0045 | +1.40 | +0.79 | 95° | misdirected |

"Aim angle" is between the wing's growth direction and pressing straight into the
local ear surface. **Above 90° the wing is growing away from the surface it is
supposed to press.** Three of the seven are *blocked* — the tip shows a gap while
the wing is already overlapping 2.1–2.5 mm further down its span, so it bottoms
out before the tip can land. For those, a longer wing is strictly worse.

### What wing change fixes the most ears

Ears retained (tip in the spring band) after extending the wing by *L* along its
growth axis and splaying it by *θ* about the root:

| ext \ splay | −15° | −10° | −5° | 0° | +5° | +10° | +15° |
|---|---|---|---|---|---|---|---|
| **−3.0 mm** | 47 | **50** | 49 | 47 | 41 | 36 | 34 |
| **−2.5 mm** | 49 | 46 | 44 | 46 | 44 | 41 | 37 |
| **−2.0 mm** | 42 | 49 | 46 | 47 | 50 | 43 | 36 |
| −1.5 mm | 35 | 40 | 44 | 41 | 38 | 42 | 37 |
| −1.0 mm | 31 | 42 | 42 | 38 | 39 | 33 | 34 |
| −0.5 mm | 22 | 37 | 36 | 32 | 34 | 28 | 34 |
| **0 (as built)** | 17 | 30 | 27 | **31** | 30 | 33 | 33 |
| +1.0 mm | 19 | 26 | 22 | 16 | 22 | 29 | 27 |
| +1.5 mm | 16 | 18 | 20 | 17 | 18 | 22 | 20 |

**Shorten the wing by 2.5–3.0 mm and splay it −10°.** That takes retention from
31/107 to **50/107** — a 61 % improvement — and the whole upper half of the table
confirms that *every* amount of lengthening makes it worse. The surface is broad
and flat around the optimum (46–50 ears across a ±1 mm, ±10° neighbourhood), so
the change is not knife-edged.

Two caveats. The pose is frozen, so a shorter wing would also re-seat — this is a
lower bound. And retention is the noisiest axis in the suite (61/107 seed
agreement), so treat 31 → 50 as a direction and a magnitude, not a precise count.

---

## Visual verification

`viz_scene.py` exports one seated ear for eyeball checking, because every number
in this report comes from a signed-distance field and a pose matrix, and both of
this project's worst bugs — the nozzle-frame error and the faceplate-into-the-
skull roll — were invisible in the metrics until someone reasoned about the
geometry.

`cad/iem/viz/seated_scene.glb` (2.76 MB, 141 k triangles) holds **sonicom/P0023**
— protrusion 11.15 mm against a population median of 11.26, seal passing at the
2.5 mm budget, overall grade *marginal* — with the ear in tan, core and faceplate
in silver, jacket and wing in blue, skirt/carrier in orange and the nozzle insert
in grey. Parts carry exactly the transforms the try-on scores them with,
including the nozzle-local → assembly cant, so a wrong-looking render means wrong
numbers. `seated_scene_meta.json` carries the ear id, its full metric row, the
seal-compliance row, the landmarks and the 4×4 seating transform.

Sanity check on the export, which is also a check on the pipeline: distance from
each part to the ear surface comes out as skirt 0.04 mm (it is the sealing part),
jacket/wing 0.03 mm (it presses), core 1.53 mm and faceplate 2.47 mm (rigid, and
standing clear of flesh). That is precisely what the metrics claim.

No PNG renders: neither `pyglet` nor `pyrender` is installed, so headless
rasterisation was not trivially available and was skipped as instructed. The GLB
opens in any glTF viewer.

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

## Earlier revisions (condensed)

The full v1/v2 tables are in git history (`fe609b7`, `01ae359`, `13e75b3`,
`5ff456f`). What each step established:

| step | change | effect |
|---|---|---|
| v1 (`387b8eb` → `8db64d7`) | first 107-ear run | Found the body stacked along the nozzle axis: 32.6 mm long against a 10–18 mm concha. Protrusion failed 104/107. Recommended canting the nozzle. |
| v2 (`d740ddc`) | nozzle cant 45° | Protrusion median 14.19 → 12.88 mm; seal fails 64 → 44 (rigid metric); clearance passes 76 → 101. A **worst-case** fix: ≤18 mm protrusion went 74 % → 93 %, but the median barely moved. |
| v2 rescore (`13e75b3`) | protrusion thresholds ≤2/≤5 → ≤10/≤14 mm | The strict pair was saturated and had been *hiding* the cant's benefit (it read 103 → 107 fails; recalibrated, 56 → 38). Shortening requirement fell from ~12 mm to ~4–6 mm. |
| v2 seal rescore (`5ff456f`) | compliance-aware seal | Overturned "seal is the binding constraint" — that was an artefact of scoring a compliant skirt as rigid. Seal passes 10 → 67. Identified the intertragic notch as 77–88 % of residual seal failures. |
| v3 (`1ae768a`) | corner roll + 0.65 mm body trim + 1 mm notch sector | **The corner roll was the win**: protrusion median 12.88 → 10.52 mm and its fails 38 → 17, from removing the protruding *corner* rather than shrinking the body. Overall fails 45 → 27. The 1 mm notch sector was worth +2 ears (net-neutral). |

Two corrections made along the way, both of which moved the numbers materially:
the canted carrier was being read in the **wrong frame** (`carrier.stl` is written
in the nozzle-local frame; only `assembly.stl` applies `nozzle_T`), and the
seating search had an **unconstrained roll** that seated 23/107 ears with the
faceplate pointing into the skull until `c_face` was added.

---

## Compliance-aware seal scoring

`tryon.py`'s seal axis is **wrong about the physics in two ways**, and both make
it pessimistic. It treats the skirt as a rigid ring — scoring the fraction of the
Ø19 rim lying within 1.5 mm of flesh — when the real part is a 0.35 mm
Shore-A-10/15 silicone flare that folds and drapes over 2–4 mm, carried on a
mag-float carrier with **1.5 mm of axial travel** to seat it. And it scores
*coverage*, when what seals is a **closed loop**: a rim touching over 90 % of its
perimeter with one continuous 36° hole leaks, while one touching 88 % broken into
a dozen 3° specks does not.

`seal_compliance.py` re-scores the seatings `align_ear.py` produced, with:

- the rim sampled at **360 points** (1°; an 18° gap criterion needs finer
  resolution than the 72 points `tryon.py` uses);
- a **deformation budget** *B*: a sample is sealed where its signed distance
  ≤ *B*. Negative means the lip is already pressed into flesh; positive means the
  silicone must span that gap and can, within budget. (One-sided by design —
  a lip pressed into flesh seals. Over-burial is bounded by the clearance axis
  and by `c_soft` in the seating cost, not here.)
- the carrier's **1.5 mm axial travel** modelled explicitly as a 1-D search along
  the nozzle axis — a seating-*depth* search inside the seal metric, not a re-run
  of the 6-DOF pose optimiser;
- a **continuity** criterion: sealed over ≥ 95 % of the perimeter **and** largest
  single unsealed arc ≤ 18°. Both, because either alone is gameable as above;
- **the as-built rim**, sampled off `carrier.stl` as the max-radius locus per
  azimuth (`--rim mesh`, the default). Since v3 the rim is no longer a circle —
  the notch sector flares it — and an analytic ring samples straight past the
  feature it is meant to measure. `--rim circle` restores the old ring, which is
  how the notch contribution above was isolated.

### Seal pass counts (v4)

| skirt budget | model | seal pass | fail | pass rate | v3 |
|---|---|---|---|---|---|
| 1.5 mm | rigid | 39 | 68 | 36 % | 45 |
| **2.5 mm** | **conservative** | **67** | 40 | **63 %** | 76 |
| **4.0 mm** | **optimistic** | **98** | 9 | **92 %** | 103 |

Scored on the **as-built** (flared) rim. With a plain Ø19 circle on the identical
seatings the same three rows read 41 / 75 / 102 — which is the flare penalty
documented in [the verdict](#the-notch-flare-is-harmful--and-that-is-my-error).
`tryon.py`'s rigid coverage metric passes 13 on v4.

### Where the gain comes from (v4)

| skirt budget | pass, drape only | pass, drape + 1.5 mm travel | travel contributes |
|---|---|---|---|
| 1.5 mm | 21 | 39 | +18 |
| 2.5 mm | 58 | 67 | +9 |
| 4.0 mm | 92 | 98 | +6 |

The mag-float carrier's 1.5 mm of axial travel is worth +9 ears at the
conservative budget and +18 at the rigid one — a large share of the seal, and the
cheaper of the two mechanisms to guarantee. **Protect it in the tolerance
stack.** Median travel actually used: 1.40 mm at the 1.5 mm budget, 1.00 mm at
2.5 mm, 0 mm at 4.0 mm.

| skirt budget | median sealed arc | median worst gap |
|---|---|---|
| 1.5 mm | 92 % | 23° |
| 2.5 mm | 98 % | 7° |
| 4.0 mm | 100 % | 0° |

### The hole is at the intertragic notch

The intertragic notch — the soft gap between tragus and antitragus at the
inferior-anterior margin of the concha — is where the skirt has no cartilage wall
to land on. Classifying each failing ear by where its largest unsealed arc sits
(notch sector = ±45° about the inferior-anterior direction, i.e. **25 % of the
perimeter**):

| skirt budget | failing ears | gap at the notch | elsewhere | notch share | enrichment vs uniform |
|---|---|---|---|---|---|
| 1.5 mm | 62 | **51** | 11 | 82 % | **3.3×** |
| 2.5 mm | 31 | **25** | 6 | 81 % | **3.2×** |
| 4.0 mm | 4 | **3** | 1 | 75 % | **3.0×** |

**This holds in v3 exactly as in v2.** The notch occupies a quarter of the rim
but takes three-quarters to four-fifths of all seal failures. The v3 notch sector
did not change this, because — as shown above — it delivered 1.0 mm of realised
reach against a 2.0 mm median requirement.

The 4 ears still leaking at the optimistic 4.0 mm budget, 3 of them at the notch:

| ear | sealed arc | worst gap | at notch |
|---|---|---|---|
| sonicom P0040 | 94 % | 21° | no |
| hutubs pp89 | 92 % | 28° | yes |
| hutubs pp10 | 91 % | 34° | yes |
| synthetic xs_shallow | 83 % | 60° | yes |

`xs_shallow` — the 7 × 4.5 mm aperture on an 8 mm concha — remains a different
regime: a 60° hole is a sixth of the rim open, and no plausible skirt closes it.
That corner needs a smaller rim, not more compliance.

---

## Priorities after v4

Overall (seed 0): **7 pass / 64 marginal / 36 fail**; the remaining failures are
driven by protrusion (24), seal (9) and retention (7).

1. **Revert the 3.25 mm notch flare** to a plain Ø19 rim. It is net-negative at
   every budget on both seeds, and it breaks ears that previously sealed
   perfectly. This is a regression to undo, not a parameter to tune. Expected
   recovery: ~8 seal passes at the conservative budget.
2. **Shorten the wing 2.5–3.0 mm and splay it −10°.** Retention 31 → 50 ears in
   the spring band. Every amount of *lengthening* makes it worse; only one ear in
   107 is genuinely short. This is the largest single available gain.
3. **Pursue notch sealing as compliance, not geometry** — thinner wall, longer
   free length, or lower durometer over the inferior 90°, so the lip drapes
   toward the cartilage flanking the notch rather than reaching into the gap
   between. The deformation-budget sweep remains the right target; it was always
   a compliance target and was mis-translated into geometry at v3.
4. **Protrusion: the corner roll approach, not the body.** v3 established the
   corner roll buys several times more per mm than trimming the body, and the
   frozen-pose shrink estimate says a further 2 mm of body is worth only 0.46 mm.
   With the core now restored, protrusion sits at 11.26 mm median / 24 fails; if
   that is not acceptable, extend the roll.
5. **Watch the jacket.** Mean clearance 3.01 → 3.94 → 4.27 → 4.63 mm over three
   revisions. It costs no grades yet, but it is the surface the wing reacts
   against and it is drifting monotonically the wrong way.

## Pipeline changes in v4

- **`retention_analysis.py`** (new): buckets each ear as retained / short /
  misdirected / blocked / overpressed against the spring's working range, reports
  the aim angle between the wing's growth direction and the local surface normal,
  and sweeps wing extension × splay to find the change that retains the most
  ears. Takes `--tryon-csv` to separate the graded failures from the wider
  population.
- **`viz_scene.py`** (new): exports one seated ear + IEM as a coloured GLB plus a
  metadata JSON, for human verification of the poses the metrics are computed
  from. Includes a dependency-free vertex-clustering decimator, since the raw
  parts are ~1.25 M triangles and no quadric simplifier is installed.
- **`seal_compliance.py`**: gained `--json-dir`, so a scoring pass can be pointed
  at a saved copy of the seatings instead of racing a concurrent `--reseat`.

Carried forward: the as-built rim sampling (`--rim mesh`, essential now that the
rim is non-circular), the nozzle-frame fix, `c_face`, `c_soft`, `c_prot`,
clearance graded on `shell`, seeded sampling, Powell over the best 4 starts.

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

The cost surface is rough and Powell finds a local minimum, so results move
between sampling seeds. v4 under `--field-seed 0` vs `7`:

| metric | median \|Δ\| | p90 \|Δ\| |
|---|---|---|
| `rim_cover` | 0.10 | 0.31 |
| `rim_gap` | 0.77 mm | 2.28 mm |
| `wing_tip` | 0.77 mm | 2.60 mm |
| `jacket_mean` | 0.82 mm | 3.06 mm |
| `protrusion` | 1.32 mm | 5.41 mm |
| `hard_min` | 0.73 mm | 2.14 mm |

| axis | same grade | fail count, seed 0 → 7 |
|---|---|---|
| clearance | 101 / 107 | 0 → 0 |
| protrusion | 68 / 107 | 24 → 18 |
| seal (compliance) | 67 / 107 | 9 → 6 |
| retention | 61 / 107 | 7 → 5 |
| **overall** | **61 / 107** | **7/64/36 → 17/62/28** |

**v4 is the noisiest revision so far and its headline pair should not be
over-read.** Seal passes at the 2.5 mm budget move 67 → 79 between seeds, and
overall passes 7 → 17. That is why the conclusions in this report rest on
**paired, same-pose comparisons** wherever possible — flared rim versus circle
rim on identical seatings, and the wing sweep against each ear's own baseline.
Those are immune to seed noise, and both replicate on both seeds.

Directionally robust: clearance (fully clear under both seeds), the flare
penalty, the retention taxonomy's shape, and the notch attribution — which is a
geometric classification of *where* a gap sits rather than a threshold
comparison, and the most stable result in the suite.

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

.venv/bin/python generate.py --all
.venv/bin/python align_ear.py --reseat --jobs 7        # fresh seating, ~12 min
.venv/bin/python tryon.py                              # protrusion/retention/clearance
.venv/bin/python seal_compliance.py                    # compliance seal, as-built rim
.venv/bin/python seal_compliance.py --rim circle       # same, plain ring: flare delta
.venv/bin/python retention_analysis.py --tryon-csv ears/aligned/tryon.csv
.venv/bin/python shrink_estimate.py --case 0,0 --case 2,2
.venv/bin/python viz_scene.py                          # viz/seated_scene.glb
```

Run reseats **one at a time** — concurrent reseats interleave their poses in
`ears/aligned/*.json`. `--json-dir` on `seal_compliance.py`,
`retention_analysis.py` and `shrink_estimate.py` reads a saved copy of the
seatings if a reseat may be running. `--field-seed` varies the sampling for a
stability check; `--cant` / `--stl-dir` score a non-default build; `--qc-png`
writes a per-ear depth map with the detected aperture and tragus marked;
`--landmarks` overrides a bad pick by hand.
