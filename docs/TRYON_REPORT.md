# Magneto IEM — virtual fit validation ("try-on") · v3

Seats the generated IEM assembly in 107 real and synthetic ears and scores the
fit. **v3** scores `cad/iem/generate.py` at commit **`1ae768a`** — posterior-
inferior corner roll, body trim, and the intertragic-notch skirt sector — with a
**fresh 6-DOF seating optimisation on all 107 ears**, against the v2 build
(`d740ddc`, canted nozzle) it replaces.

Scoring throughout: protrusion thresholds **≤10 mm pass / ≤14 mm marginal**
([calibration](#protrusion-threshold-calibration)) and the **compliance-aware
seal** ([method](#compliance-aware-seal-scoring)).

## Verdict

**v3 is a clear step forward — overall failures fall from 45 to 27, and the
corner roll did the work. But the intertragic-notch sector, the change aimed
squarely at the seal, delivered almost nothing: it is aimed correctly and built
about 1 mm short of what the ears need. And the 8 mm driver option is not
worth taking on protrusion grounds.**

| | v2 (`d740ddc`) | v3 (`1ae768a`) |
|---|---|---|
| overall pass / marginal / fail | 10 / 52 / 45 | **21 / 59 / 27** |
| protrusion median | 12.88 mm | **10.52 mm** |
| protrusion p90 / max | 17.83 / 20.26 mm | **15.37 / 19.41 mm** |
| seal pass @ 2.5 mm budget | 67 | **76** |
| seal pass @ 4.0 mm budget | 99 | **103** |
| axes driving the fails | protrusion 38, seal 8, retention 6 | protrusion 17, retention 8, seal 4 |

Read the pass column with care: **the fail count is robust, the pass count is
not.** Under a second sampling seed v3 reads 13 / 68 / 26 — failures move by one
ear, passes by eight, because a large marginal population sits near the
boundaries. The 45 → 27 improvement is far outside that noise; "21 pass" is not.

Four findings:

1. **The corner roll is the win.** Protrusion median drops 2.36 mm (12.88 →
   10.52) and its failures fall 38 → 17, even though the body trim itself
   delivered only 0.65 mm. Removing the *actual protruding corner* was worth
   several times more than shrinking the body — as the v2 sensitivity analysis
   predicted, because the protruding point was a corner diagonal to all three
   axes.
2. **The notch sector did essentially nothing: +2 ears.** Scoring the same v3
   seatings with the as-built rim versus a plain circle isolates the flare's
   contribution: 74 → 76 passes at the conservative budget, 102 → 103 at the
   optimistic. The v2 → v3 seal gain (67 → 76) is almost entirely *reseating*
   against the new body, not the notch feature.
   [Why, and what to do](#did-the-notch-sector-work-no-and-the-reason-is-specific).
3. **Retention is now the second-largest problem** (8 fails, ahead of seal's 4)
   and it is the noisiest axis. It has been drifting the wrong way as the body
   moves under it; it needs its own pass.
4. **A 2 mm smaller shell buys almost nothing.** Median protrusion 10.52 →
   10.06 mm, pass rate 43 % → 47 %.
   [Driver-downsizing estimate](#would-a-smaller-driver-help-barely).

---

## v2 → v3, side by side

Both scored identically; v3 has its own fresh seating optimisation.

### Overall (compliance seal + recalibrated protrusion)

| | pass | marginal | fail |
|---|---|---|---|
| v2 (`d740ddc`) | 10 | 52 | 45 |
| **v3 (`1ae768a`)** | **21** | **59** | **27** |
| v3, second seed | 13 | 68 | 26 |

### Per axis

| axis | v2 pass / marg / fail | v3 pass / marg / fail | change |
|---|---|---|---|
| protrusion (≤10 / ≤14 mm) | 21 / 48 / 38 | **46 / 44 / 17** | **−21 fails** |
| seal (compliance, 2.5 / 4.0 mm) | 67 / 32 / 8 | **76 / 27 / 4** | **−4 fails** |
| clearance | 101 / 5 / 1 | **105 / 2 / 0** | −1 fail |
| retention | 48 / 53 / 6 | 55 / 44 / **8** | +2 fails |

### Metric medians

| metric | v2 | v3 |
|---|---|---|
| faceplate past tragus | 12.88 mm | **10.52 mm** |
| skirt rim contact coverage (rigid) | 0.54 | **0.61** |
| skirt rim max gap (rigid) | 2.76 mm | **2.96 mm** |
| wing tip signed distance | −1.16 mm | −0.73 mm |
| jacket mean clearance | 3.94 mm | 4.27 mm |
| worst rigid interference | −0.27 mm | **+0.26 mm** |

Clearance is now fully clear — the median ear has *no* rigid interference at all
(+0.26 mm) and nothing fails. The jacket continues to drift off the concha floor
(3.94 → 4.27 mm), the same slow regression flagged in v2; it is not yet costing
grades but it is the thing to watch.

### Worst 5 ears (v3)

| ear | failing axes | cover | max gap | wing tip | protrusion |
|---|---|---|---|---|---|
| hutubs pp63 | protrusion + retention + seal | 0.49 | 4.04 | +1.45 | 17.51 |
| synthetic xl_shallow | retention + seal | 0.39 | 1.47 | +3.46 | 10.18 |
| hutubs pp21 | protrusion + seal | 0.46 | 2.27 | −1.38 | 17.12 |
| sonicom P0016 | protrusion + retention | 0.83 | 1.31 | +2.09 | 15.81 |
| hutubs pp47 | retention + seal | 0.36 | 0.22 | +1.08 | 6.03 |

Only one ear now fails three axes (v2's worst failed two, v1's four — but v1/v2
were scored differently). `xl_shallow` and `pp47` fail on retention and seal at
*low* protrusion, which is the signature of the remaining work: the shell now
fits, and what is left is holding it in and closing the lip.

---

## Did the notch sector work? No — and the reason is specific

The v3 skirt carries a 90° inferior sector with +1.75 mm of specified reach and a
0.22 mm wall. Scoring the **same v3 seatings** with the as-built rim sampled off
`carrier.stl` versus a plain Ø19 circle isolates exactly what that feature buys:

| skirt budget | circle rim (no flare) | as-built rim (with flare) | flare contributes |
|---|---|---|---|
| 1.5 mm | 44 | 45 | **+1** |
| 2.5 mm | 74 | 76 | **+2** |
| 4.0 mm | 102 | 103 | **+1** |

**Two ears.** Two candidate explanations, and the data separates them cleanly.

**It is not aimed wrong.** Measuring the built flare's angular position against
the anatomical notch direction at every seated pose: median deviation **18°**,
p90 36°, and only **2 of 107** ears put the flare peak outside the ±45° notch
sector. The generator put the feature where the anatomy is.

**It is built too short.** The flare's radial extension measured off the mesh is
**1.00 mm** (rim radius 9.52 → 10.52 mm), against a 1.75 mm specification — the
difference is the 35° cone slant and marching-cubes rounding of the peak. Now
compare that with what the ears actually require. Sweeping the deformation
budget per ear to find the smallest reach that seals it:

| additional reach required | mm |
|---|---|
| p25 | 0.9 |
| **median** | **2.0** |
| p75 | 2.7 |
| p90 | 3.4 |
| max | 5.2 |

| reach available | ears sealed |
|---|---|
| 1.0 mm (as built) | 28 / 107 (26 %) |
| 2.0 mm | 59 / 107 (55 %) |
| 2.5 mm | 76 / 107 (71 %) |
| **3.0 mm** | **90 / 107 (84 %)** |
| **3.5 mm** | **100 / 107 (93 %)** |
| 4.0 mm | 103 / 107 (96 %) |

The median ear needs **2.0 mm** and the p90 ear **3.4 mm**; the sector delivers
**1.0 mm**. That is why it bought two ears — it is on target but roughly a third
of the required depth. And the notch remains where the failures are: **82 % /
81 % / 75 %** of failing ears at the three budgets have their largest gap in the
notch sector, which is a quarter of the perimeter (3.0–3.3× enrichment).

**Recommendation: take the notch sector from ~1 mm of realised radial reach to
3.0–3.5 mm** — measured on the mesh, not specified on the slant, because the two
differ by nearly 2× here. That is worth roughly **+14 to +24 ears** on the seal
axis by the table above, and it is the single highest-value change on the board.
Verify the built extension by measuring `carrier.stl` rather than trusting the
parameter.

---

## Would a smaller driver help? Barely

A first-order estimate for the 8 mm driver decision, at frozen pose
(`shrink_estimate.py`): displace the faceplate inward by the shrink and re-take
the protrusion maximum. Y shrinks symmetrically (half per side); Z comes off the
+Z/faceplate side, since the −Z jacket face still has to bed on the concha floor.

| shell shrink | protrusion median | p90 | pass ≤10 mm | ≤14 mm | fail >14 mm |
|---|---|---|---|---|---|
| none (as built) | 10.52 mm | 15.37 | 43 % | 84 % | 16 % |
| −1 mm Y, −1 mm Z | 10.18 mm | 14.93 | 48 % | 85 % | 15 % |
| **−2 mm Y, −2 mm Z** | **10.06 mm** | **14.79** | **47 %** | **87 %** | **13 %** |
| −3 mm Y, −3 mm Z | 10.26 mm | 14.46 | 48 % | 90 % | 10 % |

**A 2 mm shrink moves median protrusion 0.46 mm and the pass rate 43 % → 47 %.**
The curve is flat and even non-monotone, for the same reason the corner roll
worked: protrusion is a *max over a corner*, so shrinking the body's faces slides
the maximum onto a neighbouring point instead of removing it. The gain shows up
mostly in the tail (fails 16 % → 13 %).

**On protrusion grounds alone, the 8 mm driver is not justified.** Note the
caveat cuts the other way too: the pose is frozen at the full-size shell's
seating, so this is a lower bound — a genuinely smaller shell would also re-seat
deeper. If the driver is attractive for acoustic or packaging reasons, this
analysis does not argue against it; it just says protrusion will not pay for it.
Getting the same 3 % of tail from the notch/retention work is far cheaper.

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

### Seal pass counts (v3)

| skirt budget | model | seal pass | fail | pass rate | v2 |
|---|---|---|---|---|---|
| 1.5 mm | rigid | 45 | 62 | 42 % | 45 |
| **2.5 mm** | **conservative** | **76** | 31 | **71 %** | 67 |
| **4.0 mm** | **optimistic** | **103** | 4 | **96 %** | 99 |

For comparison, `tryon.py`'s rigid coverage metric passes 18 on v3 (10 on v2).
Note the 1.5 mm row already passes 45 — most of that gap against the rigid metric
is the *criterion* (continuity instead of a 0.75 coverage fraction), before any
extra compliance is granted.

### Where the gain comes from (v3)

| skirt budget | pass, drape only | pass, drape + 1.5 mm travel | travel contributes |
|---|---|---|---|
| 1.5 mm | 32 | 45 | +13 |
| 2.5 mm | 62 | 76 | +14 |
| 4.0 mm | 99 | 103 | +4 |

The mag-float carrier's 1.5 mm of axial travel is still worth ~14 ears on its own
at the conservative budget — comparable to the silicone's contribution and the
cheaper of the two to guarantee. **Protect it in the tolerance stack.** Median
travel actually used falls from 1.20 mm at the 1.5 mm budget to 0.30 mm at
2.5 mm and 0 mm at 4.0 mm.

| skirt budget | median sealed arc | median worst gap |
|---|---|---|
| 1.5 mm | 92 % | 25° |
| 2.5 mm | 99 % | 2° |
| 4.0 mm | 100 % | 0° |

At the conservative budget the median v3 ear is sealed over 99 % of its rim with
a 2° worst gap. The seal problem is not population-wide — it is a minority of
ears with one large hole each, and that hole is at the notch.

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

## Priorities after v3

Overall: **21 pass / 59 marginal / 27 fail**, with the remaining failures driven
by protrusion (17), retention (8) and seal (4).

1. **Deepen the notch sector to 3.0–3.5 mm of *realised* radial reach**
   (currently 1.0 mm). Worth ~+14 to +24 ears on the seal axis — the largest
   single gain available, and the feature is already in the right place. Measure
   the built flare on `carrier.stl`; the specified 1.75 mm became 1.00 mm on the
   mesh, so specification alone is not evidence.
2. **Give retention its own pass.** It is now the second-largest failure source
   (8) and has drifted the wrong way across every revision as the body moved
   under it (v2 6 → v3 8 fails). Median `wing_tip` is −0.73 mm against a −0.5 to
   −2.0 mm target, so the median ear is *in* band — the failures are spread, and
   it is also the noisiest axis (58/107 seed agreement). Diagnose before
   redesigning.
3. **Keep chipping at protrusion via the corner, not the body.** The corner roll
   bought 2.36 mm of median where a 0.65 mm body trim bought little, and the
   shrink estimate says a further 2 mm of body is worth only 0.46 mm. If more is
   needed, extend the roll.
4. **Watch the jacket.** Mean clearance has drifted 3.01 → 3.94 → 4.27 mm across
   v2's cant and v3's trim. It costs no grades yet, but the gyroid skin is
   progressively less bedded on the concha floor, which is what the wing reacts
   against — a plausible contributor to (2).
5. **Do not buy the 8 mm driver for protrusion.** See above; it is worth ~4
   percentage points of pass rate.

## Pipeline changes in v3

One change, and it matters for reading the notch result:

**The rim is no longer a circle, so it is no longer sampled as one.** The
intertragic-notch sector flares the skirt ~1 mm further out over the inferior
90°. `seal_compliance.py` previously sampled an analytic Ø19 ring, which walks
straight past that flare and would have scored the feature as worth exactly
nothing by construction. It now samples the **as-built rim** off `carrier.stl` —
the max-radius locus per azimuth, in the nozzle-local frame, canted into place —
so it tracks whatever the generator actually produces. `--rim circle` keeps the
old ring, and the difference between the two is precisely the flare's
contribution, which is how the +2-ear figure above was measured.

Also added `shrink_estimate.py` (frozen-pose shell-shrink estimator) and a
`--json-dir` flag on it, so an analysis can be pointed at a saved copy of the
seatings rather than racing a concurrent `--reseat` over `ears/aligned/`.

Corrections from earlier revisions still in force: the nozzle-frame fix, `c_face`
(faceplate must face out of the head), `c_soft`, `c_prot`, clearance graded on
`shell` rather than `rigid`, seeded sampling, Powell refinement over the best 4
starts. See [earlier revisions](#earlier-revisions-condensed).

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
between sampling seeds. Re-seating v3 under `--field-seed 7`:

| metric | median \|Δ\| | p90 \|Δ\| |
|---|---|---|
| `rim_cover` | 0.07 | 0.23 |
| `rim_gap` | 0.55 mm | 1.72 mm |
| `wing_tip` | 0.76 mm | 3.14 mm |
| `jacket_mean` | 0.77 mm | 3.20 mm |
| `protrusion` | 1.14 mm | 4.95 mm |
| `hard_min` | 0.51 mm | 2.09 mm |

| axis | same grade | fail count, seed 0 → 7 |
|---|---|---|
| clearance | 99 / 107 | 0 → 0 |
| seal (compliance) | 86 / 107 | 4 → 2 |
| protrusion | 73 / 107 | 17 → 21 |
| retention | 58 / 107 | 8 → 7 |
| **overall** | **73 / 107** | **21/59/27 → 13/68/26** |

Seal pass counts across seeds: 45 → 51 at the 1.5 mm budget, **76 → 79** at
2.5 mm, 103 → 105 at 4.0 mm.

**Read the fail counts, not the pass counts.** Total failures are stable
(27 → 26) and every axis's fail count moves by ≤4, so the v2 → v3 improvement
(45 → 27) and the seal-budget curve are far larger than the noise. The *pass*
count is not stable (21 → 13): a large marginal population sits near the
boundaries and small pose changes tip ears across. Retention remains the noisiest
axis (58/107 agreement) and its counts are directional only. The notch
attribution — which is a geometric classification of where a gap sits, not a
threshold comparison — is the most robust result here.

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

# v3 as scored here
.venv/bin/python generate.py --all
.venv/bin/python align_ear.py --reseat --jobs 7      # fresh seating, ~10 min
.venv/bin/python tryon.py                            # protrusion, retention, clearance
.venv/bin/python seal_compliance.py                  # compliance seal, as-built rim
.venv/bin/python seal_compliance.py --rim circle     # same, analytic ring
.venv/bin/python shrink_estimate.py --case 0,0 --case 2,2
```

Run reseats **one at a time** — concurrent reseats interleave their poses in
`ears/aligned/*.json`. Pass `--json-dir` to `shrink_estimate.py` to read a saved
copy if a reseat may be running. `--field-seed` varies the sampling for a
stability check; `--cant` / `--stl-dir` score a non-default build (e.g.
`generate.py --all --cant 0 --out stl_cant0 --ear right`); `--qc-png` writes a
per-ear depth map with the detected aperture and tragus marked; `--landmarks`
overrides a bad pick by hand.
