# Magneto IEM — virtual fit validation ("try-on") · v5 (short-list) + enclosure audit

**Two things in this revision, and the first one matters more than the second.**

1. **A scoring hole was found by human review and is now closed.** The seal metric
   never checked that the sealed contact loop *surrounds the canal entrance*. On
   the v4 seatings the canal aperture sat **outside** the rim ring on more than
   half of the audited ears — a rim laid flat on the concha floor scored a perfect
   continuous seal while sealing nothing. A bore-over-canal constraint is now part
   of the seating cost. **Seal figures in v1–v4 are inflated and should not be
   quoted.**
2. **v5 geometry** (`dbc3b7b`: notch reverted to a plain Ø19 rim with a
   compliance-only sector, wing shortened to a 10.40 mm free span with −10° splay,
   k = 0.294 N/mm) was scored on a **13-ear short list**, not the full 107.

> ### Scope: this is a SHORT-LIST run (n = 13)
> Per instruction, v5 was not run on all 107 ears. The list is
> **sonicom/P0023** (the human-reviewed ear), the **4 synthetic envelope
> corners**, 2 median-protrusion ears (pp66, pp82), 2 median-seal-behaviour ears
> (pp67, P0003), 2 small-aperture-tail ears (pp69, P0016) and 2 large-aperture-tail
> ears (pp49, pp9). Population medians on the full 107 are protrusion 11.26 mm and
> basin-inscribed radius 7.3 mm; the list brackets both.
>
> **What n = 13 can support:** the enclosure audit (a geometric property, and the
> effect is large — 38 % → 85 %), the direction and rough size of the seal change,
> and the qualitative retention finding. **What it cannot support:** any P/M/F
> count as a population estimate, or a 1–2 ear difference anywhere. Differences of
> 3 ears out of 13 are inside the seed noise this suite already shows.

## Verdict

**Do not freeze for the variant-B comparison yet.** Three reasons, in order:

1. **The seal baseline just moved.** With enclosure enforced, seal at the
   conservative budget goes from 6/13 to **11/13** on identical geometry — but
   that is a *different measurement*, not an improvement in the part. Every prior
   seal number was computed without checking that the loop contained the canal.
   The design's real seal performance is currently unknown at population scale.
2. **The predicted retention gain did not materialise.** The v5 wing was expected
   to take overpressed 47 → ~28 and retained 31 → ~50 (full-107 terms). On the
   short list overpressed did fall, but retention did **not** rise: the shortening
   overshot, converting over-press into *no contact*. `wing_tip` median moved
   −1.09 → −0.20 mm and the count of ears where the wing does not touch at all
   rose 3 → 5 (seed 0) and → 7 (seed 7).
3. **Protrusion regressed under the new constraint**, median 11.26 → 13.40 mm on
   the list, because holding the bore over the canal pushes the body out. That
   trade is real and was previously hidden.

A full-107 run with the enclosure constraint is the minimum before freezing.

---

## The enclosure audit

Human review of `viz/seated_scene.glb` reported the nozzle looked like it was
sealing against flesh *near* the aperture rather than *around* it. That is
exactly the failure `seal_compliance.py` could not see: it asks only whether the
rim forms a continuous contact loop against flesh, never *where*. The concha
frame seeds the rim centre on the aperture and the Powell search was then free to
walk it away — and nothing in the cost objected, so it did.

`seal_enclosure_audit.py` adds three per-ear checks:

- **(a) aim** — angle between the nozzle axis and the rim-centre→aperture
  direction (0° = bore points straight at the canal), plus the angle to the
  inward concha normal for context;
- **(b) offset** — in-plane distance from the rim centre to the aperture,
  comparable against the rim radius;
- **(c) enclosure** — the aperture must project *inside* the rim ring **and** the
  contact loop must be continuous around it.

### Result (13 ears, 2.5 mm budget)

| | v4 seatings | v5, no constraint | **v5 + constraint** |
|---|---|---|---|
| (a) aim, median | 76.8° | 71.5° | **52.9°** |
| (b) lateral offset, median | 10.2 mm | 5.7 mm | **4.6 mm** |
| (b) lateral offset, max | 16.6 mm | 12.8 mm | **7.3 mm** |
| aperture inside the rim ring | 6/13 (46 %) | 9/13 (69 %) | **13/13 (100 %)** |
| contact loop continuous | 7/13 | 6/13 | 11/13 |
| **(c) loop actually surrounds the canal** | **6/13 (46 %)** | **5/13 (38 %)** | **11/13 (85 %)** |
| sealed but NOT enclosing (false pass) | 1/13 | 1/13 | **0/13** |

On v4 the median lateral offset (10.2 mm) was **larger than the rim radius
itself** (9.88 mm) — the sealed ring was typically centred a full radius away
from the canal it was supposed to surround. This is a genuine scoring hole, and
it inflated every seal number this report has published.

### Was P0023 specifically mis-aimed?

**Partly — the reviewer's instinct was right, though not in the way the audit
first framed it.** On the v4 seating that produced the GLB, P0023's aperture does
lie inside the rim ring (lateral offset 4.55 mm against a 9.88 mm radius), so the
loop *does* enclose the canal and it is not a false pass. But the bore is 4.5 mm
off-axis and the aperture sits **1.44 mm behind the rim plane** (`axial` −1.44) —
the rim has been pushed slightly past the canal mouth into the surrounding flesh,
with the nozzle pointing across the aperture rather than down it (aim 107.6°).
So: sealing *around* the canal, but not aimed *into* it. Under the constraint,
P0023's aim improves to 94.7° and axial goes positive.

### The constraint

Added to `seating_cost` as `c_aim`: 3 mm of lateral slack is free (real anatomy
offsets the bore from the bowl centre), quadratic beyond, weight 0.40. Toggle it
off with `align_ear.py --no-aim` to reproduce the unconstrained numbers.

---

## v5 short-list results (n = 13)

### Overall and seal

| run | overall P/M/F | seal @1.5 / 2.5 / 4.0 | protrusion median |
|---|---|---|---|
| v4, same 13 ears | 2 / 6 / 5 | 5 / 7 / 11 | 11.26 mm |
| v5 no constraint, seed 0 | 0 / 6 / 7 | 6 / 6 / 11 | 12.70 mm |
| v5 no constraint, seed 7 | 1 / 8 / 4 | 6 / 8 / 11 | 11.59 mm |
| **v5 + constraint, seed 0** | **2 / 4 / 7** | **9 / 11 / 12** | **13.40 mm** |
| **v5 + constraint, seed 7** | **1 / 6 / 6** | **8 / 9 / 12** | **12.57 mm** |

Seal at the conservative budget nearly doubles under the constraint (6 → 11 at
seed 0, 8 → 9 at seed 7) — because a bore held over the canal puts the rim in the
concha bowl, which has a continuous wall to seal against, instead of skidding
onto the floor. Protrusion pays for it (+0.7 to +2.2 mm of median). Overall P/M/F
is flat-to-worse, and at n = 13 with this suite's seed noise **none of the
overall counts is a meaningful difference.**

### Notch compliance sector — sensitivity line (clearly labelled)

The compliance-only sector cannot be scored kinematically as extra reach; its
benefit appears as a **larger achievable budget over the inferior 90°**. Scored
three ways on the constrained v5 seatings:

| budget model | ears sealed (of 13) |
|---|---|
| uniform 2.5 mm (conservative) | 11 |
| **2.5 mm base + 4.0 mm in the notch sector** | **12** |
| uniform 4.0 mm (optimistic, sensitivity bound) | 12 |

The sector captures the **entire** benefit of a uniform 4.0 mm budget while only
relaxing the inferior 90°, which is the right shape for the design. Note its
value depends on the base: with the unconstrained seatings the same three rows
read 6 / 10 / 11, so the sector is worth far more when the underlying seal is
poor. Treat the 4.0 mm column as a sensitivity bound, not a prediction.

### Retention — the shortening overshot

Same-pose paired comparison (identical v4 seatings, only the wing swapped), plus
v5 at its own seatings:

| | retained | short | misdirected | blocked | overpressed |
|---|---|---|---|---|---|
| v4 poses + v4 wing | **5** | 0 | 0 | 3 | **5** |
| v4 poses + v5 wing (paired) | 2 | 2 | 2 | 1 | 6 |
| v5 poses + v5 wing | 4 | 0 | 3 | 2 | **4** |

`wing_tip` median and contact state:

| run | wing_tip median | in band [−1.5, 0] | over-pressed < −1.5 | **not touching > 0** |
|---|---|---|---|---|
| v4, same 13 | −1.09 mm | 5 | 5 | 3 |
| v5 + constraint, seed 0 | −0.20 mm | 4 | 4 | 5 |
| v5 + constraint, seed 7 | +0.57 mm | 4 | 2 | **7** |

**Over-press did improve — but into the wrong place.** The wing went from
pressing too hard to not reaching at all; ears in the retaining band did not
increase. The v4 prediction (retained 31 → ~50 on 107) is not supported.

That prediction came from a crude proxy and I should flag why it misled: the v4
sweep *translated the wing-tip sample points* along the growth axis and rotated
them about the root. A real redesign changes the blade's shape, where its distal
surface lands, and which vertices are "the tip" — so the sweep was directionally
useful (shorten, do not lengthen) but not quantitatively predictive.
**Recommendation: split the difference — restore roughly 1–1.5 mm of the
2.4 mm that was removed**, and re-measure on the full 107 rather than a sweep.

---

## The contact contract

The graded axes say how *well* the fit scores. The contract says whether each
part is doing **the job it was designed for** — a different question, and the one
that catches a part touching flesh it should never touch. Each rule is an
explicit design intent, not a derived metric, and `tryon.py` prints the table on
every run (automatically for a single ear, or with `--contract`).

| part | intent | rule |
|---|---|---|
| skirt land | **MUST TOUCH** | continuous contact band that **encloses** the canal aperture; enclosure is hard |
| wing / rail pad | **MUST TOUCH** | antihelix, tip within the spring's working range, −1.5 to 0 mm |
| jacket ear-face | **MUST REST** | in contact with the concha floor, never more than 2.5 mm into flesh |
| nozzle + insert | **MUST NOT TOUCH** | recessed inside the skirt; any contact loads the canal wall directly |
| core / faceplate | **MUST NOT TOUCH** | load reaches the ear only through the jacket and wing |
| cable exit | **MUST CLEAR** | the ear — checkable since `17810fd` added the boot |
| plunger pads | **MUST TOUCH** | every pad reaches its site, within cam preset + spring travel |
| **STABILITY** | **MUST RESIST** | quasi-static force/moment balance against the spec loads — retention *under load*, not just contact |

### The stability row

Contact says each part touches what it should. It says nothing about whether the
assembly **stays put when something pulls on it**, which is a different failure
and the one a wearer notices. `stability.py` is a rigid-body force/moment balance
— screw theory with Coulomb friction cones, one LP per sampled load direction —
not FEA.

**Loads to resist**

| | load | direction |
|---|---|---|
| (a) skirt preload reaction | 0.31 N | outward along the nozzle axis, **always on**, never scaled |
| (b) cable tug | 0.50 N | worst direction in a 45° downward-backward cone |
| (c) inertial | 3 g × 8 g = 0.24 N | worst direction over the sphere |

**Resistance.** Each contact pushes along the ear's outward surface normal (flesh
pushes, never pulls) with N ≥ 0 and friction |t| ≤ μN, μ = 0.40. The cone is
linearised to an 8-sided inscribed pyramid (conservative). Normal budgets are
**caps**, set by what presses each contact: skirt ≤ 0.31 N (its compression),
wing ≤ k × interference with k = 0.294 N/mm, jacket a free reaction, plungers
0.18–0.49 N each *when the build has any*. Geometric interlock needs no special
term — it falls out of contact normals opposing the escape direction. The score
is the largest load scale s\* that stays feasible; s\* ≥ 1 resists the spec.

Two modelling notes worth arguing with. **μ = 0.40** for silicone/Ti on skin is a
deliberately dry, conservative pick — published skin-on-elastomer values span
~0.3–1.0, and raising it makes everything easier. And a first pass wrote the
preloads as *equalities* rather than caps, which over-constrains the balance and
returned margin 0.00 on every pose; caps are correct, because a preload bounds
how hard a contact *can* push, it does not force it to be loaded.

### P0023 on `17810fd` (tripod plungers + cable boot) — 5 pass, 3 fail

Freshly re-seated on the plunger build. Every row is now evaluable: the cable
boot exists, so `cable exit` is checkable for the first time.

| part | intent | value | | detail |
|---|---|---|---|---|
| skirt land | MUST TOUCH | 92 % / 17° | **FAIL** | band 92 % closed against a 95 % rule; worst gap 17° is *inside* the 18° tolerance. Aperture 3.1 mm off centre vs a 9.4 mm rim, so the loop does enclose the canal |
| wing / rail pad | MUST TOUCH | −0.58 mm | PASS | **see caveat below** |
| jacket ear-face | MUST REST | −1.22 mm | PASS | resting, within 2.5 mm |
| nozzle + insert | MUST NOT TOUCH | +0.29 mm | PASS | recessed, clear |
| core / faceplate | MUST NOT TOUCH | +0.63 mm | PASS | no direct flesh load |
| **plunger pads** | MUST TOUCH | **−14.64 mm** | **FAIL** | cymba −0.46, antihelix_undercut +0.36, **antitragus −14.64 NEVER REACHES**; window [−3.00, +3.75] |
| **STABILITY** | MUST RESIST | **0.04×** | **FAIL** | capacity 0.27 N vs demand 1.05 N; friction budget 0.41 N; 8 of 48 contacts interlock |
| cable exit | MUST CLEAR | +5.61 mm | **PASS** | boot clears the ear comfortably |

**Caveat on the wing row.** `wing_style` is now `"plungers"`, so there is no wing.
That row is reading the 40 most-superior vertices of `jacket_wing`, which is now
just the jacket's upper edge — it passes, but it is not measuring a retention
feature any more. It should be retired or repointed for this build.

### Does the 20.3 mm pad reach overshoot P0023? No — and that is not the problem

Pad-tip reach from the core centre is 20.29 / 20.42 / 20.17 mm (cymba /
antihelix_undercut / antitragus). Measured against P0023's actual surface at the
seated pose:

| pad | overshoot into flesh | verdict |
|---|---|---|
| cymba | **−0.46 mm** (0.46 short) | on target; trivially dialled in |
| antihelix_undercut | **+0.36 mm** proud | on target; trivially dialled back |
| antitragus | **−14.64 mm** (14.6 short) | **unreachable** |

**Two of the three pads land within half a millimetre of the ear — the reach is
well chosen, and the cam's 3 mm is far more adjustment than they need.** The
failure is not overshoot; it is that the *antitragus* pad misses this ear by
14.6 mm, which is nearly 5× the cam range and cannot be dialled out. The seated
GLB shows the same thing independently: the cymba and antihelix pads sit 0.02 and
0.04 mm off the ear, while every antitragus part is 7.8–10 mm away.

So the tripod is effectively a **bipod** on this ear, which is also why stability
barely moved: capacity rose only 0.20 → 0.27 N against 1.05 N of demand, with 8
of 48 contacts interlocking. Two working pads at the low end of the force band
add ~0.36 N of normal, and the antitragus contribution is simply absent.

The antitragus aim is (0.20, −0.96, −0.20), i.e. almost straight inferior. On
P0023 the surface in that direction falls away far faster than the other two
sites. Worth checking whether that is P0023-specific or general before re-aiming
— it is the obvious next single-ear question.

## Deriving the third plunger aim from data

`plunger_aim_search.py` sweeps candidate aims on a 10° grid and, for each ear,
casts a ray from where the generator would place the boss (`core_c` → core
ellipsoid → jacket outer surface) to find how far the cartilage actually is.
Geometry queries only — no re-seating; each ear keeps the pose it had. Run over
the 13-ear short list plus P0023 on its current `17810fd` pose.

Two windows are reported because they differ: the **as-built** leg spans
8.40–12.15 mm from the jacket surface (11.40 mm stack, −3.0 cam, +0.75 travel),
while the stated **target** working range is 1.5–6.0 mm, which describes a
considerably shorter leg.

### (a) The best aims

| aim | in target | in as-built | interlocking | median | spread p10–p90 | misses |
|---|---|---|---|---|---|---|
| `antitragus` (+0.20, −0.96, −0.20) **as built** | 1/13 | **0/13** | 0/13 | 23.5 mm | — | **10/13** |
| best inferior, target window (−0.09, −0.50, −0.86) | **6**/13 | — | 1/13 | 5.4 mm | 1.6–10.1 | 4/13 |
| **best anterior (+0.82, −0.17, −0.54)** | — | **7**/13 | 0/13 | 9.7 mm | 1.3–15.8 | — |
| (0, 0, −1) — straight at the concha floor | 8/13 | — | 1/13 | 3.2 mm | 0.7–10.9 | — |

The nominal winner on raw count, (0, 0, −1), is **degenerate**: −Z is medial, so
that leg presses on the concha floor where the jacket already rests. It earns no
interlock and adds nothing. Ranking is therefore restricted to directions with a
real inferior component, and reported separately for the anterior region.

**The as-built antitragus aim is bad across the population, not just on P0023:
it misses 10 of 13 ears entirely and lands in the usable window on none.**

### (b) Is there a good common inferior direction? No

The best inferior aim reaches only **6 of 13** ears, misses 4 outright, and
spreads 1.6–13.2 mm — about 12 mm of variation across a 3.75 mm adjustment
range. Interlock is essentially absent (0–1 of 13) for *every* inferior
direction tried. **The lower ear is too variable to serve with a fixed leg.**

The suggested tragus-inside-face aim ≈ (−0.8, −0.3, +0.5) **misses all 13 ears**,
and (+0.8, −0.3, +0.5) misses 11 — because in this design frame **+Z is lateral,
i.e. out of the ear**, so any aim with a positive Z component leaves the concha
before it can hit anything. The idea is sound with the sign corrected: the
anterior-medial direction **(+0.82, −0.17, −0.54)** aims at the tragus inner wall
and reaches **7 of 13** ears inside the as-built window, against the antitragus
leg's 0 — the best result found, and it needs no change to leg length.

### (c) Per-ear spread, best inferior aim (−0.09, −0.50, −0.86)

| ear | distance | in window |
|---|---|---|
| hutubs/pp82 | 1.61 mm | yes |
| sonicom/P0023 | 1.65 mm | yes |
| hutubs/pp66 | 2.27 mm | yes |
| sonicom/P0016 | 4.45 mm | yes |
| hutubs/pp9 | 5.41 mm | yes |
| hutubs/pp49 | 5.50 mm | yes |
| hutubs/pp69 | 7.25 mm | no |
| synthetic/xl_deep | 9.35 mm | no |
| synthetic/xl_shallow | 13.15 mm | no |
| pp67, P0003, xs_deep, xs_shallow | no hit | no |

Median 5.41 mm, p10–p90 1.6–10.1 mm. Recommendation: **drop the antitragus leg
and move leg 3 anterior to ≈ (+0.82, −0.17, −0.54)**, accepting that no fixed
third aim serves the whole population.

### Does the cymba pad get interlock credit? Yes structurally — but it rarely earns it

**The model does credit interlock.** `stability.py` builds each contact wrench
from the ear's own outward surface normal, so a pad seated under an overhang has
n · pull_out < 0 and resists pull-out *directly*, with no friction required. No
special term is needed and none is missing.

**But the cymba pad does not land on the overhang.** Measured at the pad tip
across 13 ears, n · pull_out is **+0.41 median** — the surface faces *outward* and
would help eject the shell — and only **4 of 13** ears give it interlock at all.

This is a genuine underestimate rather than a true negative, because the overhang
is present in the data: within 8 mm of the cymba pad the ear meshes carry a
**median 16.9 %** of area with n · pull_out < −0.2 (up to 40.9 %; 12.6 % over the
whole patch). The cartilage lip is there — the pad tip is simply beside it.

Two concrete sources of under-credit, both fixable:

1. **Tip-only sampling.** `iem_points` samples the pad as a disc at its tip. A pad
   tucked under a lip contacts on its shoulder, and that contact is never
   sampled, so its interlock cannot be counted. Sample the pad's full surface.
2. **Normal smoothing.** `EarField` returns a nearest-sample face normal, which
   rounds off sharp lips and understates the inward-facing component.

So the stability margins reported for the cymba leg are conservative. The
actionable move is not a model change but a geometry one: shift the cymba pad a
few millimetres to land **on** the overhang that the scans already show, then
re-measure — the interlock credit will follow automatically.

### Auto-export

Every single-ear run rewrites `viz/seated_scene.glb` and its metadata for that
ear, so each iteration can be eyeballed. `--viz-ear dataset/ear_id` forces a
specific ear on a multi-ear run.

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

> **Caveat added in v5.** Everything in this section scores whether the rim forms
> a continuous contact loop — it does **not** check that the loop surrounds the
> canal. On the v4 seatings that was false on more than half the audited ears, so
> **the v4 seal counts below are inflated.** See
> [the enclosure audit](#the-enclosure-audit). They are kept as the record of what
> the v4 run produced, not as a statement of the design's seal performance.

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

## Priorities after v5

1. **Re-run the full 107 with the enclosure constraint.** Every population figure
   in this report predates it, and the seal numbers in particular are not
   trustworthy without it. This is the gate on freezing for variant B.
2. **Give back 1–1.5 mm of wing.** The 2.4 mm shortening overshot: over-press
   fell but the wing now misses entirely on 5–7 of 13 ears. Also revisit the −10°
   splay — 3 of 13 are now *misdirected* (aim ≥ 60° off the surface normal),
   which was 0 before.
3. **Keep the compliance-only notch sector.** It captures the full benefit of a
   uniform 4.0 mm budget while relaxing only the inferior 90°, and it does not
   carry the geometric harm the v4 flare did.
4. **Expect a protrusion/seal trade and decide it deliberately.** Holding the
   bore over the canal costs 0.7–2.2 mm of median protrusion. That trade was
   previously hidden by an unconstrained seal metric; it now has to be made on
   purpose.
5. **Watch the jacket.** Mean clearance has risen monotonically across every
   revision (3.01 → 3.94 → 4.27 → 4.63 mm).

## Pipeline changes in v5

- **`seal_enclosure_audit.py`** (new): the aim / offset / enclosure checks above.
- **`c_aim` in `seating_cost`**: the bore-over-canal constraint. 3 mm of lateral
  slack free, quadratic beyond, weight 0.40. `--no-aim` disables it.
- **`--json-dir`** on `align_ear.py --reseat`, `tryon.py` and
  `seal_compliance.py`, so a subset of ears can be reseated and scored in place
  without touching the other 94 — which is what made a 13-ear run practical.

Carried forward: as-built rim sampling, the nozzle-frame fix, `c_face`, `c_soft`,
`c_prot`, clearance graded on `shell`, seeded sampling, Powell over 4 starts.

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

- **v5 is a 13-ear short-list run.** Its P/M/F counts are not population
  estimates; see the scope box at the top. The full-107 figures elsewhere in this
  report are from v4 and predate the enclosure constraint.
- **Seal figures from v1–v4 are inflated** — they did not require the sealed loop
  to contain the canal entrance.
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
.venv/bin/python seal_enclosure_audit.py               # aim / offset / enclosure
.venv/bin/python viz_scene.py                          # viz/seated_scene.glb

# short-list workflow (reseat + score a subset in place)
.venv/bin/python align_ear.py --reseat --json-dir DIR --jobs 7   # add --no-aim to
.venv/bin/python tryon.py            --json-dir DIR              # drop the bore-
.venv/bin/python seal_compliance.py  --json-dir DIR              # over-canal term
```

Run reseats **one at a time** — concurrent reseats interleave their poses in
`ears/aligned/*.json`. `--json-dir` on `seal_compliance.py`,
`retention_analysis.py` and `shrink_estimate.py` reads a saved copy of the
seatings if a reseat may be running. `--field-seed` varies the sampling for a
stability check; `--cant` / `--stl-dir` score a non-default build; `--qc-png`
writes a per-ear depth map with the detected aperture and tragus marked;
`--landmarks` overrides a bad pick by hand.
