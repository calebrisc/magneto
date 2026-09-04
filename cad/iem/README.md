# Magneto IEM — parametric STL generator (v0)

`generate.py` builds every part of the Magneto in-ear monitor from one `PARAMS`
dict. Geometry is expressed as signed-distance fields, sampled on a regular grid
and polygonised with marching cubes. That is slower than a B-rep kernel but it
is **watertight by construction** and it lets the graded gyroid lattice and the
solid Ti bodies live in the same expression tree.

```bash
python3 -m venv .venv
.venv/bin/pip install numpy scikit-image trimesh rtree   # or -r requirements.txt

.venv/bin/python generate.py --all                       # ~18 s, both ears
.venv/bin/python generate.py --part core --part faceplate
.venv/bin/python generate.py --all --magnet-preset n52_long
.venv/bin/python generate.py --part jacket_wing --voxel 0.07   # print-quality lattice
.venv/bin/python generate.py --list
```

STLs land in `stl/right/` and `stl/left/`. **`stl/` is gitignored** — the full
set is ~488 MB (the lattice alone is 1.0 M triangles) and regenerates in 18 s,
so it is a build product, not a source artifact.

---

## 1. Coordinate convention

| Axis | Direction |
|---|---|
| origin | centre of the **nozzle base** — where the nozzle stub leaves the core |
| **+X** | nozzle axis, pointing **into** the ear |
| **+Z** | toward the **faceplate**, outward from the head |
| **+Y** | up, toward the antihelix / cymba concha |

Units are mm throughout. The **right** ear is the master; the left is produced by
mirroring across the sagittal plane, implemented as `scale(-1, 1, 1)` on X with
face reversal.

The **nozzle axis is canted** relative to the body: `nozzle_cant_deg` (default
**45°**) rotates it about +Y, so +X → −Z. See §2b.

---

---

## 2. Design decisions carried in from the brief

* Right-ear master, L/R mirror variants only. **One body.** S/M/L exists on two
  parts only: the bell-tip lip (9/1) and the clamp pad extension (9/4, sized by
  fit, not by ear size — see `docs/CLAMP_SIZE_BANDS.md`).
* The seal happens at the **canal opening**. Nothing enters the canal — the
  skirt rim is the most-distal feature of the whole assembly (see §5).
* Everything is **magnetically removable**: jacket→core (3 magnets + 2 pins),
  faceplate→core (2 magnets), tip carrier→nozzle (mag-float pair + bayonet).
* Minimum printable wall 0.20 mm, minimum gyroid cell 1.00 mm — enforced, with
  a printed warning when a parameter set violates either.

---

## 2b. The nozzle cant

`docs/TRYON_REPORT.md` (107 ears, commit `fe609b7`) failed **104 / 107 on
protrusion**. The mechanism was arithmetic, and the report nailed it to 0.02 mm:

> the rigid body behind the seal plane ran x = +4.65 → −18.75, i.e. **23.4 mm
> along the nozzle axis**; the seated nozzle sits a median 55° off the concha
> normal, so 23.4 × cos 55° = **13.6 mm** of faceplate protrusion — against a
> measured median of 13.58 mm.

The cause was the design frame's own assumption, already flagged in v0: it
treated the nozzle axis and the faceplate normal as **orthogonal**, when in a real
ear the canal axis and the concha-floor normal are **30–60° apart**. The body
therefore stacked straight down the canal instead of lying in the concha bowl.

**Fix: cant the nozzle, keep the body.** `nozzle_cant_deg` rotates the nozzle axis
about the design +Y, taking +X toward −Z:

```
axis  â = (cos c, 0, −sin c)     c = nozzle_cant_deg, default 45°
```

−Z is chosen (not +Z) because it keeps the jacket's ear-facing −Z surface pointed
into the concha floor: −Ẑ · â = +sin c.

**What follows the canted axis:** the core's nozzle stub and its bayonet lugs, the
nose cone, the sound bore, the nozzle inserts, the carrier, the skirt and the seal
plane. The magnet stack and the bayonet are unchanged — they are all measured
along that axis, so their arithmetic is identical.

**What does not move:** the core ellipsoid, the faceplate, the jacket and the
wing keep their orientation to the concha exactly as before.

The axis passes through the **core centre**, so the nozzle root is buried the full
ellipsoid support distance (6.09 mm at 45°) and the base lands where the axis
exits the shell:

| cant | nozzle base (x, z) | root buried | rigid stack along the axis |
|---|---|---|---|
| 0° | (0.00, 0.00) | 8.50 mm | **23.40 mm** — reproduces the report exactly |
| 30° | (−2.43, −3.50) | 7.00 mm | 21.1 mm |
| **45°** (default) | **(−4.19, −4.31)** | **6.09 mm** | **17.96 mm** |
| 60° | (−5.77, −4.73) | 5.47 mm | ~16 mm |

`--cant 0` reproduces the old collinear geometry bit-for-bit, and the generator
measures the stack off the actual meshes on every run — at cant 0 it prints
`−18.75 .. 4.65 = 23.40 mm`, which is the report's number to the centimetre of a
millimetre. That agreement is the cross-check that the two models describe the
same part.

**The axial number understates the benefit.** Canting removes 5.4 mm from the
stack measured *along the nozzle axis*, but the protrusion that failed 104 ears is
measured along the *concha normal*, and the cant's real work is letting the
optimiser lay the 17 × 14 mm body flat in the bowl instead of down the canal. The
report's estimate for that is 23.4 × (cos 55° − cos 85°) ≈ 11.4 mm. **Only a
try-on re-run can settle it** — the generator has no ear.

### The sound bore now curves

A canted nozzle cannot share a straight bore with a driver that stays in the body
plane, so the bore is built as a **swept sphere along a polyline**: straight along
â out through the stub, one circular arc of radius `bore_arc_r` (9 mm) turning
through the cant angle, then straight into the front volume under the driver.
Diameter is a constant Ø`nozzle_bore` by construction, the turn is exactly the
cant (≤ 45° at the default), and the arc bulges *away* from the shell's lower
surface rather than toward it. The driver pocket and the driver itself do not
move — they stay in the body plane, firing −Z into the front volume, and the arc
picks the sound up directly beneath them.

### Knock-on changes

* **Jacket magnet and pin stations moved.** The old station at (cx + 5.6, 0) sat
  1.4 mm from the canted nozzle root — inside the nose boss. They are now
  magnets at (cx − 6.0, 0) and (cx + 0.5, ±4.8), pins at (cx − 3.0, ±4.0), all
  ≥ 6 mm from the nozzle exit in XY.
* **The jacket now offsets from the ellipsoid alone** (`core_body`) rather than
  from ellipsoid ∪ nose, with an explicit clearance hole cut where the canted
  nozzle passes through the shell. Otherwise the jacket wrapped the nozzle root.
* **The faceplate got smaller** (331 → 304 mm³): the canted nose no longer
  reaches the +Z cap, so the cap is now a clean ellipsoid section.
* **The front vent is anchored to the bore path** rather than to a fixed x, so it
  always opens into the front volume whatever the cant.
* **Nozzle inserts, carrier and moulds keep their own axis-aligned frame.** Their
  STLs are unchanged by the cant — they are bodies of revolution about â — and
  the assembly applies the cant transform when it places them.

### Still open

The report also asked for body shortening on top of the cant ("the remainder has
to come from body length"). That is **not** done here: the core is still
17 × 14 × 10 mm, set by the 12 mm driver carrier. Doing it means a smaller driver
or a two-piece shell.

---

## 2c. Protrusion work: corner roll, body trim, notch sector

Three changes answering `docs/TRYON_REPORT.md` v2 + the recalibration (`13e75b3`)
and the seal rescore (`5ff456f`). All three are parametric and all three are
measured and printed on every run.

### Corner roll — `corner_chamfer` (4.0 mm requested)

The worst-protruding point is a **corner**, median (−11.5, −4.1, +3.7): the −X/−Y/+Z
octant of the faceplate, diagonal to all axes, with protrusion sensitivity
0.65 / 0.49 / 0.41 mm per mm on +Y / +Z / +X. The roll cuts perpendicular to that
diagonal (`corner_dir`), which is the most material-efficient direction available,
with a `corner_roll` = 1.5 mm blend. Depth ramps with Z: below the parting plane it
is limited to what the closed internals allow, above it the full 4 mm.

**The driver pocket is the wall.** Ø12.2 needs 6.1 + 0.9 mm to the shell in every
XY direction, and `core_ry` is exactly 7.0 — zero slack. The generator computes the
allowance and reports it: **0.12 mm below the parting plane**. So the 4 mm roll
lives entirely in the faceplate, and the worst point simply relocates to the core
rim at z = 1.0. The faceplate magnet stations were moved to a diagonal pair,
(cx ± 6.5, ∓2.9), because the old axial −X station sat inside the roll.

### Body trim — `body_trim_mm` (**default 0.00 — protrusion ACCEPTED**)

**Settled 2026-08-31: protrusion is accepted as-is.** No 8 mm driver, no further
body shrink. `body_trim_mm` defaults to 0 and the core stays 17.00 × 14.00 ×
10.00 mm with a 550.8 mm³ acoustic void. The solver below is kept because the
finding it produced is the reason the decision went that way.

Asked for 5 mm of protrusion reduction, distributed across X/Y/Z by
`body_trim_w` and clamped to what the internals allow (X and Y analytically; Z by
bisecting `solve_body_trim()` against a *measured* 0.5 cc void), the answer was:

| axis | requested | feasible | what binds |
|---|---|---|---|
| X | 2.47 mm | **0.00** | faceplate magnet rim, already at 0.02 mm |
| Y | 3.91 mm | **0.00** | driver pocket wall is exactly `body_trim_min_wall` |
| Z | 2.95 mm | 1.34 | acoustic void hits the 0.5 cc floor |

**0.65 mm of the 5.00 mm asked for**, and combined with the corner roll it moved
the worst point 10.47 → 10.20 mm — 0.27 mm, measured on the meshes. The report's
worked example (6 X + 4 Y + 4 Z) would take the shell to 11 × 10 × 6 mm around a
14 mm-minimum driver pocket. 4–6 mm was never reachable without a driver change,
which is what made ACCEPT the right call.

The margin report still runs every build, so a future parameter change cannot
quietly eat a wall. Two rims stay flagged and are **pre-existing**, not caused by
any trim: a Ø2 magnet pocket needs 2 mm of rim and the annulus between the Ø12.2
driver pocket and the shell is 1.17–2.23 mm. That mounting wants a rethink on its
own schedule.

### Intertragic-notch sector — compliance, not reach

v3 gave the inferior sector 3.25 mm of extra radial flare. **v4 (`2d2c491`)
showed that is net negative: it breaks more seals than it fixes.** A radial
extension commits the lip *toward the notch opening*, where there is no flesh to
seal against — the v3 number was a deformation budget, not a reach.

**Reverted: the rim is a plain Ø19.0 mm circle all the way round** and the
4.5 mm land is unchanged. `notch_sector_ext` is 0 and should stay there; the code
path is kept only so the experiment is reproducible.

The sector now gets **compliance** instead, so the lip can drape onto the
cartilage flanking the notch rather than reaching across it (`notch_compliance`):

| | sector | elsewhere |
|---|---|---|
| land wall | **0.22 mm** | 0.25 mm |
| hinge wall | **0.15 mm** | 0.20 mm |
| hinge free length | **1.60 mm** | 0.60 mm |
| rim | Ø19.0 mm | Ø19.0 mm |
| land width | 4.50 mm | 4.50 mm |

Sector is `notch_sector_deg` 90° centred `notch_sector_center_deg` 180° from
+Y_local (inferior), with 20° smoothstep blends. Try-on v3 measured the aim 18°
off, inside that blend.

Both numbers are **measured, not asserted**. `measure_skirt_wall()` marches
inward along the cone normal at mid-land and returns the interval where the
built field is negative; every `--all` run prints
`MEASURED land wall: 0.221 mm in the sector vs 0.250 mm outside`, and
`measure_notch_reach()` confirms the sector reach is `+0.00 mm`. Mould draw is
still re-checked each run: `non-monotone rim steps 0 / 0 -> demoulds`.

> Lesson worth keeping: v3's realised reach was 1.00 mm against a 1.75 mm spec
> because of a trig error (a radial offset shifts a cone's signed distance by
> *d·cos 35°*, not *d·sin 35°*). That was fixed and verified — and then the
> feature turned out to be the wrong idea. Measuring the number did not save the
> design decision; it only made the decision legible.

---

## 3. Parts

| STL | Material / process | What it is |
|---|---|---|
| `core.stl` | Ti, LPBF | sealed shell, driver seat, nozzle stub, vents, socket + bone-sensor pockets |
| `faceplate.stl` | Ti, LPBF | +Z cap, held by 2 magnets, hollowed to add rear volume |
| `jacket_wing.stl` | Ti, LPBF (or SLS nylon for mules) | graded-gyroid ear-facing skin + antihelix wing, one piece |
| `nozzle_insert_short/med/long.stl` | Ti, LPBF or machined | tube that bayonets onto the core stub; carries the **fixed** ring magnet |
| `carrier.stl` | silicone, cast | mag-float tip carrier + sealing skirt, single material for prototypes |
| `carrier_mold_a.stl` / `_b.stl` | resin | two-part mould, alignment pins, pour spout, vent |
| `carrier_mold_core.stl` | resin | removable core rod: defines the bore, the L-slots and the magnet seat |
| `driver_carrier.stl` | resin | press-fit ring for the 10 mm dynamic driver, with a rear vent notch |
| `damper_jig.stl` | resin | punch die + magazine + plunger for Ø4 mm damper discs |
| `assembly.stl` | — | everything positioned, magnets omitted, visual check only |

---

## 4. Magnets

Numbers come from `docs/MAGFLOAT_MAGNETS.md`, **"Asymmetric pair (as-built)"**
(magpylib 5.0.1 + magpylib-force, Maxwell-stress solver). The moving ring has to
slide over the Ø5 nozzle tube, so its ID is pinned at 5.4 mm and the pair can
never be identical — every preset below is asymmetric. Select with
`--magnet-preset`.

| preset | fixed ring | moving ring | rest gap | F over 1.5 mm travel | ratio | moving mass |
|---|---|---|---|---|---|---|
| **`asym_as_built`** (default) | 7 × 3 × 1.5 | **9 × 5.4 × 2.0** | **2.25 mm** | 0.307 → 0.214 → 0.155 N | **1.98** | 0.489 g |
| `asym_light` | 7 × 3 × 1.5 | 9 × 5.4 × 1.5 | 1.25 mm | 0.405 → 0.288 → 0.199 N | 2.04 | 0.366 g |
| `asym_8mm` | 7 × 3 × 1.5 | 8 × 5.4 × 2.0 | 1.89 mm | 0.382 → 0.230 → 0.151 N | 2.53 | 0.328 g |
| `n52_long` / `n52_clean_bore` | (symmetric, first study) | — | 6.4 / 6.8 mm | — | 1.69 / 1.61 | — |

All bonded NdFeB, Br 0.65 T. `asym_light` is the compact alternative: nearly the
same regulation (2.04 vs 1.98) at half the rest gap and 25 % less moving-ring
mass — pick it if envelope or carrier inertia matters more than the last bit of
flatness. It changes only the carrier's counterbore depth and the bayonet
station, both of which are derived, so `--magnet-preset asym_light` regenerates a
consistent insert + carrier + mould set on its own.

The two symmetric `n52_*` presets are kept for reference only: their moving ring
has ID 3 mm, which cannot slide on a 5 mm tube, and the generator warns.

**Both magnets live in the nozzle insert / carrier. Nothing magnetic is recessed
into the core.** The fixed ring sits in a Ø7.1 counterbore at the +X face of the
insert's flange. The moving ring is **encapsulated in the carrier**: 0.5 mm of
silicone (`magnet_encap`) on each axial face, 0.70 mm on the OD
(`carrier_wall_od`, checked at runtime), and its ID flush with the sliding bore —
a 9 × 5.4 ring over a Ø5.2 bore leaves 0.1 mm, which will not cast, so the
mould core's post locates the ring during the pour instead. The carrier body was
widened to **Ø10.5 mm** to carry the 9 mm ring; that is also the skirt root
diameter, so the skirt's 35° flare now runs 5.82 mm instead of 7.61 mm to reach
the same Ø19 rim.

A counterbore in the carrier's core end swallows the air gap: the carrier's face
sits at x = 4.65 mm, its counterbore floor at x = 6.65 mm, and the moving ring's
−X face at x = 7.15 mm — 2.25 mm in front of the fixed ring at x = 4.90 mm.

One geometric conflict remains, flagged at runtime: the **fixed ring's ID is
3 mm, below the 4 mm nozzle bore**, so the acoustic path necks **Ø4 → Ø3 → Ø4
over 1.5 mm**. Acceptable (short and smooth), and unavoidable at this ring size.

Stray field from the fixed ring is 3.30 mT at 10 mm behind its back face, i.e.
roughly at x = −5.9 mm — inside the front volume, ~5 mm from the driver's magnet.
Worth a bench check on driver THD before committing.

## 5. The protrusion budget (why the numbers came out where they did)

The seal-at-the-opening rule means nothing may sit distal to the skirt rim. The
skirt therefore flares **forward** (a trumpet, not a normal backward-flaring
eartip), and the rim lands in the same plane as the carrier's distal face:

```
x = 0.00   core face / nozzle base
x = 0..3   core nozzle stub, Ø5, two external bayonet lugs
x = 3.20   insert socket ends
x = 4.90   insert magnet flange face  <- FIXED ring's +X face
x = 4.65   carrier core face (its Ø10.5 collar overlaps the flange by 0.25)
x = 6.65   carrier counterbore floor  (0.5 mm encapsulation web)
x = 7.15   MOVING ring -X face        <- 2.25 mm gap
x = 7.83   skirt root, r = 5.25
x = 9.20   moving ring +X face
x = 9.55   bayonet lug station (derived: ring end + 0.35)
x = 12.90  far end of the L-slot float pocket (1.5 mm travel)
x = 13.65  carrier distal face == skirt rim, Ø19.0     <- seal plane
x = 15.15  same, at full 1.5 mm float
```

**Total protrusion from the core face: 13.65 mm at rest, 15.15 mm at full float**,
and the seal plane is the frontmost feature — nothing enters the canal. Against
the `EAR_ANTHROPOMETRY.md` envelope (concha depth 8–18 mm) that fits the upper
~60 % of the range at rest. `asym_light` is 1.0 mm shorter in the counterbore but the same overall, since
`carrier_len` dominates. The `n52_*` presets add ~4 mm and are not recommended.

The `short` / `med` / `long` inserts differ only in tube length past the magnet
flange (**7 / 9 / 11 mm**, raised to reach the new lug station). The lug station
is derived from the magnet stack, so it is the same for all three inserts within
a preset and **one carrier fits all three**; the longer tubes exist to support the carrier
further out and to give deep conchas more retention, and they do protrude past
the seal plane, so `short` is the default in the assembly.

---

### The skirt's contact land

`MECH_VALIDATION.md` JOB 2 found the force is fine (0.155–0.307 N) but the
geometry never committed to spreading it: a plain 35° cone with a uniform
0.35 mm wall lands on a curved aperture as a **line** contact, giving
17.9–34.1 kPa — 4–8× the 4.27 kPa capillary-closing flag. The required minimum
band was **≥ 4.0 mm** of slant width.

The skirt is now built as an **exact 35° outer cone with the wall cut from the
inside**, so the outer face is a true conical band (Ø19.0 rim unchanged, and the
rolled rim lip is a torus tangent to that cone so the diameter is exact):

| slant station | wall | what it is |
|---|---|---|
| 0 → 1.91 mm | `skirt_wall_neck` 0.40 mm | structural neck into the carrier body |
| 1.91 → 2.71 mm | `skirt_wall_hinge` 0.20 mm | compliance groove — lets the land rock and bed flat instead of digging in at its leading edge |
| **2.71 → 7.41 mm** | `skirt_wall_land` 0.25 mm | **the contact land, 4.50 mm of slant width, Ø13.8 → Ø19.0 mm** |

Cone-normal pressure at F_max = 0.307 N over that land, printed on every run:

| contact Ø | 10 mm | 13 mm | 16 mm | 19 mm |
|---|---|---|---|---|
| **kPa** | **3.79** | **2.91** | **2.37** | **1.99** |
| verdict | borderline | borderline | borderline | comfortable |

Every diameter now clears the 4.27 kPa ischaemia flag (against 17.9–34.1 kPa for
line contact). The shear-derated 2.15 kPa stretch target is met only at Ø19; that
is the FEA's own conclusion — at Ø10 the 19 mm / 35° funnel cannot reach it
without a lower preload or a wider flare. Note also FEA §2.8: real pressure peaks
at the leading edge of contact at plausibly 1.5–3× the average, which is what the
0.20 mm compliance groove is there to blunt.

## 6. The wing mechanism: three radial mag-plungers

`docs/MECH_VALIDATION.md` (`8eb6ac1`) retired the Ti spring wing. The compliant
gyroid sheet passed the force target only by being ~14 mm long, which reinstated
exactly the overpressure the v4 shortening was added to fix — it traded a force
problem for a fit problem. The replacement decouples them: a **cam preset absorbs
per-ear geometry**, and the magnets only regulate force over a small window they
can comfortably cover. It also removes fatigue from the design entirely.

The plan then moved from one 14 × 6 mm rail to **three independent radial
plungers**, one per contact site, each with its own aim, cam and stops.

`wing_style="plungers"` (default) builds them; `wing_style="gyroid"` keeps the old
sheet as a legacy option, unchanged.

### Sites (`plunger_aims`, parametric)

| site | aim (x, y, z) | cam | boss | pad | reach | clearance to the nozzle stack |
|---|---|---|---|---|---|---|
| cymba | (+0.20, +0.98, −0.10) | 4 × 3.0 mm | 5.50 mm | Ø10.5 | 20.21 mm | +0.62 mm |
| antihelix undercut | (−0.45, +0.85, −0.28) | 4 × 3.0 mm | 5.50 mm | Ø9.0 | 20.42 mm | +3.43 mm |
| **tragus_inner** *(disabled)* | (+0.82, −0.17, −0.54) | **6 × 4.5 mm** | 6.90 mm | Ø9.0 | 22.32 mm | **−6.53 mm — CLASH** |

Each aim is normalised at load and the boss base is found by ray-casting the aim
to the point where the core's outward normal *is* the aim, then offsetting by
`clearance + jacket_thick`. A new aim places its own boss — no hand-set
coordinates. Clearance to the nozzle/insert/carrier/skirt stack is computed every
run from the real stack profile (`nozzle_stack_profile()`), with the plunger
modelled as two capsule sections because the boss is narrower than the pad.

### Per-leg enable — the shipped build is TWO legs

`plunger_enable` turns each leg on or off. **tragus_inner ships OFF**, so the
default build is the two-leg + interlock variant, ready for a stability test.
All three are still defined, reported and clearance-checked every run; the
disabled one just is not built.

> ### tragus_inner cannot be reached from this shell — and no base fixes it
>
> The extended cam **does** fit its local depth budget: 6 detents over 4.5 mm
> needs `cam_h` 4.70 mm and a 6.90 mm boss, which the leg carries fine.
>
> The aim is the problem. **(+0.82, −0.17, −0.54) sits 15.2° off the canted nozzle
> axis**, so the leg drives down the same corridor the nozzle already occupies:
> **−6.53 mm of interference**, not a near miss.
>
> Since the base is a free variable, `leg3_feasibility()` sweeps bases over the
> anterior / anterior-superior jacket surface, re-aims each at the same tragus-wall
> target, and scores clearance against the real stack profile. Result, printed
> every run: **0 of 1431 candidate bases clear 0.80 mm; the best is −6.04 mm.**
> Relaxing the filters does not rescue it either — at 20 000 directions with *no*
> aim or length limit at all, i.e. every point on the jacket, the best is still
> **−2.27 mm**.
>
> **The root cause, which the sweep also prints:** the tragus-wall target sits
> 8.01 mm off the nozzle axis at station 14.74 mm — 1.1 mm past the skirt rim
> plane and *inside* the Ø19 rim's 9.50 mm radius. **The tragus inner wall is in
> the skirt's shadow.** Any line from the shell to it has to pass through the
> skirt. This is not a boss-placement problem and cannot be solved by moving the
> leg.
>
> The levers that would actually work, none chosen unilaterally: shrink
> `skirt_max_dia`, reduce `nozzle_cant_deg` so the stack swings off that corridor,
> or accept the tragus inner wall as unreachable and let the two-leg variant carry
> retention.

### The cymba lip bias

`cymba_lip_bias` = **7°** rotates the cymba aim toward +Y, in the plane it shares
with +Y, so the pad lands *under* the cymba's overhanging lip rather than on it.
That leg also gets `cymba_pad_extra` = **1.5 mm** of extra pad diameter (Ø9.0 →
Ø10.5) on a **0.40 mm rolled shoulder** (`plunger_pad_roll`), so the extended
shoulder can tuck under the lip and interlock instead of digging in. It is a
separate STL, `plunger_pad_cymba`; the assembly picks it for that site only.

The bias narrowed the cymba boss's clearance to the nozzle insert's socket flange
to **+0.62 mm** — legal, but the tightest joint on the shell and worth watching if
`plunger_boss_od` grows.

### The plunger stack

Axial stations, measured along the aim from the boss mount face at s = 0:

```
s = -4.35 .. -1.15   cam preset ring (3.20 mm, 4 detents over 3.0 mm)
s = -1.00 ..  0.00   FIXED 5x2.5x1 N35 ring, seated on the selected cam step
s =  0.00 ..  2.75   air gap  <-- rest gap
s =  2.75 ..  3.75   MOVING 5x2.5x1 N35 ring, in the foot
s =  3.75 ..  4.55   Ti plate over the moving ring (0.80 mm)
s =  4.55 ..  5.90   silicone pad, 1.00 mm + 0.35 mm rocker crown
```

| | |
|---|---|
| **depth stack** (fixed ring back → moving ring face) | **4.75 mm** |
| **dynamic travel** | **±0.75 mm**, hard stops at s = 2.00 and 3.50 |
| force per pair (from MECH_VALIDATION §5.2) | 0.493 N at −0.75 → 0.284 N at rest → 0.180 N at +0.75, ratio 2.73 |
| cam preset | 4 detents over 3.0 mm; **6 over 4.5 mm** on `plunger_cam_ext_sites` |
| guide pin | Ø2.0 through both ring IDs, 0.30 mm polymer sleeve in the jacket bore |

The inward stop is a 0.8 mm-wall skirt on the foot that bottoms on the boss face;
`--all` measures it on the built mesh and prints
`stop skirt reaches s = 2.00 mm -> 0.75 mm of inward travel`. The outward stop is
the pin head catching in the boss bore. The pin runs through the 2.5 mm ring IDs,
so magnets and guidance are coaxial — that is why the ring geometry was chosen.

### Parts

`plunger_foot`, `plunger_pad`, `plunger_pin`, `plunger_cam` are each **one STL,
printed three times**, plus two per-site variants — `plunger_pad_cymba` (wider
shoulder) and `plunger_cam_ext` (6 × 4.5 mm); they are bodies of revolution so their standalone STLs are
axis-aligned about +X and the assembly transforms three copies onto each site.
The jacket carries the three bosses (Ø9.6, 5.5 mm tall) with the cam bore, the
sleeved pin bore and four detent notches.

> **Mass warning, printed every run:** jacket + 3 bosses is 1206 mm³ = **5.34 g of
> Ti**, and the bosses dominate it. `MECH_VALIDATION` budgeted 971 mg/side for the
> whole mechanism. If mass matters — and on an IEM it does — thin
> `plunger_boss_od` or lattice the bosses. This is the biggest open item on the
> mechanism.

Also open: pad tip reach is ~20.2–20.4 mm from the core centre at rest, which is
generous against a concha depth envelope of 8–18 mm. The 3 mm cam takes up
per-ear variation, but the aims and `plunger_boss_h` want a try-on pass before
anyone cuts metal.

### Did the bosses disturb the skirt rim?

No. The v6 seal rescore flagged P0023's skirt band at 92 % against the 95 % rule,
and the plunger bosses were the obvious suspect. They are not the cause: the
carrier's rim cross-section was sampled at 721 azimuths and compared against the
pre-plunger commit (`dbc3b7b`), and the rim radius is **bit-identical — max,
min and RMS delta all 0.0000 mm**, a flat 9.482 mm all the way round. The bosses
live on `jacket_wing`; `carrier_field` shares no parameter with them. The 92 % is
a seating/pose effect, not geometry.

## 6b. Cable exit boot

`cable_exit` — `"up_back"` (default), `"back"` or `"none"`. A tapered strain-relief
stub grown off the 2-pin socket opening so the try-on contract can score cable
clearance: 6.0 mm long, Ø5.5 → Ø3.5, Ø1.8 cable bore, raked 35° up (+Y) and back
(−X). It is smooth-unioned into the shell and the bore is cut through into the
socket pocket.

## 7. Every parameter

### Process limits

| param | default | why |
|---|---|---|
| `min_wall` | 0.20 mm | thinnest wall LPBF Ti / high-res resin will hold; violating it prints a warning |
| `min_cell` | 1.00 mm | below this a gyroid traps unfused powder; violating it prints a warning |
| `clearance` | 0.15 mm | jacket-to-core gap. Smaller and thermal/print variation binds; larger and the magnets feel loose |
| `press_clearance` | 0.15 mm | slip fit on the insert-over-stub joint. It is a bayonet, not a press, so it wants clearance |

### Magnets

| param | default | why |
|---|---|---|
| `magnet_preset` | `bonded_compact` | the only combination in the magpylib grid that clears all four force bounds at a compact gap |
| `magnet_pocket_clear` | 0.05 mm | glue gap around the ring; too much and the axial position (hence the force) drifts |
| `magnet_encap` | 0.50 mm | silicone over the moving ring's axial faces. Drives the counterbore depth; a preset whose rest gap is under 0.8 mm will warn |

### Core shell

| param | default | why |
|---|---|---|
| `core_rx` | 8.5 mm | half-length. Driver pocket needs 6.1 mm of it; the remaining 2.4 mm is the front volume |
| `core_ry` | 7.0 mm | **raised from the brief's 6.0.** A 12 mm driver-carrier ring needs 6.1 mm of half-width; 0.9 mm of wall on top gives 7.0. A 12 mm carrier simply does not fit a 12 mm-wide shell |
| `core_rz` | 5.0 mm | as briefed. 10 mm total sits inside the 8–18 mm concha-depth envelope |
| `core_wall` | 1.20 mm | LPBF Ti pressure wall with margin for post-polish |
| `faceplate_z` | 1.00 mm | parting plane. Lower and the faceplate gets thick; higher and the rim annulus outside the driver pocket gets too narrow for the Ø2 magnets |
| `cavity_cap_z` | −0.60 mm | caps the acoustic cavity 0.6 mm below the magnet pockets so the rim stays solid |
| `nose_cone_r0/r1/x0` | 5.0 / 3.2 / −6.0 mm | the cone that blends the ellipsoid into the nozzle stub; smooth-unioned with k = 1.2 mm for a fillet |

Resulting shell: **17 × 14 × 10 mm**, versus the brief's ~16 × 12 × 10. Both
increases are forced by the 12 mm driver carrier, not stylistic.

### Driver

| param | default | why |
|---|---|---|
| `driver_dia` | 10.0 mm | the specified dynamic driver |
| `driver_carrier_id` | 10.0 mm | slip fit on the driver can |
| `driver_carrier_od` | 12.0 mm | 1 mm of ring wall each side |
| `driver_carrier_h` | 3.0 mm | as briefed |
| `driver_pocket_clear` | 0.20 mm | → Ø12.2 pocket in the core. The brief's "10.2 mm pocket" is the driver's *acoustic* aperture; the seat has to clear the carrier's 12 mm OD |
| `driver_pocket_depth` | 3.0 mm | carrier sits flush with the faceplate parting plane |

### Nozzle stack

| param | default | why |
|---|---|---|
| `nozzle_cant_deg` | **45.0°** | rotation of the nozzle axis about +Y relative to core/faceplate/jacket, so the body lies in the concha plane. See §2b. `--cant 0` restores the old collinear geometry |
| `bore_arc_r` | 9.00 mm | centreline radius of the curved sound bore. Larger = gentler turn but the arc reaches further back into the shell |
| `bore_arc_end` | 1.20 mm | how far behind the nozzle base the arc ends, i.e. how much straight bore the stub gets |
| `bore_run_in` | 2.50 mm | straight bore from the end of the arc into the front volume |
| `nozzle_bore` | 4.00 mm | per `EAR_ANTHROPOMETRY.md`: 4 mm clears the 4.5 mm small-aperture tail with the skirt around it. Industry-standard 5.5–6.5 mm nozzles are why small-canal users get hurt |
| `stub_od` | 5.00 mm | 0.5 mm of wall on the bore |
| `stub_len` | 3.00 mm | enough bayonet engagement to take the skirt's side load |
| `lug_h` / `lug_w` | 0.60 / 1.50 mm | as briefed; two lugs at 180° |
| `socket_od` | 8.00 mm | insert socket wall = (8 − 5.15)/2 = 1.42 mm, leaving 0.67 mm behind a 0.75 mm-deep L-slot |
| `insert_od` | 5.00 mm | the carrier's sliding surface |
| `insert_tube_lengths` | 7 / 9 / 11 mm | tube beyond the magnet flange — see §5 |
| `damper_dia` / `damper_recess` | 4.00 / 0.30 mm | standard acoustic damper disc, recessed at the ear end |

### Mag-float carrier

| param | default | why |
|---|---|---|
| `carrier_bore` | 5.20 mm | 0.20 mm on the Ø5 tube; in cast silicone this is a working sliding fit, not a press. Also pins the moving ring's ID at 5.4 mm |
| `carrier_od` | **10.50 mm** | raised from 8.0: a 9 mm moving ring plus 0.75 mm of silicone each side. Also the skirt root diameter |
| `carrier_len` | **9.00 mm** | raised from 8.0: the encapsulated 2 mm ring plus lug width plus 1.5 mm travel plus margins does not fit in 8 mm. Checked at runtime |
| `carrier_travel` | 1.50 mm | the travel the force study was run over; enforced by the L-slot pocket length (lug width + travel) |
| `skirt_flare_deg` | 35.0° | as briefed |
| `skirt_wall` | 0.35 mm | legacy nominal; superseded by the three-zone profile below |
| `skirt_land_w` | **4.50 mm** | slant width of the conical contact land. FEA minimum is 4.0 mm at every contact diameter from Ø10 to Ø19 |
| `skirt_wall_neck` | 0.40 mm | wall at the root — carries the preload into the carrier body |
| `skirt_wall_hinge` | 0.20 mm | compliance groove behind the land, so the land rocks to match the ear instead of digging in |
| `skirt_hinge_w` | 0.60 mm | slant width of that groove |
| `skirt_wall_land` | 0.25 mm | wall through the land — thin enough to bed flat over its full width |
| `skirt_max_dia` | 19.0 mm | covers the design envelope's 4.5–14 mm aperture width with room to conform |

The skirt is modelled as **part of the carrier STL** — a single material for the
first prototypes. The production version casts the skirt in **Shore A 10–15**
over a firmer (Shore A 40–50) carrier body, which needs a two-shot mould or an
overmould step; the mould in this repo is the single-shot prototype version.

### Jacket skin (fine gyroid — unchanged, structural-rigid by design)

| param | default | why |
|---|---|---|
| `jacket_thick` | 1.60 mm | lattice + skin. Thick enough for two gyroid cells across |
| `jacket_x_clip` | −2.00 mm | the jacket stops here so it never fouls the nozzle nose |
| `gyroid_cell` | 1.20 mm | above `min_cell`, small enough for two cells through the jacket |
| `wall_face` | 0.20 mm | soft against the ear at the outer/ear face |
| `wall_root` | 0.40 mm | stiff at the rim, where the jacket takes the magnet load |
| `grade_len` | 6.00 mm | distance from the root corner over which the wall grades root → face |
| `skin_t` | 0.60 mm | solid membrane on the ear-facing surface so the lattice never touches skin |
| `perf_dia` / `perf_pitch` | 0.40 / 1.50 mm | sweat perforations through the membrane, on a square grid |
| `solid_root` | 1.00 mm | solid Ti collar before the lattice starts |

FEA measured this region at k = 185 300 N/mm over a 9 mm² patch. That is
intentional — the jacket is the rigid backbone, the wing carries all compliance.

### Wing (macro gyroid shell — see §6)

| param | default | why |
|---|---|---|
| `gyroid_cell_wing` | **12.00 mm** | top of the briefed 8–12 mm range: 1–2 cells across the envelope, so the gyroid is a single doubly-curved sheet rather than a lattice. Bigger cell → less sheet in any cut → softer |
| `wing_wall_root` / `wing_wall_tip` | **0.22 / 0.20 mm** | lightly graded; 0.20 mm is the process floor and k ∝ t³, so there is no room below |
| `wing_edge_wall` | 0.40 mm | rolled/thickened rim on every exposed sheet edge — no knife edges at the tip |
| `wing_edge_band` | 0.70 mm | distance over which the wall ramps up into that rolled edge |
| `wing_root_solid` | 1.20 mm | solid Ti transition into the jacket rim (brief: ≥ 1 mm) |
| `wing_thick` | 7.00 mm | envelope across the press direction (in XY) |
| `wing_width` | 5.00 mm | envelope depth in Z, into the concha |
| `wing_anchor_w` | **2.40 mm** | Z-width at the foot. The primary softening lever: necking 5.0 → 2.4 mm took k from 0.32 to 0.25 N/mm |
| `wing_anchor_len` | 7.00 mm | distance over which the foot opens back out to the full 5 mm |
| `wing_edge_round` | 0.85 | fraction of the half-section used as a corner radius, so the deep edge is a rounded stadium, not a flat ceiling |
| `wing_taper_deg` | 40.0° | overhang of the deep-edge taper, under the 45° self-support limit |
| `wing_taper` | 1.60 mm | depth over which that taper acts |
| `wing_rise` | 11.0 mm | tip lands this far above the core rim; sets the free span (13.2 mm) and hence k ∝ 1/L³ |
| `wing_back_deg` | 30.0° | tip angled toward −X, under the antihelix |
| `wing_root_dx` | 1.00 mm | root offset from the core centre in X |
| `wing_len` | 14.0 mm | nominal centreline length; the Bezier from rise + back angle lands close to this |
| `shell_chi` | 0.40 | ⟨cos²θ⟩ sheet-orientation factor in the stiffness estimate. The single biggest uncertainty — see §6 |

The wing centreline is a quadratic Bezier in XY; the envelope is that curve offset
by `wing_thick/2`, extruded in Z with the necked foot and a rounded deep edge.

### Jacket / core interface

| param | default | why |
|---|---|---|
| `gasket_w` / `gasket_d` | 0.60 / 0.40 mm | groove around the z = 0 parting line for a cast silicone gasket ring (sweat + grit seal) |
| `jmag_dia` / `jmag_depth` | 2.00 / 1.00 mm | three Ø2 × 1 magnets, drilled along the local surface normal on the −Z hemisphere |
| `pin_dia` / `pin_depth` | 1.00 / 1.50 mm | two locating pins. Pins take the shear, magnets only take the normal load |

Magnet stations are at (cx ± 5.6, y = 0) and (cx, y = +4.4); pins at
(cx ± 2.6, y = −4.0). The core gets holes at +0.06 mm, the jacket gets pins at
−0.06 mm.

### Electronics

| param | default | why |
|---|---|---|
| `socket_w/h/d` | 5.60 / 2.80 / 6.00 mm | standard 2-pin socket body, at the −X rim, top |
| `socket_z` | 1.20 mm | keeps the pocket above the driver pocket |
| `bone_w/h/d` | 4.00 / 3.00 / 1.50 mm | bone-conduction sensor, in the −Y (tragus) flank |
| `wire_dia` | 1.00 mm | channel from the bone pocket to the socket pocket |
| `vent_dia` | 0.80 mm | front vent (front volume → −Y/−Z exterior) and rear vent (rear cavity → −Y/−Z exterior) |

The socket pocket is wrapped in a 0.7 mm sleeve that is unioned into the shell
**before** the pocket is subtracted, so it never opens into the acoustic cavity.
The bone-sensor pocket sits in a rounded boss on the −Y flank, smooth-unioned so
there is real wall behind it.

### Mould

| param | default | why |
|---|---|---|
| `mold_wall` | 4.5 mm | block wall around the cavity |
| `mold_pin_dia` / `mold_pin_len` | 3.0 / 3.0 mm | four alignment pins on the +Y half, holes (+0.12 mm) on the −Y half |
| `mold_spout_dia` | 3.0 mm | pour spout into the ring end from +Z |
| `mold_vent_dia` | 1.2 mm | vent at the skirt rim, the last place to fill |

The mould splits on the y = 0 plane. `carrier_mold_core.stl` is a third piece: a
Ø5.2 rod with the L-slot ridges and the magnet-seat boss, plus a Ø8 grip. Silicone
demoulds off the L-slot undercuts without tearing. Assemble core → half A →
half B, pour at the ring end, vent at the rim.

### Meshing

| param | default | why |
|---|---|---|
| `voxel` | `None` | override with `--voxel`; otherwise derived from the budget |
| `budget` | 3.0e6 | voxels for solid parts. Gives ~0.08–0.19 mm depending on bbox — 0.8 mm vents land on 4–10 samples |
| `budget_lattice` | 6.0e6 | voxels for `jacket_wing`. ~0.098 mm, i.e. **2 samples across a 0.20 mm wall** |

That last number is the one real quality compromise. Marching cubes still emits a
watertight surface at 2 samples/wall, but the lattice walls come out rounded and
slightly under-thick. The generator prints a note saying so. **Before sending the
jacket to a printer, re-run `--part jacket_wing --voxel 0.07`** (≈ 3 samples/wall,
~17 M voxels, a couple of minutes, a much larger STL).

---

## 8. Verification

`generate.py` checks `trimesh.is_watertight` on every mesh and prints volume,
bounding box, triangle count, voxel size and build time, then **re-reads every
exported binary STL** and reports boundary edges — because a slicer sees the file,
not the in-memory mesh. `--all` exits non-zero if any exported surface is open.

Two marching-cubes artefacts are cleaned automatically (`drop_specks`): positive
components under 0.02 mm³ (sub-voxel islands that would be loose particles) and
enclosed voids under 1.0 mm³ (trapped-powder pockets at gyroid/skin junctions —
dropping the shell fills the pocket). The jacket sheds ~160 of these, 0.47 mm³
total. Dropping whole closed components keeps the surface closed.

`jacket_wing` re-reads as **0 holes, 2 pinch edges** — two places where the gyroid
sheet touches itself along an edge shared by four triangles. `is_watertight`
demands exactly two faces per edge so it reports False, but the surface is closed
and slicers handle it. Holes are what break a print; pinches do not. That is why
the tool reports both.

One thing worth knowing if you touch `polygonise()`: **do not** call
`merge_vertices()` on the marching-cubes output. skimage already returns a
manifold, index-shared surface; welding across a 0.2 mm lattice wall creates
non-manifold edges and loses watertightness. That bug cost an afternoon.

---

## 9. `fit_check.py`

```bash
.venv/bin/python fit_check.py --ear subject_042.stl --transform "0,0,0,0,-25,10"
```

Loads a SONICOM / HUTUBS-style ear mesh and reports, for the **wing tip** (20
most-distal vertices) and the **skirt rim** (72 points around the sealing lip),
the min/mean/max distance to the nearest ear surface — plus signed distance and
worst interference if the ear mesh is closed. Uses `trimesh.proximity`.

Alignment is **not** automated and is the hard part; `--transform` takes either
`tx,ty,tz,rx,ry,rz` (mm, degrees, XYZ intrinsic) or a path to a 4×4 matrix file.
The docstring describes the three-landmark method for deriving it. Without a
transform the tool runs but says so, loudly.

---

## 10. Known gaps / v1 list

1. **The wing's stiffness is an estimate, not a measurement.** k = 0.251 N/mm
   comes from a shell-bending model whose orientation factor χ = 0.40 is the
   dominant uncertainty (χ = 0.3–0.5 spans k = 0.19–0.31). Sibling FEA to confirm.
2. **13.8 % of the wing needs support rim-down** — the price of a 12 mm cell.
3. **Moving-ring ID is not encapsulated** — 9 × 5.4 over a Ø5.2 bore leaves
   0.1 mm, so the ring's bore face is exposed and located by the mould core post.
   Watch it for delamination on the first cast.
4. **Body shortening** — the try-on report wants it on top of the 45° cant. The
   core is still 17 × 14 × 10 mm, set by the 12 mm driver carrier.
5. **Two-shot carrier mould** — the current mould is single-shot; production
   wants a Shore A 10–15 skirt over a Shore A 40–50 body.
6. **13.65 mm protrusion** is at the upper end of what a shallow concha will
   take. The lever is `carrier_len` and the magnet gap, in that order —
   `asym_light` does not help because `carrier_len` dominates.
7. **Core dome supports** — the −Z pole is a true horizontal overhang.
8. **No front/rear volume tuning** — the acoustic cavity is whatever the shell
   leaves. Once a driver is in hand, tune `cavity_cap_z` and the vent diameters
   against a measured response.
