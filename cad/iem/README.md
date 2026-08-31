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

* Right-ear master, L/R mirror variants only; **one size**, no S/M/L.
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

## 6. The wing: a macro-scale gyroid shell

`docs/MECH_VALIDATION.md` (FEA, commit `658bf9d`) refuted the original wing: a
1.2 mm-cell graded gyroid inside a 4 × 7 mm blade is **43–85 % dense** at printable
walls and measured **k = 5 134 N/mm** against a target of ~0.3 N/mm — a rigid
wedge, not a spring. A 1.2 mm cell simply cannot be open at a 0.20 mm minimum
wall (ρ ≈ 3.09 t/a ≥ 0.51 by construction).

**The fix is scale, not density.** The wing is now **1–2 gyroid unit cells at a
12 mm cell** — at that size the gyroid is no longer a lattice, it is a single
doubly-curved 0.2 mm Ti sheet that flexes like a leaf spring. The jacket's own
fine 1.2 mm gyroid skin is unchanged; it is structural-rigid by design.

| | value |
|---|---|
| unit cell (`gyroid_cell_wing`) | **12.0 mm** |
| sheet wall (`wing_wall_root` → `wing_wall_tip`) | **0.22 → 0.20 mm** |
| envelope (free span × press direction × depth into concha) | **13.2 × 7.0 × 5.0 mm** |
| anchored foot (`wing_anchor_w` → `wing_width` over `wing_anchor_len`) | 2.4 → 5.0 mm over 7.0 mm |
| relative density, nominal 3.09·t/a | **5.4 %** |
| relative density, **measured from the SDF** incl. rolled edges | **5.2 %** |
| solid Ti transition into the jacket rim (`wing_root_solid`) | 1.2 mm |
| rolled rim on every exposed sheet edge (`wing_edge_wall`) | 0.40 mm over a 0.70 mm band |

### Stiffness — shell-bending estimate

The wing is a shell, not a continuum, so a homogenised-modulus beam model is the
wrong tool (it gives 20–120 N/mm and is meaningless at 1–2 cells: periodic
homogenisation suppresses exactly the global inextensional modes that make a
single sheet compliant). The model used instead treats the sheet as a set of
plate strips:

> **k = 1 / ∫₀^L (L−s)² / EI(s) ds**,  with  **EI(s) = D · ℓ(s) · χ**
> **D = E t(s)³ / 12(1−ν²)** — plate rigidity of the Ti sheet
> **ℓ(s) = L_A · A_cut(s)** — sheet chord length in the cut, from the stereological
> identity L_A = (π/4)·S_V with S_V = 3.09/a for a sheet gyroid
> **χ = `shell_chi` = 0.40** — mean ⟨cos²θ⟩ orientation factor, because the sheet
> meets any cut plane at ~45–55°, so only a fraction of each strip resists
> press-direction bending

E = 110 GPa, ν = 0.31 (Ti-6Al-4V). Integrated per station over the free span:

| | value |
|---|---|
| **tip stiffness k** | **0.251 N/mm**  (target 0.15–0.35) |
| **F at 1.0 mm** | **0.251 N**  (target 0.2–0.4 N) |
| **F at 1.5 mm** | **0.376 N**  (target 0.2–0.4 N) |
| sheet chord at the foot / at the tip | 3.4 / 7.1 mm |

Sensitivity: χ is the weak link. χ = 0.30 → k = 0.19 N/mm; χ = 0.50 → k = 0.31.
The whole band still lands inside the target, which is why these parameters were
chosen, but **this is an estimate, not a measurement — the sibling FEA agent has
to confirm it.** The generator prints the number on every run so a parameter
change cannot silently drift out of band.

The three tuning levers, in order of authority:

1. `wing_anchor_w` — how much of the sheet is anchored at the foot. Necking the
   foot from 5.0 to 2.4 mm over the first 7 mm is what took k from 0.32 to 0.25.
2. `gyroid_cell_wing` — bigger cell, less sheet in any cut, softer. 12 mm is the
   top of the briefed 8–12 mm range and is already used.
3. `wing_wall_*` — k ∝ t³, so this is the strongest lever, but 0.20 mm is the
   process floor and there is nowhere left to go.

### Printability, rim-down

| surface | worst | p99 (area-weighted) | area over 45° |
|---|---|---|---|
| **as-built wing sheet (y > rim)** | 89.8° | **82.9°** | **13.8 %** |
| whole jacket + wing | 90.0° | 83.9° | 17.3 % |
| wing *solid envelope* (a bounding shape, not the part) | 89.8° | 87.3° | 35.9 % |
| best of 42 sampled build directions, on the envelope | 76.4° | 45.5° | 0.5 % |

**This is the honest cost of going macro-scale.** A fine gyroid is self-supporting
because every cell wall turns over within ~1 mm; a 12 mm-cell sheet has long,
shallow runs and its saddle regions are locally horizontal. 13.8 % of the wing's
area needs support rim-down. Either accept supports on the wing (it is a
throwaway surface, the ear side is the jacket skin), or tilt the plate — the
scan's best direction drops the envelope's over-45° area to 0.5 %.

Minimum wall is enforced: the generator warns if `wall_face` < `min_wall` or
`gyroid_cell` < `min_cell`. Every exposed sheet edge is rolled to 0.40 mm, so
there are no knife edges at the tip or along the deep edge, and the wing
cross-section is a rounded stadium (`wing_edge_round`) rather than a box.

The core's −Z dome still has a genuine horizontal pole and needs supports or a
tilted plate; that is normal for an ear shell and is not addressed here.

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
