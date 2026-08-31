# Mechanical validation — Ti gyroid wing, and the magnet contact-pressure budget

Jobs 1 and 2 were run 2026-08-30 against `cad/iem/generate.py` at commit `387b8eb`.
Section 3 re-verifies the **redesigned** wing at commit `8db64d7` and was run 2026-08-31.

1. **Does the graded gyroid wing behave as a 0.2–0.4 N plateau spring over 0.5–2.5 mm of compression?**
   Answer: **no, and not by a small margin.** It is ~1.3 × 10⁴ times too stiff. It is
   not a spring at all; it is a rigid wedge.
2. **Is the mag-float skirt's preload spread over enough contact band to be comfortable
   for all-day wear?** Answer: **only if the contact band is at least ~4 mm wide on the
   funnel.** The current design does not specify a band width, and if the ear lands on the
   cone as a line contact the pressure is 4–8× over the ischaemia flag.
3. **Which wing concept ships?** **The mag-plunger** (§5): 5×2.5×1 mm N35 rings, two
   pairs, 2.75 mm rest gap, 0.49→0.18 N over ±0.75 mm in 4.75 mm of depth, no fatigue
   risk. The Ti sheet was *also* made to pass (§4: k = 0.281 N/mm, 376 MPa at 1.5 mm) but
   only at a 13.99 mm span, which gives back the entire v4 shortening that fixed
   overpressure — it trades a force problem for a fit problem.
4. **Does the redesigned macro-gyroid wing hit 0.15–0.35 N/mm?** Answer: **not yet — it is
   3–7× too stiff at k = 1.06 N/mm**, and the generator's own estimate (0.251 N/mm) is
   optimistic by 4.2× for a structural reason. But the response is a smooth, stable,
   near-linear spring with no snap-through, and it is ~1250–4800× softer than the old block, so
   the architecture is right and only the magnitude is wrong. See §3.

Scripts: `cad/iem/fea/`. Every piece of geometry is imported from `generate.py` —
the same SDF functions that write the STLs — so nothing here can drift from the
printed part. Raw solver output: `fea/rve_results.json`, `fea/wing_results.json`,
`fea/pressure_results.json`.

---

# JOB 1 — Wing and jacket stiffness FEA

## 1.0 Method, and what is honestly *not* in it

**Solver.** [scikit-fem](https://github.com/kinnala/scikit-fem) 11.0.0, linear
elasticity, trilinear hexahedral elements, 2×2×2 Gauss (full integration — no
reduced integration, so no hourglassing). Linear systems solved by smoothed-aggregation
AMG-preconditioned CG (pyamg 5.3.0) with the six rigid-body modes supplied as the
near-nullspace; every solve converged to a relative residual < 1 × 10⁻⁷.

**Material.** Ti-6Al-4V: E = 110 GPa, ν = 0.31, σ_y = 900 MPa. (LPBF as-printed is
usually quoted 950–1100 MPa; 900 MPa is the conservative wrought-annealed figure, as
specified.) Flag level for fatigue relevance: 450 MPa = σ_y/2.

**Meshing.** Structured "voxel" hex meshes: the SDF is sampled on a regular grid and
elements whose centroid is inside are kept, retaining only the largest *face-connected*
component (6-neighbour, so the mesh cannot hang together through a single shared node
or edge, which would be a mechanism). This is the standard micro-FE approach for TPMS
and trabecular bone. It is mildly **stiff-biased** and it produces staircase stress
spikes — both stated where they matter. gmsh tetrahedral meshing is also implemented
(`_common.tet_from_sdf`) but was not needed; P1 tets are markedly worse in bending than
trilinear hexes, which matters for a cantilever.

**Two-scale, and why.** The wing lattice cannot be volume-meshed whole. The generator
emits a 1.2 mm gyroid cell with a 0.20–0.40 mm wall inside a ~12 × 14 × 7 mm blade.
Resolving that wall needs ≈ 0.05 mm elements, i.e. O(10⁷) elements over the blade —
far past what a Python solver carries. So:

- **§1.2** characterises the lattice by *direct 3-D micro-FE on the exact
  `generate.gyroid` field* — real lattice FEA, multi-cell coupons, no homogenisation
  theory in the geometry;
- **§1.3** runs the macro blade as a continuum whose element-by-element modulus comes
  from §1.2, mapped through **the generator's own grading law** (same `root_dist`,
  same `grade_len`, same `solid_root` collar).

**Nonlinearity — what is and is not included.**

| | in the model? |
|---|---|
| Contact nonlinearity (rigid frictionless platen, growing contact patch) | **yes** — active-set node-to-plane, mesh pre-rotated so the platen normal is a coordinate axis, so each contact constraint is a single Dirichlet dof and the two tangential dofs stay free. Nodes whose reaction goes tensile are released and the step re-solved. |
| Geometric nonlinearity (large displacement / rotation) | **no — and it provably cannot matter.** See below. |
| Plasticity | **no.** Reported instead as "the linear result is fiction past here". |

*Why skipping geometric nonlinearity is legitimate here, not a cop-out:* the whole point
of a geometrically nonlinear run would be to catch large-rotation softening. But the
wing reaches the entire 0.2–0.4 N design band at **39–78 nanometres** of platen travel
(§1.3), i.e. a tip rotation of order 10⁻⁵ rad over an 11.55 mm span. Second-order terms
are ~10⁻¹⁰ of first-order. A corotational or updated-Lagrangian solve would return the
same number to ten significant figures. The displacements at which geometric
nonlinearity *would* bite (0.5–2.5 mm) are displacements the part cannot reach without
being destroyed first — the linear stress there is 3.5–310× yield. **This is the honest
simpler version the brief allowed for, and this is the reason.**

Buckling is treated in §1.5 by the argument that actually settles it, rather than by an
eigenvalue that would come back at O(10⁴) N.

## 1.1 What the generator actually builds

Straight out of `PARAMS` and `G` at commit `387b8eb`:

| | value |
|---|---|
| Wing blade | 4.00 mm thick (XY, the bending direction) × 7.00 mm wide (Z, into the concha) |
| Deep-edge taper | 2.60 mm of depth at 40° (self-supporting) |
| Centreline | quadratic Bézier, (−7.50, 5.55) → (−7.50, 13.35) → (−13.27, 17.15) |
| Root plane | y = y_root = 7.15 mm (core_ry + clearance) |
| **Free span, root → tip (chord)** | **11.55 mm** (arc length beyond the root: 11.94 mm) |
| Gyroid | sheet gyroid, cell 1.20 mm, wall graded 0.40 mm (root) → 0.20 mm (face) over grade_len = 6.00 mm |
| Solid collar | fully dense Ti within 1.00 mm of the root |
| Meshed solid volume | 276.5 mm³ (free span only) |

**The first red flag is in the density.** Measured by direct integration of the
generator's own `gyroid()` field:

| wall (mm) | 0.20 | 0.25 | 0.30 | 0.35 | 0.40 |
|---|---|---|---|---|---|
| relative density ρ*/ρs | 0.426 | 0.530 | 0.664 | 0.726 | 0.852 |

This is not a lattice in the compliant sense. At the *thinnest* point of the grade it is
43% dense; at the root it is 85% dense — essentially solid metal with decorative voids.
A 1.2 mm cell simply cannot be made open at a 0.20 mm minimum printable wall: sheet-gyroid
relative density scales as ≈ 3.09 · t/a, so ρ* ≥ 0.51 for t = 0.2, a = 1.2 by construction,
before the metric-field correction.

## 1.2 Lattice characterisation — direct micro-FE on the real gyroid

Uniaxial compression of an N-cell coupon of the exact `generate.gyroid` field. Bottom
face rollered (u_z = 0, tangentially free), top face given u_z = −δ (tangentially free),
lateral faces traction-free — a free-lateral-expansion uniaxial-stress test, so
E_eff = (F/A_gross)/(δ/L) directly. Coupon 2 × 2 × 3 cells, element size = wall/4.

| wall (mm) | ρ*/ρs | elements | dof | **E_eff (GPa)** | E*/Es | K_t (peak wall stress / apparent stress) |
|---|---|---|---|---|---|---|
| 0.20 | 0.426 | 70 648 | 314 301 | **17.98** | 0.163 | 8.05 |
| 0.25 | 0.530 | 45 967 | 202 305 | **26.00** | 0.236 | 6.07 |
| 0.30 | 0.664 | 32 632 | 135 933 | **41.54** | 0.378 | 4.10 |
| 0.35 | 0.726 | 22 242 | 92 289 | **52.10** | 0.474 | 3.57 |
| 0.40 | 0.852 | 17 664 | 69 375 | **78.01** | 0.709 | 2.66 |

**Convergence** (governing thinnest wall, 0.20 mm):

| element size | wall/2 | wall/3 | wall/4 | wall/5 |
|---|---|---|---|---|
| E_eff (GPa), 1×1×2 coupon | 16.54 | 14.44 | 15.02 | 14.49 |

| coupon size | 1×1×2 | 2×2×3 | 3×3×3 |
|---|---|---|---|
| E_eff (GPa) | 15.02 | 17.98 | 19.02 |

Mesh convergence is within ~4% by wall/4. Coupon-size convergence is upward (a small
coupon is softened by its free surfaces) and has not fully plateaued at 3×3×3 — the bulk
value is probably ~19–21 GPa. **The 17.98 GPa used downstream is therefore, if anything,
an under-estimate of the stiffness**, which makes the conclusion below stronger, not weaker.
The wing is only 3.3 cells thick, so the free-surface-softened value is arguably the right
one for it anyway.

Fitted scaling law across the five points:

> **E*/Es = 0.942 · (ρ*/ρs)^2.111**   (all five points within ±5.3%)

An exponent of 2.11 is high for a sheet gyroid (1.4–1.6 is typical at low density) and
reflects that at ρ* = 0.43–0.85 the "walls" are thick blocks, not membranes.

## 1.3 Wing force–displacement

Free span outboard of y = y_root; root face encastre; rigid frictionless flat platen
advancing along the in-plane normal to the root→tip chord, n̂ = (+0.866, +0.500, 0).
That is the blade's **weak** bending axis (I = 37.3 mm⁴ vs 114 mm⁴ the other way), i.e.
the most favourable direction for a soft response. Mesh 47 418 hexes / 159 558 dof at
h = 0.18 mm. Graded modulus 18.0 GPa (face) → 110 GPa (solid root collar).

### The elastic range — this is the part that is physically real

| platen travel δ | force F | secant k | contact nodes | max von Mises |
|---|---|---|---|---|
| 0.01 µm | 0.0513 N | 5134 N/mm | 24 | 0.23 MPa |
| 0.03 µm | 0.154 N | 5134 N/mm | 24 | 0.68 MPa |
| **0.10 µm** | **0.513 N** | 5134 N/mm | 24 | 2.27 MPa |
| 0.30 µm | 1.54 N | 5134 N/mm | 24 | 6.80 MPa |
| 1.00 µm | 5.13 N | 5134 N/mm | 24 | 22.7 MPa |

Perfectly linear (constant contact patch) through this range. **k = 5134 N/mm.**

Where the design band actually falls:

| target | platen travel needed | max von Mises there |
|---|---|---|
| 0.20 N | **39 nm** | 0.88 MPa |
| 0.30 N | **58 nm** | 1.33 MPa |
| 0.40 N | **78 nm** | 1.77 MPa |

**The wing traverses the entire intended force band in 39 nanometres — about 1/300th of
typical LPBF as-printed surface roughness (Ra 10–20 µm). The design intent is a spring; the
part is, to any measurement you could make on it, rigid.**

### The requested 0.25 mm sweep to 2.5 mm

Reported as asked. **Everything from the first row onward is linear-elastic fiction:**
the von Mises stress at δ = 0.25 mm is already 3.5× yield, and by 2.5 mm it is ~310× yield.
No such force is reachable; the part (and the ear) would be destroyed first. The numbers are
here to show the *scale* of the mismatch, not as a prediction.

| δ (mm) | F (N) | secant k (N/mm) | contact nodes | max von Mises (MPa) | > 450 MPa? | > 900 MPa? |
|---|---|---|---|---|---|---|
| 0.25 | 1 791 | 7 163 | 114 | 3 120 | **flag** | **flag** |
| 0.50 | 4 256 | 8 512 | 168 | 4 345 | **flag** | **flag** |
| 0.75 | 7 752 | 10 340 | 233 | 7 033 | **flag** | **flag** |
| 1.00 | 12 520 | 12 520 | 318 | 10 440 | **flag** | **flag** |
| 1.25 | 19 710 | 15 770 | 409 | 14 500 | **flag** | **flag** |
| 1.50 | 30 430 | 20 280 | 545 | 19 030 | **flag** | **flag** |
| 1.75 | 47 760 | 27 290 | 738 | 28 210 | **flag** | **flag** |
| 2.00 | 85 760 | 42 880 | 920 | 87 120 | **flag** | **flag** |
| 2.25 | 131 200 | 58 320 | 1 001 | 185 600 | **flag** | **flag** |
| 2.50 | 182 400 | 72 950 | 1 119 | 282 700 | **flag** | **flag** |

The curve **stiffens monotonically**, it does not plateau — because as the flat platen
advances into the curved blade the contact patch grows (24 → 700+ nodes), and a growing
contact patch on a near-rigid body is a hardening contact, the opposite of the intended
behaviour.

### Stress flag, stated plainly

- **At the design-target force (0.2–0.4 N) the peak von Mises is 0.9–1.8 MPa — 0.2% of the
  450 MPa flag.** The wing is nowhere near stressed. It simply does not move.
- **At the design-target displacement (0.5–2.5 mm) the peak von Mises is 4 345–283 000 MPa,
  i.e. 5–314× yield.** Extrapolating the linear (constant-contact) regime, first yield at
  900 MPa arrives at δ ≈ 40 µm (F ≈ 200 N) and the 450 MPa fatigue flag at δ ≈ 20 µm
  (F ≈ 100 N).
- Both K_t from §1.2 (2.7–8.1) and the voxel staircase inflate the *peak*; neither changes
  the conclusion by anything like the four orders of magnitude that would be needed.

### Mesh convergence and analytic cross-check

| element size h | 0.30 mm | 0.24 mm | 0.18 mm | 0.15 mm |
|---|---|---|---|---|
| hexes | 10 324 | 20 602 | 47 418 | 82 743 |
| meshed volume (mm³) | 278.75 | 284.80 | 276.54 | 279.26 |
| **k in the elastic regime (N/mm)** | — | **5 983** | **5 134** | **5 144** |
| k at δ = 0.25 mm (N/mm) | 5 991 | 7 707 | 7 163 | 7 403 |

The elastic-regime stiffness — the number the whole verdict rests on — agrees to **0.2%**
between h = 0.18 mm and h = 0.15 mm. Converged.

Analytic cross-check, tip-loaded cantilever, 4 × 7 mm section, L = 11.55 mm,
E_eff ≈ 20 GPa: Euler–Bernoulli k = 1 455 N/mm; Timoshenko (which matters — L/t = 2.9,
this is a stubby block, not a beam) k = 1 330 N/mm (shear adds 9% compliance). The FEA figure (5 134 N/mm) is
*higher* because a flat platen against a curved blade contacts a patch part-way along
the span, not a point at the tip: the lever arm is shorter and local indentation
stiffness adds in series. The two agree to within the factor you would expect from that
difference in load introduction, which is the point of the cross-check.

## 1.4 Verdict on JOB 1's question

> *verify the graded gyroid wing produces a plateau-ish force of 0.2–0.4 N when compressed
> 0.5–2.5 mm*

**Not verified — refuted.**

| | intended | as designed |
|---|---|---|
| stiffness | ~0.1–0.4 N/mm | **5 134 N/mm** |
| travel to cover 0.2 → 0.4 N | 2.0 mm | **39 nm** |
| force at 1 mm | 0.3–0.4 N | **12 520 N** |
| shape of the curve | plateau (softening) | **hardening** |

Discrepancy: **1.3 × 10⁴ to 5 × 10⁴×**.

## 1.5 Is there a buckling / softening plateau anywhere?

No — and this is worth stating precisely, because gyroid lattices *do* famously have a
crush plateau. It is the wrong kind of plateau.

Gibson–Ashby collapse estimates at the measured relative densities
(σ_plastic ≈ 0.3 σ_y ρ^1.5, σ_elastic-buckling ≈ 0.05 Es ρ³):

| wall (mm) | ρ* | plastic collapse | elastic cell-wall buckling | governs |
|---|---|---|---|---|
| 0.20 | 0.426 | 75.0 MPa | 424.8 MPa | plastic |
| 0.25 | 0.530 | 104.2 MPa | 819.2 MPa | plastic |
| 0.30 | 0.664 | 146.1 MPa | 1 609 MPa | plastic |
| 0.35 | 0.726 | 167.1 MPa | 2 106 MPa | plastic |
| 0.40 | 0.852 | 212 MPa | 3 400 MPa | plastic |

Plastic collapse precedes elastic buckling by 5–16× at every grade. So:

1. **There is no elastic buckling plateau.** The cells yield long before they buckle. At
   ρ* ≥ 0.43 they are too stocky to buckle.
2. **There is a plastic collapse plateau — at ~75 MPa over the 4 × 7 mm section ≈ 2 100 N,
   and it is permanent.** That is an energy-absorbing crush plateau (the thing gyroids are
   used for in crash structures), 5 000–10 000× above the target force, and it destroys the
   part. It is not a return spring.
3. **Macro buckling is not the loading mode.** The blade is loaded transversely, not as a
   column, and at L/t = 2.9 it is not slender anyway.

**No mechanism in this geometry produces a reversible 0.2–0.4 N plateau.**

## 1.6 Jacket skin region

Through-thickness compression of a real 3 × 3 mm patch of the jacket wall, taken from the
full `part_jacket_wing` SDF at the core's bottom pole (the flattest part of the shell, and
clear of every magnet pocket and locating pin), so the patch carries the real perforated
0.60 mm skin plus the graded gyroid behind it. Because the shell is a curved offset of the
core, the two faces are selected on the core's own signed distance, not on a z-plane: the
core-facing face is encastre (the core behind it is 1.20 mm solid Ti, ~70× stiffer than the
shell), the ear-facing skin gets a uniform normal displacement with the tangential dofs free.

| element size h | 0.09 mm | 0.07 mm | 0.06 mm |
|---|---|---|---|
| hexes | 13 534 | 29 632 | 46 618 |
| dof | 59 523 | 121 251 | 183 924 |
| solid fill of the 1.6 mm wall | 0.685 | 0.706 | 0.699 |
| **through-thickness k over the 9 mm² patch (N/mm)** | **188 300** | **188 600** | **185 300** |

Converged to within 2%. Note this one is **direct micro-FE**: the mesh resolves the real
0.20–0.22 mm gyroid walls and the real 0.40 mm perforations, so every element is solid Ti at
110 GPa — applying the homogenised modulus here as well would knock the stiffness down for a
porosity the mesh already represents. The measured 0.699 solid fill also cross-validates §1.2
(0.6 mm skin at ~94% after perforations + 1.0 mm of lattice at ρ* ≈ 0.46 predicts 0.64–0.70).

**k = 185 300 N/mm over 9 mm², i.e. 20 590 N/mm per mm² of contact.**

Against the wing's 5 134 N/mm over its ~0.7 mm² initial contact patch (≈ 7 330 N/mm per mm²):

| basis of comparison | ratio jacket : wing | meets the >5× intent? |
|---|---|---|
| each part loaded over its own natural contact patch (absolute k) | **36×** | yes |
| like-for-like, per unit of contact area | **≈ 3×** | **no, just misses** |

The honest answer is that the ratio is not well defined, because the two parts are in
different modes — the jacket is in through-thickness compression against a rigid core, the
wing is a stubby cantilever in bending — and the answer swings between 3× and 36× depending
on how you normalise.

> *is it >5× stiffer than the wing, as intended?*

**On absolute stiffness yes, comfortably; on a like-for-like per-area basis it just misses.
And either way the comparison is not the reassuring one it was meant to be.** The intent behind "jacket ≫ wing" is that the jacket is the rigid backbone and the wing
carries all the compliance. In fact **both** are rigid: the wing is 5 134 N/mm and is
supposed to be ~0.3 N/mm. Whether the ratio is 3× or 36× means nothing, because the
denominator is four orders of magnitude too big. **This test passes and the design still
fails.**

## 1.7 What would actually hit the target

Not asked for, but it falls straight out of the same numbers and is the useful part.

**It is the *form*, not the density.** Extrapolating the fitted law
E*/Es = 0.942 ρ^2.111, the effective modulus that a 4 × 7 × 11.55 mm block would need in
order to give k ≈ 0.4 N/mm is **≈ 1.9 MPa** — which requires ρ* ≈ 0.006, i.e. a gyroid cell
of ~108 mm at the 0.20 mm minimum wall. Physically impossible, and not close. Even the
loosest lattice that fits inside a 4 mm blade (2.0 mm cell, 2 cells across, 0.20 mm wall,
ρ* ≈ 0.31, E* ≈ 8.9 GPa) still lands at ~1 830 N/mm — 4 600× too stiff. **A 4 × 7 mm block
of any printable metal, of any density, at this span, cannot be a 0.3 N spring.** 1.9 MPa is
elastomer territory (silicone is 1–10 MPa), which is the other way to read the result.

**The form that does work is a slender leaf, not a block.** For a Ti-6Al-4V cantilever at the
same L = 11.55 mm:

| leaf thickness t | width b for k ≈ 0.4 N/mm | bending stress at 0.4 N |
|---|---|---|
| 0.20 mm (= min_wall) | 2.80 mm | 247 MPa |
| 0.25 mm | 1.43 mm | 158 MPa |
| 0.30 mm | 0.83 mm | 110 MPa |

All comfortably under the 450 MPa fatigue flag. Note L/t = 58 at t = 0.20 mm — slender
enough that large-deflection behaviour is real, which is where a genuine plateau can come
from. A straight leaf would still reach ~1.0 N at 2.5 mm, so hitting a *flat* 0.2–0.4 N band
over 2 mm wants a deliberately softening element — a pre-curved (post-buckling) leaf, a
constant-force flexure, or a rolling cam — not a straight cantilever and certainly not a
lattice block.

## 1.8 Limitations, stated plainly

- **Two-scale, not direct.** The lattice was characterised on 1–3-cell coupons and carried
  into the macro blade as a graded isotropic continuum. A sheet gyroid is mildly
  cubic-anisotropic and the grading gradient is not represented at cell scale. Given the
  conclusion is a factor of 10⁴, an anisotropy correction of tens of percent is irrelevant.
- **Voxel meshes are stiff-biased** and produce staircase stress spikes. Peak von Mises
  values are over-estimates; stiffness values are mild over-estimates. Convergence is shown.
- **No geometric nonlinearity, no plasticity.** Justified in §1.0 for the elastic range and
  flagged as fiction outside it.
- **Root is idealised as encastre.** The real wing is fused to the jacket shell, which has
  some compliance, so the real k is somewhat *lower* than 5 134 N/mm. Not by 10⁴.
- **No ear.** The platen is rigid. A real ear is soft (cartilage/skin, effective modulus of
  order 10⁻¹–10¹ MPa) — which is exactly the point: with a wing this stiff, **all** of the
  compliance in the wing/ear pair lives in the ear, and the contact force is set by
  insertion geometry, not by design (see §2.5).
- **Load direction.** Only the weak in-plane bending axis was run — the softest case. Any
  other direction is stiffer.
- **σ_y = 900 MPa** as specified; LPBF as-printed is typically 950–1100 MPa. Using 1000 MPa
  moves the yield-onset displacement from ~40 µm to ~44 µm.

---

# JOB 2 — Magnet contact-pressure budget

No tissue FEA. This is a load-spreading budget: take the preload the mag-float study
actually produced, spread it over the skirt's contact band, compare against published
sustained-skin-load thresholds. Script: `fea/magnet_pressure.py`.

## 2.1 Thresholds used, and how much to trust them

| threshold | value | what it actually is | source | trust |
|---|---|---|---|---|
| **Capillary closing pressure** | **4.27 kPa (32 mmHg)** | Mean *intravascular* pressure in the arteriolar limb of nailfold capillary loops (n = 125). | Landis EM, *Heart* 1930;15:209–228. Confirmed by modern servo-null: 37.7 ± 3.7 mmHg arteriolar (Shore AC, *Br J Clin Pharmacol* 2000;50:501–513). | High as a measurement, **low as a design limit** |
| **Shear-derated sustained target** | **2.15 kPa (16 mmHg)** | Landis ÷ 2. With ~9.8 kPa of shear present, the pressure needed to occlude blood flow halves. | Bennett L, Kavner D, Lee BK, Trainor FA, *Arch Phys Med Rehabil* 1979;60(7):309–314 (n = 4, thenar eminence). | Canonical but small-n; use as a safety factor |
| **Measured injurious device pressure** | **6.3–12.3 kPa (47.6 ± 29 to 91.9 ± 42.4 mmHg)** | NIV mask on the nasal bridge, I-Scan mapping, n = 20. Discomfort correlated with measured pressure; authors note it exceeded capillary closing pressure and sometimes diastolic BP. | Brill A-K et al., *ERJ Open Res* 2018;4(2):00168-2017 | High — closest good analogue (device on thin skin over cartilage) |
| Duration inflection | ~1–2 h | Damage is deformation-driven at short durations, ischaemia-driven at long ones; the critical load falls sharply around 1–2 h. | Kosiak 1959/1961; Gefen A et al., *J Biomech* 2008; Gefen, *Ostomy Wound Manage* 2009;55(7) & 55(9) | High for the shape, weak for exact numbers |
| Momentary pressure-discomfort threshold, concha | ≫ 100 kPa | Algometer PDT mapped at 6 concha points, n = 80; **tragus is the most sensitive site**. | Yuan X, Wang Z, Feng F et al., "Measurement of pressure discomfort threshold in auricular concha for in-ear wearables design," *Applied Ergonomics* 2023;113:104078 | Citation verified; **the per-site kPa map is paywalled and was not verified** |

**Important caveats, so these numbers are not misused.**

- 32 mmHg is an *intravascular* pressure in a fingertip capillary, not an externally applied
  interface pressure and not a comfort threshold. Fletcher (*Wounds UK* 2023;19(3):52–61) is
  quotable: *"it is not a universal discriminator for harm."* Wheelchair-seating reviews note
  that in practice **no** cushion achieves < 32 mmHg at the ischium. **Treat 4.27 kPa as an
  ischaemia flag, not a pass/fail.**
- **PDT/PPT is not a wear-time limit.** It is a momentary pain threshold from a seconds-long
  ramp and is strongly probe-area dependent — 1–2 orders of magnitude above the sustained
  threshold. (The widely-circulated 564 kPa ear-region figure from Shah & Luximon 2021 is
  inflated by a 3 mm probe; the 133–193 kPa concha range attributed to Yuan 2023 is
  unverified secondary. Neither is used here.)
- **Tissue tolerates far more when loading is intermittent.** Transtibial prosthetic sockets
  routinely see 121 ± 32 kPa static / 254 ± 61 kPa dynamic without injury (Al-Fakih et al.,
  *Sensors* 2016;16(7):1119). An IEM taken out every couple of hours is in a much kinder
  regime than a continuously-worn NIV mask. Hourly removal is genuinely protective;
  30-second adjustments are not.
- **The ear is one of the two most common medical-device-related pressure injury sites**
  (with the feet) — oxygen tubing over the helix is the classic cause. Gefen et al., "SECURE
  prevention," *J Wound Care* 2020;29(Sup2a):S1–S52 is the consensus document.

**Verdict bands used below** (deliberately conservative, aimed at all-day wear):

| | |
|---|---|
| ≤ 2.15 kPa | **comfortable** — safe even with shear |
| 2.15 – 4.27 kPa | **borderline** — under the ischaemia flag but with no shear margin |
| 4.27 – 6.34 kPa | **too much** for sustained wear |
| > 6.34 kPa | **too much (≫)** — in or above the band measured to injure skin over the nasal bridge |

## 2.2 Skirt geometry and preload

From `generate.py` (`G.__init__`, preset `asym_as_built`) and `docs/MAGFLOAT_MAGNETS.md`:

| | |
|---|---|
| Rim outer diameter | 19.00 mm (mid-wall radius 9.325 mm) |
| Root diameter | 10.50 mm (= carrier OD) |
| Flare half-angle | 35.0° from the axis |
| Axial run | 5.82 mm (x = 7.83 → 13.65) |
| **Available slant length** | **7.10 mm** — the widest band the funnel can physically offer |
| Skirt wall | 0.35 mm |
| Preload (fixed 7×3×1.5 / moving 9×5.4×2.0, bonded NdFeB, rest gap 2.25 mm) | **0.307 N max / 0.214 N rest / 0.155 N min** |

**Two pressures are reported, and the difference matters.**

- **Nominal** P = F / (π · D · w) — the figure the brief asked for: axial force over the
  wetted band area.
- **Cone-normal** P = F / (π · D · w · sin 35°) — the physically correct one. For a cone
  whose meridian makes 35° with the axis, only sin 35° = 0.574 of the surface normal is
  axial, so the true contact pressure is **1.743×** the nominal. All verdicts below are
  taken on the cone-normal figure, which is the conservative choice.

## 2.3 Contact pressure — F = 0.307 N (max, full compression)

Each cell: *nominal kPa | cone-normal kPa | verdict*.

| D_contact | w = 0.5 mm | w = 1.0 mm | w = 2.0 mm | w = 3.0 mm |
|---|---|---|---|---|
| **10 mm** | 19.54 \| **34.07** \| too much ≫ | 9.77 \| **17.04** \| too much ≫ | 4.89 \| **8.52** \| too much ≫ | 3.26 \| **5.68** \| too much |
| **13 mm** | 15.03 \| **26.21** \| too much ≫ | 7.52 \| **13.11** \| too much ≫ | 3.76 \| **6.55** \| too much ≫ | 2.51 \| **4.37** \| too much |
| **16 mm** | 12.22 \| **21.30** \| too much ≫ | 6.11 \| **10.65** \| too much ≫ | 3.05 \| **5.32** \| too much | 2.04 \| **3.55** \| borderline |
| **19 mm** | 10.29 \| **17.93** \| too much ≫ | 5.14 \| **8.97** \| too much ≫ | 2.57 \| **4.48** \| too much | 1.71 \| **2.99** \| borderline |

## 2.4 Contact pressure — rest (0.214 N) and minimum (0.155 N)

**F = 0.214 N (rest):**

| D_contact | w = 0.5 mm | w = 1.0 mm | w = 2.0 mm | w = 3.0 mm |
|---|---|---|---|---|
| **10 mm** | 13.62 \| **23.75** \| ≫ | 6.81 \| **11.88** \| ≫ | 3.41 \| **5.94** \| too much | 2.27 \| **3.96** \| borderline |
| **13 mm** | 10.48 \| **18.27** \| ≫ | 5.24 \| **9.14** \| ≫ | 2.62 \| **4.57** \| too much | 1.75 \| **3.05** \| borderline |
| **16 mm** | 8.51 \| **14.85** \| ≫ | 4.26 \| **7.42** \| ≫ | 2.13 \| **3.71** \| borderline | 1.42 \| **2.47** \| borderline |
| **19 mm** | 7.17 \| **12.50** \| ≫ | 3.59 \| **6.25** \| too much | 1.79 \| **3.13** \| borderline | 1.20 \| **2.08** \| **comfortable** |

**F = 0.155 N (min, far travel):**

| D_contact | w = 0.5 mm | w = 1.0 mm | w = 2.0 mm | w = 3.0 mm |
|---|---|---|---|---|
| **10 mm** | 9.87 \| **17.20** \| ≫ | 4.93 \| **8.60** \| ≫ | 2.47 \| **4.30** \| too much | 1.64 \| **2.87** \| borderline |
| **13 mm** | 7.59 \| **13.23** \| ≫ | 3.80 \| **6.62** \| ≫ | 1.90 \| **3.31** \| borderline | 1.27 \| **2.21** \| borderline |
| **16 mm** | 6.17 \| **10.75** \| ≫ | 3.08 \| **5.38** \| too much | 1.54 \| **2.69** \| borderline | 1.03 \| **1.79** \| **comfortable** |
| **19 mm** | 5.19 \| **9.05** \| ≫ | 2.60 \| **4.53** \| too much | 1.30 \| **2.26** \| borderline | 0.87 \| **1.51** \| **comfortable** |

## 2.5 Wing contact pressure, for comparison

The brief asked for the wing's pressure at its intended 0.3–0.4 N plateau over an 8–12 mm
contact patch. The blade is 7.0 mm wide in Z, so "width 7 mm" is full-face contact and
"width 1 mm" is an edge line contact.

| F | patch 8 mm × 1 mm | × 2 mm | × 4 mm | × 7 mm |
|---|---|---|---|---|
| 0.30 N | 37.50 kPa ≫ | 18.75 kPa ≫ | 9.38 kPa ≫ | 5.36 kPa too much |
| 0.40 N | 50.00 kPa ≫ | 25.00 kPa ≫ | 12.50 kPa ≫ | 7.14 kPa ≫ |

| F | patch 12 mm × 1 mm | × 2 mm | × 4 mm | × 7 mm |
|---|---|---|---|---|
| 0.30 N | 25.00 kPa ≫ | 12.50 kPa ≫ | 6.25 kPa too much | **3.57 kPa borderline** |
| 0.40 N | 33.33 kPa ≫ | 16.67 kPa ≫ | 8.33 kPa ≫ | 4.76 kPa too much |

**Only one case clears** — 0.30 N spread over the full 12 × 7 mm face, at 3.57 kPa, and only
as *borderline*. Everything narrower is over.

**And there is a much bigger problem, which JOB 1 exposes.** A 0.3 N wing pressure assumes
the wing is a 0.3 N spring. It is not (§1.4): it is a rigid wedge at 5 134 N/mm. **A rigid
wedge has no design contact pressure at all** — the force is whatever the ear's own
deformation generates against a fixed geometry, i.e. it is set by the insertion depth and
the individual's concha size, and it varies across the population by however much the
anthropometry varies. Over the design envelope in `docs/EAR_ANTHROPOMETRY.md` (concha width
10–24 mm, a 2.4:1 range), a fixed rigid wing that is comfortable on a large concha will be a
crushing interference fit on a small one. Every number in this section is contingent on the
wing first being made compliant.

## 2.6 Minimum contact band width

Slant band width needed to hold the cone-normal pressure at or under each threshold, at
**F_max = 0.307 N**:

| D_contact | for 2.15 kPa (shear-derated) | for 4.27 kPa (capillary) | for 6.34 kPa (NIV low) | available on the funnel |
|---|---|---|---|---|
| 10 mm | **7.92 mm** — exceeds the funnel | 3.99 mm | 2.69 mm | 7.10 mm |
| 13 mm | 6.10 mm | 3.07 mm | 2.07 mm | 7.10 mm |
| 16 mm | 4.95 mm | 2.49 mm | 1.68 mm | 7.10 mm |
| 19 mm | 4.17 mm | 2.10 mm | 1.41 mm | 7.10 mm |

**Answers to "the minimum contact band width required to stay under threshold at max force":**

- **Against the capillary-closing flag (4.27 kPa), across every contact diameter from 10 to
  19 mm: w ≥ 4.0 mm.** Well inside the 7.10 mm the funnel offers.
- **Against the shear-derated target (2.15 kPa): w ≥ 7.9 mm at Ø10 mm, which this skirt
  cannot deliver** (7.10 mm available). It is reachable at Ø13 mm and above (6.1 / 5.0 /
  4.2 mm). So the small-ear case — which is exactly the population tail the anthropometry
  doc flags — is the one that cannot be brought under the shear-derated target with the
  present 19 mm / 35° funnel.

## 2.7 Verdict on JOB 2

**The force is fine. The geometry does not currently commit to spreading it.**

0.155–0.307 N is a small force; over 4 mm of band at any of these diameters it is
comfortable-to-borderline. But `generate.py` specifies the skirt as a continuous 35° cone
with a 0.35 mm wall and **does not define a contact land at all** — the band width is
whatever the ear happens to produce. If it lands as a line contact (w ≈ 0.5 mm, which is
what a thin cone against a curved ear canal aperture will naturally do), the pressure is
**17.9–34.1 kPa — 4 to 8× the ischaemia flag, and 1.5–5× the band measured to injure skin
over the nasal bridge.**

Recommendations, in order:

1. **Add an explicit compliant contact land of ≥ 4 mm slant width** to the skirt, at the
   diameters a small ear will touch (Ø10–13 mm). This is the single change that moves every
   cell of §2.3 from "too much" to "borderline or better".
2. **Soften the land, don't just widen it.** A 0.35 mm silicone wall on a 35° cone is stiff
   in hoop; a thinner, corrugated, or bellows land distributes contact along its width
   instead of peaking at the leading edge (see the limitation below).
3. **Design to 4.27 kPa, log the shear-derated 2.15 kPa as the stretch target**, and accept
   that at Ø10 mm the stretch target needs either a larger flare or a lower preload.
4. Consider the **9×5.4×1.5 alternative** in `MAGFLOAT_MAGNETS.md` (0.405 → 0.288 → 0.199 N)
   only if the band width goes up with it — its higher max force makes the pressure problem
   worse, not better.

## 2.8 Limitations, stated plainly

- **P = F/A assumes uniform pressure over the band. It will not be uniform.** A thin
  silicone cone pressed against tissue peaks at the leading edge of contact; the real peak
  is plausibly 1.5–3× the average. That makes every "borderline" here optimistic.
- **No tissue model**, as scoped. Tissue is nonlinear, viscoelastic and it relaxes — the
  pressure an hour after insertion is lower than at insertion. Not modelled.
- **The contact diameter is an input, not a prediction.** Which diameter a given ear touches
  depends on aperture size and angle (azimuth 5–53°, elevation −70° to +50° per
  `EAR_ANTHROPOMETRY.md`); 10–19 mm brackets it, it does not resolve it.
- **Shear is applied as a blanket ÷2 factor** from a 4-subject 1979 study. It is the standard
  citation and it is weak. Treat 2.15 kPa as a stretch target, not a measurement.
- **Duration is not modelled.** All thresholds here are "sustained". An IEM removed every
  1–2 hours sits in a much more forgiving regime; the pressure–time literature (Kosiak,
  Reswick & Rogers, Gefen's sigmoid critique) says so but does not give a number this
  analysis could apply.
- **Yuan 2023's per-site concha PDT map was not obtained** (paywalled). It is the single most
  on-point reference for this device and is worth buying before the next round.

---

# Macro-gyroid wing (as-redesigned)

Re-run 2026-08-31 against `cad/iem/generate.py` at commit `8db64d7`, in which the wing
was rebuilt as a compliant shell in response to §1.4/§1.7. **The wing changed again in
`dbc3b7b` while this was solving — see the box below.** Script
`cad/iem/fea/wing2_stiffness.py`; run log `fea/wing2_run.log`, results
`fea/wing2_results.json`.

**Headline: the redesign works as a concept and is the right shape of thing — but it is
still about 4× stiffer than the generator's own estimate and about 3–7× stiffer than the
0.15–0.35 N/mm target. k = 1.06 ± 0.03 N/mm, F(1.0 mm) = 1.09 N.**

> ### ⚠ Geometry has moved since this run
>
> This section was solved against commit `8db64d7`. Commit `dbc3b7b` ("shorten + splay the
> wing") then changed it again: **free span 13.161 → 10.398 mm** and the root wall
> 0.22 → 0.20 mm. Everything below still describes `8db64d7`.
>
> The generator's own estimate for the two geometries is 0.2506 → **0.2941 N/mm**. Using
> that model as a transfer function (same model both sides, so the *ratio* is far more
> trustworthy than either absolute), the FEA number scales to
>
> > **k ≈ 1.06 × (0.2941 / 0.2506) ≈ 1.24 N/mm** for the current wing,
>
> which is **3.5–8.3×** the 0.15–0.35 N/mm target — i.e. the shortening made it *stiffer*,
> and the gap slightly worse. Shortening the span was the opposite of recommendation 1 in
> §3.4 (k ∝ L⁻³). That projection is an estimate, not a solve; re-run
> `fea/wing2_stiffness.py` (~90 min) to confirm before acting on it.

## 3.1 What changed, and why the FEA had to change with it

| | old (§1.1–1.6) | redesigned |
|---|---|---|
| structure | 4 × 7 mm blade filled with a 1.2 mm gyroid | 7 × 5 mm envelope holding a **12 mm** gyroid — one to two cells, i.e. a doubly-curved sheet |
| wall | 0.20 → 0.40 mm graded | 0.22 → 0.20 mm, rolled to 0.40 mm on every exposed edge |
| relative density | 0.43 – 0.85 | **0.054** |
| root | 1.0 mm solid collar | 1.2 mm solid plug into the jacket rim |
| free span | 11.55 mm | **13.16 mm** |
| solid volume | 276.5 mm³ | **36.4 mm³** |

Four method changes follow, and each one moves the answer materially:

**(a) The mesh has to resolve the wall.** A 0.20 mm sheet spanning 13 mm is an
aspect ratio of 65. Voxel hexes were run at h = 0.12 / 0.10 / 0.08 mm — 1.7 / 2.0 / 2.5
elements through the wall.

**(b) Direct solver, not AMG.** Smoothed-aggregation *and* rootnode AMG, both with the
six rigid-body modes as near-nullspace, **stall** on this shell — 2000 CG iterations to a
relative residual of only 5.6 × 10⁻⁷ (446 s), because a thin shell's near-nullspace is
dominated by bending modes that rigid-body vectors do not span. SuperLU factorises the
same 100k-dof system in **17 s**. `_common.solve_fixed` now prefers direct up to 260k dof.

**(c) Geometric nonlinearity is now in the model.** For the old block the entire force
band lived inside 78 nm, so large-displacement terms were provably irrelevant. Here the
tip moves 1–2 mm on a 0.20 mm sheet — 5–10× the wall thickness — so a total-Lagrangian
**St Venant-Kirchhoff** solve was implemented (`set_config` / `_residual` / `_tangent`,
analytic material tangent A_iJkL = δ_ik S_JL + F_iM C_MJNL F_kN). *Verification: at u = 0
the assembled nonlinear tangent reproduces the linear stiffness matrix to
5.8 × 10⁻¹⁶ relative — machine precision.* Newton converges quadratically, 4–6 iterations
per step.

**(d) The load must not clamp the cross-section.** This one is worth stating loudly.
Prescribing a uniform normal displacement on every node of the tip patch — the obvious
"rigid platen" boundary condition — returns **k = 20.8 N/mm, twenty times the true
value**, because for an *open* thin section cross-section distortion is most of the
compliance and clamping the load patch removes it. The sweep instead uses **scalar
displacement control**: a distributed tip load pattern c, with the work-conjugate mean tip
deflection c·u driven to the target through a bordered system, so the load factor λ *is*
the platen force. This loads the tip without stiffening it, and it walks through limit
points, so a snap-through would appear as a negative tangent rather than as divergence.
The linear branch of the same routine reproduces the independent force-controlled
stiffness to 5 significant figures.

## 3.2 Linear tip stiffness, and mesh convergence

| element size h | 0.12 mm | 0.10 mm | 0.08 mm |
|---|---|---|---|
| elements through the 0.20 mm wall | 1.7 | 2.0 | 2.5 |
| hexes | 21 240 | 36 398 | 71 575 |
| meshed volume (mm³) | 36.70 | 36.40 | 36.65 |
| **k (N/mm)** | **1.0771** | **1.0159** | **1.0808** |

> **k = 1.058 ± 0.030 N/mm (spread 6.1%)**

Convergence is *non-monotone*, and that is expected rather than alarming: voxelising a
0.20 mm sheet makes the represented wall snap between 2 and 3 voxels as the grid moves
relative to the surface, and bending stiffness goes as t³. The meshed volume is stable to
0.8%, so the mass is right; it is the local t³ that jitters. The ±3% spread, not a
Richardson extrapolation, is the honest error bar. (h = 0.08 costs ~25 min of SDF
sampling and is run with `--fine`.)

### Why the generator's estimate is low: it is a lower bound by construction

Two independent bounds bracket the FEA, and they explain the 4× gap:

| model | k (N/mm) | what it assumes |
|---|---|---|
| generator `wing_report`, D·chord·χ with χ = 0.40 | **0.251** | the sheet contributes only its own local plate rigidity E t³/12(1−ν²) — **lower bound**, it ignores that material offset from the neutral axis carries bending |
| **FEA (this work)** | **1.06** | — |
| rigid-section, ∫(z−z̄)²dA measured off the mesh | **21.1** | plane sections remain plane — **upper bound**, it ignores cross-section distortion of an open section |

The FEA sits between them, 4.2× above the lower bound and 20× below the upper. So the
generator's model is not "wrong by a fudge factor" — it is structurally the wrong model,
and no choice of χ fixes it: even χ = 1.0 gives only 0.627 N/mm. Tuning χ up to ~1.7 would
reproduce the FEA number for *this* geometry, but it would not generalise, because what
the FEA is actually capturing is the offset second moment partially resisted by section
distortion, which is not proportional to D·chord at all.

## 3.3 Force–displacement, 0 → 2 mm

Displacement-controlled, h = 0.12 mm. Stress is the von Mises of the **Cauchy** stress.
Two columns are given: `vM max` is the single worst quadrature point, which on a voxel
mesh sits on a staircase corner and is a known over-estimate; `vM p99.9` is the 99.9th
percentile and is the number to design to.

| δ (mm) | F linear (N) | **F nonlinear (N)** | k_secant | k_tangent | vM p99.9 (MPa) | vM max (MPa) | > 450 MPa? |
|---|---|---|---|---|---|---|---|
| 0.25 | 0.2693 | **0.2695** | 1.078 | — | 99.7 | 183.7 | no |
| 0.50 | 0.5386 | **0.5400** | 1.080 | 1.082 | 199.6 | 364.8 | no |
| 0.75 | 0.8078 | **0.8123** | 1.083 | 1.089 | 299.6 | 543.3 | no |
| **1.00** | 1.0771 | **1.0872** | 1.087 | 1.100 | **400.0** | 719.4 | no (89% of flag) |
| 1.25 | 1.3464 | **1.3655** | 1.092 | 1.113 | 500.7 | 893.2 | **flag** |
| **1.50** | 1.6157 | **1.6481** | 1.099 | 1.130 | **601.5** | 1064.7 | **flag** |
| 1.75 | 1.8849 | **1.9358** | 1.106 | 1.151 | 703.4 | 1234.1 | **flag** |
| 2.00 | 2.1542 | **2.2297** | 1.115 | 1.175 | 806.5 | 1401.5 | **flag** |

**Behaviour: a smooth, mildly stiffening spring. No snap-through, no buckling, no
softening.** The tangent stiffness rises monotonically from 1.077 to 1.175 N/mm — a 9%
membrane-stiffening effect over the full 2 mm, and geometric nonlinearity adds only
0.1–3.5% to the force. The tangent never goes negative and Newton never needed a limit
point handled, so there is no instability anywhere in 0–2 mm. That is a genuinely good
result for the *character* of the spring; the problem is purely its magnitude.

## 3.4 Verdict

| | target | generator estimate | **FEA** |
|---|---|---|---|
| stiffness k | 0.15 – 0.35 N/mm | 0.251 N/mm | **1.06 ± 0.03 N/mm** |
| F at 1.0 mm | 0.15 – 0.35 N | 0.25 N | **1.09 N** |
| F at 1.5 mm | 0.23 – 0.53 N | 0.38 N | **1.65 N** |
| character | plateau-ish | linear | smooth, mildly stiffening |

**Verdict: FAIL on magnitude, PASS on character.** k is **3.0× the top of the target band
and 7.1× the bottom**. The generator's estimate is optimistic by 4.2×, for the structural
reason in §3.2 — this is a model error, not a tuning error.

The redesign is nonetheless a huge improvement on the old block — ~1250× softer than the
old blade's analytic weak-axis stiffness (1330 N/mm), ~4800× softer than the 5134 N/mm
§1.3 reported (see §3.6 on why those differ) — and is now the right *kind* of structure: a thin-shell spring with a smooth, stable,
near-linear response. It needs one more factor of 3–7 in compliance, not another
architecture.

### Stress margin

At the **target** force (0.35 N) the wing is at δ ≈ 0.33 mm and vM p99.9 ≈ 130 MPa — a
3.5× margin on the 450 MPa fatigue flag and 6.9× on yield. Comfortable.

At the **actual** force it delivers over its working travel, the margin is thin: the
450 MPa flag is crossed at **δ ≈ 1.12 mm** and yield (900 MPa, on p99.9) at δ ≈ 2.2 mm.
So as built the wing has adequate static margin to 2 mm but only ~1.1 mm of
fatigue-safe travel — less than the 1.5 mm the mag-float travel implies. **Softening it to
target fixes the stress margin at the same time**, because both scale with the same
sheet thickness.

### What to change

Bending stiffness of a shell goes as t³ at fixed geometry, so the cheapest lever is wall
thickness — but 0.20 mm is already `min_wall`. The remaining levers, in order of
leverage:

1. **Lengthen the free span.** k ∝ L⁻³ for tip loading: 13.16 → 17.5 mm (+33%) buys the
   full 2.4× needed to reach the middle of the band, with no change in wall or process.
2. **Narrow the load-bearing chord** — reduce `wing_thick` (7 mm) or open the macro cell
   further (12 → 16–18 mm) so fewer sheet chords cross the section. k is roughly linear
   in chord, so 7 → 3 mm would give ~2.3×.
3. **Neck the anchor harder.** `wing_anchor_w` = 2.40 mm over `wing_anchor_len` = 7 mm
   already softens the foot; deepening that is the least disruptive knob.
4. Do **not** reach for a plateau mechanism yet. The response is smooth and stable; get
   the magnitude right first, then decide whether a flat band is worth a pre-curved
   (post-buckling) element.

## 3.5 A hard flat platen is a different, worse load case

The sweep above loads the tip, matching the generator's tip-stiffness model. A rigid flat
platen — the limiting case of a *hard* ear — behaves differently, and worth recording:

| δ (mm) | F (N) | k (N/mm) | contact nodes | vM p99.9 (MPa) | vM max (MPa) |
|---|---|---|---|---|---|
| 0.05 | 0.2257 | 4.514 | 6 | 47.1 | 106.4 |
| 0.10 | 0.4509 | 4.509 | 2 | 94.0 | 209.0 |
| 0.25 | 1.1273 | 4.509 | 2 | 235.1 | 522.6 |

Linear elasticity with the active-set frictionless contact of §1.0, h = 0.12 mm. Two
things differ from the tip-load case:

* **the platen does not touch the tip.** The outermost point of the wing along the press
  normal sits at **s/L_free = 0.633** — mid-span. The lever arm is ~⅔ of the span, and
  k ∝ L⁻³, so that alone accounts for most of the 4.2× rise to **4.51 N/mm**;
* **contact is 2–6 nodes**, i.e. 0.03–0.09 mm². The platen is locally dimpling the sheet,
  not bending the wing. The contact set *shrinks* with stroke (6 → 2) rather than
  spreading, which is the signature of local indentation rather than global bending.

A real ear is soft (cartilage/skin, order 10⁻¹–10¹ MPa) and will wrap and spread the load
rather than point-load a 0.2 mm sheet, so the tip-load case is the better design model.
But the platen case flags a real manufacturing/robustness concern: **this shell has very
little local indentation resistance**, and anything hard and small pressing on it (an
insertion tool, a drop onto a hard edge, a fingernail) loads a few mm² of 0.2 mm Ti.
The 0.40 mm rolled edges help; the sheet faces are unprotected.

## 3.6 Correction to §1.3 (old wing)

While building this model I found a frame bug in the original `run_wing`: it computed the
rotated node coordinates `p_rot = R @ m.p` to choose contact nodes and prescribed values,
but assembled K on the **unrotated** mesh and constrained the global-z dof. The press
therefore acted along global z, not along the stated n̂ = (0.866, 0.500, 0). Fixed in
`_common.rotate_mesh`, which rotates the mesh itself (isotropic elasticity is
frame-indifferent, so this is exact) and is now used by both wing scripts.

**This does not change any conclusion in §1.** The number reported there — 5134 N/mm —
was measured in an oblique direction rather than the weak bending axis, so it is not the
weak-axis stiffness it was labelled as. What carries §1's conclusion is the independent
analytic bound: a 4 × 7 mm section on an 11.55 mm span has a Timoshenko weak-axis
stiffness of **1330 N/mm**, which is the *softest* the old block could possibly have been
and is still ~4400× the 0.3 N/mm target. §1.4's "1.3 × 10⁴ to 5 × 10⁴×" should be read as
"~4 × 10³ to 5 × 10⁴×, depending on load direction". The old block was superseded by
commit `8db64d7` and was not re-run, because `generate.py` no longer builds it.

## 3.7 Limitations

- **Voxel hexes on a thin shell.** 1.7–2.5 elements through the wall; trilinear hexes lock
  in bending (biasing stiffness **up**) while thickness quantisation jitters it ±3%. A
  mid-surface shell-element model would be the better tool and would likely land somewhat
  *below* 1.06 N/mm — i.e. the true value is probably at the low end of, or just under,
  the quoted band. It would not close a 3–7× gap.
- **St Venant-Kirchhoff, no plasticity.** Valid at these strains (<1%); the material is
  elastic everywhere below the yield flag, and rows past δ ≈ 2.2 mm would need plasticity.
- **Encastre at y = y_root.** The real root plug merges into the jacket rim, which has
  some compliance, so the real k is somewhat *lower* than reported.
- **One load direction** (the in-plane normal to the root→tip chord, the design press
  direction). Loads out of that plane were not swept.
- **No ear.** Rigid platen or distributed tip load; no tissue model.
- **h = 0.08 mesh point** comes from a standalone run; the default `main()` sweeps
  h = 0.12 and 0.10 and needs `--fine` for the third.


---

# Ti-sheet wing: parameters that pass (superseded, kept for the record)

Run 2026-08-31. Script `cad/iem/fea/wing2_stiffness.py` (`screen()` / `run_verify()`),
raw output `fea/wing2_candidates.json`, `fea/sweep_*.json`.

**The user has since chosen the mag-plunger wing (§5), so this route is not being built.
It is recorded because it did reach the target and the knob calibration is reusable.**

## 4.1 Baseline confirmed on `dbc3b7b`

FEA on the current geometry gives **k = 1.4008 N/mm** (L_free = 10.40 mm), against the
1.24 N/mm projected in the §3 box — the scale-up was the right direction but 13% low.
vM p99.9 at 1.5 mm = 822 MPa. Fails both criteria.

## 4.2 Knob calibration (FEA-measured, not assumed)

| knob | effect | verdict |
|---|---|---|
| **free span L** | k ∝ L⁻²·⁹⁵, σ ∝ L⁻¹·⁶² (fitted over L = 10.40–14.35 mm) | **the only strong knob, and it helps both** |
| **`wing_edge_band`** 0.70 → 0.25 mm | k × 0.60, σ × 0.81 at fixed span | **second real knob** — shrinks the stiff 0.40 mm rolled-edge zone without thinning the edge itself |
| `wing_root_solid` 1.2 → 0.6 mm | k × 0.89, σ × 0.91 | mild, free |
| `gyroid_cell_wing` 12 → 15 mm | k × 1.06 (*stiffer*); ≥18 mm **disconnects the sheet** (voxel volume halves) | **unusable** — cell must stay 12.0 |
| `wing_thick` 7 → 5.5 mm | k × 0.80 but σ × 1.10; 4.5 mm degenerates | bad trade |
| `wing_anchor_w` 1.4 → 0.7 mm | k × 0.85, σ flat | weak |
| **relief slits** (prototyped in the FEA script, not in `generate.py`) | n = 1: k × 0.85 but σ **up** 4%. n ≥ 2: the ribbons disconnect — 39–55% of voxels dropped | **rejected** |

Slits fail because slitting a plate into side-by-side ribbons does not reduce bending
stiffness about the same axis (n ribbons of width b/n still sum to bt³/12); the only gain
is losing the plate factor and the anticlastic constraint, and that is swamped by the
stress concentration at each new edge.

## 4.3 The passing parameter sets

Both hold `gyroid_cell_wing` = 12.0 and the wall at the 0.20 mm print floor.

| | **candidate B (recommended)** | candidate A (gentler edge ramp) |
|---|---|---|
| `wing_shorten` | 0.0 | 0.0 |
| `wing_rise` | 11.7 | 11.7 |
| `wing_edge_band` | **0.20** | **0.25** |
| `wing_root_solid` | 0.6 | 0.6 |
| resulting L_free | 13.99 mm | 13.99 mm |
| **k (3-mesh converged)** | **0.2809 ± 0.0081 N/mm** | **0.3040 ± 0.0111 N/mm** |
| F at 1.0 mm | 0.2689 N | 0.2878 N |
| F at 1.5 mm | 0.4045 N | 0.4331 N |
| **vM p99.9 at 1.5 mm** | **375.9 MPa** | **387.1 MPa** |
| vM max at 1.5 mm (staircase) | 806.5 MPa | 849.4 MPa |
| **verdict** | **PASS** k and stress | **PASS** k and stress |

Mesh convergence for B: k = 0.2848 / 0.2696 / 0.2881 N/mm at h = 0.12 / 0.10 / 0.08 mm
(1.7 / 2.0 / 2.5 elements through the 0.20 mm wall), spread 6.6%. Response is smooth and
near-linear — k_tangent drifts only 0.2686 → 0.2776 N/mm over the full 2 mm, no
snap-through.

**The cost:** the span has to go back to 13.99 mm, i.e. the v4 shortening (`wing_shorten`
= 2.75 mm, added because 44% of ears were overpressed) has to be given back in full. That
is the constraint that binds this route — a compliant Ti sheet needs length, and the ear
does not have it to give. It is also why the plunger won.

---

# Mag-plunger wing (the chosen design)

User-proposed replacement for the Ti spring: a rigid contact rail (~14 × 6 mm, soft pad on
top) riding on guide pins from the jacket, pushed toward the antihelix by repelling
bonded/sintered NdFeB pairs. A cam preset takes up coarse per-ear engagement, so the
magnets only cover dynamic travel about a rest point.

Script `cad/iem/fea/mag_plunger.py`, results `fea/mag_plunger_results.json`.
Method is identical to `docs/MAGFLOAT_MAGNETS.md` (magpylib 5.0.1 + magpylib-force 0.3.1,
Maxwell-stress/meshed-dipole direct force solve, 200-cell meshing, rings as
`CylinderSegment` φ = 0–360, discs as `Cylinder`, like poles facing).

**Method check:** reproduces the published 7×3×1.5 mm N52 row of `MAGFLOAT_MAGNETS.md` to
within **0.04%** at every gap (5.625 / 3.328 / 1.490 / 0.821 N at 0.5 / 1.0 / 2.0 / 3.0 mm).
*(magpylib 5 is SI — lengths in metres. Feeding it millimetres scales every force by 1e6;
this check is what caught that.)*

## 5.1 The ±1.5 mm requirement does not fit in 5 mm of depth

The force band 0.15–0.5 N permits a max/min **ratio of only 3.33**. Over ±1.5 mm the face
gap swings from (rest − 1.5) to (rest + 1.5), so flatness needs rest ≫ 1.5 mm. But the
stack depth is 2t + rest, and the budget is 5 mm.

None of the four listed pucks fits. Closest miss, 2 pairs, at the deepest in-budget rest
gap (2.95 mm):

| puck | material | ratio achieved | F near → far |
|---|---|---|---|
| 3×1 | bonded / N35 | 11.15 | 0.151 → 0.014 N |
| 4×1 | bonded / N35 | 8.12 | 0.279 → 0.034 N |
| 4×2×1 | bonded / N35 | 8.25 | 0.141 → 0.017 N |
| **5×2.5×1** | bonded / N35 | **7.54** | 0.242 → 0.032 N |

Bigger magnets flatten the curve (their field decays over a longer length scale) but cost
depth. Ratio over ±1.5 mm, N35 (needs ≤ 3.33):

| puck | r=2 | r=3 | r=4 | r=5 | r=6 | r=7 | max in-budget rest gap |
|---|---|---|---|---|---|---|---|
| 5×2.5×1 | 16.49 | 7.34 | 5.04 | 4.15 | 3.65 | 3.29 | 3.0 mm |
| 8×4×1.5 | 9.60 | 5.56 | 3.94 | 3.19 | 2.81 | 2.58 | 2.0 mm |
| 10×5×2 | 7.84 | 4.92 | 3.62 | 2.95 | 2.57 | 2.35 | 1.0 mm |
| 14×7×3 | 5.29 | 3.88 | 3.12 | 2.64 | 2.33 | 2.12 | — |

**The binding constraint is depth, not magnet strength.** The minimum-depth design that
does hold the band over ±1.5 mm is **8×4×1.5 mm N35, one pair, rest 4.95 mm** —
0.499 → 0.262 → 0.155 N, ratio 3.22, 848 mg — but it needs **7.95 mm of stack depth,
1.6× the budget.**

## 5.2 Recommended spec (fits the 5 mm budget, at ±0.75 mm travel)

> ### **5 × 2.5 × 1 mm N35 sintered NdFeB rings — TWO pairs — rest gap 2.75 mm**

| | |
|---|---|
| Travel covered | **±0.75 mm** dynamic (1.5 mm total), not ±1.5 mm |
| Force, both pairs summed | **0.493 N** at −0.75 mm → **0.284 N** at rest → **0.180 N** at +0.75 mm |
| Ratio over travel | **2.73** (band allows 3.33) |
| Axial stiffness at rest | 0.186 N/mm |
| **Stack depth** | **4.75 mm** (2 × 1.0 mm pucks + 2.75 mm gap) — fits the 5.0 mm budget |
| Guide pins | 2.0 mm OD × 5.0 mm engaged → **L/D = 2.5** (≥ 2 satisfied) |
| Rocking restoring stiffness | **4.66 N·mm/rad** = 0.081 N·mm/deg, two plungers 10 mm apart |
| Magnet mass | 442 mg (4 pucks) |
| Total added mass | **971 mg/side** — magnets 442 + Ti rail 14×6×0.8 298 + 2 pins 139 + pad 92 |
| Stray field at the bone sensor, 8 mm | **3.94 mT** |

Single-pair alternative if depth or mass is tighter: **4×1 mm N35, one pair, rest 2.15 mm**
— 0.497 → 0.156 N, ratio 3.18, 4.15 mm deep, 188 mg. It has **no** restoring torque about
the rail axis, so the guide pins carry all the rocking moment; the two-pair layout is worth
the 254 mg.

**Two pairs vs one.** Two plungers at the rail ends self-level the rail: rock it by θ and
the compressed plunger pushes harder while the extended one pushes less, giving a
torsional restoring stiffness of ½·k·s² = 4.66 N·mm/rad. One pair gives none, and the pins
take the whole moment as friction — which is exactly the stiction that would ruin the
force window.

**Stray field.** 3.94 mT at 8 mm is the single-puck on-axis bound and is conservative: the
two pucks of a pair are anti-parallel, so at 8 mm they behave as a quadrupole and cancel
faster than this. For scale, `MAGFLOAT_MAGNETS.md` accepted 3.30 mT at 10 mm for the
nozzle float. Worth a bench check against the bone sensor, not a blocker.

## 5.3 Plunger vs Ti sheet

| | **mag-plunger (5×2.5×1 N35, 2 pairs)** | Ti sheet candidate B |
|---|---|---|
| Force window | 0.493 → 0.180 N over ±0.75 mm | 0.067 → 0.543 N over 0.25–2.0 mm |
| Usable dynamic travel | ±0.75 mm (cam takes the rest) | 2.0 mm continuous |
| Travel tolerance | **cam preset absorbs per-ear variation** — the magnets never see it | none: the ear's position *is* the deflection |
| Force at the far end of travel | 0.180 N, still in band | 0.543 N at 2 mm, over band |
| Fatigue / yield risk | **none — no strained part** | 450 MPa flag at 1.79 mm; yield ~2.4 mm |
| Span demanded of the concha | rail sits in the 5 mm depth, no reach needed | **13.99 mm** — gives back the entire v4 shortening |
| Part count | rail + pad + 2 pins + 4 magnets = **8** | 1 (printed with the jacket) |
| Added mass | 971 mg/side | ~110 mg (25 mm³ Ti) |
| Stray field at 8 mm | 3.94 mT | none |
| Print risk | none (all parts are simple) | 0.20 mm wall at the LPBF floor, 0.20 mm edge-roll ramp |

**Why the plunger wins on the thing that matters.** The Ti sheet passes the force target
only by being 14 mm long, which reinstates the overpressure that the v4 shortening was
added to fix — it trades a force problem for a fit problem. The plunger decouples the two:
the cam absorbs per-ear geometry, and the magnets only regulate force over a small window
they can comfortably cover. It also removes fatigue from the design entirely.

**What it costs:** 8 parts instead of 1, ~860 mg more mass per side, a stray field near the
bone sensor, and a hard constraint that dynamic travel stay within ±0.75 mm — the cam has
to be good. If ±1.5 mm turns out to be genuinely needed, the depth budget must rise to
~8 mm (8×4×1.5 N35, one pair, rest 4.95 mm).

## 5.4 Limitations

- **No tissue or friction model.** Guide-pin friction and pad hysteresis are not in these
  numbers, and stiction is the most likely way the force window degrades in practice.
  L/D = 2.5 is a rule of thumb for no-cock, not a computed result.

  > **Cross-reference — where the tissue side is now covered.**
  > `docs/TRYON_REPORT.md` carries a quasi-static rigid-body stability check
  > (`cad/iem/stability.py`) that puts the plunger forces from §5 against a real
  > ear: 0.31 N skirt reaction, a 0.5 N cable tug and 3 g inertial, resisted by
  > contact normals and Coulomb friction. Its verdict on P0023 turns almost
  > entirely on the **skin friction coefficient**, which neither document
  > measures — the design fails at μ = 0.4 and passes from μ ≈ 0.7 upward. That
  > makes a pad-on-skin friction measurement the highest-value bench test for the
  > first build; see the sensitivity grid in TRYON_REPORT. The plunger force band
  > used there (0.18–0.49 N per pad) is §5's output.
- **Rail and pin masses are estimates** from nominal dimensions (0.8 mm Ti rail, 2.0 mm
  pins, 1 mm silicone pad), not from a CAD solid.
- **Stray field is the on-axis single-puck value.** The pair cancellation and the sensor's
  actual off-axis position are not modelled; treat 3.94 mT as an upper bound.
- **The cam is unmodelled.** Its preset resolution sets how much of the ±0.75 mm window is
  left for dynamic travel, and nothing here checks that budget.
- Magnet grades assume Br = 1.20 T (N35) at room temperature; no thermal derating.

---

# Files

| | |
|---|---|
| `cad/iem/fea/_common.py` | material, SDF→mesh (voxel-hex and gmsh-tet), graded-modulus assembly, AMG solver, von Mises recovery |
| `cad/iem/fea/rve_homogenise.py` | JOB 1a — direct micro-FE on the real gyroid, §1.2 |
| `cad/iem/fea/wing_stiffness.py` | JOB 1b — wing contact FEA + jacket patch + collapse estimates, §1.3–1.6 |
| `cad/iem/fea/magnet_pressure.py` | JOB 2 — pressure budget, §2.2–2.6 |
| `cad/iem/fea/wing2_stiffness.py` | §3 — redesigned macro-gyroid wing: St Venant-Kirchhoff nonlinear solid, scalar displacement control, section bounds |
| `cad/iem/fea/wing2_run.log` | §3 — raw solver log for the 2026-08-31 run |
| `cad/iem/fea/mag_plunger.py` | §5 — mag-plunger sizing (magpylib), the chosen design |
| `cad/iem/fea/mag_plunger_results.json` | §5 — sweep, closest-miss and recommended spec |
| `cad/iem/fea/*.json` | raw solver output for every table above |

Reproduce: `cd cad/iem && .venv/bin/python fea/rve_homogenise.py && .venv/bin/python fea/wing_stiffness.py && .venv/bin/python fea/magnet_pressure.py && .venv/bin/python fea/wing2_stiffness.py`
(`rve_homogenise.py` must run first; `wing_stiffness.py` reads its JSON. `wing2_stiffness.py`
is independent and takes ~90 min — add `--fine` for the h = 0.08 mesh point. It checkpoints
`wing2_results.json` after every block. Note `wing_stiffness.py` now builds the *current*
wing envelope, so it no longer reproduces §1; see §3.6.)
