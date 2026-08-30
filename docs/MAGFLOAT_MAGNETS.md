# Ring Magnet Axial Repulsion — Force vs Gap

Sizing study for a coaxial, like-poles-facing NdFeB ring magnet pair used as a
contactless preload/return spring (working travel 1.5 mm, target preload 0.20–0.40 N,
allowable range 0.15–0.5 N across travel).

## Method

- Library: `magpylib` 5.0.1 + `magpylib-force` 0.3.1 (Maxwell-stress-tensor /
  meshed-dipole force solver — direct force computation, not a manual dU/dz
  virtual-work finite difference).
- Each ring modeled as a `CylinderSegment` (φ = 0–360°, i.e. a full annulus),
  axial magnetization `M = Br/μ0` with **Br = 1.45 T** (N52), giving
  M ≈ 1.1539×10⁶ A/m.
- Source ring magnetized +z, target ring magnetized −z, so the two faces across
  the gap are both effective North poles → axial repulsion. Force reported is
  `getFT(source, target, anchor=target_position)[0][2]` (axial force component
  on the target, N).
- Each ring meshed into 200 dipole cells (`obj.meshing = 200`) for the target
  side of the force calc. Mesh-convergence check on 7×3×1.5 at gap 4.6 mm:
  50/100/200/400/800 cells → 0.390/0.399/0.394/0.395/0.395 N (converged to
  <1% by 200 cells — used throughout).
- Mass = π/4·(OD²−ID²)·t · density, density = 7.5 g/cc (N52).

## Force (N) vs face-to-face gap

| Ring OD×ID×t (mm) | 0.5mm | 0.75mm | 1.0mm | 1.25mm | 1.5mm | 2.0mm | 2.5mm | 3.0mm | Mass (g) |
|---|---|---|---|---|---|---|---|---|---|
| 6×3×1   | 2.919 | 2.036 | 1.476 | 1.105 | 0.850 | 0.540 | 0.370 | 0.268 | 0.159 |
| 6×3×1.5 | 4.094 | 2.955 | 2.212 | 1.705 | 1.348 | 0.898 | 0.638 | 0.474 | 0.239 |
| 6×3×2   | 4.900 | 3.639 | 2.794 | 2.204 | 1.780 | 1.227 | 0.895 | 0.678 | 0.318 |
| 7×3×1.5 | 5.625 | 4.269 | 3.328 | 2.652 | 2.153 | 1.490 | 1.085 | 0.821 | 0.353 |
| 8×4×1.5 | 6.827 | 5.161 | 3.999 | 3.161 | 2.543 | 1.724 | 1.232 | 0.920 | 0.424 |
| 5×2.5×1 | 2.171 | 1.450 | 1.022 | 0.751 | 0.571 | 0.360 | 0.246 | 0.177 | 0.110 |

At every gap in the 0.5–3.0 mm range requested, all six candidates produce far
more than the 0.20–0.40 N target (2–7 N at 0.5 mm, still 0.18–0.9 N at 3.0 mm).
To hit the target preload band, the rest gap must be pushed out to roughly
2.5–7 mm depending on ring size — a follow-up scan (0.5–8.0 mm, 0.25 mm step,
then 0.05 mm resolution search over rest position) was run to find, for each
candidate, the rest gap whose ±0.75 mm window (1.5 mm total travel) keeps
force inside [0.15, 0.5] N with the rest-point force inside [0.20, 0.40] N,
and minimizes the max/min force ratio over that window.

## Best rest-gap solution per candidate (window = rest ± 0.75 mm)

| Ring | Rest gap | F(rest−0.75) | F(rest) | F(rest+0.75) | Ratio max/min |
|---|---|---|---|---|---|
| 6×3×1   | 3.30 mm | 0.359 N | 0.226 N | 0.153 N | 2.35 |
| 6×3×1.5 | 4.75 mm | 0.286 N | 0.205 N | 0.150 N | 1.90 |
| 6×3×2   | 5.80 mm | 0.269 N | 0.201 N | 0.152 N | 1.76 |
| 7×3×1.5 | 6.40 mm | 0.263 N | 0.200 N | 0.155 N | **1.69** |
| 8×4×1.5 | 6.80 mm | 0.259 N | 0.203 N | 0.161 N | **1.61** |
| 5×2.5×1 | 2.50 mm | 0.448 N | 0.246 N | 0.152 N | 2.95 |

Larger rings need a larger rest gap to bleed the force down into the target
band — but at that larger gap the fixed 1.5 mm travel is a smaller fraction of
the gap, so the force curve is flatter and the ratio is lower. Smaller/thinner
rings (6×3×1, 5×2.5×1) hit the band at a much more compact 2.5–3.3 mm gap but
swing 2.3–3.0× over the travel.

## Recommendation

**Best 2 candidates (lowest force ratio over travel):**

1. **8×4×1.5 mm** — rest gap **6.8 mm**, F = 0.259 → 0.203 → 0.161 N over the
   1.5 mm travel band, ratio **1.61**, mass 0.424 g/ring.
2. **7×3×1.5 mm** — rest gap **6.4 mm**, F = 0.263 → 0.200 → 0.155 N, ratio
   **1.69**, mass 0.353 g/ring.

Both need a fairly large (6.4–6.8 mm) air gap relative to their OD, which
flattens the force-vs-travel curve at the cost of a bigger cavity than the
originally-scoped 0.5–3 mm gap range. If a compact gap is the harder
constraint, **6×3×1 mm** (rest gap 2.9–3.3 mm, ratio ~2.4–2.6, mass 0.159 g)
is the best small-package fallback, at the cost of a wider force swing.

---

## Compact-gap options (rest gap 1.0–2.5 mm)

Follow-up: same method (magpylib 5.0.1 + magpylib-force 0.3.1, Maxwell-stress
solver, 200-cell mesh, like-poles-facing coaxial rings), re-run for weaker
materials so the 0.20–0.40 N preload lands at a compact 1–2.5 mm rest gap
instead of the 6+ mm gap the N52 study required. Gaps swept 0.75–3.0 mm in
0.25 mm steps (7 ring sizes × 3 materials). Window = rest ± 0.75 mm (1.5 mm
travel), same acceptance rule as the N52 study (window ⊆ [0.15, 0.5] N, rest
∈ [0.20, 0.40] N), searched both symmetric and offset window placements.

**Materials:** sintered ferrite Y30/C8 (Br = 0.40 T, ρ ≈ 5.0 g/cc), bonded
NdFeB (Br = 0.65 T, ρ ≈ 6.0 g/cc), N35 sintered NdFeB (Br = 1.20 T, ρ = 7.5 g/cc).
**Ring sizes:** 6×3×1, 6×3×1.5, 7×3×1.5, 8×4×1.5, 5×2.5×1, 4×2×0.5, 4×2×1 mm.

### Result: only one combination clears every bound in this grid

| Material | Ring | Rest gap | F(rest−0.75) | F(rest) | F(rest+0.75) | Ratio | Mass |
|---|---|---|---|---|---|---|---|
| **Bonded NdFeB (0.65 T)** | **7×3×1.5 mm** | **2.25 mm** | 0.433 N (@1.5mm) | 0.254 N | 0.165 N (@3.0mm) | **2.62** | 0.283 g |

No ferrite ring and no N35-sintered ring in the requested size list satisfies
all four bounds simultaneously inside the 0.75–3.0 mm / 1.5 mm-travel grid:

- **Ferrite (0.40 T) — too weak.** Best candidate 8×4×1.5 at the smallest
  feasible rest gap (1.5 mm, set by the 0.75 mm gap floor + 0.75 mm half-window):
  F(0.75mm)=0.393 N, F(1.5mm)=0.194 N, F(2.25mm)=0.110 N. Misses on both ends —
  rest force 6 mN under the 0.20 N floor, far-travel force 40 mN under the
  0.15 N floor. Would need either a slightly smaller minimum gap (~0.6 mm) or
  a bigger/thicker ring than tested.
- **N35 sintered (1.20 T) — too strong.** Closest candidate 5×2.5×1 at the
  largest feasible rest gap (2.25 mm): F(1.5mm)=0.391 N, F(2.25mm)=0.202 N,
  F(3.0mm)=0.121 N. Rest point lands fine in-band, but far-travel force is
  29 mN under the 0.15 N floor — this ring/material combo decays too fast
  over 1.5 mm at gaps this small. A thinner/smaller ring than 4×2×0.5 (weaker
  still) or a shorter working travel would be needed.

### Recommendation

**Bonded NdFeB, 7×3×1.5 mm ring, rest gap ≈ 2.25 mm** is the single valid
compact-gap solution: preload 0.254 N at rest, swinging 0.433 N → 0.165 N
over the 1.5 mm travel band (ratio 2.62), 0.283 g/ring. It trades regulation
flatness for a ~3× smaller rest gap than the N52 pick (2.25 mm vs 6.8 mm).

### Stray field 10 mm behind the fixed magnet (driver-interaction proxy)

On-axis B-field, 10 mm behind the source ring's back face (opposite side from
the working gap), computed with `magpy.getB` for the source ring alone:

| Pick | Ring / Br / rest gap | Stray field @ 10 mm behind |
|---|---|---|
| Best N52 pick | 8×4×1.5 mm, Br=1.45 T, rest 6.8 mm | **8.26 mT** |
| Best ferrite pick (closest-miss, see above) | 8×4×1.5 mm, Br=0.40 T, rest 1.5 mm | **2.28 mT** |
| Compact bonded-NdFeB pick (valid solution) | 7×3×1.5 mm, Br=0.65 T, rest 2.25 mm | 3.30 mT |

The N52 pick leaks ~3.6× more stray field at 10 mm than the ferrite pick,
despite sitting much farther back (14.9 mm vs 12.25 mm from the ring), simply
because N52 is ~3.6× stronger in remanence. Ferrite is the gentlest on a
nearby sensor/driver; the compact bonded-NdFeB solution sits in between.

---

## Asymmetric pair (as-built)

The moving ring must slide over a 5 mm tube, so the pair is no longer
identical: fixed ring stays **7×3×1.5 mm** (bonded NdFeB study winner), the
moving ring's ID is fixed at **5.4 mm** (clearance over the 5 mm tube) with
OD ∈ {7, 8, 9} mm and t ∈ {1.5, 2.0} mm — 6 moving-ring geometries. Same
method as above (magpylib 5.0.1 + magpylib-force 0.3.1, 200-cell mesh,
like-poles-facing, asymmetric center-to-center separation = gap + t_fixed/2 +
t_moving/2), gaps swept 0.5–3.0 mm in 0.25 mm steps, window = rest ± 0.75 mm.

**Bonded NdFeB (Br = 0.65 T) reaches the floor fine** — 4 of 6 moving-ring
geometries produce a valid window; the thin-wall 7×5.4 mm moving rings
(0.8 mm wall) do not, they're too weak throughout the swept range. N35
sintered (Br = 1.2 T) was also run across all 6 geometries as a check: **no
combination fits** — force is 2–4× too high everywhere in 0.5–3.0 mm, so N35
is not needed here (bonded NdFeB already clears the target).

| Fixed / Moving | Rest gap | F(rest−0.75) | F(rest) | F(rest+0.75) | Ratio | Moving mass | Fixed mass |
|---|---|---|---|---|---|---|---|
| 7×3×1.5 / 7×5.4×1.5 | — | — | — | — | no valid window (too weak) | 0.140 g | 0.283 g |
| 7×3×1.5 / 7×5.4×2.0 | — | — | — | — | no valid window (too weak) | 0.187 g | 0.283 g |
| 7×3×1.5 / 8×5.4×1.5 | 1.57 mm | 0.428 N | 0.242 N | 0.150 N | 2.85 | 0.246 g | 0.283 g |
| 7×3×1.5 / 8×5.4×2.0 | 1.89 mm | 0.382 N | 0.230 N | 0.151 N | 2.53 | 0.328 g | 0.283 g |
| 7×3×1.5 / 9×5.4×1.5 | 1.25 mm | 0.405 N | 0.288 N | 0.199 N | 2.04 | 0.366 g | 0.283 g |
| **7×3×1.5 / 9×5.4×2.0** | **2.25 mm** | **0.307 N** | **0.214 N** | **0.155 N** | **1.98** | 0.489 g | 0.283 g |

(All values at exact swept grid points — 0.50/1.25/1.50/2.25/3.00 mm — no
interpolation.)

### Winner

**Fixed 7×3×1.5 mm vs moving 9×5.4×2.0 mm, bonded NdFeB, rest gap 2.25 mm** —
lowest ratio (1.98): 0.307 N → 0.214 N → 0.155 N across the 1.5 mm travel
band (window 1.5–3.0 mm), rest preload 0.214 N (inside 0.20–0.40 N).

**Lighter/more-compact alternative:** 7×3×1.5 mm vs moving **9×5.4×1.5 mm**,
rest gap **1.25 mm** (window 0.50–2.00 mm): 0.405 → 0.288 → 0.199 N, ratio
2.04 — nearly the same regulation at a much smaller rest gap and 25% less
moving-ring mass (0.366 g vs 0.489 g), useful if carrier mass or envelope
matters more than the last bit of flatness.
