# Clamp pad S/M/L bands vs one universal stroke (2026-09-04)

## Verdict (2026-09-04)

1. **Sizing fixes reach, not retention.** Three pad lengths seat the cymba pad on 96/102 ears
   vs 51/102 with one universal 6 mm stroke. But the along-aim reach spans 25 mm and has **no
   correlation with ear size** (r = -0.04 vs basin radius, -0.19 vs crus distance), so a wearer
   cannot pick S/M/L by how big their ear is. Sizes would be picked by trying the three pads.
2. **A single pad does not hold, at any size, in either load case.** Cone (cable down/back) and
   sphere (any yank) both fail on ~97% of ears at mu 0.6. Pad area (60-200 mm2) and skirt
   push-out (0.31-0.15 N) barely move the margin. The binding constraint is **moment balance**: one
   interlock point in the cymba is a lever, and a pull on the nozzle rotates the body about it.
   The two-pad build (cymba + antihelix) reaches 0.57x on the same ear where one pad gives 0.00x.
3. **The 9/1 clamp is not what this harness tested.** The clamp's fixed jaw is the body edge tucked
   under the crus; that is the second reaction that balances the lever. In the cached seated poses
   (plunger-era optimiser) the jacket touches the concha floor at 0-4 points and is not under the
   crus, so the model cannot credit the jaw. Before anything else: re-seat by construction (body
   edge under the crus, 8/31 item), then rerun `size_bands.py --stability` and `crus_bands.py`.
4. **Decision this does support:** one body, pad extension in three lengths chosen by fit, stroke
   ~7-8 mm per size (5-6 gaps), sunk pocket if the S size is to fit. Crus-travel bands from
   `crus_bands.py`: S 4-9, M 10-16, L 17-22 mm (M and L stacks fit; S needs a recess).
5. **Open question the FEA cannot answer:** whether the crus jaw + one pad gives >= 1x. That is the
   first thing the Form 3 fit shell should tell us on Cale's ear.

102 real ears (P0023 excluded), cymba site only, aim [0.2, 0.97, -0.1]. `reach` = along-aim distance from the jacket face to the undercut patch under the cymba lip; it is the length the pad extension must span. Script: `cad/iem/size_bands.py`.

reach over 96 ears with an undercut: min -0.05  p10 2.75  median 7.28  p90 16.47  max 24.92 mm  (span 24.97); no undercut found on 6 ears.

## Windows

| config | band (reach, mm) | compacted | stroke | gaps | stack length | fits in compacted? | ears | seated |
|---|---|---|---|---|---|---|---|---|
| universal | all | 4.50 | 6.00 | 4 | 6.80 | NO | 102 | 51 (50%) |
| S | -0.05–8.10 | -1.05 | 10.14 | 7 | 9.80 | NO | 60 | 60 (100%) |
| M | 8.35–16.24 | 7.35 | 9.90 | 7 | 9.80 | NO | 26 | 26 (100%) |
| L | 16.69–24.92 | 15.69 | 10.24 | 7 | 9.80 | yes | 10 | 10 (100%) |

Banded: 96/102 ears seated with the right size vs 51/102 universal.

## Stability, mu = 0.6

Two load cases. **cone** = cable hanging down and back (the 8/31 matrix convention). **sphere** = every pull direction including a straight outward yank along the nozzle axis, the case a hook has to survive. Margin >= 1x passes.

| config | ears | seated | cone pass @0.5 N | cone pass @0.2 N | cone median 0.5 N | sphere pass @0.5 N | sphere pass @0.2 N | median straight pull-out capacity |
|---|---|---|---|---|---|---|---|---|
| universal | 102 | 51 | 3 (3%) | 8 (8%) | 0.00x | 1 (1%) | 3 (3%) | 0.20 N |
| S | 60 | 60 | 2 (3%) | 8 (13%) | 0.00x | 1 (2%) | 4 (7%) | 0.31 N |
| M | 26 | 26 | 1 (4%) | 1 (4%) | 0.00x | 1 (4%) | 1 (4%) | 0.21 N |
| L | 10 | 10 | 0 (0%) | 0 (0%) | 0.01x | 0 (0%) | 0 (0%) | 0.17 N |
| banded | 96 | 96 | 3 (3%) | 9 (9%) | 0.00x | 2 (2%) | 5 (5%) | 0.27 N |

## Can a wearer pick the size by ear size?

- **ear scale (basin inscribed radius)**: r = -0.04 with reach (n=96); picking S/M/L by terciles of it lands 30% of ears in their correct reach band.
- **crus helicis distance**: r = -0.19 with reach (n=37); picking S/M/L by terciles of it lands 19% of ears in their correct reach band.

Reach is a body-to-lip distance in the seated pose, so it depends on how the body sits, not only on how big the ear is. If the correlation is weak, size selection has to be by fit (try the three pads), not by a tape measure.


# Crus-pinch clamp (the 9/1 design), sized S/M/L by fit

41 ears with a detected crus (of 102); undercut found at the crus on 36. `travel` = distance from the jacket face (on the line core-centre -> patch) to the crus overhang patch. Script: `cad/iem/crus_bands.py`.

travel: min 4.04  p10 9.68  median 14.62  p90 20.14  max 22.52 mm (span 18.48)

| size | band (travel, mm) | compacted | stroke | gaps | stack length | fits (recess 0 / 3 mm) | ears | pass @0.5 N | pass @0.2 N | median margin 0.5 N |
|---|---|---|---|---|---|---|---|---|---|---|
| S | 4.04–9.00 | 3.04 | 6.95 | 5 | 7.80 | no / no | 4 | 0 (0%) | 0 (0%) | 0.03x |
| M | 10.37–16.37 | 9.37 | 8.00 | 6 | 8.80 | yes / yes | 22 | 1 (5%) | 2 (9%) | 0.00x |
| L | 17.42–22.52 | 16.42 | 7.10 | 5 | 7.80 | yes / yes | 10 | 0 (0%) | 0 (0%) | 0.12x |

All 36 ears, pad seated by construction: pass 1 @0.5 N, 2 @0.2 N; median straight pull-out capacity 0.30 N vs demand ~1.05 N (0.31 skirt + 0.5 tug + 0.24 inertial); median interlocking contacts 14.

Universal-vs-banded is not the question here: with the pad placed on the overhang by construction every ear is 'seated', so the pass rate is the ceiling for a single crus pad. Sizing only decides whether a physical stack can deliver the pad to that overhang (the fits column).


## What the pad has to be (sweep over pad cap and skirt push-out)

Pad placed on the cymba overhang, banded sizes (every ear seated), mu 0.6, 0.5 N cable tug. Pad normal-force cap = 4 kPa pain-onset pressure x pad area; skirt = sustained push-out of the seal.

| pad area | skirt push-out | ears | cone pass | sphere pass | cone median | sphere median |
|---|---|---|---|---|---|---|
| 60 mm² (cap 0.24 N) | 0.31 N | 96 | 3 (3%) | 2 (2%) | 0.00x | 0.00x |
| 120 mm² (cap 0.48 N) | 0.31 N | 96 | 4 (4%) | 2 (2%) | 0.01x | 0.00x |
| 120 mm² (cap 0.48 N) | 0.15 N | 96 | 2 (2%) | 2 (2%) | 0.02x | 0.01x |
| 200 mm² (cap 0.80 N) | 0.15 N | 96 | 2 (2%) | 2 (2%) | 0.02x | 0.01x |