// ============================================================
// Magneto shell — 1978 Corvette C3 (Silver Anniversary) study
// v0 shape study: ergonomic mouse body wearing C3 design language.
//
// True C3: 4483 x 1753 x 1219 mm, wheelbase 2489.
// At mouse length 120 mm, scale = 1:37.36 →
//   true-scale width 46.9 (mouse needs ~66: +40% "widebody"),
//   true-scale height 32.6 (kept close: hump 37),
//   wheelbase 66.6 → axles at x≈27 / x≈93 (arch positions kept true).
//
// Construction: pairwise-hulled ellipsoid stations, floor-clipped.
// Front = car nose = button end. X: 0 (nose) → 120 (tail).
// ============================================================

part = "mouse"; // "mouse" (solid, printable-ish) | preview handled by colors

/* BODY station tables: [x, half_width, height] — plan widths unchanged
   from v0 (the top-down outline stays); heights scaled up to car-slab
   proportions so the body reads as a body, not a dome skirt. */
stations = [
    [  2, 15, 17],   // nose tip
    [ 12, 24, 21],   // over pop-up headlight zone
    [ 27, 30, 26],   // front fender peak (front axle)
    [ 45, 28, 28],   // hood / cowl — coke-bottle begins
    [ 60, 27, 29],   // waist / door tops
    [ 70, 28, 29.5], // deck line under the cabin
    [ 78, 30, 30],
    [ 86, 32, 30],
    [ 93, 33, 29],   // rear haunch (rear axle) — widest point
    [108, 28, 26],   // rear deck
    [118, 22, 18],   // Kamm tail
];

/* CABIN (greenhouse): narrower than the body, sits inboard on the deck —
   the step where it meets the body is the segregation line. Peak stays 37. */
cabin_stations = [
    [ 56, 16, 30.5],  // windshield base
    [ 64, 18, 34.5],
    [ 72, 19.5, 36.5],
    [ 78, 20, 37],    // roof peak = palm hump (unchanged)
    [ 86, 19, 35.5],
    [ 95, 17, 32],    // fastback glass
    [104, 14, 28.5],  // glass meets rear deck
];

floor_clip = 400;
$fn = 96;

module station(s) {
    translate([s[0], 0, 0])
        scale([6, s[1], s[2]])
            sphere(r = 1, $fn = 72);
}

module body_loft() {
    for (i = [0 : len(stations) - 2])
        hull() { station(stations[i]); station(stations[i+1]); }
}

module cabin_loft() {
    for (i = [0 : len(cabin_stations) - 2])
        hull() { station(cabin_stations[i]); station(cabin_stations[i+1]); }
}

module shell_solid() {
    difference() {
        intersection() {
            union() { body_loft(); cabin_loft(); }
            translate([-10, -floor_clip/2, 0])
                cube([160, floor_clip, 100]);   // flat floor
        }

        // door-panel scallop = thumb / ring-finger groove (both sides)
        for (sy = [-1, 1])
            translate([62, sy * 33, 14])
                scale([30, 8, 9]) sphere(r = 1, $fn = 72);

        // wheel-arch cutouts in the skirt (front + rear axle, both sides)
        for (ax = [27, 93], sy = [-1, 1])
            translate([ax, sy * 34, -7])
                rotate([90, 0, 0])
                    cylinder(r = 13, h = 12, center = true, $fn = 72);

        // circular click-zone rings, engraved (inductive zones beneath)
        for (sy = [-1, 1])
            translate([30, sy * 13, 22])
                ring_engrave(r = 9);

        // four round taillights, engraved into the Kamm tail
        for (sy = [-17, -8, 8, 17])
            translate([121, sy, 9])
                rotate([0, 90, 0])
                    cylinder(r = 3.2, h = 6, center = true, $fn = 48);

        // scroll slot at the cowl / windshield base
        translate([51, 0, 30]) cube([16, 4.5, 20], center = true);
    }
}

/* engraved circle outline that follows the hood surface (deep thin tube) */
module ring_engrave(r) {
    rotate_extrude($fn = 72)
        translate([r, 0]) square([0.8, 60], center = false);
}

/* ---- Silver Anniversary two-tone via z-split, preview colors ---- */
beltline = 14;   // silver over charcoal split height
module mouse() {
    // tinted-glass cabin — visually segregated from the body
    color(c = [0.16, 0.17, 0.19]) render()
        intersection() { shell_solid(); cabin_loft(); }
    // silver upper body (cabin excluded so the glass break line reads)
    color(c = [0.80, 0.82, 0.85]) render()
        difference() {
            intersection() { shell_solid(); slab(beltline + 0.5, 60); }
            cabin_loft();
        }
    color(c = [0.70, 0.16, 0.20]) render()
        intersection() { shell_solid(); slab(beltline - 0.5, beltline + 0.5); } // red pinstripe
    color(c = [0.23, 0.24, 0.26]) render()
        intersection() { shell_solid(); slab(-1, beltline - 0.5); }
    // scroll wheel (cosmetic, for the study)
    color(c = [0.35, 0.35, 0.37])
        translate([51, 0, 24])
            rotate([90, 0, 0])
                cylinder(r = 7, h = 3.6, center = true, $fn = 64);
}
module slab(z0, z1) { translate([-20, -60, z0]) cube([170, 120, z1 - z0]); }

if (part == "mouse") mouse();
