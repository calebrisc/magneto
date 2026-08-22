// ============================================================
// Magneto M1 bench sled ("the puck")
// Holds one NVE AAT00x-10E-EVB01 breakout VERTICALLY (sensor
// end down) over a magnetic strip taped to the desk.
//
// Why vertical: field above a striped pattern rotates in the
// vertical plane along the travel axis; the AAT is an in-plane
// angle sensor, so its chip plane must stand in that plane.
// That's why NVE put the chip at the end of the board.
//
// Parts (select with `part`, or use the render script):
//   puck   — the sled body
//   wedge  — friction wedge that clamps the board in the slot
//   shims  — gap-calibration plates (0.5/1/2/3 mm, notch-coded)
//   coupon — small test block: print FIRST to verify slot fit
//
// Gap-setting procedure (M1 task #1):
//   1. Stack shims on the strip to the gap you want.
//   2. Drop board through slot until its bottom edge rests on
//      the stack. Chip faces the "F" (front) wall relief.
//   3. Push wedge down behind the board until snug.
//   4. Slide the shims out. Gap = stack height. Measure signal.
// ============================================================

part = "puck"; // "puck" | "wedge" | "shims" | "coupon"

/* ---- measured / datasheet inputs (mm) ---- */
board_w   = 9.8;   // EVB01 width  (spans travel direction X)
board_t   = 1.6;   // PCB thickness
chip_clear= 1.35;  // relief depth for the DFN chip + solder (0.8 pkg + margin)

/* ---- fit clearances — tune after printing the coupon ---- */
slot_clear_x = 0.6;  // total X clearance around board width
board_gap_y  = 0.5;  // Y clearance in board zone (board drops freely)
wedge_gap_y  = 2.5;  // Y space behind board for the wedge

/* ---- puck body ---- */
body_x = 55;  body_y = 36;  body_z = 14;  corner_r = 4;
chan_w = 16;  chan_h = 2;    // strip channel in the bottom, runs along X
snug_top = 6;                // snug slot: z = chan_h .. snug_top
well_x = 17; well_y = 12.5;  // open well above the slot (header/wire room)
hole_d = 4.2;                // zip-tie / tow holes
eps = 0.1;

slot_x = board_w + slot_clear_x;                 // 10.4
slot_y = board_t + board_gap_y + wedge_gap_y;    // 4.6
slot_y0 = -(board_t + board_gap_y)/2 - 0.0;      // front face of slot
// board zone: y in [slot_y0, slot_y0+board_t+board_gap_y]; wedge behind it

$fn = 64;

module rounded_slab(x, y, z, r) {
    linear_extrude(z)
        offset(r) offset(-r)
            square([x, y], center = true);
}

module puck() {
    difference() {
        rounded_slab(body_x, body_y, body_z, corner_r);

        // strip channel through the bottom, full length
        translate([-body_x/2 - eps, -chan_w/2, -eps])
            cube([body_x + 2*eps, chan_w, chan_h + eps]);

        // snug slot (board + wedge), connects channel to well
        translate([-slot_x/2, slot_y0, chan_h - eps])
            cube([slot_x, slot_y, snug_top - chan_h + 2*eps]);

        // chip relief groove in the front wall (chip faces here)
        translate([-2.3, slot_y0 - chip_clear, chan_h - eps])
            cube([4.6, chip_clear + eps, snug_top - chan_h + 2*eps]);

        // open well above the slot
        translate([-well_x/2, -well_y/2 + 1, snug_top])
            cube([well_x, well_y, body_z - snug_top + eps]);

        // zip-tie / tow holes, one per corner
        for (sx = [-1, 1], sy = [-1, 1])
            translate([sx*21, sy*12.5, -eps])
                cylinder(d = hole_d, h = body_z + 2*eps);

        // "F" marks the front (chip-facing) wall, engraved on top
        translate([0, -14.5, body_z - 0.6])
            linear_extrude(0.7)
                text("F", size = 6, halign = "center", valign = "center");

        // project label, engraved on top rear
        translate([0, 13, body_z - 0.6])
            linear_extrude(0.7)
                text("MAGNETO M1", size = 4.5, halign = "center", valign = "center");
    }
}

/* ---- wedge: drops in behind the board, taper self-locks ---- */
wedge_h = 24;
module wedge() {
    // flat face bears on the board's back; taper rides the rear wall
    hull() {
        cube([board_w - 0.4, 1.5, eps]);                    // thin end (bottom)
        translate([0, 0, wedge_h]) cube([board_w - 0.4, 2.9, eps]); // thick end
    }
    translate([-2.3, -1, wedge_h])   // T-handle to push/pull
        cube([board_w + 4.2, 5, 5]);
}

/* ---- gap shims: notches on the edge = thickness / 0.5 ---- */
shim_x = 22; shim_y = 14;
module shim(t, notches) {
    difference() {
        cube([shim_x, shim_y, t]);
        for (i = [0 : notches - 1])
            translate([4 + i*4.5, shim_y, -eps])
                cylinder(r = 1.5, h = t + 2*eps);
    }
}
module shims() {
    ts = [0.5, 1.0, 2.0, 3.0];
    for (i = [0 : len(ts) - 1])
        translate([0, i * (shim_y + 4), 0])
            shim(ts[i], i == 0 ? 1 : ts[i]/0.5);
}

/* ---- fit coupon: slot cross-section only, print this first ---- */
module coupon() {
    difference() {
        translate([-13, -10, 0]) cube([26, 20, 8]);
        translate([-slot_x/2, slot_y0, -eps])
            cube([slot_x, slot_y, 8 + 2*eps]);
        translate([-2.3, slot_y0 - chip_clear, -eps])
            cube([4.6, chip_clear + eps, 8 + 2*eps]);
    }
}

if (part == "puck")   puck();
if (part == "wedge")  wedge();
if (part == "shims")  shims();
if (part == "coupon") coupon();
