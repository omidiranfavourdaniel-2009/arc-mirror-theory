# What Arc-Mirror is, in plain language

A parabolic dish is the shape you want for a solar concentrator. An ordinary two-axis CNC lathe does not cut a parabola as one native move. It does cut circles (G02/G03).

Arc-Mirror replaces the parabola with a short stack of circular arcs. Each arc has its centre on the optical axis. At one chosen point on the segment — the left end, or the midpoint — the arc matches the parabola in height and slope. Then the circle continues to the far end of the segment.

At that far end the height and slope are not perfect. Those jumps get smaller as you add segments. That is the whole idea: a finishing pass the lathe already knows how to run, with a form error designed to sit under ordinary machine positioning error.

## Two knobs

- **Segmentation:** even radial steps, or extra arcs where the curvature mismatch is worse (adaptive).
- **Matching point:** left endpoint (canonical, simplest proof) or midpoint (smaller leading constants: about 4× in sag, 2× in slope).

This repository defaults to **adaptive + left-endpoint**. Midpoint is the stronger shop choice at the same segment count. Both are in the generator and the live dashboard.

## What the numbers mean

Surface and slope figures in the demo are **simulations of the commanded surface**. They are not a measured blank. Heat, tool deflection, and axis error on a real lathe will dominate the mathematical arc mismatch on a 200 mm dish.

## Next hardware step

Cut one 200 mm blank on a two-axis lathe from the generated G-code. Measure form. Then talk about optical performance.

The long write-up stays off this repository on purpose.
