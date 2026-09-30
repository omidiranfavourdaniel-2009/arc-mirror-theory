# Arc-Mirror Theory

**A manufacturing-first way to cut a solar dish on an ordinary CNC lathe.**

The finishing pass is a short list of circular interpolations (G02/G03). The machine already knows how to cut a circle. It does not natively cut a parabola as one block.

This repository is the computational side of that idea: construct the arcs, check the surface, write tool-compensated G-code, and export a watertight STL shell. Optical numbers here are **geometric simulations of the commanded surface**. No machined part is claimed.

Author: **Omidiran Favour Daniel**  
Federal University of Petroleum Resources Effurun (FUPRE), Delta State, Nigeria  
`omidiranfavourdaniel@gmail.com`

---

## Why this exists

Precision aspheres usually want specialist machines. A two-axis lathe already interpolates lines and circles. Arc-Mirror builds a rotationally symmetric reflector whose finishing path is those circles, with the centre of each arc locked to the optical axis.

Two independent choices:

| Segmentation | Matching point |
|---|---|
| uniform | left endpoint |
| adaptive (curvature-mismatch packing) | segment midpoint (smaller leading constants) |

Default demo: **adaptive + left-endpoint**.

The journal is private. Do not upload it here.

---

## Quick start

```bash
pip install -r requirements.txt
python arc_mirror_omega.py
```

That run validates the construction, prints the four family members, and writes G-code plus a self-certified STL.

Open `index.html` in a browser for the dashboard.

Live demo: _https://arc-mirror-theory.netlify.app_

Preprint: _https://doi.org/10.5281/zenodo.22690971_

---

## Files

- `arc_mirror_omega.py` — generator, validation, G-code, STL
- `index.html` — interactive dashboard
- `sample_gcode.nc` — example G-code
- `requirements.txt` — numpy
- `LICENSE`

---

## What the script will not do

- It will not invent a measured form error. CNC positioning and thermal sag are larger than the mathematical arc mismatch on a 200 mm f/1.5 blank.
- It will not quote a market price per square metre. Use shop quotations.
- Tool compensation for a **concave** optical surface uses `R_path = R - R_tool`. A convex mandrel flips the sign.

---

## Status

| Item | State |
|---|---|
| Construction + error analysis | written (kept private) |
| Simulator, G-code, STL | this repo |
| Measured prototype | not yet |

Next hardware step: a controlled 200 mm blank on a two-axis lathe, then measure form and spot — do not treat the nanometre simulation as the manufactured optic.

---

## License

Copyright © 2026 Omidiran Favour Daniel. See `LICENSE`.
