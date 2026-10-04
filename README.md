# Arc-Mirror

**A manufacturing-first way to cut a solar dish on an ordinary CNC lathe.**

The finishing pass is a short list of circular interpolations (G02/G03). The machine already knows how to cut a circle; Arc-Mirror constructs a rotationally symmetric reflector from circular arcs whose centres are locked to the optical axis.

This repository is the computational side of that idea: construct the arcs, check the commanded surface, generate tool-compensated G-code, and export a watertight STL shell. Optical numbers here are **geometric/computational results for the commanded surface**. No machined part is claimed.

Author: **Omidiran Favour Daniel**  
Federal University of Petroleum Resources Effurun (FUPRE), Delta State, Nigeria  
`omidiranfavourdaniel@gmail.com`

---

## Why this exists

Precision aspheres can require specialist manufacturing processes. Arc-Mirror explores a different manufacturing route: use the circular interpolation already available on a conventional CNC lathe, while controlling the approximation error between the circular arcs and the target parabola.

Two independent choices define the segmentation:

| Segmentation | Matching point |
|---|---|
| uniform | left endpoint |
| uniform | segment midpoint |
| adaptive | left endpoint |
| adaptive | segment midpoint |

The adaptive midpoint construction is the selected 2 m case in the current technical work.

---

## Current 2 m design case

The repository supports the **2 m aperture, f = 3 m (f/D = 1.5)** solar-dish case used in the IAS Global Energy Hackathon 2026 technical work.

| Quantity | Current result/status |
|---|---|
| Aperture | 2 m |
| Focal length | 3 m |
| Segmentation | adaptive midpoint |
| Number of rings | 98 |
| Maximum commanded surface error | 15.94 nm |
| Maximum commanded slope error | 0.01215 mrad |
| Aperture area | 3.142 m² |
| Estimated receiver power | 2.1–2.6 kW* |
| 95% focal-spot diameter | ~28 mm† |
| Physical 2 m dish | not yet built/measured |

\* Estimate from 900 W/m² direct irradiance and a stated 75–92% optical-efficiency assumption; it is not measured receiver output.

† Monte Carlo ray-trace result for the commanded geometry; it is a simulation result, not a measured focal spot.

The 98-ring result is a **numerically verified commanded geometry**, not a claim about the accuracy of a finished machined optic.

---

## Quick start

Create a Python environment, install the dependencies, then run the supplied generator/validation entry point:

```bash
pip install -r requirements.txt
python arc_mirror_generator.py
```

Run the repository validation:

```bash
python run_validation.py
```

Other entry points:

```bash
python ring_calculator.py
python ray_trace.py
python thermal_expansion_worstcase.py
```

For the interactive browser demonstration, open `index.html` if it is included in your checkout.

---

## What is in the repository

- `arc_mirror_generator.py` — Arc-Mirror construction, validation, G-code generation, and STL export.
- `ring_calculator.py` — ring-count/error calculation entry point.
- `ray_trace.py` — optical validation/reference entry point.
- `thermal_expansion_worstcase.py` — first-order thermal-expansion sensitivity entry point.
- `index.html` — interactive dashboard, when included.
- `requirements.txt` — Python dependencies.
- `LICENSE.txt` — software license.

---

## Manufacturing concept

The finishing geometry is expressed as circular interpolation blocks (G02/G03). For a concave optical surface, the supplied generator applies tool-radius compensation using:

```text
R_path = R - R_tool
```

The generated toolpath is a computational output and still requires **CAM verification, machine-specific confirmation, fixturing, and physical metrology** before machining.

The repository does not claim that the generated G-code is ready to run on an arbitrary CNC machine without those checks.

---

## What the script will not do

- It will not invent a measured form error. The nanometre-scale values describe the commanded mathematical geometry; machining accuracy, thermal drift, fixturing, tool condition, and shell deformation still have to be measured.
- It will not claim that the 2 m dish has been manufactured. The physical demonstrator remains a validation step.
- It will not turn the 2.1–2.6 kW estimate into measured solar-thermal output. Receiver power, steam production, focal spot, and sterilisation performance require physical testing.
- It will not quote an unsupported market price or percentage cost reduction. Actual manufacturing cost should be established from shop quotations and the first controlled build.

---

## Reproducibility and evidence status

| Result | Status |
|---|---|
| Arc construction and error analysis | computationally implemented |
| 98-ring adaptive-midpoint geometry | numerically verified |
| 15.94 nm maximum surface error | numerically verified |
| 0.01215 mrad maximum slope error | numerically verified |
| ~28 mm 95% focal spot | Monte Carlo simulation |
| 2.1–2.6 kW receiver power | estimate from stated assumptions |
| Physical 2 m dish | not yet built/measured |
| Measured receiver power and steam output | not yet measured |

The repository is intended to make the computational claims inspectable and reproducible; it is not a substitute for machining and experimental validation.

---

## Next hardware step

The next gate is physical validation: verify the toolpath in CAM, confirm machine swing/stock/fixturing, machine a controlled demonstrator, measure the reflector form and focal spot, and compare measured receiver performance with the computational estimate.

---

## License

Copyright © 2026 Omidiran Favour Daniel. See `LICENSE.txt`.
