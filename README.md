# Arc-Mirror Theory – Computational Implementation

Python implementation of **Arc-Mirror Theory**: a constructive framework for approximating parabolic mirrors using piecewise circular arcs with controlled discontinuities.

This repository contains the computational tools that realise the mathematical theory, including surface construction, error analysis, optical validation, and manufacturing output (G-code and STL).

## Overview

Arc-Mirror Theory approximates a parabolic (or general C²) mirror profile with a sequence of circular arcs. Each arc is constrained to have its centre on the optical axis and is matched to the target surface in both position and slope at its left endpoint. The resulting controlled discontinuities vanish as the number of segments increases.

The implementation supports:

- Exact tangential-arc construction (Theorem 1)
- Curvature-mismatch evaluation (Theorem 2)
- Controlled discontinuity analysis
- Uniform and adaptive segmentation
- Surface and slope error metrics
- Monte-Carlo ray tracing with finite solar disk
- Tool-radius compensated G-code generation (native G02/G03)
- Watertight STL shell export

## Key Features

| Feature                        | Description                                              |
|--------------------------------|----------------------------------------------------------|
| Tangential Arc Construction    | Axis-constrained circular arcs matched in position & slope |
| Adaptive Segmentation          | Curvature-informed node distribution                     |
| Error Analysis                 | Surface, slope, and reflected-angle metrics              |
| Ray Tracing                    | Monte-Carlo simulation with solar-disk sampling          |
| Manufacturing Output           | Tool-compensated G-code and watertight STL               |
| Validation Suite               | Self-consistency tests against theoretical predictions   |

## Repository Contents

- Core `ArcMirror` class implementing the full construction pipeline
- Adaptive and uniform segmentation strategies
- Optical performance evaluation and ray tracing
- CNC toolpath generation with correct tool-radius compensation
- STL export with topology certification
- Example scripts and numerical validation cases

## Theoretical Background

The mathematical foundation is developed in the companion paper:

> Omidiran Favour Daniel, *Arc-Mirror Theory: Controlled-Discontinuity Circular Arc Approximation of Parabolic Mirrors with Error Estimates, Adaptive Segmentation, and Manufacturing Implementation*, 2026.

A complete formula reference and construction manual is also available.

## Notes

- All reported optical figures are geometric simulations of the commanded surface.
- No machined part or interferogram is included in this release.
- The adaptive segmentation is a practical curvature-informed heuristic.

## Citation

If you use this implementation, please cite the associated paper and this repository.

## License

All rights reserved.  
This code is provided for research and evaluation purposes. Redistribution or commercial use requires explicit permission from the author.
