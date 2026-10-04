#!/usr/bin/env python3
"""Reproduce the Arc-Mirror 2 m Monte-Carlo focal-spot result."""
from pathlib import Path
import sys

ROOT = Path(__file__).resolve().parent
sys.path.insert(0, str(ROOT / "generator"))
from arc_mirror_omega import ArcMirror  # noqa: E402

m = ArcMirror(diameter=2000.0, f_number=1.5, n_rings=98,
              adaptive=True, matching="midpoint")
spot = m.ray_trace(n_rays=5000, seed=0)

print("Arc-Mirror Monte-Carlo ray trace")
print("Design: D=2000 mm, f=3000 mm, 98 adaptive-midpoint rings")
print(f"Rays: {spot['n_hits']}")
print(f"RMS radius: {spot['rms_radius']:.3f} mm")
print(f"95% focal-spot diameter: {spot['spot_diameter_95']:.3f} mm")
print(f"Maximum radius: {spot['max_radius']:.3f} mm")
print("Status: simulation only; no physical focal-spot measurement.")
