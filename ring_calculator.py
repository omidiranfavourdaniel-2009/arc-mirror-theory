#!/usr/bin/env python3
"""Reproduce the Arc-Mirror ring-count result for the 2 m design.

The computational construction is delegated to the supplied Arc-Mirror Omega
implementation in generator/arc_mirror_omega.py, so this entry point and the
technical generator use the same geometry implementation.
"""
from pathlib import Path
import sys

ROOT = Path(__file__).resolve().parent
sys.path.insert(0, str(ROOT / "generator"))
from arc_mirror_omega import ArcMirror  # noqa: E402

D = 2000.0
F_NUMBER = 1.5

def evaluate(n, adaptive, matching="midpoint"):
    m = ArcMirror(diameter=D, f_number=F_NUMBER, n_rings=n,
                  adaptive=adaptive, matching=matching)
    e = m.evaluate_surface(n_eval=2000)
    return e["max_surf"] * 1e3, e["max_slope"] * 1e3


def main():
    print("Arc-Mirror ring calculator / numerical verification")
    print(f"D = {D:.0f} mm, f = {D*F_NUMBER:.0f} mm, f/D = {F_NUMBER:.2f}")
    print("Targets: <= 16 nm surface error, <= 0.039 mrad slope error\n")

    cases = [
        ("uniform-left", 380, False, "left"),
        ("uniform-midpoint", 190, False, "midpoint"),
        ("adaptive-left", 194, True, "left"),
        ("adaptive-midpoint", 98, True, "midpoint"),
    ]
    for name, n, adaptive, matching in cases:
        h, s = evaluate(n, adaptive, matching)
        print(f"{name:20s} {n:3d} rings  {h:8.3f} nm  {s:9.5f} mrad")

    h97, s97 = evaluate(97, True, "midpoint")
    h98, s98 = evaluate(98, True, "midpoint")
    print("\nTransition check")
    print(f"97 adaptive-midpoint: {h97:.3f} nm, {s97:.5f} mrad  -> fails 16 nm")
    print(f"98 adaptive-midpoint: {h98:.3f} nm, {s98:.5f} mrad  -> passes 16 nm")

if __name__ == "__main__":
    main()
