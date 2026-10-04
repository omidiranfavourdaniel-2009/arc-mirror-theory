#!/usr/bin/env python3
"""Repository-level checks against the audited 2 m design."""
from pathlib import Path
import sys

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "generator"))
from arc_mirror_omega import ArcMirror  # noqa: E402


def eval_case(n, adaptive=True):
    m = ArcMirror(diameter=2000.0, f_number=1.5, n_rings=n,
                  adaptive=adaptive, matching="midpoint")
    e = m.evaluate_surface(n_eval=2000)
    return e["max_surf"] * 1e3, e["max_slope"] * 1e3


def main():
    h97, s97 = eval_case(97)
    h98, s98 = eval_case(98)
    hu190, su190 = eval_case(190, adaptive=False)

    print("Arc-Mirror repository validation")
    print(f"97 adaptive-midpoint : {h97:.3f} nm, {s97:.5f} mrad")
    print(f"98 adaptive-midpoint : {h98:.3f} nm, {s98:.5f} mrad")
    print(f"190 uniform-midpoint: {hu190:.3f} nm, {su190:.5f} mrad")

    assert h97 > 16.0
    assert h98 <= 16.0
    assert s98 <= 0.039
    assert abs(h98 - 15.9413) < 0.02
    assert abs(s98 - 0.0121501) < 0.00002
    print("PASS: audited 2 m headline result reproduced.")

if __name__ == "__main__":
    main()
