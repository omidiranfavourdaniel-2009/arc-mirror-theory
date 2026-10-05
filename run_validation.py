#!/usr/bin/env python3
"""Self-contained validation checks for the audited 2 m Arc-Mirror design.

This file intentionally does not import Arc-Mirror repository modules.  It
reconstructs only the geometry and numerical surface/slope evaluation needed
to verify the published 2 m headline result.
"""

import numpy as np


# ---------------------------------------------------------------------------
# Design constants used by the repository's audited 2 m case
# ---------------------------------------------------------------------------
DIAMETER_MM = 2000.0
F_NUMBER = 1.5
FOCAL_LENGTH_MM = DIAMETER_MM * F_NUMBER
SURFACE_TARGET_NM = 16.0
SLOPE_TARGET_MRAD = 0.039


def regularisation_epsilon(f_mm: float, r_max_mm: float, n_rings: int) -> float:
    """Return the principled adaptive regularisation used by Arc-Mirror."""
    s_max = 2.0 * (
        (4.0 * f_mm**2 + r_max_mm**2) ** 0.25
        - (4.0 * f_mm**2) ** 0.25
    )
    return s_max / (2.0 * np.sqrt(2.0) * f_mm**1.5 * n_rings)


def build_radii(n_rings: int, adaptive: bool) -> np.ndarray:
    """Build the same uniform/adaptive radial partition used for validation."""
    r_max = DIAMETER_MM / 2.0

    if not adaptive:
        return np.linspace(0.0, r_max, n_rings + 1)

    # Adaptive construction: equidistribute sqrt(Delta-kappa + epsilon).
    n_fine = 10_000
    r_fine = np.linspace(1e-6, r_max, n_fine)
    f = FOCAL_LENGTH_MM
    eps = regularisation_epsilon(f, r_max, n_rings)

    delta_kappa = r_fine**2 / (4.0 * f**2 + r_fine**2) ** 1.5
    density = np.sqrt(delta_kappa + eps)

    cumulative = np.zeros_like(r_fine)
    cumulative[1:] = np.cumsum(
        0.5 * (density[:-1] + density[1:]) * np.diff(r_fine)
    )
    cumulative /= cumulative[-1]

    targets = np.linspace(0.0, 1.0, n_rings + 1)
    radii = np.interp(targets, cumulative, r_fine)
    radii[0] = 0.0
    return radii


def evaluate_case(n_rings: int, adaptive: bool = True, n_eval: int = 2000):
    """Evaluate maximum surface and slope error for the midpoint family."""
    f = FOCAL_LENGTH_MM
    radii = build_radii(n_rings, adaptive)

    # Construct each axis-centred circular arc, matched at the segment
    # midpoint.  For y=r^2/(4f), the tangent-circle centre is
    # R=sqrt(r_m^2 + (2f)^2), cy=y_m+2f.
    segments = []
    for i in range(n_rings):
        r_start = radii[i]
        r_end = radii[i + 1]
        r_match = 0.5 * (r_start + r_end)
        y_match = r_match**2 / (4.0 * f)
        radius = np.sqrt(r_match**2 + (2.0 * f) ** 2)
        cy = y_match + 2.0 * f

        if r_end >= radius:
            raise ValueError(
                f"Ring {i + 1}: endpoint {r_end:.6g} mm reaches/exceeds "
                f"arc radius {radius:.6g} mm."
            )
        segments.append((r_start, r_end, radius, cy))

    r_eval = np.linspace(0.0, DIAMETER_MM / 2.0, n_eval)
    y_eval = np.empty_like(r_eval)
    slope_eval = np.empty_like(r_eval)
    seg_idx = 0

    for i, r in enumerate(r_eval):
        while seg_idx < len(segments) - 1 and r > segments[seg_idx][1]:
            seg_idx += 1
        _, _, radius, cy = segments[seg_idx]
        root = np.sqrt(max(0.0, radius**2 - r**2))
        y_eval[i] = cy - root
        slope_eval[i] = r / np.sqrt(max(1e-12, radius**2 - r**2))

    y_true = r_eval**2 / (4.0 * f)
    slope_true = r_eval / (2.0 * f)

    # Surface error is converted from mm to nm; slope from radians to mrad.
    max_surface_nm = float(np.max(np.abs(y_eval - y_true)) * 1e6)
    max_slope_mrad = float(np.max(np.abs(slope_eval - slope_true)) * 1e3)
    return max_surface_nm, max_slope_mrad


def main() -> None:
    h97, s97 = evaluate_case(97, adaptive=True)
    h98, s98 = evaluate_case(98, adaptive=True)
    h190, s190 = evaluate_case(190, adaptive=False)

    print("Arc-Mirror repository validation (self-contained)")
    print(f"97 adaptive-midpoint : {h97:.3f} nm, {s97:.5f} mrad")
    print(f"98 adaptive-midpoint : {h98:.3f} nm, {s98:.5f} mrad")
    print(f"190 uniform-midpoint: {h190:.3f} nm, {s190:.5f} mrad")

    # Boundary case: 97 rings must fail the 16 nm surface target.
    assert h97 > SURFACE_TARGET_NM

    # Selected design: 98 adaptive-midpoint rings must pass both stated
    # headline targets and reproduce the audited numerical result.
    assert h98 <= SURFACE_TARGET_NM
    assert s98 <= SLOPE_TARGET_MRAD
    assert abs(h98 - 15.9413) < 0.02
    assert abs(s98 - 0.0121501) < 0.00002

    # Independent uniform comparison retained as a reference check.
    assert h190 <= SURFACE_TARGET_NM

    print("PASS: audited 2 m headline result reproduced without ArcMirror imports.")


if __name__ == "__main__":
    main()
