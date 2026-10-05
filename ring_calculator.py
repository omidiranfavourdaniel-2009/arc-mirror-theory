#!/usr/bin/env python3
"""
Arc-Mirror Theory - Ring Calculator
===================================

Finds the smallest number of rings (circular arcs) that keeps a parabolic
Arc-Mirror within a target surface error and a target slope error, for each
of the four constructions in the family:

    uniform-left      uniform-midpoint
    adaptive-left     adaptive-midpoint

Method (no fixed multipliers anywhere)
--------------------------------------
1. Leading-order estimate from the paper's formulas. With R = D/2 and
   lambda(r) = r^2 / (8 f^3), the worst-segment errors are

       uniform,  left      eps_h = R^4 / (16 f^3 n^2)    eps_s = R^3 / (8 f^3 n)
       uniform,  midpoint  eps_h = R^4 / (64 f^3 n^2)    eps_s = R^3 / (16 f^3 n)
       adaptive, left      eps_h = R^4 / (64 f^3 n^2)    eps_s = R^3 / (16 f^3 n)
       adaptive, midpoint  eps_h = R^4 / (256 f^3 n^2)   eps_s = R^3 / (32 f^3 n)

   (Theorem 4 and Section 5.7; adaptive rows equidistribute sqrt(lambda).)
   Inverting each formula for both targets gives a starting estimate.

2. Direct numerical verification. For a candidate n the partition is built
   using the same uniform or regularised sqrt(delta-kappa) construction used by
   the generator, every ring's axis-constrained arc is constructed from
   Theorem 1, and the maximum sampled surface and slope deviations are evaluated
   against the parabola. The smallest n that meets BOTH targets is returned.

All figures describe the commanded geometry. They are simulations, not
measurements of a machined part.
"""

import math
import sys

import numpy as np

# Benchmark from the paper: D = 200 mm, f/1.5, 120 uniform left-endpoint rings.
BENCHMARK_SURFACE_NM = 16.0
BENCHMARK_SLOPE_MRAD = 0.039

ARCMIN_PER_MRAD = (180.0 / math.pi) * 60.0 / 1000.0

FAMILIES = ("left", "midpoint")
SEGMENTATIONS = ("uniform", "adaptive")


# ----------------------------------------------------------------------
# Partition and numerical error evaluation
# ----------------------------------------------------------------------

def build_nodes(R, f, n, adaptive):
    """Radial node positions 0 = r_0 < ... < r_n = R."""
    if not adaptive:
        return np.linspace(0.0, R, n + 1)

    # Regularised sqrt(delta-kappa) equidistribution (Theorem 6, Section 5.3).
    s_max = 2.0 * ((4.0 * f ** 2 + R ** 2) ** 0.25 - (4.0 * f ** 2) ** 0.25)
    eps = s_max / (2.0 * math.sqrt(2.0) * f ** 1.5 * n)

    r = np.linspace(0.0, R, 40_001)
    dk = r ** 2 / (4.0 * f ** 2 + r ** 2) ** 1.5
    rho = np.sqrt(dk + eps)
    cum = np.zeros_like(r)
    cum[1:] = np.cumsum(0.5 * (rho[:-1] + rho[1:]) * np.diff(r))
    cum /= cum[-1]
    nodes = np.interp(np.linspace(0.0, 1.0, n + 1), cum, r)
    nodes[0] = 0.0
    nodes[-1] = R
    return nodes


def exact_errors(R, f, nodes, matching):
    """
    Worst surface error (nm) and worst slope error (mrad) of the arc
    construction on the given nodes. Returns None when an arc cannot cover
    its ring (the full-coverage condition of Theorem 3 fails).
    """
    r0 = nodes[:-1]
    r1 = nodes[1:]
    r_match = r0 if matching == "left" else 0.5 * (r0 + r1)

    # Theorem 1: axis-constrained arc tangent to the parabola at r_match.
    # (At r_match = 0 this is the osculating limit R = 2f, c_y = 2f.)
    arc_R = np.sqrt(r_match ** 2 + (2.0 * f) ** 2)
    arc_cy = r_match ** 2 / (4.0 * f) + 2.0 * f

    if np.any(r1 >= arc_R):
        return None

    t = np.linspace(0.0, 1.0, 9)
    r = r0[:, None] + (r1 - r0)[:, None] * t[None, :]
    root = np.sqrt(arc_R[:, None] ** 2 - r ** 2)

    y_arc = arc_cy[:, None] - root
    m_arc = r / root
    y_true = r ** 2 / (4.0 * f)
    m_true = r / (2.0 * f)

    surface_nm = float(np.max(np.abs(y_arc - y_true))) * 1e6
    slope_mrad = float(np.max(np.abs(m_arc - m_true))) * 1e3
    return surface_nm, slope_mrad


# ----------------------------------------------------------------------
# Leading-order estimates from the paper's formulas
# ----------------------------------------------------------------------

# (surface constant, slope constant): eps_h = a * R^4/(f^3 n^2), eps_s = b * R^3/(f^3 n)
_CONSTANTS = {
    ("uniform", "left"): (1.0 / 16.0, 1.0 / 8.0),
    ("uniform", "midpoint"): (1.0 / 64.0, 1.0 / 16.0),
    ("adaptive", "left"): (1.0 / 64.0, 1.0 / 16.0),
    ("adaptive", "midpoint"): (1.0 / 256.0, 1.0 / 32.0),
}


def formula_estimate(R, f, seg, matching, target_nm, target_mrad):
    """Leading-order n for the surface target, the slope target, and both."""
    a, b = _CONSTANTS[(seg, matching)]
    eps_h = target_nm * 1e-6      # mm
    eps_s = target_mrad * 1e-3    # rad
    n_h = math.sqrt(a * R ** 4 / (f ** 3 * eps_h))
    n_s = b * R ** 3 / (f ** 3 * eps_s)
    return n_h, n_s, max(n_h, n_s)


# ----------------------------------------------------------------------
# Search for the smallest n meeting all targets
# ----------------------------------------------------------------------

def evaluate(R, f, n, seg, matching):
    nodes = build_nodes(R, f, n, seg == "adaptive")
    return exact_errors(R, f, nodes, matching)


def meets(result, target_nm, target_mrad):
    if result is None:
        return False
    surface_nm, slope_mrad = result
    return surface_nm <= target_nm and slope_mrad <= target_mrad


def smallest_n(R, f, seg, matching, target_nm, target_mrad, n_limit=2_000_000):
    """Smallest n whose exact errors meet every target, or None."""
    est = formula_estimate(R, f, seg, matching, target_nm, target_mrad)[2]

    def ok(n):
        return meets(evaluate(R, f, n, seg, matching), target_nm, target_mrad)

    hi = max(2, int(math.ceil(est * 1.25)))
    while not ok(hi):
        hi *= 2
        if hi > n_limit:
            return None

    lo = 1
    while lo < hi:                      # smallest passing n, assuming monotone
        mid = (lo + hi) // 2
        if ok(mid):
            hi = mid
        else:
            lo = mid + 1
    n = lo

    # Guard against non-monotone behaviour near the answer.
    while n > 1 and ok(n - 1):
        n -= 1
    k = 1
    while k <= 12:
        if not ok(n + k):
            n = n + k + 1
            k = 1
            continue
        k += 1
    return n


# ----------------------------------------------------------------------
# Report
# ----------------------------------------------------------------------

def analyse(diameter, focal_length, target_nm, target_mrad):
    R = diameter / 2.0
    f = focal_length
    rows = {}
    for seg in SEGMENTATIONS:
        for matching in FAMILIES:
            n_h, n_s, n_est = formula_estimate(R, f, seg, matching,
                                               target_nm, target_mrad)
            n = smallest_n(R, f, seg, matching, target_nm, target_mrad)
            result = evaluate(R, f, n, seg, matching) if n else None
            rows[(seg, matching)] = {
                "n_surface_est": n_h, "n_slope_est": n_s, "n_est": n_est,
                "n": n, "result": result,
            }
    return rows


def print_report(diameter, focal_length, matching, target_nm, target_mrad,
                 rows):
    line = "=" * 58
    thin = "-" * 58
    print()
    print(line)
    print("  Arc-Mirror Theory - Ring Calculator")
    print(line)
    print(f"  Mirror   : D = {diameter:g} mm, f = {focal_length:g} mm")
    print(f"             (f/{focal_length / diameter:.2f})")
    print(f"  Targets  : surface <= {target_nm:g} nm")
    print(f"             slope   <= {target_mrad:g} mrad")
    print(f"  Matching : {matching}")
    print(thin)
    print("  Smallest ring count meeting both targets")
    print("  (verified by direct numerical evaluation of the arcs)")
    print()
    print("  construction        est.  rings  surf(nm) slope(mrad)")
    for seg in SEGMENTATIONS:
        for fam in FAMILIES:
            row = rows[(seg, fam)]
            name = f"{seg}-{fam}"
            mark = "*" if fam == matching else " "
            if row["n"] is None:
                print(f" {mark}{name:<18} {row['n_est']:5.0f}  not reachable")
                continue
            s_nm, s_mrad = row["result"]
            print(f" {mark}{name:<18} {row['n_est']:5.0f} {row['n']:6d} "
                  f"{s_nm:8.2f} {s_mrad:10.4f}")
    print("  * = selected family, est. = leading-order formula")
    print(thin)

    uni = rows[("uniform", matching)]
    ada = rows[("adaptive", matching)]
    if uni["n"] and ada["n"]:
        if uni["n_surface_est"] >= uni["n_slope_est"]:
            binding = "surface error"
        else:
            binding = "slope error"
        saving = uni["n"] - ada["n"]
        pct = 100.0 * saving / uni["n"]
        print("  Recommended number of rings:")
        print(f"    Uniform  : {uni['n']} rings")
        print(f"    Adaptive : {ada['n']} rings")
        print(f"  Adaptive uses {saving} fewer rings ({pct:.0f}% fewer)")
        print(f"  Binding target: {binding}")
        s_nm, s_mrad = ada["result"]
        print()
        print(f"  Adaptive: {s_nm:.2f} nm, {s_mrad:.4f} mrad,")
        print(f"  {2.0 * s_mrad * ARCMIN_PER_MRAD:.3f} arcmin reflected-angle proxy")
    else:
        print("  No ring count within the search limit meets")
        print("  every target.")
    print(line)
    print("  Simulated geometry; not a measurement.")
    print(line)
    print()


# ----------------------------------------------------------------------
# Input
# ----------------------------------------------------------------------

def _ask_float(prompt, default=None):
    while True:
        raw = input(prompt).strip()
        if raw == "":
            if default is not None:
                return default
            print("  Please enter a number.")
            continue
        try:
            value = float(raw)
        except ValueError:
            print("  Please enter a number.")
            continue
        if value <= 0:
            print("  Value must be positive.")
            continue
        return value


def _ask_choice(prompt, options, default):
    while True:
        raw = input(prompt).strip()
        if raw == "":
            return default
        if raw in options:
            return raw
        print(f"  Please enter one of: {', '.join(options)}")


def main():
    print()
    print("=" * 58)
    print("  Arc-Mirror Theory - Ring Calculator")
    print("=" * 58)
    print()
    diameter = _ask_float("Enter mirror diameter (mm): ")
    focal_length = _ask_float(
        f"Enter focal length (mm) [Enter = {1.5 * diameter:g}]: ",
        default=1.5 * diameter)

    print()
    print("Matching family:")
    print("  1. left-endpoint   (canonical)")
    print("  2. midpoint")
    choice = _ask_choice("Enter choice [1/2, Enter = 2]: ", ("1", "2"), "2")
    matching = "left" if choice == "1" else "midpoint"

    print()
    print("Targets:")
    print(f"  1. Benchmark: {BENCHMARK_SURFACE_NM:g} nm surface, "
          f"{BENCHMARK_SLOPE_MRAD:g} mrad slope  [recommended]")
    print("  2. Custom surface and slope targets")
    choice = _ask_choice("Enter choice [1/2, Enter = 1]: ", ("1", "2"), "1")
    if choice == "1":
        target_nm = BENCHMARK_SURFACE_NM
        target_mrad = BENCHMARK_SLOPE_MRAD
    else:
        target_nm = _ask_float("  Maximum surface error (nm): ")
        target_mrad = _ask_float("  Maximum slope error (mrad): ")

    print()
    print("Calculating...")
    rows = analyse(diameter, focal_length, target_nm, target_mrad)
    print_report(diameter, focal_length, matching, target_nm, target_mrad,
                 rows)


if __name__ == "__main__":
    try:
        main()
    except (KeyboardInterrupt, EOFError):
        print()
        sys.exit(1)
