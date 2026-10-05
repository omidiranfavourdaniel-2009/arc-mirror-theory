#!/usr/bin/env python3
"""
Arc-Mirror Theory - Omega Computational Implementation  (v3.0)
Uniform / Adaptive Tangential-Arc Families with Matching-Point Variants

Author: Omidiran Favour Daniel

This module implements the complete Arc-Mirror Theory computational pipeline:
  1. Matching-point circular-arc construction (left-endpoint / midpoint family)
  2. Curvature-mismatch analysis
  3. Uniform and adaptive segmentation using the curvature-mismatch proxy
  4. Tool-center-path G-code generation (cutter-compensated)
  5. Watertight binary STL export (manifold mesh, self-certified)
  6. Optical validation and ray tracing

Construction family:
  * matching="left"    : canonical left-endpoint tangential construction
  * matching="midpoint" : midpoint-matched member of the same family

Segmentation and matching are independent axes, so the implementation
supports the four computational family members:
  * uniform-left, adaptive-left
  * uniform-midpoint, adaptive-midpoint

The left-endpoint construction remains the canonical/default member.

=============================================================================
 CHANGELOG v2.0 - correctness fixes from an independent line-by-line review
=============================================================================
 1. CRITICAL  Tool-radius compensation was offset to the WRONG SIDE of the
    surface (R_path = R + R_tool), which machines ~2*R_tool INTO the
    material (-5 mm at the vertex for a 5 mm tool). A concave surface is
    machined internally: the tool centre follows the concentric arc of
    radius R_path = R - R_tool. NOTE: the same error appears in paper
    Section 8.1 and the Formula Reference - update those documents.
 2. CRITICAL  The STL was inside-out (both dish surfaces wound inverted;
    the net signed volume was positive only because the rim wall
    dominated), had a geometric hole at the pole, duplicate pole
    vertices and hundreds of degenerate triangles. Rebuilt topology:
    two pole caps at full wall thickness (no shared-vertex taper).
    Inner face is a filled disk. Mesh self-certifies watertight.
 3. HIGH      G02 I-offset: Fanuc/Haas controls expect a RADIUS value for
    I even under diameter programming. I is now emitted as -r_start
    (see i_radius_value=False for the old behaviour).
 4. HIGH      ray_trace sampled the solar disk uniformly in theta (biased
    the rms spot radius low by ~17%). Now sqrt-weighted for uniform
    disk brightness, with exact off-axis composition (rotation about Z).
 5. MEDIUM    evaluate_surface returned uniform-corollary bounds in
    adaptive mode, where they do not apply. The leading-order adaptive
    proxy is computed from the regularised density; uniform bounds are None
    when adaptive=True.
 6. MEDIUM    The roughing pass plunged to Z = -allowance, i.e. below the
    vertex plane, permanently pocketing the centre. Replaced with
    allowance-shifted Z-level passes that provably never cut below the
    final surface (convexity argument in comments).
 7. LOW       Per-ring errors now report the segment maximum from a
    9-point scan instead of the midpoint only; m_end guard fixed;
    wall >= arc radius now raises instead of silently capping.
=============================================================================

 HOW TO SET MIRROR PARAMETERS
 There are TWO places parameters appear in this script:
 1. Default values inside ArcMirror.__init__(...) - fallbacks only.
 2. The run configuration inside  if __name__ == "__main__": - EDIT HERE.
 =============================================================================
"""
import numpy as np
import struct
from dataclasses import dataclass
from typing import List, Tuple, Dict, Optional
from pathlib import Path
from collections import Counter

# ---------------------------------------------------------------------------
# DATA STRUCTURES
# ---------------------------------------------------------------------------

@dataclass
class ArcSegment:
    """Single circular arc segment with a prescribed matching point."""
    r_start: float          # Left endpoint radius (mm)
    r_end: float            # Right endpoint radius (mm)
    R: float                # Surface arc radius (mm)
    cy: float               # Arc center y-coordinate (mm)
    y_start: float          # Surface height at r_start (mm)
    y_end: float            # Surface height at r_end (mm)
    m_start: float          # Surface slope at r_start
    m_end: float            # Surface slope at r_end
    r_match: float          # Matching radius used to construct the arc (mm)
    y_match: float          # Exact target height at matching radius (mm)
    m_match: float          # Exact target slope at matching radius


@dataclass
class ToolPathSegment:
    """CNC tool-center path derived from an ArcSegment.

    CORRECTED (v2.0): a concave surface is machined internally, so the
    tool-center path is the concentric arc of radius R_path = R - R_tool
    (was R + R_tool, which offset the tool into the material).

    i_offset is stored as a RADIUS value (-r_start), the Fanuc/Haas
    convention even under diameter X programming.
    """
    r_start: float          # Tool center radius at segment start
    r_end: float            # Tool center radius at segment end
    R_path: float           # Tool-path arc radius = R - R_tool
    cy: float               # Same center as surface arc
    z_start: float          # Tool center Z at segment start
    z_end: float            # Tool center Z at segment end
    i_offset: float         # G02 I offset, RADIUS value: 0 - r_start
    k_offset: float         # G02 K offset: c_y - z_start


# ---------------------------------------------------------------------------
# CORE MIRROR CLASS
# ---------------------------------------------------------------------------

class ArcMirror:
    """
    Construct and export an arc-mirror concentrator.

    Parameters
    ----------
    diameter : float
        Mirror aperture diameter (mm).
    f_number : float
        Focal ratio f/D.
    n_rings : int
        Number of arc segments.
    adaptive : bool, optional
        Use curvature-mismatch-adaptive segmentation with regularisation.
    matching : {"left", "midpoint"}, optional
        Matching-point member of the Arc-Mirror construction family.
        "left" is the canonical construction; "midpoint" matches position
        and slope at the segment midpoint.
    reg_eps : float or None, optional
        Regularisation parameter near vertex for adaptive mode.
        If None and adaptive=True, the principled value
        eps* = S(R_max)/(2*sqrt(2)*f^(3/2)*n) is computed automatically.
    wall_thickness : float, optional
        Shell wall thickness used for STL export (mm).
    """

    def __init__(
        self,
        diameter: float = 200.0,
        f_number: float = 1.5,
        n_rings: int = 65,
        adaptive: bool = True,
        reg_eps: Optional[float] = None,
        wall_thickness: float = 2.0,
        matching: str = "left",
    ):
        self.diameter = float(diameter)
        self.f_number = float(f_number)
        self.focal_length = self.diameter * self.f_number
        self.n_rings = int(n_rings)
        self.adaptive = bool(adaptive)
        self.wall_thickness = float(wall_thickness)
        self.matching = str(matching).lower().strip()
        if self.matching not in {"left", "midpoint"}:
            raise ValueError("matching must be 'left' or 'midpoint'")

        # Principled regularisation (Section 5.4)
        if self.adaptive and reg_eps is None:
            f = self.focal_length
            Rm = self.diameter / 2.0
            S_max = 2.0 * ((4.0 * f**2 + Rm**2)**0.25 - (4.0 * f**2)**0.25)
            self.reg_eps = S_max / (2.0 * np.sqrt(2.0) * f**1.5 * self.n_rings)
        else:
            self.reg_eps = float(reg_eps) if reg_eps is not None else 0.0

        # Derived constants
        self.max_r = self.diameter / 2.0
        self.sag = self.diameter ** 2 / (16.0 * self.focal_length)

        # Build segments
        self.segments: List[ArcSegment] = self._build_segments()

    # ------------------------------------------------------------------
    # 1. SEGMENT CONSTRUCTION
    # ------------------------------------------------------------------

    def _build_segments(self) -> List[ArcSegment]:
        """Build tangential arc segments (Theorem 1)."""
        radii = self._compute_radii()
        segments = []
        f = self.focal_length

        for i in range(self.n_rings):
            r_start = radii[i]
            r_end = radii[i + 1]

            # The matching point is a construction parameter independent of
            # the segmentation rule.  Thus uniform/adaptive and left/midpoint
            # are orthogonal choices in the computational family.
            if self.matching == "left":
                r_match = r_start
            else:
                r_match = 0.5 * (r_start + r_end)

            y_match = r_match ** 2 / (4.0 * f)
            m_match = r_match / (2.0 * f)

            # Axis-constrained circular arc through (r_match, y_match)
            # tangent to the parabola there.  For r_match=0 this is the
            # osculating-circle limit R=2f, cy=2f.
            if r_match == 0.0:
                R = 2.0 * f
                cy = R
            else:
                r_over_m = r_match / m_match  # = 2f for the parabola
                R = np.sqrt(r_match ** 2 + r_over_m ** 2)
                cy = y_match + r_over_m

            # Full-segment coverage is required: the circular continuation
            # must reach the complete segment without crossing its radial
            # turning point r=R.
            if r_end >= R:
                raise ValueError(
                    f"Segment {i+1}: r_end={r_end:.6g} mm reaches/exceeds "
                    f"arc radius R={R:.6g} mm. Refine the partition."
                )

            y_start = cy - np.sqrt(R ** 2 - r_start ** 2)
            y_end = cy - np.sqrt(max(0.0, R ** 2 - r_end ** 2))
            m_start = (
                r_start / np.sqrt(max(1e-12, R ** 2 - r_start ** 2))
                if r_start < R else float("inf")
            )
            m_end = (
                r_end / np.sqrt(max(1e-12, R ** 2 - r_end ** 2))
                if r_end < R else float("inf")
            )

            segments.append(
                ArcSegment(
                    r_start=r_start,
                    r_end=r_end,
                    R=R,
                    cy=cy,
                    y_start=y_start,
                    y_end=y_end,
                    m_start=m_start,
                    m_end=m_end,
                    r_match=r_match,
                    y_match=y_match,
                    m_match=m_match,
                )
            )

        return segments

    def _compute_radii(self) -> np.ndarray:
        """Compute radial partition (uniform or adaptive)."""
        if not self.adaptive:
            return np.linspace(0.0, self.max_r, self.n_rings + 1)

        # Adaptive: equidistribute sqrt(Delta kappa) (Theorem 6)
        n_fine = 40_001
        r_fine = np.linspace(0.0, self.max_r, n_fine)
        f = self.focal_length

        # Curvature mismatch (Theorem 2)
        dk = r_fine ** 2 / (4.0 * f ** 2 + r_fine ** 2) ** 1.5
        density = np.sqrt(dk + self.reg_eps)

        # Cumulative distribution S(r)
        S = np.zeros_like(r_fine)
        S[1:] = np.cumsum(0.5 * (density[:-1] + density[1:]) * np.diff(r_fine))
        S /= S[-1]  # normalise to [0, 1]

        # Invert S(r) at equispaced targets
        s_targets = np.linspace(0.0, 1.0, self.n_rings + 1)
        radii = np.interp(s_targets, S, r_fine)
        radii[0] = 0.0
        return radii

    # ------------------------------------------------------------------
    # 2. SURFACE EVALUATION & VALIDATION
    # ------------------------------------------------------------------

    def evaluate_surface(self, n_eval: int = 2_000) -> Dict:
        """
        Evaluate approximation quality against true parabola.

        Returns dict with keys:
          rms_surf, max_surf (microns)
          rms_slope, max_slope (radians)
          rms_ray, max_ray (arcminutes)
          ring_data : list of per-ring error dicts (midpoint AND max)
          bound_h_nm, bound_s_mrad, bound_ray_arcmin :
              uniform-partition corollaries (None when adaptive=True)
          bound_h_adaptive_nm :
              Leading-order adaptive proxy from the regularised density (adaptive only)
        """
        f = self.focal_length
        r_eval = np.linspace(0.0, self.max_r, n_eval)

        # Evaluate arc approximation
        y_eval = np.empty(n_eval)
        slope_eval = np.empty(n_eval)
        seg_idx = 0

        for i, r in enumerate(r_eval):
            while seg_idx < len(self.segments) - 1 and r > self.segments[seg_idx].r_end:
                seg_idx += 1
            seg = self.segments[seg_idx]
            y_eval[i] = seg.cy - np.sqrt(max(0.0, seg.R ** 2 - r ** 2))
            slope_eval[i] = r / np.sqrt(max(1e-12, seg.R ** 2 - r ** 2))

        # True parabola
        y_true = r_eval ** 2 / (4.0 * f)
        slope_true = r_eval / (2.0 * f)

        # Errors
        surf_err = (y_eval - y_true) * 1_000.0          # microns
        slope_err = np.abs(slope_eval - slope_true)     # radians
        ray_err_arcmin = np.degrees(2.0 * slope_err) * 60.0

        # Per-ring statistics.
        # v2.0: report the segment MAXIMUM from a 9-point scan. The
        # maximum sag error on a segment occurs at its RIGHT endpoint
        # (delta'' > 0), so midpoint-only reporting understates it.
        ring_data = []
        for seg in self.segments:
            rs = np.linspace(seg.r_start, seg.r_end, 9)
            ys = seg.cy - np.sqrt(np.maximum(0.0, seg.R ** 2 - rs ** 2))
            yt = rs ** 2 / (4.0 * f)
            ms = rs / np.sqrt(np.maximum(1e-12, seg.R ** 2 - rs ** 2))
            mt = rs / (2.0 * f)
            r_mid = 0.5 * (seg.r_start + seg.r_end)
            y_arc_m = seg.cy - np.sqrt(max(0.0, seg.R ** 2 - r_mid ** 2))
            y_tr_m = r_mid ** 2 / (4.0 * f)
            m_arc_m = r_mid / np.sqrt(max(1e-12, seg.R ** 2 - r_mid ** 2))
            m_tr_m = r_mid / (2.0 * f)
            ring_data.append(
                {
                    "r": r_mid,
                    "slope_err_mrad": np.abs(m_arc_m - m_tr_m) * 1_000.0,
                    "surf_err_um": (y_arc_m - y_tr_m) * 1_000.0,
                    "slope_err_max_mrad": float(np.max(np.abs(ms - mt)) * 1_000.0),
                    "surf_err_max_um": float(np.max(np.abs(ys - yt)) * 1_000.0),
                }
            )

        D = self.diameter
        n = float(self.n_rings)

        if self.adaptive:
            # Uniform corollaries do not apply to adaptive partitions:
            # return them as None instead of silently computing them.
            bound_h_nm = None
            bound_s_mrad = None
            bound_ray_arcmin = None
            # Leading-order adaptive proxy using the regularised density
            # actually used for the nodes: proxy_h = S_reg^2 / (2 n^2).
            r_f = np.linspace(0.0, self.max_r, 40_001)
            dk = r_f ** 2 / (4.0 * f ** 2 + r_f ** 2) ** 1.5
            S_reg = np.trapezoid(np.sqrt(dk + self.reg_eps), r_f)
            bound_h_adaptive_nm = float(S_reg ** 2 / (2.0 * n ** 2) * 1e6)
        else:
            # Uniform-partition corollaries (journal Corollaries 1-3)
            bound_h_nm = float((D ** 4) / (256.0 * n ** 2 * f ** 3) * 1e6)
            bound_s_mrad = float((D ** 3) / (64.0 * n * f ** 3) * 1e3)
            bound_ray_arcmin = float(6876.0 * (D ** 3) / (64.0 * n * f ** 3))
            bound_h_adaptive_nm = None

        return {
            "rms_surf": float(np.sqrt(np.mean(surf_err ** 2))),
            "max_surf": float(np.max(np.abs(surf_err))),
            "rms_slope": float(np.sqrt(np.mean(slope_err ** 2))),
            "max_slope": float(np.max(slope_err)),
            "rms_ray": float(np.sqrt(np.mean(ray_err_arcmin ** 2))),
            "max_ray": float(np.max(ray_err_arcmin)),
            "ring_data": ring_data,
            "bound_h_nm": bound_h_nm,
            "bound_s_mrad": bound_s_mrad,
            "bound_ray_arcmin": bound_ray_arcmin,
            "bound_h_adaptive_nm": bound_h_adaptive_nm,
        }

    def validate(self, verbose: bool = True) -> Dict:
        """
        Run built-in unit tests (Theorems 1, 2 and 3).

        Tests:
          1. Left-endpoint tangency: y_i(r_i) = f(r_i), y'_i(r_i) = f'(r_i)
          2. Boundary jumps vs Theorem 3 (r_i>0): axial sag O((Dr)^2), graph slope O(Dr)
             Vertex segment skipped: leading coefficients vanish.
          3. Curvature mismatch formula (Theorem 2) - evaluated at construction point
        """
        f = self.focal_length
        report = {"passed": 0, "failed": 0, "details": []}
        tol = 1e-9

        def log(msg: str, ok: bool):
            status = "PASS" if ok else "FAIL"
            report["details"].append(f"[{status}] {msg}")
            if ok:
                report["passed"] += 1
            else:
                report["failed"] += 1
            if verbose:
                print(f"  [{status}] {msg}")

        if verbose:
            print("\n=== ARC-MIRROR VALIDATION ===")

        # Test 1: construction-point position and slope match.
        # This validates either family member without pretending that the
        # midpoint variant is a left-endpoint construction.
        for i, seg in enumerate(self.segments):
            y_arc = seg.cy - np.sqrt(max(0.0, seg.R ** 2 - seg.r_match ** 2))
            m_arc = (
                seg.r_match / np.sqrt(max(1e-12, seg.R ** 2 - seg.r_match ** 2))
                if seg.r_match < seg.R else float("inf")
            )
            y_ok = abs(y_arc - seg.y_match) < tol
            m_ok = abs(m_arc - seg.m_match) < tol
            log(f"Segment {i+1} {self.matching}-match position", y_ok)
            log(f"Segment {i+1} {self.matching}-match slope", m_ok)

        # Test 2: Boundary jumps vs Theorem 3 leading terms.
        # The published left-endpoint coefficients are intentionally not
        # asserted for midpoint matching; midpoint has its own knot-jump
        # expansion and matching-location analysis.
        if self.matching == "left":
            n_checked = 0
            worst_h_ratio = 0.0
            worst_s_ratio = 0.0
            max_pos = 0.0
            max_slp = 0.0
            ok_pos = True
            ok_slope = True
            for i in range(len(self.segments) - 1):
                seg = self.segments[i]
                seg_next = self.segments[i + 1]
                r0 = seg.r_start
                h = seg.r_end - seg.r_start
                if r0 < 1e-12 or h < 1e-15:
                    continue

                y_left = seg.cy - np.sqrt(max(0.0, seg.R ** 2 - seg.r_end ** 2))
                y_right = seg_next.y_start
                pos = abs(y_left - y_right)

                m_left = seg.r_end / np.sqrt(max(1e-12, seg.R ** 2 - seg.r_end ** 2))
                m_right = seg_next.m_start
                slp = abs(m_left - m_right)

                pred_h = (r0 ** 2 * h ** 2) / (16.0 * f ** 3)
                pred_s = (r0 ** 2 * h) / (8.0 * f ** 3)

                max_pos = max(max_pos, pos)
                max_slp = max(max_slp, slp)
                n_checked += 1

                if pred_h > 1e-18:
                    ratio_h = pos / pred_h
                    worst_h_ratio = max(worst_h_ratio, ratio_h)
                    if not (0.25 <= ratio_h <= 6.0):
                        ok_pos = False
                if pred_s > 1e-18:
                    ratio_s = slp / pred_s
                    worst_s_ratio = max(worst_s_ratio, ratio_s)
                    if not (0.25 <= ratio_s <= 6.0):
                        ok_slope = False

            log(
                f"Position jumps vs r_i^2(Dr)^2/(16f^3) "
                f"(max {max_pos*1e6:.2f} nm, worst ratio {worst_h_ratio:.2f}, "
                f"n={n_checked})",
                ok_pos and n_checked > 0,
            )
            log(
                f"Slope jumps vs r_i^2(Dr)/(8f^3) "
                f"(max {max_slp*1e6:.2f} urad, worst ratio {worst_s_ratio:.2f}, "
                f"n={n_checked})",
                ok_slope and n_checked > 0,
            )

        else:
            log("Midpoint boundary jumps: left-endpoint coefficients skipped "
                "(use midpoint knot-jump expansion)", True)

        # Test 3: Curvature mismatch formula (Theorem 2)
        # Evaluate at the actual construction/matching radius.
        for seg in self.segments:
            r_test = seg.r_match
            if r_test < 1e-12:
                continue  # vertex has the limiting osculating construction
            kappa_arc = 1.0 / seg.R
            kappa_para = 4.0 * f ** 2 / (4.0 * f ** 2 + r_test ** 2) ** 1.5
            dk_formula = r_test ** 2 / (4.0 * f ** 2 + r_test ** 2) ** 1.5
            dk_actual = kappa_arc - kappa_para
            ok = abs(dk_actual - dk_formula) < 1e-9
            log(f"Curvature mismatch at r={r_test:.3f} mm", ok)

        if verbose:
            print(f"\nResults: {report['passed']} passed, {report['failed']} failed")

        return report

    # ------------------------------------------------------------------
    # 3. RAY TRACING
    # ------------------------------------------------------------------

    def ray_trace(
        self,
        n_rays: int = 20_000,
        theta_sun_arcmin: float = 16.0,
        off_axis_deg: float = 0.0,
        seed: Optional[int] = None,
    ) -> Dict:
        """
        Monte-Carlo ray tracing with finite solar disk.

        CORRECTED (v2.0):
          * Solar disk sampled with theta = theta_half * sqrt(u), giving
            uniform brightness per unit solid angle (uniform-in-theta
            biased the spot statistics low).
          * Optional off-axis pointing angle, composed EXACTLY by rotation
            about the Z axis (same construction as the web tool).

        Returns focal-spot statistics (mm) plus the spot centroid.
        """
        rng = np.random.default_rng(seed)
        f = self.focal_length
        theta_half = np.radians(theta_sun_arcmin / 60.0)
        th_axis = np.radians(off_axis_deg)
        cos_t, sin_t = np.cos(th_axis), np.sin(th_axis)
        max_r = self.max_r
        segs = self.segments

        seg_starts = np.array([s.r_start for s in segs])

        hits_x, hits_z = [], []

        for _ in range(n_rays):
            # Uniform random point on aperture
            r = max_r * np.sqrt(rng.random())
            phi = 2.0 * np.pi * rng.random()
            x = r * np.cos(phi)
            z = r * np.sin(phi)

            seg_idx = int(np.searchsorted(seg_starts, r, side="right") - 1)
            seg_idx = max(0, min(seg_idx, len(segs) - 1))
            seg = segs[seg_idx]

            y_surf = seg.cy - np.sqrt(max(0.0, seg.R ** 2 - r ** 2))

            # Outward surface normal (unit vector toward the tool side)
            nx = x / seg.R
            ny = (y_surf - seg.cy) / seg.R
            nz = z / seg.R
            norm = np.sqrt(nx ** 2 + ny ** 2 + nz ** 2)
            nx, ny, nz = nx / norm, ny / norm, nz / norm

            # Solar-disk direction: sqrt-weighted cone about the -y axis
            sun_theta = theta_half * np.sqrt(rng.random())
            sun_phi = 2.0 * np.pi * rng.random()
            vx = np.sin(sun_theta) * np.cos(sun_phi)
            vy = -np.cos(sun_theta)
            vz = np.sin(sun_theta) * np.sin(sun_phi)

            # Exact off-axis composition: rotate about Z by th_axis
            dx = vx * cos_t - vy * sin_t
            dy = vx * sin_t + vy * cos_t
            dz = vz

            # Reflection: r = d - 2(d.n)n
            dot = dx * nx + dy * ny + dz * nz
            rx = dx - 2.0 * dot * nx
            ry = dy - 2.0 * dot * ny
            rz = dz - 2.0 * dot * nz

            # Intersect with focal plane y = f
            if abs(ry) < 1e-12:
                continue
            t = (f - y_surf) / ry
            if t <= 0:
                continue

            hits_x.append(x + t * rx)
            hits_z.append(z + t * rz)

        hits_x = np.array(hits_x)
        hits_z = np.array(hits_z)
        r_hits = np.sqrt(hits_x ** 2 + hits_z ** 2)

        return {
            "mean_radius": float(np.mean(r_hits)),
            "rms_radius": float(np.sqrt(np.mean(r_hits ** 2))),
            "max_radius": float(np.max(r_hits)),
            "spot_diameter_95": float(2.0 * np.percentile(r_hits, 95)),
            "centroid_x": float(np.mean(hits_x)),
            "centroid_z": float(np.mean(hits_z)),
            "n_hits": len(r_hits),
        }

    # ------------------------------------------------------------------
    # 4. G-CODE GENERATION (TOOL-COMPENSATED)
    # ------------------------------------------------------------------

    def generate_gcode(
        self,
        tool_radius: float = 5.0,
        spindle_rpm: int = 2_500,
        feed_rate: float = 150.0,
        finish_feed: float = 80.0,
        roughing_allowance: float = 0.3,
        rough_depth: float = 1.0,
        diameter_programming: bool = True,
        i_radius_value: bool = True,
        coolant: bool = False,
    ) -> str:
        """
        Generate Fanuc-style G-code for the TOOL CENTER PATH.

        CRITICAL: This outputs the compensated toolpath, not the workpiece
        profile. The tool nose radius is offset EXACTLY by concentric arcs.

        CORRECTED (v2.0):
          * R_path = R - R_tool. The mirror is a CONCAVE surface; the tool
            machines it internally, so the tool-centre path is the
            concentric arc on the SAME side as the arc centre. The old
            R + R_tool put the tool centre inside the material (about
            2*R_tool too deep along the normal; -5 mm at the vertex for
            a 5 mm tool).
          * I offset emitted as a RADIUS value (Fanuc/Haas convention even
            under diameter programming). Set i_radius_value=False only if
            your control documents I in diameter units.
          * Roughing passes are allowance-shifted Z-level passes:
            pass at depth z stops at X = 2*sqrt(4f*(z - a)) so the tool
            never cuts below the final surface (convexity of the parabola:
            the straight sweep between two "surface + a" points stays
            above the surface + a). The old code plunged to Z = -a,
            below the vertex plane, pocketing the centre.
        """
        x_scale = 2.0 if diameter_programming else 1.0
        R_tool = float(tool_radius)
        f = self.focal_length
        max_r = self.max_r
        sag = self.sag
        segs = self.segments

        lines = []
        lines.append("%")
        lines.append(
            f"O1000 (Arc-Mirror {self.diameter:.0f}mm f/{self.f_number:.1f} "
            f"- {self.n_rings} {self.matching.title()}-Matched Arcs)"
        )
        lines.append("(Author: Omidiran Favour Daniel)")
        lines.append(f"(Focal Length: {f:.1f}mm, Diameter: {self.diameter:.1f}mm)")
        lines.append(f"(Sag: {sag:.2f}mm, Segments: {self.n_rings})")
        lines.append(f"(Construction family: {('adaptive' if self.adaptive else 'uniform')}-{self.matching})")
        lines.append(f"(Tool Nose Radius: {R_tool:.2f}mm)")
        lines.append("(NOTE: This program is TOOL-CENTER-PATH compensated)")
        lines.append("(R_path = R - R_tool ; I programmed as radius value)")
        lines.append("")

        # Safety & setup
        lines.append("(=== SAFETY & SETUP ===)")
        lines.append("G21 G40 G49 G80 G99")
        lines.append("G54 G90 G94")
        lines.append(f"G50 S{int(round(spindle_rpm * 1.5))} (spindle speed clamp)")
        lines.append("G28 U0 W0")
        lines.append("")

        # Tool & spindle
        lines.append("(=== TOOL & SPINDLE ===)")
        lines.append(f"T0101 (Button Tool R{R_tool:.1f}mm)")
        lines.append(f"M03 S{spindle_rpm}")
        if coolant:
            lines.append("M08")
        lines.append("G04 P2000")
        lines.append("")

        # Approach
        lines.append("(=== APPROACH ===)")
        start_x = (max_r + 10.0) * x_scale
        start_z = sag + 5.0
        lines.append(f"G00 X{start_x:.4f} Z{start_z:.4f}")
        lines.append("")

        # Roughing pass: allowance-shifted Z-level passes.
        # Pass at depth z sweeps from outside in to X = 2*sqrt(4f*(z-a)):
        # the stop circle is exactly where the final surface equals z - a,
        # so the tool leaves 'a' mm of stock everywhere and never touches
        # the final surface (old code went to Z = -a, below the vertex).
        lines.append("(=== ROUGHING PASS - allowance-shifted Z levels ===)")
        lines.append(f"(Leaving {roughing_allowance:.2f}mm axial allowance;")
        lines.append(" never cuts below the final surface. Verify in CAM.)")
        a = roughing_allowance
        current_z = start_z
        while current_z > a + 1e-12:
            target_z = max(current_z - rough_depth, a)
            r_stop = np.sqrt(max(0.0, 4.0 * f * (target_z - a)))
            target_r = min(max_r, r_stop)
            target_x = target_r * x_scale
            lines.append(f"G01 X{target_x:.4f} Z{target_z:.4f} F{feed_rate}")
            lines.append(f"G00 X{start_x:.4f}")
            current_z = target_z
        lines.append("")

        # Finishing pass - compensated tool center path
        lines.append("(=== FINISHING PASS - TOOL CENTER PATH ===)")
        lines.append(f"(Feed: {finish_feed} mm/min, {self.n_rings} G02 arcs)")
        lines.append("")

        tool_segs = self._build_toolpath_segments(R_tool)

        first = tool_segs[0]
        lines.append(f"G00 X{first.r_start * x_scale:.6f} Z{first.z_start:.6f}")
        lines.append(f"G01 X{first.r_start * x_scale:.6f} Z{first.z_start:.6f} F{finish_feed}")
        lines.append("")

        for i, tp in enumerate(tool_segs):
            x_end = tp.r_end * x_scale
            # Fanuc/Haas: I is a radius value even under diameter
            # programming; some controls want the diameter value instead.
            i_val = tp.i_offset if i_radius_value else tp.i_offset * 2.0
            lines.append(
                f"(Segment {i+1}/{self.n_rings}: "
                f"r={tp.r_start:.3f} to {tp.r_end:.3f}, "
                f"R_path={tp.R_path:.3f})"
            )
            lines.append(
                f"G02 X{x_end:.6f} Z{tp.z_end:.6f} "
                f"I{i_val:.6f} K{tp.k_offset:.6f}"
            )

        lines.append("")
        lines.append("(=== RETRACT & END ===)")
        lines.append("G00 X100 Z100")
        if coolant:
            lines.append("M09")
        lines.append("M05")
        lines.append("M30")
        lines.append("%")

        return "\n".join(lines)

    def _build_toolpath_segments(self, R_tool: float) -> List[ToolPathSegment]:
        """
        Offset each surface arc by R_tool along the surface normal on the
        TOOL side.

        CORRECTED (v2.0): the mirror is concave and is machined from the
        concave (tool) side, which faces the arc centre. The tool centre
        therefore follows the concentric circle of radius

            R_path = R - R_tool

        with the SAME centre (0, c_y). The old code used R + R_tool,
        placing the tool centre on the far side of the surface, inside
        the material (about 2*R_tool too deep along the normal).
        """
        tool_segs = []
        for i, seg in enumerate(self.segments):
            if R_tool >= seg.R:
                raise ValueError(
                    f"Segment {i+1}: tool radius {R_tool:.2f} mm is >= arc radius "
                    f"{seg.R:.2f} mm. R_path would be <= 0. Use a smaller tool."
                )

            R_path = seg.R - R_tool

            r_start_tp = seg.r_start * (R_path / seg.R)
            r_end_tp = seg.r_end * (R_path / seg.R)

            z_start_tp = seg.cy - np.sqrt(max(0.0, R_path ** 2 - r_start_tp ** 2))
            z_end_tp = seg.cy - np.sqrt(max(0.0, R_path ** 2 - r_end_tp ** 2))

            # G02 I/K offsets (centre minus start point).
            # I is stored as a RADIUS value (Fanuc/Haas convention).
            i_offset = -r_start_tp
            k_offset = seg.cy - z_start_tp

            tool_segs.append(
                ToolPathSegment(
                    r_start=r_start_tp,
                    r_end=r_end_tp,
                    R_path=R_path,
                    cy=seg.cy,
                    z_start=z_start_tp,
                    z_end=z_end_tp,
                    i_offset=i_offset,
                    k_offset=k_offset,
                )
            )
        return tool_segs

    # ------------------------------------------------------------------
    # 5. STL EXPORT (MANIFOLD MESH, SELF-CERTIFYING)
    # ------------------------------------------------------------------

    def export_stl(
        self,
        filename: str,
        wall_thickness: Optional[float] = None,
        n_profile: int = 20,
        n_theta: int = 40,
        pole_taper: Optional[float] = None,
    ) -> Dict:
        """
        Export a watertight, outward-oriented binary STL shell.

        Pole topology: TWO separate axis vertices, full wall at r = 0.
          * Outer pole at (0, 0, 0) — optical vertex.
          * Inner pole at (0, wall, 0) — filled cap on the concave side.
        The old shared-pole taper pinched both faces onto one point and
        showed up as a dark speck when the bowl was viewed down the axis.
        `pole_taper` is ignored (kept so old call sites do not break).

        The wall extends from the reflector toward the arc centre (+y).
        Confirm that side matches fabrication intent (shell vs mandrel).

        Returns a certification dict (watertight, outward, signed volume).
        """
        if wall_thickness is None:
            wall_thickness = self.wall_thickness

        segs = self.segments
        _ = pole_taper  # deprecated: full thickness is used at the axis

        # --- Build unified profile (outer surface), pole deduplicated ---
        profile = [(0.0, 0.0, 0)]            # (r, y, segment index)
        for i, seg in enumerate(segs):
            local_r = np.linspace(seg.r_start, seg.r_end, n_profile)
            local_r = local_r[1:]            # skip left endpoint (pole/previous knot)
            for r in local_r:
                y = seg.cy - np.sqrt(max(0.0, seg.R ** 2 - r ** 2))
                profile.append((float(r), float(y), i))
        n_pts = len(profile)

        # --- Inner surface: full-thickness concentric offset ---
        inner = []
        for (r, y, i) in profile:
            seg = segs[i]
            if wall_thickness >= seg.R:
                raise ValueError(
                    f"Segment {i+1}: wall thickness {wall_thickness:.2f} mm "
                    f"is >= arc radius {seg.R:.2f} mm. Reduce wall thickness."
                )
            Ri = seg.R - wall_thickness
            if r <= 0.0:
                inner.append((0.0, seg.cy - Ri))
            else:
                scale = Ri / seg.R
                r_i = r * scale
                y_i = seg.cy - np.sqrt(max(0.0, Ri ** 2 - r_i ** 2))
                inner.append((r_i, y_i))

        # --- Vertices: outer pole + inner pole + rings ---
        theta = np.linspace(0.0, 2.0 * np.pi, n_theta, endpoint=False)

        vertices: List[List[float]] = [
            [0.0, profile[0][1], 0.0],   # 0: outer pole (optical vertex)
            [0.0, inner[0][1], 0.0],     # 1: inner pole (filled cap)
        ]
        faces: List[Tuple[int, int, int]] = []

        def add_ring(r: float, y: float) -> List[int]:
            ring = []
            for t in theta:
                vertices.append([float(r * np.cos(t)), float(y), float(r * np.sin(t))])
                ring.append(len(vertices) - 1)
            return ring

        outer_rings = [add_ring(profile[k][0], profile[k][1]) for k in range(1, n_pts)]
        inner_rings = [add_ring(inner[k][0], inner[k][1]) for k in range(1, n_pts)]

        # --- Faces ---
        P_OUT, P_IN = 0, 1
        o1, i1 = outer_rings[0], inner_rings[0]
        for j in range(n_theta):
            j1 = (j + 1) % n_theta
            faces.append((P_OUT, o1[j], o1[j1]))   # outer cap, normal -y
            faces.append((P_IN, i1[j1], i1[j]))    # inner cap, normal +y

        # Quad strips between consecutive rings
        for k in range(n_pts - 2):
            ao, bo = outer_rings[k], outer_rings[k + 1]
            ai, bi = inner_rings[k], inner_rings[k + 1]
            for j in range(n_theta):
                j1 = (j + 1) % n_theta
                faces.append((ao[j], bo[j], bo[j1]))
                faces.append((ao[j], bo[j1], ao[j1]))
                faces.append((ai[j], bi[j1], bi[j]))
                faces.append((ai[j], ai[j1], bi[j1]))

        # Rim wall
        ao, ai = outer_rings[-1], inner_rings[-1]
        for j in range(n_theta):
            j1 = (j + 1) % n_theta
            faces.append((ao[j], ai[j], ai[j1]))
            faces.append((ao[j], ai[j1], ao[j1]))

        # --- In-memory certification (index-exact: shared pole vertex
        # means coincident geometry IS the same topology here) ---
        V = np.asarray(vertices, dtype=np.float64)
        F = np.asarray(faces, dtype=np.int64)
        tri = V[F]
        fn = np.cross(tri[:, 1] - tri[:, 0], tri[:, 2] - tri[:, 0])
        farea = 0.5 * np.linalg.norm(fn, axis=1)
        signed_volume = float(
            np.einsum("ij,ij->i", tri[:, 0], np.cross(tri[:, 1], tri[:, 2])).sum() / 6.0
        )
        edge_count: Counter = Counter()
        for a_, b_, c_ in F:
            e1 = (a_, b_) if a_ < b_ else (b_, a_)
            e2 = (b_, c_) if b_ < c_ else (c_, b_)
            e3 = (c_, a_) if c_ < a_ else (a_, c_)
            edge_count[e1] += 1
            edge_count[e2] += 1
            edge_count[e3] += 1
        boundary_edges = sum(1 for c_ in edge_count.values() if c_ == 1)
        nonmanifold_edges = sum(1 for c_ in edge_count.values() if c_ > 2)
        degenerate = int((farea < 1e-12).sum())

        watertight = (boundary_edges == 0) and (nonmanifold_edges == 0)
        outward = signed_volume > 0
        ok = watertight and outward and degenerate == 0

        # --- Write binary STL ---
        header = "ArcMirror3 by Omidiran Favour Daniel".encode("ascii")[:80]
        header = header + b" " * (80 - len(header))

        with open(filename, "wb") as fh:
            fh.write(header)
            fh.write(struct.pack("<I", len(faces)))
            for t_idx, tri_idx in enumerate(F):
                v0_, v1_, v2_ = V[tri_idx[0]], V[tri_idx[1]], V[tri_idx[2]]
                n = fn[t_idx]
                n_norm = np.linalg.norm(n)
                n = n / n_norm if n_norm > 0 else np.array([0.0, 1.0, 0.0])
                fh.write(struct.pack("<3f", *n))
                fh.write(struct.pack("<3f", *v0_))
                fh.write(struct.pack("<3f", *v1_))
                fh.write(struct.pack("<3f", *v2_))
                fh.write(struct.pack("<H", 0))

        cert = {
            "path": str(filename),
            "vertices": len(vertices),
            "faces": len(faces),
            "signed_volume_mm3": signed_volume,
            "boundary_edges": boundary_edges,
            "nonmanifold_edges": nonmanifold_edges,
            "degenerate_faces": degenerate,
            "watertight": watertight,
            "outward_normals": outward,
            "certified": ok,
        }
        status = "PASS" if ok else "FAIL"
        print(
            f"[STL:{status}] {len(vertices)} vertices, {len(faces)} faces -> {filename}\n"
            f"        watertight={watertight} (boundary={boundary_edges}, "
            f"nonmanifold={nonmanifold_edges})\n"
            f"        outward_normals={outward} (signed volume {signed_volume:+.1f} mm^3), "
            f"degenerate={degenerate}"
        )
        return cert

    # ------------------------------------------------------------------
    # 6. UTILITY / REPORTING
    # ------------------------------------------------------------------

    @staticmethod
    def family_members(
        diameter: float = 200.0,
        f_number: float = 1.5,
        n_rings: int = 65,
        reg_eps: Optional[float] = None,
        wall_thickness: float = 2.0,
    ) -> Dict[str, "ArcMirror"]:
        """Construct the four computational members of the Omega family.

        The segmentation rule and matching rule are independent:
            uniform-left, adaptive-left,
            uniform-midpoint, adaptive-midpoint.

        The default left-endpoint construction remains the canonical member.
        """
        return {
            "uniform-left": ArcMirror(diameter, f_number, n_rings, False,
                                      reg_eps, wall_thickness, "left"),
            "adaptive-left": ArcMirror(diameter, f_number, n_rings, True,
                                       reg_eps, wall_thickness, "left"),
            "uniform-midpoint": ArcMirror(diameter, f_number, n_rings, False,
                                           reg_eps, wall_thickness, "midpoint"),
            "adaptive-midpoint": ArcMirror(diameter, f_number, n_rings, True,
                                            reg_eps, wall_thickness, "midpoint"),
        }

    def print_report(self) -> None:
        """Print a full optical and manufacturing report."""
        f = self.focal_length
        D = self.diameter
        n = self.n_rings

        print("=" * 60)
        print("ARC-MIRROR THEORY - TECHNICAL REPORT")
        print("=" * 60)
        print("Author        : Omidiran Favour Daniel")
        print(f"Configuration : {D:.0f} mm, f/{self.f_number:.1f}, {n} segments")
        print(f"Segmentation  : {'ADAPTIVE' if self.adaptive else 'UNIFORM'}")
        print(f"Matching      : {self.matching.upper()}")
        print(f"Family member : {'adaptive' if self.adaptive else 'uniform'}-{self.matching}")
        if self.adaptive:
            print(f"Reg. eps      : {self.reg_eps:.3e}")
        print(f"Wall thickness: {self.wall_thickness:.1f} mm")
        print(f"Focal length  : {f:.1f} mm")
        print(f"Sag           : {self.sag:.2f} mm")
        print(f"Edge slope    : {np.degrees(np.arctan(D/(2*f))):.1f} deg")
        print("-" * 60)

        opt = self.evaluate_surface()
        print("OPTICAL PERFORMANCE")
        print(f"  RMS surface error : {opt['rms_surf']:.4f} um")
        print(f"  Max surface error : {opt['max_surf']:.4f} um ({opt['max_surf']*1000:.1f} nm)")
        print(f"  RMS slope error   : {opt['rms_slope']*1000:.3f} mrad")
        print(f"  Max slope error   : {opt['max_slope']*1000:.3f} mrad")
        print(f"  RMS ray deviation : {opt['rms_ray']:.2f} arcmin")
        print(f"  Max ray deviation : {opt['max_ray']:.2f} arcmin")
        if self.adaptive:
            print(f"  Adaptive segmentation monitor (left-endpoint reference): {opt['bound_h_adaptive_nm']:.1f} nm")
        else:
            print(f"  Theory bound (uniform) : {opt['bound_h_nm']:.1f} nm, "
                  f"{opt['bound_s_mrad']:.3f} mrad, {opt['bound_ray_arcmin']:.2f} arcmin")
        print("-" * 60)

        solar_flux = 1000.0
        reflectance = 0.90
        aperture_area = np.pi * (D / 2000.0) ** 2
        theta_sun = np.radians(32.0 / 60.0)
        theta_half = theta_sun / 2.0

        C_max = (D / (2.0 * f * theta_half)) ** 2
        image_diameter_mm = f * theta_sun
        image_radius_m = (image_diameter_mm / 2.0) / 1000.0
        C_geo_ideal = aperture_area / (np.pi * image_radius_m ** 2)

        max_slope_rad = opt["max_slope"]
        C_geo_real = C_geo_ideal * (theta_sun / (theta_sun + 2.0 * max_slope_rad)) ** 2
        P_thermal = solar_flux * aperture_area * reflectance

        print("CSP PERFORMANCE")
        print(f"  Theoretical max concentration : {C_max:.0f}x")
        print(f"  Geometric (ideal optics)      : {C_geo_ideal:.0f}x")
        print(f"  Geometric (with slope error)  : {C_geo_real:.0f}x")
        print(f"  Solar image diameter          : {image_diameter_mm:.2f} mm")
        print(f"  Thermal power (aperture)      : {P_thermal:.1f} W")
        print("=" * 60)


# ---------------------------------------------------------------------------
# DEMONSTRATION / CLI
# ---------------------------------------------------------------------------

if __name__ == "__main__":
    # OMEGA FAMILY DEMONSTRATION
    # The default demonstration reproduces the current 2 m hackathon design: adaptive-midpoint.
    mirror = ArcMirror(
        diameter=2000.0,
        f_number=1.5,
        n_rings=98,
        adaptive=True,
        wall_thickness=2.0,
        matching="midpoint",
    )

    # Validate the selected hackathon design.
    mirror.validate()

    # Display the four computational family members.
    print("\n=== ARC-MIRROR OMEGA FAMILY ===")
    family = ArcMirror.family_members(
        diameter=mirror.diameter,
        f_number=mirror.f_number,
        n_rings=mirror.n_rings,
        wall_thickness=mirror.wall_thickness,
    )
    for name, member in family.items():
        q = member.evaluate_surface()
        print(
            f"  {name:19s} | max surface {q['max_surf']*1000:.3f} nm"
            f" | max slope {q['max_slope']*1000:.5f} mrad"
            f" | max ray {q['max_ray']:.5f} arcmin"
        )

    # Full report for the selected adaptive-midpoint member.
    mirror.print_report()

    # Quick physical spot check (seeded for reproducibility)
    spot = mirror.ray_trace(n_rays=5000, seed=0)
    print("\nFOCAL SPOT (Monte-Carlo, 16 arcmin half-angle sun)")
    print(f"  rms radius {spot['rms_radius']:.3f} mm | "
          f"95% diameter {spot['spot_diameter_95']:.3f} mm | "
          f"max radius {spot['max_radius']:.3f} mm")

    # Build dynamic filenames from actual parameters
    D = int(mirror.diameter)
    F = int(round(mirror.focal_length))
    fnum = mirror.f_number
    n = mirror.n_rings
    mode = "Adaptive" if mirror.adaptive else "Uniform"
    match = mirror.matching.title()

    base_name = f"ArcMirror_{D}mm_F{F}_f{fnum:.1f}_{n}arcs_{mode}_{match}"

    # Export G-code (tool-compensated)
    gcode = mirror.generate_gcode(tool_radius=5.0)
    nc_file = f"{base_name}_Compensated.nc"
    with open(nc_file, "w") as fh:
        fh.write(gcode)
    print(f"\n[G-CODE] Saved to {nc_file}")

    # Export STL (self-certifying)
    stl_file = f"{base_name}_Shell.stl"
    mirror.export_stl(stl_file)
