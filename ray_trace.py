#!/usr/bin/env python3
"""
Arc-Mirror Theory – Omega Ray-Trace Demonstration
=========================================================
Self-contained visual demonstration of the full Omega family:
  uniform-left, adaptive-left, uniform-midpoint, adaptive-midpoint.

Shows surface profile, meridian rays, and Monte-Carlo focal spot
for any chosen family member.

Author: Omidiran Favour Daniel
"""

import numpy as np
import matplotlib.pyplot as plt
from dataclasses import dataclass
from typing import List, Dict, Optional

# ----------------------------------------------------------------------
# DATA STRUCTURES
# ----------------------------------------------------------------------

@dataclass
class ArcSegment:
    r_start: float
    r_end: float
    R: float
    cy: float
    y_start: float
    y_end: float
    m_start: float
    m_end: float
    r_match: float
    y_match: float
    m_match: float


# ----------------------------------------------------------------------
# CORE MIRROR (Omega family)
# ----------------------------------------------------------------------

class ArcMirror:
    """
    Omega-family constructor supporting:
      matching = "left" | "midpoint"
      adaptive = True | False
    """

    def __init__(
        self,
        diameter: float = 2000.0,
        f_number: float = 1.5,
        n_rings: int = 98,
        adaptive: bool = True,
        matching: str = "midpoint",
        reg_eps: Optional[float] = None,
    ):
        self.diameter = float(diameter)
        self.f_number = float(f_number)
        self.focal_length = self.diameter * self.f_number
        self.n_rings = int(n_rings)
        self.adaptive = bool(adaptive)
        self.matching = str(matching).lower().strip()
        if self.matching not in {"left", "midpoint"}:
            raise ValueError("matching must be 'left' or 'midpoint'")

        self.max_r = self.diameter / 2.0
        self.sag = self.diameter ** 2 / (16.0 * self.focal_length)

        if self.adaptive and reg_eps is None:
            f = self.focal_length
            Rm = self.max_r
            S_max = 2.0 * ((4.0 * f**2 + Rm**2)**0.25 - (4.0 * f**2)**0.25)
            self.reg_eps = S_max / (2.0 * np.sqrt(2.0) * f**1.5 * self.n_rings)
        else:
            self.reg_eps = float(reg_eps) if reg_eps is not None else 0.0

        self.segments = self._build_segments()

    def _compute_radii(self) -> np.ndarray:
        if not self.adaptive:
            return np.linspace(0.0, self.max_r, self.n_rings + 1)

        n_fine = 10_000
        r_fine = np.linspace(1e-6, self.max_r, n_fine)
        f = self.focal_length
        dk = r_fine ** 2 / (4.0 * f**2 + r_fine**2) ** 1.5
        density = np.sqrt(dk + self.reg_eps)
        S = np.zeros_like(r_fine)
        S[1:] = np.cumsum(0.5 * (density[:-1] + density[1:]) * np.diff(r_fine))
        S /= S[-1]
        s_targets = np.linspace(0.0, 1.0, self.n_rings + 1)
        radii = np.interp(s_targets, S, r_fine)
        radii[0] = 0.0
        return radii

    def _build_segments(self) -> List[ArcSegment]:
        radii = self._compute_radii()
        segments = []
        f = self.focal_length

        for i in range(self.n_rings):
            r_start = radii[i]
            r_end = radii[i + 1]

            if self.matching == "left":
                r_match = r_start
            else:
                r_match = 0.5 * (r_start + r_end)

            y_match = r_match ** 2 / (4.0 * f)
            m_match = r_match / (2.0 * f) if r_match > 1e-15 else 0.0

            if r_match < 1e-15:
                R = 2.0 * f
                cy = R
            else:
                r_over_m = r_match / m_match
                R = np.sqrt(r_match ** 2 + r_over_m ** 2)
                cy = y_match + r_over_m

            if r_end >= R:
                raise ValueError(
                    f"Segment {i+1}: r_end={r_end:.6g} reaches/exceeds R={R:.6g}"
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
                    r_start=r_start, r_end=r_end, R=R, cy=cy,
                    y_start=y_start, y_end=y_end,
                    m_start=m_start, m_end=m_end,
                    r_match=r_match, y_match=y_match, m_match=m_match,
                )
            )
        return segments

    def evaluate_surface(self, n_eval: int = 2000) -> Dict:
        f = self.focal_length
        r_eval = np.linspace(0.0, self.max_r, n_eval)
        y_eval = np.empty_like(r_eval)
        slope_eval = np.empty_like(r_eval)
        seg_idx = 0
        for i, r in enumerate(r_eval):
            while seg_idx < len(self.segments) - 1 and r > self.segments[seg_idx].r_end:
                seg_idx += 1
            seg = self.segments[seg_idx]
            y_eval[i] = seg.cy - np.sqrt(max(0.0, seg.R ** 2 - r ** 2))
            slope_eval[i] = r / np.sqrt(max(1e-12, seg.R ** 2 - r ** 2))

        y_true = r_eval ** 2 / (4.0 * f)
        slope_true = r_eval / (2.0 * f)
        surf_err = (y_eval - y_true) * 1_000.0          # microns
        slope_err = np.abs(slope_eval - slope_true)     # radians
        ray_err = np.degrees(2.0 * slope_err) * 60.0    # arcmin

        return {
            "max_surf_nm": float(np.max(np.abs(surf_err)) * 1000.0),
            "max_slope_mrad": float(np.max(slope_err) * 1000.0),
            "max_ray_arcmin": float(np.max(ray_err)),
            "rms_surf_nm": float(np.sqrt(np.mean(surf_err**2)) * 1000.0),
            "rms_slope_mrad": float(np.sqrt(np.mean(slope_err**2)) * 1000.0),
            "rms_ray_arcmin": float(np.sqrt(np.mean(ray_err**2))),
        }

    @staticmethod
    def family_members(diameter=200.0, f_number=1.5, n_rings=65):
        return {
            "uniform-left": ArcMirror(diameter, f_number, n_rings, False, "left"),
            "adaptive-left": ArcMirror(diameter, f_number, n_rings, True, "left"),
            "uniform-midpoint": ArcMirror(diameter, f_number, n_rings, False, "midpoint"),
            "adaptive-midpoint": ArcMirror(diameter, f_number, n_rings, True, "midpoint"),
        }


# ----------------------------------------------------------------------
# RAY COLLECTION
# ----------------------------------------------------------------------

def collect_meridian_rays(mirror: ArcMirror, n_rays: int = 13):
    segs = mirror.segments
    starts = np.array([s.r_start for s in segs])
    f = mirror.focal_length
    xs = np.linspace(-mirror.max_r * 0.96, mirror.max_r * 0.96, n_rays)
    traces = []
    for x in xs:
        r = abs(float(x))
        idx = int(np.searchsorted(starts, r, side="right") - 1)
        idx = max(0, min(idx, len(segs) - 1))
        seg = segs[idx]
        y = seg.cy - np.sqrt(max(0.0, seg.R**2 - r**2))
        nx = x / seg.R
        ny = (y - seg.cy) / seg.R
        nrm = np.hypot(nx, ny)
        nx, ny = nx / nrm, ny / nrm
        dx, dy = 0.0, -1.0
        dot = dx * nx + dy * ny
        rx = dx - 2.0 * dot * nx
        ry = dy - 2.0 * dot * ny
        if abs(ry) < 1e-12:
            continue
        t = (f - y) / ry
        xf, yf = x + t * rx, y + t * ry
        traces.append(((x, y + 40.0), (x, y), (xf, yf)))
    return traces


def collect_spot(mirror: ArcMirror, n_rays: int = 8000, sun_half_angle_arcmin: float = 16.0, seed: int = 1):
    rng = np.random.default_rng(seed)
    f = mirror.focal_length
    theta_half = np.radians(sun_half_angle_arcmin / 60.0)
    segs = mirror.segments
    starts = np.array([s.r_start for s in segs])
    hx, hz = [], []
    for _ in range(n_rays):
        r = mirror.max_r * np.sqrt(rng.random())
        phi = 2.0 * np.pi * rng.random()
        x = r * np.cos(phi)
        z = r * np.sin(phi)
        idx = int(np.searchsorted(starts, r, side="right") - 1)
        idx = max(0, min(idx, len(segs) - 1))
        seg = segs[idx]
        y = seg.cy - np.sqrt(max(0.0, seg.R**2 - r**2))
        nx = x / seg.R
        ny = (y - seg.cy) / seg.R
        nz = z / seg.R
        nrm = np.sqrt(nx*nx + ny*ny + nz*nz)
        nx, ny, nz = nx/nrm, ny/nrm, nz/nrm
        # Uniform-brightness solar disk: sqrt(area) sampling makes rays uniform
        # over the disk area. 16 arcmin is the solar half-angle, not full diameter.
        theta = theta_half * np.sqrt(rng.random())
        sp = rng.uniform(0.0, 2.0 * np.pi)
        dx = np.sin(theta) * np.cos(sp)
        dy = -np.cos(theta)
        dz = np.sin(theta) * np.sin(sp)
        dnm = np.sqrt(dx*dx + dy*dy + dz*dz)
        dx, dy, dz = dx/dnm, dy/dnm, dz/dnm
        dot = dx*nx + dy*ny + dz*nz
        rx = dx - 2.0*dot*nx
        ry = dy - 2.0*dot*ny
        rz = dz - 2.0*dot*nz
        if abs(ry) < 1e-12:
            continue
        t = (f - y) / ry
        if t <= 0:
            continue
        hx.append(x + t*rx)
        hz.append(z + t*rz)
    return np.array(hx), np.array(hz)


def surface_profile(mirror: ArcMirror, n: int = 400):
    r = np.linspace(-mirror.max_r, mirror.max_r, n)
    y = np.empty_like(r)
    starts = np.array([s.r_start for s in mirror.segments])
    for i, ri in enumerate(r):
        idx = int(np.searchsorted(starts, abs(ri), side="right") - 1)
        idx = max(0, min(idx, len(mirror.segments) - 1))
        seg = mirror.segments[idx]
        y[i] = seg.cy - np.sqrt(max(0.0, seg.R**2 - ri**2))
    return r, y


# ----------------------------------------------------------------------
# MAIN DEMO
# ----------------------------------------------------------------------

def main():
    print("=" * 70)
    print("ARC-MIRROR OMEGA – RAY-TRACE DEMONSTRATION")
    print("=" * 70)

    # Choose the member to visualise (change these three lines as needed)
    # Default to the current hackathon design case.
    # These are commanded-geometry ray-trace settings, not measured hardware.
    diameter = 2000.0
    f_number = 1.5
    n_rings = 98
    matching = "midpoint"      # "left" or "midpoint"
    adaptive = True

    mirror = ArcMirror(
        diameter=diameter,
        f_number=f_number,
        n_rings=n_rings,
        adaptive=adaptive,
        matching=matching,
    )

    print(f"Configuration : {diameter:.0f} mm, f/{f_number:.1f}, {n_rings} segments")
    print(f"Family member : {'adaptive' if adaptive else 'uniform'}-{matching}")
    print()

    # Optical metrics
    opt = mirror.evaluate_surface()
    print("OPTICAL PERFORMANCE")
    print(f"  Max surface error : {opt['max_surf_nm']:.1f} nm")
    print(f"  Max slope error   : {opt['max_slope_mrad']:.3f} mrad")
    print(f"  Max ray deviation : {opt['max_ray_arcmin']:.2f} arcmin")
    print()

    # Family comparison
    print("=== OMEGA FAMILY (same n) ===")
    family = ArcMirror.family_members(diameter, f_number, n_rings)
    for name, member in family.items():
        q = member.evaluate_surface()
        print(
            f"  {name:19s} | max surface {q['max_surf_nm']:7.1f} nm"
            f" | max slope {q['max_slope_mrad']:7.4f} mrad"
            f" | max ray {q['max_ray_arcmin']:6.3f} arcmin"
        )
    print()

    # Collect rays and spot
    print("Ray tracing...")
    traces = collect_meridian_rays(mirror, n_rays=13)
    hx, hz = collect_spot(mirror, n_rays=8000)
    rh = np.hypot(hx, hz)
    r_prof, y_prof = surface_profile(mirror)

    # Figure
    fig = plt.figure(figsize=(12, 5))
    gs = fig.add_gridspec(2, 2, width_ratios=(1.2, 1), height_ratios=(1, 1))

    # (Top left) Meridian rays
    ax1 = fig.add_subplot(gs[0, 0])
    ax1.plot(r_prof, y_prof, "b-", lw=1.2, alpha=0.8, label="Arc surface")
    for (x1, y1), (x2, y2), (xf, yf) in traces:
        ax1.plot([x1, x2], [y1, y2], "cyan", lw=0.8, alpha=0.5)
        ax1.plot([x2, xf], [y2, yf], "lime", lw=1.2, alpha=0.8)
    ax1.axhline(y=mirror.focal_length, color="yellow", ls="--", lw=1, alpha=0.7, label="Focal plane")
    ax1.set_title(f"Meridian rays\n({'adaptive' if adaptive else 'uniform'}-{matching})")
    ax1.set_xlabel("x (mm)")
    ax1.set_ylabel("y (mm)")
    ax1.legend(loc="upper right", fontsize=8)
    ax1.set_aspect("equal")
    ax1.grid(True, alpha=0.2)

    # (Bottom left) Surface vs true parabola
    ax2 = fig.add_subplot(gs[1, 0])
    ax2.plot(r_prof, y_prof, "b-", lw=1.2, label="Arc surface")
    r_true = np.linspace(-mirror.max_r, mirror.max_r, 200)
    y_true = r_true**2 / (4.0 * mirror.focal_length)
    ax2.plot(r_true, y_true, "r--", lw=1.0, alpha=0.6, label="True parabola")
    ax2.set_title("Surface profile vs true parabola")
    ax2.set_xlabel("x (mm)")
    ax2.set_ylabel("y (mm)")
    ax2.legend(fontsize=8)
    ax2.grid(True, alpha=0.2)

    # (Right) Focal spot
    ax3 = fig.add_subplot(gs[:, 1])
    ax3.scatter(hx, hz, s=1.5, c="yellow", alpha=0.6, edgecolors="none")
    ax3.set_title(
        f"Focal spot (8000 rays)\n"
        f"95% diameter = {2.0 * np.percentile(rh, 95):.2f} mm"
    )
    ax3.set_xlabel("x (mm)")
    ax3.set_ylabel("z (mm)")
    ax3.set_aspect("equal")
    ax3.grid(True, alpha=0.2)

    plt.tight_layout()
    outfile = f"ray_trace_omega_{'adaptive' if adaptive else 'uniform'}_{matching}.png"
    plt.savefig(outfile, dpi=150, bbox_inches="tight")
    plt.close(fig)
    print(f"Figure saved as {outfile}")
    print("=" * 70)


if __name__ == "__main__":
    main()
