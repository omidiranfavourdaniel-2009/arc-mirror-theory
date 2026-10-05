#!/usr/bin/env python3
"""First-order thermal-expansion sensitivity reference.

This is deliberately labelled a sensitivity model, not physical validation.
The final 2 m optical baseline is 0.01215 mrad; the 0.037 mrad value belongs
to the separate 200 mm reference case and is not used here.
==================================================
Arc-Mirror Theory – Thermal Analysis
==================================================
Worst-case steady-sun thermal analysis for a scaled parabolic dish.

Computes:
  - Front-to-back temperature gradients
  - Thermal-bending figure & slope errors
  - Focal-length drift
  - Concentration retention under solar loading

Generates:
  1. thermal_omega_plots.png   (6-panel dashboard)
  2. thermal_omega_report.txt  (full text summary)

Author: Omidiran Favour Daniel
"""

import sys
import numpy as np
from scipy.integrate import solve_ivp
import matplotlib.pyplot as plt
from dataclasses import dataclass, field
from typing import List, Dict

# ---------------------------------------------------------------------------
# MATERIAL & REFERENCE CONSTANTS
# ---------------------------------------------------------------------------
RHO_AL = 2700.0            # kg/m³
C_AL = 896.0               # J/(kg·K)
K_AL = 167.0               # W/(m·K)
ALPHA_AL = 23.6e-6         # 1/K  (CTE)
SIGMA = 5.670374419e-8     # Stefan-Boltzmann

# Reference floors from Arc-Mirror Theory (200 mm, f/1.5, adaptive-left n=120)
ARC_THEORY_FLOOR_NM = 16.0
ARC_THEORY_SLOPE_RAD = 0.037e-3
CNC_POSITIONING_NM = 1000.0


# ---------------------------------------------------------------------------
# DATA CLASSES
# ---------------------------------------------------------------------------

@dataclass
class Geometry:
    diameter_mm: float = 2000.0
    f_number: float = 1.5
    wall_thickness_mm: float = 20.0

    def __post_init__(self):
        self.D = self.diameter_mm / 1000.0
        self.f = self.D * self.f_number
        self.t = self.wall_thickness_mm / 1000.0
        self.A_aperture = np.pi * (self.D / 2.0) ** 2
        self.A_surface_total = 2.0 * self.A_aperture
        self.mass = RHO_AL * self.A_surface_total * self.t


@dataclass
class Environment:
    label: str
    T_amb_C: float
    I_solar: float
    h_front: float
    h_back: float
    reflectance: float
    eps_front: float = 0.08
    eps_back: float = 0.08

    def __post_init__(self):
        self.T_amb_K = self.T_amb_C + 273.15
        self.q_abs = self.I_solar * (1.0 - self.reflectance)


# ---------------------------------------------------------------------------
# LAYERED FINITE-DIFFERENCE WALL MODEL
# ---------------------------------------------------------------------------

class LayeredWallModel:
    def __init__(self, geom: Geometry, env: Environment, n_nodes: int = 9):
        self.geom = geom
        self.env = env
        self.n = n_nodes
        self.dx = geom.t / (n_nodes - 1)

    def rhs(self, t, T):
        geom, env = self.geom, self.env
        k, dx, rc = K_AL, self.dx, RHO_AL * C_AL
        dTdt = np.zeros(self.n)

        # Front boundary
        q_cond = k * (T[1] - T[0]) / dx
        q_net = (env.q_abs
                 - env.h_front * (T[0] - env.T_amb_K)
                 - env.eps_front * SIGMA * (T[0]**4 - env.T_amb_K**4))
        dTdt[0] = (q_cond + q_net) / (rc * dx / 2.0)

        # Interior
        for i in range(1, self.n - 1):
            dTdt[i] = k * (T[i+1] - 2.0*T[i] + T[i-1]) / dx**2 / rc

        # Back boundary
        q_cond_b = k * (T[self.n-2] - T[self.n-1]) / dx
        q_net_b = (-env.h_back * (T[self.n-1] - env.T_amb_K)
                   - env.eps_back * SIGMA * (T[self.n-1]**4 - env.T_amb_K**4))
        dTdt[self.n-1] = (q_cond_b + q_net_b) / (rc * dx / 2.0)
        return dTdt

    def simulate(self, duration_s: float, n_eval: int = 600):
        T0 = np.full(self.n, self.env.T_amb_K)
        sol = solve_ivp(self.rhs, (0.0, duration_s), T0,
                        t_eval=np.linspace(0.0, duration_s, n_eval),
                        method="BDF", rtol=1e-9, atol=1e-7)
        return sol.t, sol.y - 273.15

    def steady_state(self, settle_s: float = 7200.0):
        T0 = np.full(self.n, self.env.T_amb_K)
        sol = solve_ivp(self.rhs, (0.0, settle_s), T0, method="BDF",
                        rtol=1e-10, atol=1e-8)
        return sol.y[:, -1] - 273.15


# ---------------------------------------------------------------------------
# METRICS
# ---------------------------------------------------------------------------

def gradient_metrics(geom: Geometry, T_nodes_C):
    dT_grad = T_nodes_C[0] - T_nodes_C[-1]
    dk_thermal = ALPHA_AL * dT_grad / geom.t
    r_max = geom.D / 2.0
    dy_rim_nm = 0.5 * dk_thermal * r_max**2 * 1e9
    slope_thermal_urad = dk_thermal * r_max * 1e6
    return dT_grad, dk_thermal, dy_rim_nm, slope_thermal_urad


def focal_shift_um(geom: Geometry, dT_bulk: float) -> float:
    return geom.f * ALPHA_AL * dT_bulk * 1e6


def concentration_with_error(geom: Geometry, extra_slope_rad: float,
                             base_slope_rad: float = ARC_THEORY_SLOPE_RAD):
    theta_sun = np.radians(32.0 / 60.0)
    image_diameter_mm = geom.f * theta_sun * 1000.0
    image_radius_m = (image_diameter_mm / 2.0) / 1000.0
    C_geo_ideal = geom.A_aperture / (np.pi * image_radius_m**2)
    total_slope = base_slope_rad + extra_slope_rad
    C_geo_real = C_geo_ideal * (theta_sun / (theta_sun + 2.0 * total_slope))**2
    ray_dev_arcmin = np.degrees(2.0 * total_slope) * 60.0
    return C_geo_ideal, C_geo_real, ray_dev_arcmin


# ---------------------------------------------------------------------------
# DUAL LOGGER
# ---------------------------------------------------------------------------

class DualLogger:
    def __init__(self, filename: str):
        self.terminal = sys.stdout
        self.log = open(filename, "w")

    def write(self, message: str):
        self.terminal.write(message)
        self.log.write(message)

    def flush(self):
        self.terminal.flush()
        self.log.flush()


# ---------------------------------------------------------------------------
# MAIN ANALYSIS
# ---------------------------------------------------------------------------

def main():
    sys.stdout = DualLogger("thermal_omega_report.txt")

    geom = Geometry(diameter_mm=2000.0, f_number=1.5, wall_thickness_mm=20.0)

    scenarios = [
        Environment("Nominal (clean, ventilated, mild sun)",
                    T_amb_C=30, I_solar=1000, h_front=10, h_back=10, reflectance=0.90),
        Environment("Hot & Still (dusty, low wind, hot afternoon)",
                    T_amb_C=38, I_solar=1050, h_front=4, h_back=4, reflectance=0.86),
        Environment("Worst Case (trapped face layer + ventilated back)",
                    T_amb_C=42, I_solar=1100, h_front=2.0, h_back=15.0, reflectance=0.80),
    ]

    bulk_drift_1k_nm = geom.f * ALPHA_AL * 1.0 * 1e9

    print("=" * 82)
    print("ARC-MIRROR OMEGA – STEADY-SUN THERMAL ANALYSIS")
    print(f"Geometry : Ø{geom.D:.1f} m | wall = {geom.wall_thickness_mm:.0f} mm | f = {geom.f:.1f} m")
    print("=" * 82)
    print(f"{'Scenario':48s} {'ΔT_grad':>9s} {'Sag err':>9s} {'Δf':>8s} {'C_geo':>10s}")
    print(f"{'':48s} {'(mK)':>9s} {'(nm)':>9s} {'(µm)':>8s} {'(× ideal)':>10s}")
    print("-" * 82)

    results = []
    for env in scenarios:
        model = LayeredWallModel(geom, env, n_nodes=9)
        T_ss = model.steady_state()
        dT_grad, dk, dy_rim_nm, slope_urad = gradient_metrics(geom, T_ss)
        dT_bulk = np.mean(T_ss) - env.T_amb_C
        dfocus = focal_shift_um(geom, dT_bulk)
        C_ideal, C_real, ray_arcmin = concentration_with_error(
            geom, extra_slope_rad=slope_urad * 1e-6)

        results.append(dict(
            env=env, T_ss=T_ss, dT_grad=dT_grad, dk=dk,
            dy_rim_nm=dy_rim_nm, slope_urad=slope_urad,
            dT_bulk=dT_bulk, dfocus=dfocus,
            C_ideal=C_ideal, C_real=C_real, ray_arcmin=ray_arcmin,
        ))

        print(f"{env.label:48s} {dT_grad*1000:9.3f} {dy_rim_nm:9.3f} "
              f"{dfocus:8.2f} {C_real/C_ideal:9.2%}")

    print("-" * 82)
    worst = results[-1]
    print("\n" + "=" * 30 + " DETAILED BREAKDOWN (WORST CASE) " + "=" * 30)
    print(f"  Front-face steady temperature : {worst['T_ss'][0]:.3f} °C "
          f"({worst['T_ss'][0]-worst['env'].T_amb_C:+.2f} K above ambient)")
    print(f"  Back-face steady temperature  : {worst['T_ss'][-1]:.3f} °C")
    print(f"  Front-to-back gradient ΔT     : {worst['dT_grad']*1000:.3f} mK")
    print(f"  Thermal bending curvature     : {worst['dk']:.3e} 1/m")
    print(f"  Extra rim sag error           : {worst['dy_rim_nm']:.2f} nm "
          f"({worst['dy_rim_nm']/ARC_THEORY_FLOOR_NM:.1f}× arc-theory floor)")
    print(f"  Extra slope error             : {worst['slope_urad']:.3f} µrad")
    print(f"  Focal-length shift (defocus)  : {worst['dfocus']:.2f} µm")
    print(f"  Combined ray deviation        : {worst['ray_arcmin']:.3f} arcmin")
    print(f"  Geometric concentration       : {worst['C_real']:.0f}× "
          f"(ideal {worst['C_ideal']:.0f}×, {worst['C_real']/worst['C_ideal']:.1%} retained)")
    print("=" * 82 + "\n")

    # Transient series
    DURATION = 8 * 3600.0
    t_nom, T_nom = LayeredWallModel(geom, scenarios[0], n_nodes=9).simulate(DURATION)
    t_wc, T_wc = LayeredWallModel(geom, scenarios[-1], n_nodes=9).simulate(DURATION)

    grad_nom = T_nom[0] - T_nom[-1]
    grad_wc = T_wc[0] - T_wc[-1]
    dy_nom = 0.5 * (ALPHA_AL * grad_nom / geom.t) * (geom.D / 2)**2 * 1e9
    dy_wc = 0.5 * (ALPHA_AL * grad_wc / geom.t) * (geom.D / 2)**2 * 1e9
    dfocus_nom = focal_shift_um(geom, np.mean(T_nom, axis=0) - scenarios[0].T_amb_C)
    dfocus_wc = focal_shift_um(geom, np.mean(T_wc, axis=0) - scenarios[-1].T_amb_C)

    # 2-D parameter sweeps
    T_amb_range = np.linspace(25, 45, 12)
    I_range = np.linspace(700, 1150, 12)
    sweep_configs = [
        ("Insulated back (gradient collapses)", dict(h_front=4, h_back=0.6, eps_back=0.02)),
        ("Trapped front layer + ventilated back", dict(h_front=2, h_back=15, eps_back=0.08)),
    ]
    sweep_results = {}
    for name, cfg in sweep_configs:
        grid = np.zeros((len(T_amb_range), len(I_range)))
        for i, Ta in enumerate(T_amb_range):
            for j, I in enumerate(I_range):
                env = Environment(name, T_amb_C=Ta, I_solar=I, reflectance=0.85, **cfg)
                model = LayeredWallModel(geom, env, n_nodes=5)
                T_ss = model.steady_state()
                _, _, dy_rim_nm, _ = gradient_metrics(geom, T_ss)
                grid[i, j] = dy_rim_nm
        sweep_results[name] = grid

    # Restore stdout
    sys.stdout = sys.stdout.terminal
    print("[REPORT] Text summary saved to thermal_omega_report.txt")

    # ------------------------------------------------------------------
    # PLOTTING DASHBOARD
    # ------------------------------------------------------------------
    fig, axes = plt.subplots(2, 3, figsize=(16, 9))

    # (a) Transient temperatures
    ax = axes[0, 0]
    ax.plot(t_nom / 60, T_nom[0], color="#2e86ab", lw=2, label="Nominal — front")
    ax.plot(t_nom / 60, T_nom[-1], color="#2e86ab", lw=1, ls="--", label="Nominal — back")
    ax.plot(t_wc / 60, T_wc[0], color="#d1495b", lw=2, label="Worst case — front")
    ax.plot(t_wc / 60, T_wc[-1], color="#d1495b", lw=1, ls="--", label="Worst case — back")
    ax.set_xscale("log")
    ax.set_xlabel("Time (minutes, log scale)")
    ax.set_ylabel("Temperature (°C)")
    ax.set_title("(a) Front/back temperature approach to steady state")
    ax.legend(fontsize=8)
    ax.grid(alpha=0.3)

    # (b) Sag error
    ax = axes[0, 1]
    ax.plot(t_nom / 3600, dy_nom, color="#2e86ab", lw=2, label="Nominal")
    ax.plot(t_wc / 3600, dy_wc, color="#d1495b", lw=2, label="Worst case")
    ax.axhline(ARC_THEORY_FLOOR_NM, color="gray", ls="--", lw=1,
               label=f"Arc floor ({ARC_THEORY_FLOOR_NM:.0f} nm)")
    ax.set_xlabel("Time (hours)")
    ax.set_ylabel("Gradient-bending sag error (nm)")
    ax.set_title("(b) Figure error from thermal gradient (8 h sun)")
    ax.legend(fontsize=8)
    ax.grid(alpha=0.3)

    # (c) Defocus
    ax = axes[0, 2]
    ax.plot(t_nom / 3600, dfocus_nom, color="#2e86ab", lw=2, label="Nominal")
    ax.plot(t_wc / 3600, dfocus_wc, color="#d1495b", lw=2, label="Worst case")
    ax.set_xlabel("Time (hours)")
    ax.set_ylabel("Focal shift Δf (µm)")
    ax.set_title("(c) Defocus over 8 h steady sun")
    ax.legend(fontsize=8)
    ax.grid(alpha=0.3)

    # (d)–(e) Heatmaps
    for ax, name in zip([axes[1, 0], axes[1, 1]],
                         ["Insulated back (gradient collapses)",
                          "Trapped front layer + ventilated back"]):
        grid = sweep_results[name]
        im = ax.pcolormesh(I_range, T_amb_range, grid, shading="auto", cmap="inferno")
        fig.colorbar(im, ax=ax, label="Sag error (nm)")
        ax.set_xlabel("Irradiance (W/m²)")
        ax.set_ylabel("Ambient temp (°C)")
        ax.set_title(f"(d/e) {name}", fontsize=9)

    # (f) Error hierarchy
    ax = axes[1, 2]
    labels = ["Gradient\nfigure err.\n(worst case)", "Arc-approx.\nfloor\n(theory)",
              "CNC\npositioning", "Bulk thermal\ndrift (1 K)"]
    values = [worst["dy_rim_nm"], ARC_THEORY_FLOOR_NM, CNC_POSITIONING_NM, bulk_drift_1k_nm]
    colors = ["#d1495b", "#2e86ab", "#f18f01", "#8d99ae"]
    ax.bar(labels, values, color=colors)
    ax.set_yscale("log")
    ax.set_ylabel("Error magnitude (nm, log scale)")
    ax.set_title("(f) Thermal error vs other sources (2.0 m dish)")
    ax.grid(alpha=0.3, axis="y")

    fig.suptitle(
        f"Arc-Mirror Omega Thermal Analysis — "
        f"{geom.diameter_mm:.0f} mm, {geom.wall_thickness_mm:.0f} mm wall",
        fontsize=13,
    )
    fig.tight_layout(rect=[0, 0, 1, 0.96])
    fig.savefig("thermal_omega_plots.png", dpi=150)
    print("[PLOTS] Saved dashboard to thermal_omega_plots.png")
    print("=" * 70)


if __name__ == "__main__":
    main()
