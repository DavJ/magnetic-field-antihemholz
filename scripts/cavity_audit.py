"""Regenerate the numerical audit and its figure from the repository root.

Run: python scripts/cavity_audit.py
Optional figure: python scripts/cavity_audit.py --plot (requires .[report]).
Outputs are fixed repository paths under docs/results and docs/figures.
"""

import argparse
from math import exp, pi
from pathlib import Path
import json
import platform
import sys

import numpy as np
import scipy
from scipy.constants import c, mu_0
from scipy.optimize import brentq
from scipy.special import spherical_jn

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))

from cavity_modes import cylindrical_tm_modes, meridian_mesh, solve_modes
from cavity_profile import HISTORICAL_SOURCE, integral_profile, open_integral
from toroidal_potential import ThinToroidPotential


def summary(result):
    return {
        "radial_cells": result.mesh.radial_cells,
        "axial_cells": result.mesh.axial_cells,
        "frequencies_hz": result.frequencies_hz.tolist(),
        "max_algebraic_relative_residual": float(result.relative_residuals.max()),
        "electric_magnetic_energy_ratios": result.electric_magnetic_energy_ratios.tolist(),
    }


def run_audit():
    coil_radius, outer_radius, z0 = 0.1, 0.09, 0.1
    exact = cylindrical_tm_modes(0.1, 0.2, 6)
    exact_f = np.array([row["frequency_hz"] for row in exact])
    cylinder = []
    for nr in (8, 16, 32, 64):
        r = np.linspace(0, 0.1, nr + 1)
        result = solve_modes(meridian_mesh(r, np.full_like(r, 0.1), 2 * nr), 6)
        row = summary(result)
        row["relative_frequency_errors"] = (result.frequencies_hz / exact_f - 1).tolist()
        cylinder.append(row)
    sphere_root = brentq(lambda x: spherical_jn(1, x) + x * spherical_jn(1, x, derivative=True), 2, 3)
    sphere_exact = c * sphere_root / (2 * pi * 0.1)
    sphere = []
    for nr in (16, 32, 64, 128):
        r = np.linspace(0, 0.1 * (1 - 1 / (4 * nr ** 2)), nr + 1)
        result = solve_modes(meridian_mesh(r, np.sqrt(0.1 ** 2 - r ** 2), 2 * nr), 1)
        row = summary(result)
        row["equatorial_truncation_over_radius"] = 1 / (4 * nr ** 2)
        row["relative_frequency_error_to_full_sphere"] = float(result.frequencies_hz[0] / sphere_exact - 1)
        sphere.append(row)
    profiles = {}
    finest = None
    for name, historical in (("historical_literal", True), ("corrected", False)):
        rows = []
        for nr in (12, 24, 48, 96):
            r = np.linspace(0, outer_radius, nr + 1)
            height = integral_profile(r, z0, coil_radius, historical=historical)
            result = solve_modes(meridian_mesh(r, height, 2 * nr), 4)
            row = summary(result)
            if rows:
                row["relative_change_from_previous_mesh"] = (
                    result.frequencies_hz / np.array(rows[-1]["frequencies_hz"]) - 1
                ).tolist()
            rows.append(row)
        profiles[name] = rows
        if not historical:
            finest = result
    # Separately check integration error on a curved-wall mesh.
    r = np.linspace(0, outer_radius, 25)
    mesh = meridian_mesh(r, integral_profile(r, z0, coil_radius), 48)
    q3 = solve_modes(mesh, 4, quadrature_order=3)
    q5 = solve_modes(mesh, 4, quadrature_order=5)
    sample_r = np.array([0, 0.25, 0.5, 0.75, 0.9]) * coil_radius
    # A smooth static gauge leaves the axis target and fields unchanged.
    # Dimensionless R=Phi=1 is used only for this mathematical control.
    potential = ThinToroidPotential(1, 1)
    az = lambda r, z: potential.vector(r, z)[1]
    chi = lambda r, z: 0.2 * r ** 2 * z * exp(-(r ** 2 + z ** 2) / 2)
    target = open_integral(az, 0, 1)
    top = float(integral_profile(np.array([0, 0.9]), 1, 1)[-1])
    shifted_top = brentq(lambda z: open_integral(az, 0.9, z) + chi(0.9, z) - target, 0, 1)
    # Ordinary low-frequency reference: ideal continuous toroidal winding.
    inner, outer, height, turns, current, drive_hz, pickup_turns = 0.09, 0.11, 0.01, 100, 0.1, 10000, 1
    flux = mu_0 * turns * current * height / (2 * pi) * np.log(outer / inner)
    report = {
        "status": "numerical Maxwell baseline; not experimental validation or a physical scalar source",
        "environment": {"python": platform.python_version(), "numpy": np.__version__, "scipy": scipy.__version__},
        "historical_source": HISTORICAL_SOURCE,
        "geometry": {
            "coil_radius_m": coil_radius, "axis_half_height_m": z0, "closing_wall_radius_m": outer_radius,
            "construction": "mirror the positive integral profile about z=0 and close at r=0.09 m; empty PEC cavity",
            "source_placement": "reference flux ring is outside the closing wall; eigenproblem includes no actual source or port",
        },
        "profile_samples": {
            "r_m": sample_r.tolist(),
            "historical_literal_half_height_m": integral_profile(sample_r, z0, coil_radius, historical=True).tolist(),
            "corrected_half_height_m": integral_profile(sample_r, z0, coil_radius).tolist(),
        },
        "cylinder_analytic_modes": exact,
        "cylinder_convergence": cylinder,
        "sphere_analytic_fundamental": {"radius_m": 0.1, "ka": sphere_root, "frequency_hz": sphere_exact},
        "sphere_domain_and_mesh_convergence": sphere,
        "profile_convergence": profiles,
        "quadrature_check_relative_frequency_change_3_to_5": (q5.frequencies_hz / q3.frequencies_hz - 1).tolist(),
        "open_integral_gauge_control": {
            "chi_over_flux": "0.2*(r/R)^2*(z/R)*exp(-((r/R)^2+(z/R)^2)/2)",
            "r_over_R": 0.9, "target_over_flux": target,
            "original_Z_over_R": top, "gauge_shifted_Z_over_R": shifted_top,
            "physical_E_B_change": "zero; static smooth gauge, no hardware change",
        },
        "ordinary_induction_reference": {
            "assumptions": "air/vacuum core, ideal continuous toroidal winding; quasistatic; ignores leads and leakage",
            "inner_radius_m": inner, "outer_radius_m": outer, "axial_height_m": height,
            "drive_turns": turns, "current_peak_A": current, "frequency_hz": drive_hz,
            "pickup_turns_linking_one_core": pickup_turns,
            "core_flux_peak_Wb": float(flux),
            "pickup_emf_peak_V": float(2 * pi * drive_hz * pickup_turns * flux),
        },
    }
    return report, finest


def plot(report, result, output):
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    from matplotlib.tri import Triangulation

    plt.rcParams.update({"font.size": 10, "svg.hashsalt": "maxwell-cavity-audit"})
    fig, axes = plt.subplots(1, 2, figsize=(10, 4.3), layout="constrained")
    r = np.linspace(0, 0.09, 97)
    for historical, label, color in (
        (True, "2016 formula (k passed as m)", "#b55a30"),
        (False, "Corrected elliptic parameter", "#146e8c"),
    ):
        z = integral_profile(r, 0.1, 0.1, historical=historical)
        axes[0].plot(100 * r, 100 * z, label=label, color=color)
        axes[0].plot(100 * r, -100 * z, color=color)
        axes[0].plot([9, 9], [-100 * z[-1], 100 * z[-1]], color=color)
    axes[0].set(title="Explicit closed comparison geometries", xlabel="r (cm)", ylabel="z (cm)")
    axes[0].legend(fontsize=8, loc="upper left")
    axes[0].grid(alpha=0.2)
    nodes, cells = result.mesh.nodes, result.mesh.cells
    triangles = np.concatenate((cells[:, [0, 1, 2]], cells[:, [0, 2, 3]]))
    tri = Triangulation(100 * nodes[:, 0], 100 * nodes[:, 1], triangles)
    im = axes[1].tricontourf(tri, result.relative_h_phi(), levels=np.linspace(0, 1, 13), cmap="viridis")
    axes[1].set(
        title=f"Corrected cavity: first m=0 TM mode\n{result.frequencies_hz[0] / 1e9:.5f} GHz",
        xlabel="r (cm)", ylabel="z (cm)",
    )
    fig.colorbar(im, ax=axes[1], label="Relative H_phi amplitude")
    fig.suptitle("Empty vacuum PEC cavity; chosen side closure; no toroid, ports or losses", fontsize=11)
    output.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(output, metadata={"Date": None})
    plt.close(fig)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--plot", action="store_true")
    args = parser.parse_args()
    report, result = run_audit()
    output = ROOT / "docs/results/cavity-audit.json"
    output.parent.mkdir(parents=True, exist_ok=True)
    output.write_text(json.dumps(report, indent=2, allow_nan=False) + "\n", encoding="utf-8")
    if args.plot:
        plot(report, result, ROOT / "docs/figures/cavity-audit.svg")
    print(f"Wrote {output.relative_to(ROOT)}")
    print("Corrected cavity frequencies (GHz):", result.frequencies_hz / 1e9)


if __name__ == "__main__":
    main()
