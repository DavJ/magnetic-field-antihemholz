"""Regenerate the assumed-component transfer examples and optional figure.

Run from the checkout: python scripts/transfer_audit.py [--plot]
The report records model predictions, not experimental data.
"""

import argparse
from dataclasses import asdict, replace
import json
from pathlib import Path
import platform
import sys

import numpy as np
import scipy

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))

from toroid_transfer import RectangularToroid, SCENARIOS, phasor_json, reference_circuit


def run_audit():
    source_peak = 1.0
    frequencies = np.geomspace(100, 1e6, 81)
    points = np.array([1000, 10000, 100000, 1000000])
    cases = {name: reference_circuit(name) for name in SCENARIOS}
    cases["linked-50-ohm"] = reference_circuit(load_resistance_ohm=50)
    cases["balanced-no-feedthrough"] = reference_circuit("balanced-pair", feedthrough_capacitance_f=0)
    cases["balanced-10pF-feedthrough"] = reference_circuit("balanced-pair", feedthrough_capacitance_f=10e-12)
    results = {}
    for name, circuit in cases.items():
        sweep, sample = circuit.response(frequencies, source_peak), circuit.response(points, source_peak)
        results[name] = {
            "components": asdict(circuit), "signed_magnetic_coupling": circuit.signed_coupling,
            "output_voltage_v_real": sweep.output_voltage_v.real.tolist(),
            "output_voltage_v_imaginary": sweep.output_voltage_v.imag.tolist(),
            "primary_current_a_real": sweep.primary_current_a.real.tolist(),
            "primary_current_a_imaginary": sweep.primary_current_a.imag.tolist(),
            "maximum_relative_average_power_error": float(sweep.relative_power_errors.max()),
            "sample_points": [{
                "frequency_hz": float(f), "primary_current_a": phasor_json(i),
                "output_voltage_v": phasor_json(v), "source_average_power_w": float(power),
            } for f, i, v, power in zip(points, sample.primary_current_a, sample.output_voltage_v, sample.source_power_w)],
        }
    geometry = RectangularToroid()
    # First-order independent standard uncertainties for h, a and b, not
    # measured tolerances and not a bound on all experimental artefacts.
    sigma_dimension = 0.0001
    log_ratio = np.log(geometry.outer_radius_m / geometry.inner_radius_m)
    relative_sigma_m = np.sqrt(
        (sigma_dimension / geometry.axial_height_m) ** 2
        + (sigma_dimension / (geometry.inner_radius_m * log_ratio)) ** 2
        + (sigma_dimension / (geometry.outer_radius_m * log_ratio)) ** 2
    )
    sigma_pair_m = np.sqrt(2) * relative_sigma_m * geometry.linked_mutual_h(10)
    null = cases["balanced-pair"]
    step_m = geometry.linked_mutual_h(10) * 1e-5
    plus = replace(null, mutual_inductance_h=step_m).response(10000).output_voltage_v[0]
    minus = replace(null, mutual_inductance_h=-step_m).response(10000).output_voltage_v[0]
    voltage_sigma = abs((plus - minus) / (2 * step_m)) * sigma_pair_m
    report = {
        "status": "assumed-component quasistatic circuit predictions; not measured data or a scalar-field detection",
        "environment": {"python": platform.python_version(), "numpy": np.__version__, "scipy": scipy.__version__},
        "geometry": asdict(geometry),
        "primary_inductance_h": geometry.primary_inductance_h,
        "mutual_inductance_10_turn_pickup_h": geometry.linked_mutual_h(10),
        "source_open_circuit_peak_v": source_peak, "source_resistance_ohm": 50,
        "phasor_convention": "peak amplitudes, exp(+i*omega*t), currents into dotted winding terminals",
        "frequencies_hz": frequencies.tolist(), "cases": results,
        "illustrative_geometry_uncertainty": {
            "assumptions": "independent 0.1 mm one-sigma uncertainties in h,a,b on each toroid; exact turn counts and mu_r=1",
            "single_mutual_relative_standard_uncertainty": float(relative_sigma_m),
            "pair_net_mutual_standard_uncertainty_h": float(sigma_pair_m),
            "pair_complex_output_voltage_standard_uncertainty_v_at_10kHz": float(voltage_sigma),
            "scope": "first-order geometric contribution only; not a noise floor, total error budget or anomaly threshold",
        },
        "limitations": [
            "Continuous primary winding and confined core flux; discrete-turn and lead fields are not solved.",
            "R,L,C values beyond the ideal primary inductance and linked M are assumed, not measured.",
            "Pair examples ignore cross-core mutual inductances; C values are assumed aggregate components.",
            "The sweep is a lumped-model calculation, not a validated hardware bandwidth.",
            "No closed conducting cavity, radiofrequency feed, radiation, ferromagnetic dispersion or new scalar coupling.",
        ],
    }
    return report


def plot(report, output):
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    plt.rcParams.update({"font.size": 10, "svg.hashsalt": "toroid-transfer-audit"})
    fig, axes = plt.subplots(1, 2, figsize=(10, 4.6), layout="constrained")
    frequency = np.asarray(report["frequencies_hz"])
    groups = (
        (("linked", "One toroid, 1 Mohm load", "#176d91"),
         ("linked-50-ohm", "One toroid, 50 ohm load", "#b14e29"),
         ("sum-pair", "Pair with pickup polarity summing", "#39824c")),
        (("balanced-pair", "Balanced pair, Cf = 1 pF", "#176d91"),
         ("balanced-10pF-feedthrough", "Balanced pair, Cf = 10 pF", "#6d5099"),
         ("mismatched-pair", "1% height mismatch, Cf = 1 pF", "#b14e29")),
    )
    for ax, group in zip(axes, groups):
        for key, label, color in group:
            case = report["cases"][key]
            magnitude = np.hypot(case["output_voltage_v_real"], case["output_voltage_v_imaginary"])
            ax.loglog(frequency, magnitude * 1e6, label=label, color=color)
        ax.set(xlabel="Frequency (Hz)", ylabel="Output peak voltage (microvolt)")
        ax.grid(which="both", alpha=0.2)
        ax.legend(fontsize=8, loc="upper left")
    axes[0].set_title("Loaded ordinary-induction signal")
    axes[1].set_title("Ordinary residuals near a magnetic null")
    fig.suptitle("1 V peak Thevenin source, 50 ohm source resistance; assumed R,L,C; no cavity", fontsize=11)
    output.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(output, metadata={"Date": None})
    plt.close(fig)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--plot", action="store_true")
    args = parser.parse_args()
    report = run_audit()
    output = ROOT / "docs/results/toroid-transfer.json"
    output.parent.mkdir(parents=True, exist_ok=True)
    output.write_text(json.dumps(report, indent=2, allow_nan=False) + "\n", encoding="utf-8")
    if args.plot:
        plot(report, ROOT / "docs/figures/toroid-transfer.svg")
    print(f"Wrote {output.relative_to(ROOT)}")
    for key in ("linked", "balanced-pair", "mismatched-pair", "sum-pair"):
        sample = report["cases"][key]["sample_points"][1]
        print(key, "10 kHz output peak (V):", sample["output_voltage_v"]["peak_magnitude"])


if __name__ == "__main__":
    main()
