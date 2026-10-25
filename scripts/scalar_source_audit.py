"""Regenerate the scalar-source/detector controls, not an experimental signal."""

import json
from pathlib import Path
import platform
import sys

import numpy as np
import scipy

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))

from scalar_gate import (
    LocalizedGaugePulse, analytic_open_divergence, diagnose_potentials,
    probe_gauge_integral, retarded_path_potential,
)


def run_audit():
    time, position = 1.1, [0.3, -0.2, 0.4]
    cases = {}
    for name, closed, charges in (
        ("closed_current_path", True, True),
        ("open_path_with_endpoint_charge", False, True),
        ("incomplete_open_path_without_charge", False, False),
    ):
        potential = lambda t, x: retarded_path_potential(t, x, closed=closed, include_endpoint_charge=charges)
        cases[name] = diagnose_potentials(potential, time, position)
    pulse = LocalizedGaugePulse()
    gauge_control = diagnose_potentials(pulse.potential, 0.8, position)
    gauge_control["analytic_G"] = pulse.analytic_g(0.8, position)
    return {
        "status": "conditional analytic source/detector result with finite numerical controls; not a new UBT action or measured data",
        "environment": {"python": platform.python_version(), "numpy": np.__version__, "scipy": scipy.__version__},
        "numerical_units": "dimensionless c=mu0=epsilon0=1; probe integrals per unit transverse area, not SI voltages",
        "retarded_path_sample": {"time": time, "position": position, "cases": cases,
                                 "exact_open_path_boundary_div_A": analytic_open_divergence(time, position)},
        "time_dependent_maxwell_gauge_control": gauge_control,
        "conserved_probe_gauge_control": probe_gauge_integral(80),
        "probe_quadrature_control_order_64": probe_gauge_integral(64),
        "conditional_prediction": {
            "extra_scalar_voltage_v": 0,
            "assumptions": [
                "complete locally conserved source current, including feeds and surface/endpoint charges",
                "retarded solution, zero initial scalar data and no incoming scalar wave",
                "ordinary conserved detector, no additional theory-derived coupling",
            ],
            "frameworks": "Maxwell has no independent local G observable; the specified lambda=1 extended AB model also has no forced G here",
            "not_claimed": "No universal no-go for all UBT scalar completions, no full composite-UBT gauge reduction, no removal of homogeneous scalar solutions, no rejection of Aharonov-Bohm holonomy",
        },
        "source_snapshots": [
            {"repository": "UBT-Institute/unified-biquaternion-theory", "commit": "fc4c6834b9f052a21a4054418024764ac046fa6e",
             "path": "canonical/ACTION.en.md", "blob_sha": "5a063201f4a13f4f5cf37c77d6729b52ae4ca1a0",
             "used_for": "fundamental action remains DEFINED_FAMILY_NOT_FINALIZED"},
            {"repository": "UBT-Institute/unified-biquaternion-theory", "commit": "fc4c6834b9f052a21a4054418024764ac046fa6e",
             "path": "canonical/interactions/qed.tex", "blob_sha": "dc8683c29ca8f530b974a78eae31125bfc185230",
             "used_for": "declared effective local U(1) symmetry and conserved current only; not endorsement of the entire file"},
            {"repository": "DavJ/research", "commit": "49e78d7a389ee7d3e502a9f036b02dc7a67b8e6b",
             "path": "audits/scalar_g/README.md", "blob_sha": "1f018e262f54cf66305880ee82c77e211f09769e",
             "used_for": "existing restricted scalar-energy audit; not rerun or replaced here"},
        ],
        "verification_scope": "The continuum Green-function, integration-by-parts and Noether arguments are written analytically in the companion document. Finite tests do not prove those general statements. No Lean formalization is claimed.",
    }


def main():
    report = run_audit()
    output = ROOT / "docs/results/scalar-source-gate.json"
    output.parent.mkdir(parents=True, exist_ok=True)
    output.write_text(json.dumps(report, indent=2, allow_nan=False) + "\n", encoding="utf-8")
    print(f"Wrote {output.relative_to(ROOT)}")
    for name, case in report["retarded_path_sample"]["cases"].items():
        print(name, "G:", case["G"])
    print("Complete probe interaction change:", report["conserved_probe_gauge_control"]["complete_interaction_change"])


if __name__ == "__main__":
    main()
