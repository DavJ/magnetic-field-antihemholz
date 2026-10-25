"""Regenerate curved-Maxwell gauge controls and toroid energy comparisons."""

from dataclasses import asdict
import json
from pathlib import Path
import platform
import sys

import numpy as np
import scipy
from scipy.constants import G as NEWTON_G, c, mu_0, physical_constants
from scipy.integrate import quad

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))

from gravity_gate import (
    de_sitter_gauge_divergence, de_sitter_metric,
    diagnose_covariant_potential, energy_equivalent_scales,
    polynomial_gauge_covector,
)
from toroid_transfer import RectangularToroid, phasor_json, reference_circuit


def run_audit():
    q, hubble = [0.7, 0.3, -0.2, 0.4], 0.4
    pure = diagnose_covariant_potential(polynomial_gauge_covector, de_sitter_metric, q)
    magnetic = lambda point: np.array([0, -0.6 * point[2], 0.6 * point[1], 0])
    transformed = lambda point: magnetic(point) + polynomial_gauge_covector(point)
    before = diagnose_covariant_potential(magnetic, de_sitter_metric, q)
    after = diagnose_covariant_potential(transformed, de_sitter_metric, q)
    exact_g = de_sitter_gauge_divergence(q[0], hubble)

    toroid, circuit = RectangularToroid(), reference_circuit("linked")
    response = circuit.response(10000, source_peak_v=1)
    i1, i2 = response.primary_current_a[0], response.pickup_current_a[0]
    v1, v2 = response.primary_voltage_v[0], response.output_voltage_v[0]
    net_ampere_turns = toroid.drive_turns * i1 + 10 * i2
    mean_core = toroid.permeance_h * abs(net_ampere_turns) ** 2 / 4
    leakage_inductance = circuit.pickup_inductance_h - toroid.permeance_h * 10 ** 2
    mean_leakage = leakage_inductance * abs(i2) ** 2 / 4
    mean_magnetic = (circuit.primary_inductance_h * abs(i1) ** 2
                     + circuit.pickup_inductance_h * abs(i2) ** 2
                     + 2 * circuit.mutual_inductance_h * (i1 * i2.conjugate()).real) / 4
    mean_electric = (circuit.primary_capacitance_f * abs(v1) ** 2
                     + circuit.load_capacitance_f * abs(v2) ** 2
                     + circuit.feedthrough_capacitance_f * abs(v1 - v2) ** 2) / 4
    dc_energy = toroid.primary_inductance_h / 2
    dc_integral = quad(lambda r: (mu_0 * toroid.drive_turns / (2 * np.pi * r)) ** 2
                       / (2 * mu_0) * 2 * np.pi * r * toroid.axial_height_m,
                       toroid.inner_radius_m, toroid.outer_radius_m, epsabs=1e-20)[0]
    snapshots = [
        ("canonical/ACTION.en.md", "5a063201f4a13f4f5cf37c77d6729b52ae4ca1a0",
         "Full unique microscopic action is not finalized; no G coupling inferred from an unspecified action."),
        ("canonical/gr_closure/gr_recovery_status.yaml", "84ccee77175689b2ace46b8fd08b8f01979bcb15",
         "Local classical gravity is marked CLOSED after adoption of its postulate; full action and microscopic source reduction remain OPEN."),
        ("canonical/gr_closure/gravity_dynamical_postulate.en.md", "0f6ff8693ba4557aaa8614c0ae5853e304f430ad",
         "One-coupling split-jet Palatini/HT local Einstein-Lambda branch. Its Lambda is integration data, not an identification with electromagnetic G."),
        ("canonical/CANONICAL_DEFINITIONS.md", "1b9e057d3f631fdbcce49bee9d6576e6f15cdcd6",
         "Central tetrad metric definition only; not endorsement of unrelated older representation claims."),
        ("canonical/interactions/qed.tex", "dc8683c29ca8f530b974a78eae31125bfc185230",
         "Declared effective U(1) and Maxwell photon term only; not a derived full composite-Theta reduction."),
        ("canonical/t_munu/step3_einstein_with_matter.tex", "ad9ab0205968cccde43936dacd813da4236ac3cf",
         "Independent-metric Hilbert variation with an inserted Einstein-Hilbert action; this does not derive a physical G source from a Theta-only action."),
    ]
    return {
        "status": "Analytical Maxwell gauge obstruction with finite controls; UBT scalar-gravity coupling is not yet derived; no experimental data.",
        "environment": {"python": platform.python_version(), "numpy": np.__version__, "scipy": scipy.__version__},
        "constants_used_si": {"c_m_s": c, "newton_G_m3_kg_s2": NEWTON_G, "mu0_N_A2": mu_0,
                              "newton_G_standard_uncertainty": physical_constants["Newtonian constant of gravitation"][2]},
        "curved_gauge_control": {
            "units": "dimensionless c=mu0=1; signature (-,+,+,+); A_mu=(-phi,A)",
            "background": "de Sitter flat slicing, a(t)=exp(H*t); prescribed background, not a backreaction solve",
            "coordinates": q, "H": hubble, "background_R": 12 * hubble ** 2,
            "chi": "t^2+(x^2+y^2+z^2)/2, local gauge function",
            "analytic_G": exact_g, "pure_gauge": pure,
            "step_refinement_absolute_errors": {
                str(h): abs(diagnose_covariant_potential(polynomial_gauge_covector, de_sitter_metric, q, step=h)["G"] - exact_g)
                for h in (0.02, 0.01, 0.005)
            },
            "nonzero_maxwell_field_before": before, "same_field_after_gauge": after,
            "max_abs_stress_change": float(np.max(np.abs(np.array(after["stress_energy_cov"]) - before["stress_energy_cov"]))),
        },
        "toroid_energy_comparison_si": {
            "geometry": asdict(toroid),
            "primary_inductance_h": toroid.primary_inductance_h,
            "dc_current_a": 1,
            "dc_magnetic_energy_j": dc_energy,
            "dc_independent_field_volume_integral_j": dc_integral,
            "dc_energy_equivalent_scales": energy_equivalent_scales(dc_energy),
            "loaded_ac_example": {
                "frequency_hz": 10000, "source_peak_v": 1, "pickup_turns": 10,
                "circuit_parameters": asdict(circuit),
                "primary_current_a": phasor_json(i1), "pickup_current_a": phasor_json(i2),
                "mean_core_magnetic_energy_j": mean_core, "mean_pickup_leakage_energy_j": mean_leakage,
                "mean_total_inductive_energy_j": mean_magnetic,
                "core_plus_leakage_minus_total_j": mean_core + mean_leakage - mean_magnetic,
                "mean_explicit_capacitive_energy_j": mean_electric,
                "mean_lumped_em_energy_scales": energy_equivalent_scales(mean_magnetic + mean_electric),
            },
            "one_joule_reference": energy_equivalent_scales(1),
            "scope": "Assumed continuous air-core winding and lumped components. U/c^2 and the 1 m point-mass comparison are energy scales, not a conserved apparatus source or predicted accelerometer signal. Leads, supports, generator, heat and total energy transfer are not modeled. Internal battery-to-field energy transfer does not add total rest energy to an isolated device.",
        },
        "ubt_source_snapshot": {
            "repository": "UBT-Institute/unified-biquaternion-theory",
            "commit": "fc4c6834b9f052a21a4054418024764ac046fa6e",
            "files": [{"path": p, "blob_sha": sha, "used_for": use} for p, sha, use in snapshots],
            "independently_verified_scope": "Read pinned source content and scope statements. The GR closure proofs/verifiers were not rerun here. No UBT action or metric was modified.",
        },
        "not_claimed": [
            "No proof that all UBT completions lack an additional physical scalar.",
            "No gravitational G source or detector gain derived for this toroid.",
            "No stability endorsement or independent verification of the cosmological extended-electromagnetism proposal.",
            "No general continuum theorem proved by finite numerical samples or by Lean in this audit.",
        ],
    }


def main():
    report = run_audit()
    output = ROOT / "docs/results/scalar-gravity-audit.json"
    output.parent.mkdir(parents=True, exist_ok=True)
    output.write_text(json.dumps(report, indent=2, allow_nan=False) + "\n", encoding="utf-8")
    print(f"Wrote {output.relative_to(ROOT)}")
    curved = report["curved_gauge_control"]
    print("Curved pure-gauge G:", curved["pure_gauge"]["G"])
    print("Maxwell stress change:", curved["max_abs_stress_change"])
    ac = report["toroid_energy_comparison_si"]["loaded_ac_example"]
    print("Mean lumped EM energy, J:", ac["mean_lumped_em_energy_scales"]["energy_j"])
    print("Energy mass equivalent, kg:", ac["mean_lumped_em_energy_scales"]["energy_mass_equivalent_kg"])


if __name__ == "__main__":
    main()
