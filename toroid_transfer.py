"""Quasistatic toroid/pickup transfer with loading and capacitive feedthrough.

The finite rectangular core uses a continuous, uniformly wound toroid model.
The circuit then solves coupled inductors with backreaction, winding losses,
a Thevenin source and explicit capacitances. It is separate from the empty
GHz cavity model. Component values in reference_circuit are assumptions,
not measurements or a model of a new scalar field. Phasors use exp(+i*omega*t)
and peak amplitudes. See docs/toroid-transfer.md for the circuit and limits.
"""

import argparse
from dataclasses import asdict, dataclass
import json
from math import copysign, isfinite, log1p, pi, sqrt

import numpy as np
from scipy.constants import mu_0


def _finite(value, name, *, positive=False):
    value = float(value)
    if not isfinite(value) or (value <= 0 if positive else value < 0):
        adjective = "positive" if positive else "nonnegative"
        raise ValueError(f"{name} must be finite and {adjective}")
    return value


def _turns(value, *, signed=False):
    if (
        isinstance(value, (bool, np.bool_)) or not isinstance(value, (int, np.integer))
        or (value == 0 if signed else value <= 0)
    ):
        raise ValueError("turns must be a nonzero integer" if signed else "turns must be a positive integer")
    return int(value)


@dataclass(frozen=True)
class RectangularToroid:
    inner_radius_m: float = 0.09
    outer_radius_m: float = 0.11
    axial_height_m: float = 0.01
    drive_turns: int = 100
    relative_permeability: float = 1.0

    def __post_init__(self):
        for name in ("inner_radius_m", "outer_radius_m", "axial_height_m", "relative_permeability"):
            object.__setattr__(self, name, _finite(getattr(self, name), name, positive=True))
        object.__setattr__(self, "drive_turns", _turns(self.drive_turns))
        if self.outer_radius_m <= self.inner_radius_m:
            raise ValueError("outer_radius_m must exceed inner_radius_m")

    @property
    def permeance_h(self):
        """Geometric coefficient P: Phi=P*N*I, L=P*N^2, M=P*N*Np."""
        return mu_0 * self.relative_permeability * self.axial_height_m / (2 * pi) * log1p(
            (self.outer_radius_m - self.inner_radius_m) / self.inner_radius_m
        )

    @property
    def flux_per_ampere_wb_a(self):
        """Flux through the cross-section per primary ampere; not N*Phi."""
        return self.permeance_h * self.drive_turns

    @property
    def primary_inductance_h(self):
        return self.permeance_h * self.drive_turns ** 2

    def linked_mutual_h(self, pickup_turns=10):
        """A pickup linking the full cross-section; sign gives dot polarity."""
        return self.permeance_h * self.drive_turns * _turns(pickup_turns, signed=True)


@dataclass
class TransferResponse:
    frequencies_hz: np.ndarray
    primary_current_a: np.ndarray
    pickup_current_a: np.ndarray
    primary_voltage_v: np.ndarray
    output_voltage_v: np.ndarray
    source_current_a: np.ndarray
    source_power_w: np.ndarray
    resistor_power_w: np.ndarray
    relative_power_errors: np.ndarray


@dataclass(frozen=True)
class CoupledCircuit:
    primary_inductance_h: float
    pickup_inductance_h: float
    mutual_inductance_h: float
    primary_resistance_ohm: float = 0.55
    pickup_resistance_ohm: float = 0.2
    source_resistance_ohm: float = 50.0
    load_resistance_ohm: float | None = 1e6
    primary_capacitance_f: float = 50e-12
    load_capacitance_f: float = 100e-12
    feedthrough_capacitance_f: float = 1e-12

    def __post_init__(self):
        for name in ("primary_inductance_h", "pickup_inductance_h", "source_resistance_ohm"):
            object.__setattr__(self, name, _finite(getattr(self, name), name, positive=True))
        for name in (
            "primary_resistance_ohm", "pickup_resistance_ohm", "primary_capacitance_f",
            "load_capacitance_f", "feedthrough_capacitance_f",
        ):
            object.__setattr__(self, name, _finite(getattr(self, name), name))
        if self.load_resistance_ohm is not None:
            object.__setattr__(self, "load_resistance_ohm", _finite(
                self.load_resistance_ohm, "load_resistance_ohm", positive=True,
            ))
        mutual = float(self.mutual_inductance_h)
        bound = sqrt(self.primary_inductance_h) * sqrt(self.pickup_inductance_h)
        if not isfinite(mutual) or abs(mutual) > bound * (1 + 8 * np.finfo(float).eps):
            raise ValueError("mutual inductance requires |M| <= sqrt(L1*L2) for passive magnetic energy")
        # Clamp only last-bit roundoff at perfect coupling, preserving PSD.
        object.__setattr__(self, "mutual_inductance_h", copysign(min(abs(mutual), bound), mutual))

    @property
    def signed_coupling(self):
        return self.mutual_inductance_h / sqrt(self.primary_inductance_h) / sqrt(self.pickup_inductance_h)

    def response(self, frequencies_hz, source_peak_v=1.0):
        """Solve the two winding equations and both terminal KCL equations.

        I1,I2 enter dotted winding terminals. V1,V2 are measured from their
        dots to the common return. Cf joins V1 to V2; Cp and Rload||Cload
        go to the return. The source is Vg in series with Rs. None for Rload
        means an open resistive load; Cload may still load the pickup.
        At DC power uses V*I; sinusoidal frequencies use Re(V*I*)/2.
        """
        frequencies = np.atleast_1d(np.asarray(frequencies_hz, dtype=float))
        if frequencies.ndim != 1 or frequencies.size == 0 or not np.all(np.isfinite(frequencies)) or np.any(frequencies < 0):
            raise ValueError("frequencies_hz must be a nonempty finite nonnegative scalar or 1D array")
        source = complex(source_peak_v)
        if not isfinite(source.real) or not isfinite(source.imag):
            raise ValueError("source_peak_v must be finite")
        if np.any(frequencies == 0) and source.imag != 0:
            raise ValueError("a DC source voltage must be real")
        omega = 2 * pi * frequencies
        jw = 1j * omega
        load_g = 0 if self.load_resistance_ohm is None else 1 / self.load_resistance_ohm
        yp = jw * self.primary_capacitance_f
        yl = load_g + jw * self.load_capacitance_f
        yf = jw * self.feedthrough_capacitance_f
        matrix = np.zeros((len(frequencies), 4, 4), dtype=complex)
        matrix[:, 0, 0] = self.primary_resistance_ohm + jw * self.primary_inductance_h
        matrix[:, 1, 1] = self.pickup_resistance_ohm + jw * self.pickup_inductance_h
        matrix[:, 0, 1] = matrix[:, 1, 0] = jw * self.mutual_inductance_h
        matrix[:, 0, 2] = matrix[:, 1, 3] = -1
        matrix[:, 2, 0] = matrix[:, 3, 1] = 1
        matrix[:, 2, 2] = 1 / self.source_resistance_ohm + yp + yf
        matrix[:, 3, 3] = yl + yf
        matrix[:, 2, 3] = matrix[:, 3, 2] = -yf
        rhs = np.zeros((len(frequencies), 4, 1), dtype=complex)
        rhs[:, 2, 0] = source / self.source_resistance_ohm
        solution = np.linalg.solve(matrix, rhs)[..., 0]
        i1, i2, v1, v2 = solution.T
        # KCL avoids subtracting nearly equal source/port voltages near an
        # antiresonance, where that subtraction degrades the real-power check.
        source_current = i1 + yp * v1 + yf * (v1 - v2)
        average = np.where(frequencies == 0, 1.0, 0.5)
        source_power = average * np.real(source * np.conj(source_current))
        resistor_power = average * (
            self.source_resistance_ohm * abs(source_current) ** 2
            + self.primary_resistance_ohm * abs(i1) ** 2
            + self.pickup_resistance_ohm * abs(i2) ** 2 + load_g * abs(v2) ** 2
        )
        power_errors = abs(source_power - resistor_power) / np.maximum(
            np.maximum(abs(source_power), abs(resistor_power)), 1e-30,
        )
        return TransferResponse(frequencies, i1, i2, v1, v2, source_current, source_power, resistor_power, power_errors)


SCENARIOS = ("linked", "control", "balanced-pair", "mismatched-pair", "sum-pair")


def reference_circuit(scenario="linked", *, pickup_turns=10, load_resistance_ohm=1e6, feedthrough_capacitance_f=1e-12):
    """Illustrative air-core examples, with explicit assumed aggregate R,L,C.

    Pair cases use two series-opposed primary windings and two local pickups
    in series. Cross-core inductance and secondary cross-coupling are ignored.
    The 1% mismatch is a 1% height increase in the second toroid. Cp and Cf
    in pair cases are assumed aggregate capacitances, not series reductions.
    """
    if scenario not in SCENARIOS:
        raise ValueError(f"unknown scenario: {scenario}")
    turns = _turns(pickup_turns)
    geometry = RectangularToroid()
    l1 = geometry.primary_inductance_h
    mutual = geometry.linked_mutual_h(turns)
    l2 = 10e-6 * (turns / 10) ** 2
    r1, r2 = 0.55, 0.2 * turns / 10
    if scenario == "control":
        mutual = 0
    elif scenario in ("balanced-pair", "mismatched-pair", "sum-pair"):
        second = RectangularToroid(axial_height_m=0.0101 if scenario == "mismatched-pair" else 0.01)
        l1 += second.primary_inductance_h
        mutual += (1 if scenario == "sum-pair" else -1) * second.linked_mutual_h(turns)
        l2 *= 2
        r1 *= 2
        r2 *= 2
    return CoupledCircuit(
        l1, l2, mutual, primary_resistance_ohm=r1, pickup_resistance_ohm=r2,
        load_resistance_ohm=load_resistance_ohm, feedthrough_capacitance_f=feedthrough_capacitance_f,
    )


def phasor_json(value):
    value = complex(value)
    return {"real": value.real, "imaginary": value.imag, "peak_magnitude": abs(value), "phase_degrees": float(np.angle(value, deg=True))}


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--scenario", choices=SCENARIOS, default="linked")
    parser.add_argument("--frequency", type=float, default=10000, help="Hz")
    parser.add_argument("--source-peak-v", type=float, default=1, help="Thevenin open-circuit source peak voltage")
    parser.add_argument("--pickup-turns", type=int, default=10)
    parser.add_argument("--load-ohms", type=float, default=1e6)
    parser.add_argument("--feed-capacitance-pf", type=float, default=1)
    args = parser.parse_args(argv)
    try:
        circuit = reference_circuit(
            args.scenario, pickup_turns=args.pickup_turns, load_resistance_ohm=args.load_ohms,
            feedthrough_capacitance_f=args.feed_capacitance_pf * 1e-12,
        )
        result = circuit.response(args.frequency, args.source_peak_v)
    except (ValueError, np.linalg.LinAlgError) as exc:
        parser.error(str(exc))
    report = {
        "model": "linear quasistatic lumped toroid/pickup circuit; not the GHz cavity",
        "status": "assumed component parameters, no measured data or new scalar coupling",
        "scenario": args.scenario, "frequency_hz": args.frequency, "source_peak_v": args.source_peak_v,
        "parameters": asdict(circuit), "signed_magnetic_coupling": circuit.signed_coupling,
        "primary_current_a": phasor_json(result.primary_current_a[0]),
        "output_voltage_v": phasor_json(result.output_voltage_v[0]),
        "source_average_power_w": float(result.source_power_w[0]),
        "resistor_average_power_w": float(result.resistor_power_w[0]),
        "relative_power_error": float(result.relative_power_errors[0]),
    }
    print(json.dumps(report, indent=2, allow_nan=False))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
