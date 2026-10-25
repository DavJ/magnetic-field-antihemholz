"""Magnetic fields of ideal circular loops in cylindrical coordinates.

All positions are global coordinates; ``axial_offset`` locates a loop's center
on the z axis. B0 is the signed field at that loop's center, in the desired
field units. Lengths must use one consistent unit. Coordinates are scalars.
"""

from math import hypot, isfinite, pi, sqrt

from scipy.special import ellipe, ellipk, ellipkm1


def K(k):
    """Complete elliptic integral of the first kind, with modulus k."""
    return ellipk(k ** 2)


def E(k):
    """Complete elliptic integral of the second kind, with modulus k."""
    return ellipe(k ** 2)


def _finite(value, name):
    value = float(value)
    if not isfinite(value):
        raise ValueError(f"{name} must be finite")
    return value


class AntiHemholtz:
    """Opposite-current loops centered at +axial_offset and -axial_offset.

    The historical class spelling is retained for compatibility. B0 is the
    center field of the positive-offset loop alone, not the pair's center
    field (which is zero).
    """

    def __init__(self, B0, loop_radius, axial_offset):
        self.loop_clockwise = LoopField(B0, loop_radius, axial_offset)
        self.B0 = self.loop_clockwise.B0
        self.loop_radius = self.loop_clockwise.loop_radius
        self.axial_offset = self.loop_clockwise.axial_offset
        self.loop_anticlockwise = LoopField(-self.B0, self.loop_radius, -self.axial_offset)

    def Bz(self, r, z):
        return self.vector(r, z)[1]

    def Br(self, r, z):
        return self.vector(r, z)[0]

    def vector(self, r, z):
        br1, bz1 = self.loop_clockwise.vector(r, z)
        br2, bz2 = self.loop_anticlockwise.vector(r, z)
        return br1 + br2, bz1 + bz2

    def normal_vector(self, r, z):
        """Return the unit meridional vector perpendicular to the field.

        A direction is undefined at field nulls, including the pair's center.
        """
        br, bz = self.vector(r, z)
        norm = hypot(br, bz)
        if norm == 0 or not isfinite(norm):
            raise ValueError("normal is undefined where the magnetic field is zero or non-finite")
        return -bz / norm, br / norm


class LoopField:
    """Field of a filamentary loop of positive radius, centered on the z axis."""

    def __init__(self, B0, loop_radius, axial_offset=0):
        self.B0 = _finite(B0, "B0")
        self.loop_radius = _finite(loop_radius, "loop_radius")
        self.axial_offset = _finite(axial_offset, "axial_offset")
        if self.loop_radius <= 0:
            raise ValueError("loop_radius must be positive")

    def alpha(self, r):
        r = _finite(r, "r")
        if r < 0:
            raise ValueError("cylindrical radius r must be non-negative")
        return r / self.loop_radius

    def beta(self, z):
        return (_finite(z, "z") - self.axial_offset) / self.loop_radius

    def gama(self, r, z):
        """Return local z/r; the historical helper spelling is retained."""
        alpha = self.alpha(r)
        if alpha == 0:
            raise ValueError("z/r is undefined on the axis")
        return self.beta(z) / alpha

    def q(self, r, z):
        return (1 + self.alpha(r)) ** 2 + self.beta(z) ** 2

    def k(self, r, z):
        return sqrt(4 * self.alpha(r) / self.q(r, z))

    def Baxial(self, z):
        """Axial field at global coordinate z, including the loop offset."""
        local_z = _finite(z, "z") - self.axial_offset
        return self.B0 * (self.loop_radius / hypot(self.loop_radius, local_z)) ** 3

    def vector(self, r=0, z=0):
        """Return (Br, Bz), raising ValueError on a current-carrying wire."""
        alpha, beta = self.alpha(r), self.beta(z)
        if self.B0 == 0:
            return 0.0, 0.0
        if alpha == 0:
            return 0.0, self.Baxial(z)

        if alpha < 1e-3:
            # The elliptic expression for Br subtracts nearly equal terms
            # here. Use the current-free axis expansion instead:
            # Bz = b - r^2 b''/4 + r^4 b''''/64 + O(r^6),
            # Br = -r b'/2 + r^3 b'''/16 + O(r^5), with b = Baxial.
            inverse_distance = 1 / hypot(1, beta)
            axis_field = self.B0 * inverse_distance ** 3
            u = beta * inverse_distance
            v = alpha * inverse_distance
            bz = axis_field * (
                1 + 0.75 * v ** 2 * (1 - 5 * u ** 2)
                + 45 / 64 * v ** 4 * (1 - 14 * u ** 2 + 21 * u ** 4)
            )
            br = axis_field * v * u * (1.5 + 15 / 16 * v ** 2 * (3 - 7 * u ** 2))
            return br, bz

        # Compute 1-m directly: subtracting m from 1 loses precision near the
        # wire, and can round m to 1 even for a finite off-wire field.
        delta = (1 - alpha) ** 2 + beta ** 2
        if delta == 0:
            raise ValueError("magnetic field is singular on the ideal loop wire")
        q = (1 + alpha) ** 2 + beta ** 2
        complement = delta / q
        first = ellipkm1(complement)
        second = ellipe(1 - complement)
        factor = self.B0 / (pi * sqrt(q))
        bz = factor * (second * (1 - alpha ** 2 - beta ** 2) / delta + first)
        br = factor * beta / alpha * (
            second * (1 + alpha ** 2 + beta ** 2) / delta - first
        )
        return float(br), float(bz)

    def Bz(self, r, z):
        return self.vector(r, z)[1]

    def Br(self, r, z):
        return self.vector(r, z)[0]
