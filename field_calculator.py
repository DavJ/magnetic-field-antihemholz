from scipy.special import ellipk, ellipe, ellipkm1
from numpy import pi, sqrt, linspace, NaN
from pylab import plot, xlabel, ylabel, suptitle, legend, show


def K(k):
    return ellipk(k ** 2.0)  # Elliptic integral, first kind, as a function of k

def E(k):
    return ellipe(k ** 2.0)  # Elliptic integral, second kind, as a function of k


class AntiHemholtz:

    def __init__(self, B0, loop_radius, axial_offset):
        self.B0 = B0
        self.loop_radius = loop_radius
        self.axial_offset = axial_offset
        self.loop_clockwise = LoopField(B0, loop_radius, axial_offset)
        self.loop_anticlockwise = LoopField(-B0, loop_radius, -axial_offset)

    def Bz(self, r, z):
        return self.loop_clockwise.Bz(r, z) + self.loop_aniclockwise.Bz(r, z)

    def Br(self, r, z):
        return self.loop_clockwise.Br(r, z) + self.loop_anticlockwise.Br(r, z)

    def vector(self, r, z):
        return self.Br(r, z), self.Bz(r, z)

    def normal_vector(self, r, z):
        br, bz = self.vector(r, z)
        norm = sqrt(br**2 + bz **2)
        return bz / norm, -br / norm


class LoopField:

    def __init__(self, B0, loop_radius, axial_offset):
        self.B0 = B0
        self.loop_radius = loop_radius
        self.axial_offset = axial_offset

    def vector(self, r=0, z=0):
        return (self.Br(r,z), self.Bz(r, z))

    def alpha(self, r):
        return r / self.loop_radius

    def beta(self, z):
        return z / self.loop_radius

    def gama(self, r, z):
        return z / r

    def q(self, r, z):
        return (1 + self.alpha(r)) ** 2 + self.beta(z) ** 2

    def k(self, r, z):
        return sqrt(4 * self.alpha(r) / self.q(r, z))

        # On-Axis field = f(current and radius of loop, x of measurement point)
    def Baxial(self, z):
        if self.loop_radius == 0:
            if z == 0:
                return NaN
            else:
                return 0.0
        else:
            return (self.B0 * self.loop_radius ** 3) / (self.loop_radius ** 2 + z ** 2) ** 1.5

    # Axial field component = f(current and radius of loop, r and x of meas. point)
    def Bz(self, r, z):
        if r == 0:
            if z == 0:
                return self.B0  # central field
            else:
                return self.Baxial(z)  # axial field
        else:  # axial component, any location
            return self.B0 * \
                   (E(self.k(r, z)) * ((1.0 - self.alpha(r) ** 2 - self.beta(z) ** 2) / (self.q(r, z) - 4 * self.alpha(r))) + K(
                       self.k(r, z))) \
                   / pi / sqrt(self.q(r, z))

    # Radial field component = f(current and radius of loop, r and x of meas. point)
    def Br(self, r, z):
        if r == 0:
            return 0  # no radial component on axis!
        else:  # radial component, any location other than axis.
            return self.B0 * self.gama(z, r) * \
                   (E(self.k(r, z)) * ((1.0 + self.alpha(r) ** 2 + self.beta(z) ** 2) / (self.q(r, z) - 4 * self.alpha(r))) - K(
                       self.k(r, z))) \
                   / pi / sqrt(self.q(r, z))
