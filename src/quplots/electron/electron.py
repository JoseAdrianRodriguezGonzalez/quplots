import numpy as np
from .coordinates import Cartesian_definition, CartesianToSpherical
from .radial import radial
from .sphericalHarmonics import RealSpherical, ImaginarySpherical
from .waveFunction import normalized_wf, normalized_wf3D, normalized_wf3D_complex
from ..utils.json_utils import save_wavefunction


class electron:
    def __init__(self, n, l, m, radius):
        if n <= 0 or l < 0 or l >= n:
            raise ValueError(f"Valores inválidos para n y l: n={n}, l={l}")
        if abs(m) > l:
            raise ValueError(f"m debe estar entre -l y l: l={l}, m={m}")
        self.n = n
        self.l = l
        self.m = m
        self.radius = radius

    def get_cartesian_grid(self, limit_inferior=-10, limit_superior=10, granularity=51):
        return Cartesian_definition(limit_inferior, limit_superior, granularity)

    def get_spherical_grid(self, x, y, z):
        return CartesianToSpherical(x, y, z)

    def compute_radial(self):
        return radial(self.radius, self.n, self.l)

    def compute_real_spherical(self):
        return RealSpherical(self.l, self.m)

    def compute_imaginary_spherical(self):
        return ImaginarySpherical(self.l, self.m)

    def compute_wavefunction_2D(self, d=2000):
        return normalized_wf(d, self.n, self.l, self.m)

    def compute_wavefunction_3D(self):
        return normalized_wf3D(self.n, self.l, self.m)

    def safe_wavefunction_3D_json(self, output, cartesian_options=(-10, 10, 51)):
        x, y, z = self.get_cartesian_grid(*cartesian_options)
        psi = self.compute_wavefunction_3D()
        try:
            save_wavefunction(x, y, z, psi, output)
            print(f"file saved in {output}")
        except Exception as e:
            print(f"Error in {e}")

    def compute_wavefunction_3D_complex(self):
        return normalized_wf3D_complex(self.n, self.l, self.m)

    def getN(self):
        return self.n

    def getL(self):
        return self.l

    def getM(self):
        return self.m

    def getRadius(self):
        return self.radius

