from scipy import __version__ as scipy_version
from packaging.version import Version

if Version(scipy_version) >= Version("1.15.0"):
    from scipy.special import sph_harm_y

    def sph_harm_compat(m, n, theta, phi):
        return sph_harm_y(n, m, theta, phi)
else:
    from scipy.special import sph_harm

    def sph_harm_compat(m, n, theta, phi):
        return sph_harm(m, n, theta, phi)
