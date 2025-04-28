import numpy as np
from scipy.special import sph_harm
from .radial import radial   
from .coordinates import * 
def normalized_wf(d,n,l,m):
    r     = np.linspace(0, d, 10000)
    pr = radial(r, n, l)**2 * r**2 * (r[1]-r[0])
    max_r = r[ np.where(np.cumsum(pr) >0.95)[0][0] ]
    # Set coordinates grid to assign a certain probability to each point (x, y) in the plane
    x = y = np.linspace(-max_r, max_r, 501)
    x, y = np.meshgrid(x, y)
    r = np.sqrt((x ** 2 + y ** 2))
    # Ψnlm(r,θ,φ) = Rnl(r).Ylm(θ,φ)
    psi = radial(r, n, l) * sph_harm(m, l, 0, np.arctan(x / (y + 1e-7)))
    psi_sq = np.abs(psi)**2
    return psi_sq,max_r,n,l,m  
def normalized_wf3D(n,l,m):
    x,y,z = Cartesian_definition()
    r,phi,theta = CartesianToSpherical(x,y,z)
    psi = radial(r, n, l) * np.real(sph_harm(m, l, phi, theta))
    return psi
def normalized_wf3D_complex(n,l,m):
    x,y,z = Cartesian_definition()
    r,phi,theta = CartesianToSpherical(x,y,z)
    psi = radial(r, n, l) * sph_harm(m, l, phi, theta)*np.real(sph_harm(m, l, phi, theta))
    return psi