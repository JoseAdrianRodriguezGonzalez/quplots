import numpy as np
from scipy.special import sph_harm
def RealSpherical(l,m):
        #Creade grid of phi and theta angles for ploting surface mesh
        phi, theta = np.linspace(0, np.pi, 100), np.linspace(0, 2*np.pi, 100)
        phi, theta = np.meshgrid(phi, theta)

        #Calcualte spherical harmonic with given m and l
        Y = sph_harm(abs(m), l, theta, phi)
        if m<0:
            Y=np.sqrt(2)*(-1)**m*Y.imag
        elif m>0:
            Y=np.sqrt(2)*(-1)**m*Y.real
        R=abs(Y)
        # Let's normalize color scale
        fcolors    = Y.real
        fmax, fmin = fcolors.max(), fcolors.min()
        fcolors    = (fcolors - fmin)/(fmax - fmin)
        if m>l:
            raise ValueError("l and m must be postive, it means, that l can't be smaller than m")
        return R,phi,theta,fcolors,l,m
def ImaginarySpherical(l,m):
    #Creade grid of phi and theta angles for ploting surface mesh
    phi, theta = np.linspace(0, np.pi, 100), np.linspace(0, 2*np.pi, 100)
    phi, theta = np.meshgrid(phi, theta)

    #Calcualte spherical harmonic with given m and l
    Ylm = sph_harm(m, l, theta, phi)
    R=abs(Ylm)

    # Let's normalize color scale
    fcolors    = Ylm.imag
    fmax, fmin = fcolors.max(), fcolors.min()
    fcolors    = (fcolors - fmin)/(fmax - fmin)
    if m>l:
            raise ValueError("l and m must be postive, it means, that l can't be smaller than m")
    return R,phi,theta,fcolors,l,m