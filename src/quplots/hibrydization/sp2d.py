import numpy as np 
from ..electron.waveFunction import *
def sp2d():
    Spsi =  normalized_wf3D(3,0,0)
    Pxpsi=  normalized_wf3D(3,1,1)
    Pypsi = normalized_wf3D(3,1,-1)
    Dx2y2psi = normalized_wf3D(3,2,2)
    PSI2D_1=1/2*Spsi+1/np.sqrt(2)*Pxpsi+1/2*Dx2y2psi
    PSI2D_2=1/2*Spsi-1/np.sqrt(2)*Pxpsi+1/2*Dx2y2psi
    PSI2D_3=1/2*Spsi+1/np.sqrt(2)*Pypsi-1/2*Dx2y2psi
    PSI2D_4=1/2*Spsi-1/np.sqrt(2)*Pypsi-1/2*Dx2y2psi
    return PSI2D_1,PSI2D_2,PSI2D_3,PSI2D_4