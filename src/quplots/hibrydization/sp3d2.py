import numpy as np
from ..electron.waveFunction import *
def sp3d2():
    Spsi =  normalized_wf3D(3,0,0)
    Pxpsi=  normalized_wf3D(3,1,1)
    Pypsi = normalized_wf3D(3,1,-1)
    Pzpsi = normalized_wf3D(3,1,0)
    Dz2psi = normalized_wf3D(3,2,0)
    Dx2y2psi = normalized_wf3D(3,2,2)
    PSIP3D2_1=np.sqrt(1/6)*(Spsi+np.sqrt(3)*Pzpsi+np.sqrt(2)*Dz2psi)
    PSIP3D2_2=np.sqrt(1/6)*(Spsi+np.sqrt(3)*Pzpsi-np.sqrt(1/2)*Dz2psi+np.sqrt(3/2)*Dx2y2psi)
    PSIP3D2_3=np.sqrt(1/6)*(Spsi+np.sqrt(3)*Pzpsi-np.sqrt(1/2)*Dz2psi-np.sqrt(3/2)*Dx2y2psi)
    PSIP3D2_4=np.sqrt(1/6)*(Spsi+np.sqrt(3)*Pzpsi-np.sqrt(1/2)*Dz2psi+np.sqrt(3/2)*Dx2y2psi)
    PSIP3D2_5=np.sqrt(1/6)*(Spsi+np.sqrt(3)*Pzpsi-np.sqrt(1/2)*Dz2psi-np.sqrt(3/2)*Dx2y2psi)
    PSIP3D2_6=np.sqrt(1/6)*(Spsi-np.sqrt(3)*Pzpsi+np.sqrt(2)*Dz2psi)
    return PSIP3D2_1,PSIP3D2_2,PSIP3D2_3,PSIP3D2_4,PSIP3D2_5,PSIP3D2_6