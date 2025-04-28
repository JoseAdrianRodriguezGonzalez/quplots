import numpy as np
from ..electron.waveFunction import *
def sp3d():
        Spsi =  normalized_wf3D(3,0,0)
        Pxpsi=  normalized_wf3D(3,1,1)
        Pypsi = normalized_wf3D(3,1,-1)
        Pzpsi = normalized_wf3D(3,1,0)
        Dz2psi = normalized_wf3D(3,2,0)
        PSIP3D_1=np.sqrt(1/3)*(Spsi+np.sqrt(2)*Pzpsi)
        PSIP3D_2=1/np.sqrt(3)*(Spsi-np.sqrt(1/2)*Pxpsi+np.sqrt(3/2)*Pypsi)
        PSIP3D_3=1/np.sqrt(3)*(Spsi-np.sqrt(1/2)*Pxpsi-np.sqrt(3/2)*Pzpsi)
        PSIP3D_4=1/np.sqrt(2)*(Dz2psi+Pzpsi)
        PSIP3D_5=1/np.sqrt(2)*(-Dz2psi+Pzpsi)
        return PSIP3D_1,PSIP3D_2,PSIP3D_3,PSIP3D_4,PSIP3D_5