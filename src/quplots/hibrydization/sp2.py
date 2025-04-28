from ..electron.waveFunction import normalized_wf3D
import numpy as np
def sp2():
    Spsi = normalized_wf3D(2,0,0)
    Pxpsi= normalized_wf3D(2,1,1)
    Pzpsi= normalized_wf3D(2,1,0)
    Pypsi= normalized_wf3D(2,1,-1)
    PSI2_1=1/np.sqrt(3)*(Spsi+np.sqrt(2)*Pzpsi)
    PSI2_2=1/np.sqrt(3)*(Spsi-np.sqrt(1/2)*Pxpsi-np.sqrt(3/2)*Pzpsi)
    PSI2_3=1/np.sqrt(3)*(Spsi-np.sqrt(1/2)*Pxpsi+np.sqrt(3/2)*Pypsi)
    return PSI2_1,PSI2_2,PSI2_3