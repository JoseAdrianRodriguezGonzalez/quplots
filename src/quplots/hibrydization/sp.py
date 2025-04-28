import numpy as np
from ..electron.waveFunction import normalized_wf3D
def sp():
        Spsi = normalized_wf3D(2,0,0)
        Pzpsi= normalized_wf3D(2,1,0)
        PSI1=1/np.sqrt(2)*(Spsi-Pzpsi)
        PSI2=1/np.sqrt(2)*(Spsi+Pzpsi)
        return PSI1,PSI2