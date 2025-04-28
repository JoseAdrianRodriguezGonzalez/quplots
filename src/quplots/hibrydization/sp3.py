import numpy as np
from ..electron.waveFunction import *
def sp3():
    Pxpsi=  normalized_wf3D_complex(2,1,1)
    Spsi =  normalized_wf3D_complex(2,0,0)
    Pypsi = normalized_wf3D_complex(2,1,-1)
    Pzpsi = normalized_wf3D_complex(2,1,0)
    PSI3_1=1/2*(Spsi+Pxpsi+Pypsi+Pzpsi)
    PSI3_2=1/2*(Spsi+Pxpsi-Pypsi-Pzpsi)
    PSI3_3=1/2*(Spsi-Pxpsi-Pypsi+Pzpsi)
    PSI3_4=1/2*(Spsi-Pxpsi+Pypsi-Pzpsi)
    return PSI3_1,PSI3_2,PSI3_3,PSI3_4