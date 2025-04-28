import numpy as np
def Cartesian_definition():
    xyz = np.linspace(-10, 10, 51)
    x,y,z = np.meshgrid(xyz, xyz, xyz, sparse=False)
    return x,y,z
def CartesianToSpherical(x,y,z):
    r = np.sqrt(x**2 + y**2 + z**2)
    phi = np.arctan2(y+1e-10, x)
    theta = np.where( np.isclose(r, 0.0), np.zeros_like(r), np.arccos(z/r) )
    return r,phi,theta