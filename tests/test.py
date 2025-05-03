from quplots import electron
from quplots import hybridization
from quplots import plots

# Inserta la ruta absoluta del directorio src/ al path
p=plots()
e1=electron(2,1,0,2)
#p.plot_radial_(e1, color="red", linestyle="--", label="Orbital 2p")

#p.plot_radial_(n=3, l=0, d=12, color="red", linewidth=2)
#p.plot_wf_3d(e1)

#p.plot_spherical_imaginary(e1.compute_imaginary_spherical())
#p.plot_spherical_real(e1.compute_real_spherical())
#p.plot_wf_2d(e1.compute_wavefunction_2D())
#p.plot_wf_3d(e1.compute_wavefunction_3D())

p.plot_sp()