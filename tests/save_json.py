from quplots import electron,plots

e = electron(1, 0, 0, 1.0)
e.safe_wavefunction_3D_json("example.json")

p=plots()
p.plot_wf_3d(
    e,
    title="",
    colorscale='rdbu',
    opacity=0.7,
    surface_count=8,
    reversescale=True,
    showscale=False,
    width=700,
    height=700
)