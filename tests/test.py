from quplots import plots
p = plots()
p.plot_sp(
    colorscale='Viridis',
    reversescale=True,
    opacity=0.6,
    surface_count=5,
    lighting=dict(ambient=0.4, diffuse=0.9, specular=0.5, roughness=0.9, fresnel=0.1),
    lightposition=dict(x=100, y=200, z=0),
    showscale=False,
    name="Orbital"
)

p.plot_sp2(colorscale='Plasma',
            reversescale=False,
            opacity=0.7,
            surface_count=5)

p.plot_sp3(colorscale='Cividis',
           reversescale=False,
           opacity=0.7,
           surface_count=5)

p.plot_sp2d(colorscale='Inferno',
            reversescale=False,
            opacity=0.7,
            surface_count=5)

p.plot_sp3d(colorscale='Magma',
            reversescale=False,
            opacity=0.7,
            surface_count=5)

p.plot_sp3d2(colorscale='Turbo',
             reversescale=False,
             opacity=0.7,
             surface_count=5)