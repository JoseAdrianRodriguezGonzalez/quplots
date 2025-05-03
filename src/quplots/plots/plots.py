import numpy as np
import  matplotlib.pyplot as plt
import plotly.graph_objects as go 
from ..electron.electron import electron
from scipy.special import sph_harm, genlaguerre, factorial
from ..electron.radial import radial
from ..electron.coordinates import *
from ..hibrydization.sp import sp
from ..hibrydization.sp2 import sp2
from ..hibrydization.sp3 import sp3
from ..hibrydization.sp2d import sp2d
from ..hibrydization.sp3d import sp3d
from ..hibrydization.sp3d2 import sp3d2
from ..hibrydization.sp3d import sp3d  
class plots:
    def __init__(self, *args):
        super(plots, self).__init__(*args)
    def plot_radial_(self,elect=None,d=None,n=None,l=None,**plot_kwargs):
        if elect is not None and isinstance(elect,electron):
            R_function=elect.compute_radial()
            r=np.linspace(0,elect.getRadius(),1000)
            d=elect.getRadius()
        else:
            r = np.linspace(0, d, 1000)
            R_function = radial(r, n, l)
        plt.plot(r, r**2 * R_function**2)
        plt.title(f'$R_{{n, l}}(r)$ distancia = {d}')
        plt.xlabel(r'$r [a_0]$')
        plt.ylabel(r'$R_nl(r) r^2$')
        plt.show()
    def plot_spherical_real(self,elect=None,R=None,theta=None,phi=None,fcolors=None,l=None,m=None,**plots_kwargs):
        if elect is not None and isinstance(elect,electron):
            package=elect.compute_real_spherical()
            x=package[0]*np.sin(package[1]) * np.cos(package[2])
            y=package[0]*np.sin(package[1]) * np.sin(package[2])
            z=package[0]*np.cos(package[1])
            surfacecolor=package[3]
            l_val,m_val=package[4],package[5]
        else:
            x=R*np.sin(phi) * np.cos(theta)
            y=R*np.sin(phi) * np.sin(theta)
            z=R*np.cos(phi)
            surfacecolor=fcolors
            l_val,m_val=l,m
        surface_args = {
            "x": x,
            "y": y,
            "z": z,
            "surfacecolor": surfacecolor,
            "colorscale": plots_kwargs.pop("colorscale", "balance"),
            "showscale": plots_kwargs.pop("showscale", False)
        }
        fig=go.Figure(data=[go.Surface(**surface_args)])
        layout_args = {
            "title": plots_kwargs.pop("title", f"<b>Y<sub>{l_val},{m_val}</sub></b>"),
            "autosize": plots_kwargs.pop("autosize", False),
            "width": plots_kwargs.pop("width", 700),
            "height": plots_kwargs.pop("height", 700),
            "margin": plots_kwargs.pop("margin", dict(l=65, r=50, b=65, t=90)),
            "paper_bgcolor": plots_kwargs.pop("paper_bgcolor", "black"),
            "scene": plots_kwargs.pop("scene", dict(
                xaxis=dict(visible=False),
                yaxis=dict(visible=False),
                zaxis=dict(visible=False))),
            "font": plots_kwargs.pop("font", dict(color="white"))
        }
        fig.update_layout(**layout_args,**plots_kwargs)
        fig.show()
    def plot_spherical_imaginary(self,elect=None,R=None,theta=None,phi=None,fcolors=None,l=None,m=None,**plots_kwargs):
        if elect is not None and isinstance(elect,electron):
            package=elect.compute_imaginary_spherical()
            x=package[0]*np.sin(package[1]) * np.cos(package[2])
            y=package[0]*np.sin(package[1]) * np.sin(package[2])
            z=package[0]*np.cos(package[1])
            surfacecolor=package[3]
            l_val,m_val=package[4],package[5]
        else:
            x=R*np.sin(phi) * np.cos(theta)
            y=R*np.sin(phi) * np.sin(theta)
            z=R*np.cos(phi)
            surfacecolor=fcolors
            l_val,m_val=l,m
        surface_args = {
            "x": x,
            "y": y,
            "z": z,
            "surfacecolor": surfacecolor,
            "colorscale": plots_kwargs.pop("colorscale", "balance"),
            "showscale": plots_kwargs.pop("showscale", False)
        }
        fig=go.Figure(data=[go.Surface(**surface_args)])
        layout_args = {
            "title": plots_kwargs.pop("title", f"<b>Y<sub>{l_val},{m_val}</sub></b>"),
            "autosize": plots_kwargs.pop("autosize", False),
            "width": plots_kwargs.pop("width", 700),
            "height": plots_kwargs.pop("height", 700),
            "margin": plots_kwargs.pop("margin", dict(l=65, r=50, b=65, t=90)),
            "paper_bgcolor": plots_kwargs.pop("paper_bgcolor", "black"),
            "scene": plots_kwargs.pop("scene", dict(
                xaxis=dict(visible=False),
                yaxis=dict(visible=False),
                zaxis=dict(visible=False))),
            "font": plots_kwargs.pop("font", dict(color="white"))
        }
        fig.update_layout(**layout_args,**plots_kwargs)
        fig.show()
    def plot_wf_2d(self,elect=None,psi_sq=None,max_r=None,n=None,l=None,m=None,**plot_kwargs):
        fig, ax = plt.subplots()
        package=elect.compute_wavefunction_2D() if(elect is not None and isinstance(elect,electron)) else (psi_sq,max_r,n,l,m) 
        levels = plot_kwargs.pop("levels", 40)
        cmap = plot_kwargs.pop("cmap", "RdBu")
        title_fontsize = plot_kwargs.pop("title_fontsize", 15)
        xlabel = plot_kwargs.pop("xlabel", r"$x$")
        ylabel = plot_kwargs.pop("ylabel", r"$y$")
        title = plot_kwargs.pop("title", r"$|\psi_{{({0}, {1}, {2})}}|^2$".format(package[2], package[3], package[4]))
        ax.contour(package[0], levels=levels, cmap=cmap, extent=[-package[1], package[1], -package[1], package[1]], **plot_kwargs)
        ax.set_title(title, fontsize=title_fontsize)
        ax.set_aspect("equal")
        ax.set_xlabel(xlabel)
        ax.set_ylabel(ylabel)
        ax.set_xlim(-package[1], package[1])
        ax.set_ylim(-package[1], package[1])
        plt.show()
    def plot_wf_3d(self,elect,**kwargs_plot):
        psi= elect.compute_wavefunction_3D() if elect is not None and isinstance(elect,electron) else elect
        x,y,z = Cartesian_definition()
        isosurface_keys = {
        "x", "y", "z", "value", "isomin", "isomax", "opacity", "surface_count", 
        "caps", "colorscale", "lighting", "lightposition", "reversescale", 
        "showscale", "name", "hoverinfo"
        }
        isosurface_kwargs = {k: v for k, v in kwargs_plot.items() if k in isosurface_keys}
        layout_kwargs = {k: v for k, v in kwargs_plot.items() if k not in isosurface_keys}

        fig = go.Figure(data=go.Isosurface(
                x=x.flatten(),
                y=y.flatten(),
                z=z.flatten(),
                value=abs(psi).flatten(),
                isomin=-.75 * abs(psi).min(),
                isomax=.75 * abs(psi).max(),
                opacity=0.5,
                surface_count=6,
                caps=dict(x_show=False, y_show=False, z_show=False),
                colorscale='RdBu',
                **isosurface_kwargs
            ))
        n_val,l_val,m_val=elect.getN(),elect.getL(),elect.getM()
        fig.update_layout(
            title = f"<b>|ψ<sub>{n_val},{l_val},{m_val}</sub>(x,y,z)|</b>",
            paper_bgcolor="black",
            scene=dict(
                xaxis=dict(visible=False),
                yaxis=dict(visible=False),
                zaxis=dict(visible=False)
            ),
            font=dict(color="white"),
            **layout_kwargs
        )
        fig.show()
    def plot_sp(self,**kwargs):
        PSI1, PSI2 = sp()
        packageCartesian=Cartesian_definition()
        layout_kwargs = kwargs.get('layout', {})
        isosurface_kwargs = kwargs.get('isosurface', {})
        # Grafico
        fig = go.Figure()
        for psi in [PSI1,PSI2]:
            fig.add_trace(go.Isosurface(
                x=packageCartesian[0].flatten(),
                y=packageCartesian[1].flatten(),
                z=packageCartesian[2].flatten(),
                value=abs(psi).flatten(),
                #colorscale='RdYlBu_r',
                isomin=-.75 * abs(psi).min(),
                isomax=.75 * abs(psi).max(),
                surface_count=6,
                opacity=0.5,
                caps=dict(x_show=False, y_show=False, z_show=False),
                **isosurface_kwargs
            ))
        fig.update_layout(
            title = layout_kwargs.get("title", "<b>sp: Lineal</b>"),
            paper_bgcolor=layout_kwargs.get('paper_bgcolor', "black"),
            scene=layout_kwargs.get('scene', dict(
                xaxis=dict(visible=False),
                yaxis=dict(visible=False),
                zaxis=dict(visible=False)
            )),
            font=layout_kwargs.get('font', dict(color="white")),
            **{k: v for k, v in layout_kwargs.items() if k not in ['title', 'paper_bgcolor', 'scene', 'font']}
        )
        fig.show()
    def plot_sp2(self,**kwargs):
        PSI2_1, PSI2_2, PSI2_3 = sp2()
        # Grafico
        packageCartesian=Cartesian_definition()
        layout_kwargs = kwargs.get('layout', {})
        isosurface_kwargs = kwargs.get('isosurface', {})
        fig = go.Figure()
        for psi in [PSI2_1,PSI2_2,PSI2_3]:
            fig.add_trace(go.Isosurface(
                x=packageCartesian[0].flatten(),
                y=packageCartesian[1].flatten(),
                z=packageCartesian[2].flatten(),
                value=abs(psi).flatten(),
                #colorscale='RdYlBu_r',
                isomin=-.75 * abs(psi).min(),
                isomax=.75 * abs(psi).max(),
                surface_count=6,
                opacity=0.5,
                caps=dict(x_show=False, y_show=False, z_show=False),
                **isosurface_kwargs
            ))
        fig.update_layout(
            title = layout_kwargs.get("title", "<b>sp<sup>2</sup>: Trigonal Planar</b>"),
            paper_bgcolor=layout_kwargs.get('paper_bgcolor', "black"),
            scene=layout_kwargs.get('scene', dict(
                xaxis=dict(visible=False),
                yaxis=dict(visible=False),
                zaxis=dict(visible=False)
            )),
            font=layout_kwargs.get('font', dict(color="white")),
            **{k: v for k, v in layout_kwargs.items() if k not in ['title', 'paper_bgcolor', 'scene', 'font']}
        )
        fig.show()
    def plot_sp3(self,**kwargs):
        PSI3_1, PSI3_2, PSI3_3, PSI3_4 = sp3()
        packageCartesian=Cartesian_definition()
        layout_kwargs = kwargs.get('layout', {})
        isosurface_kwargs = kwargs.get('isosurface', {})
        fig = go.Figure()
        for psi in [PSI3_1,PSI3_2,PSI3_3,PSI3_4]:
            fig.add_trace(go.Isosurface(
                x=packageCartesian[0].flatten(),
                y=packageCartesian[1].flatten(),
                z=packageCartesian[2].flatten(),
                value=abs(psi).flatten(),
                #colorscale='RdYlBu_r',
                isomin=-.75 * abs(psi).min(),
                isomax=.75 * abs(psi).max(),
                surface_count=6,
                opacity=0.5,
                caps=dict(x_show=False, y_show=False, z_show=False),
                **isosurface_kwargs
            ))
        fig.update_layout(
            title = layout_kwargs.get("title", "<b>sp<sup>3</sup>: Tetrahedral</b>"),
            paper_bgcolor=layout_kwargs.get('paper_bgcolor', "black"),
            scene=layout_kwargs.get('scene', dict(
                xaxis=dict(visible=False),
                yaxis=dict(visible=False),
                zaxis=dict(visible=False)
            )),
            font=layout_kwargs.get('font', dict(color="white")),
            **{k: v for k, v in layout_kwargs.items() if k not in ['title', 'paper_bgcolor', 'scene', 'font']}
        )
        fig.show()
    def plot_sp2d(self,**kwargs):
        PSI2D_1,PSI2D_2,PSI2D_3,PSI2D_4 = sp2d()
        packageCartesian=Cartesian_definition()
        layout_kwargs = kwargs.get('layout', {})
        isosurface_kwargs = kwargs.get('isosurface', {})
        fig = go.Figure()
        for psi in [PSI2D_1,PSI2D_2,PSI2D_3,PSI2D_4]:
            fig.add_trace(go.Isosurface(
                x=packageCartesian[0].flatten(),
                y=packageCartesian[1].flatten(),
                z=packageCartesian[2].flatten(),
                value=abs(psi).flatten(),
                #colorscale='RdYlBu_r',
                isomin=-.75 * abs(psi).min(),
                isomax=.75 * abs(psi).max(),
                surface_count=6,
                opacity=0.5,
                caps=dict(x_show=False, y_show=False, z_show=False),
                **isosurface_kwargs
            ))
        fig.update_layout(
            title = layout_kwargs.get("title", "<b>sp<sup>2</sup>d: Square Planar</b>"),
            paper_bgcolor=layout_kwargs.get('paper_bgcolor', "black"),
            scene=layout_kwargs.get('scene', dict(
                xaxis=dict(visible=False),
                yaxis=dict(visible=False),
                zaxis=dict(visible=False)
            )),
            font=layout_kwargs.get('font', dict(color="white")),
            **{k: v for k, v in layout_kwargs.items() if k not in ['title', 'paper_bgcolor', 'scene', 'font']}
        )
        fig.show()
    def plot_sp3d(self,**kwargs): 
        PSIP3D_1, PSIP3D_2, PSIP3D_3, PSIP3D_4, PSIP3D_5 = sp3d()
        packageCartesian=Cartesian_definition()
        layout_kwargs = kwargs.get('layout', {})
        isosurface_kwargs = kwargs.get('isosurface', {})
        fig = go.Figure()
        for psi in [PSIP3D_1, PSIP3D_2, PSIP3D_3, PSIP3D_4, PSIP3D_5]:    
            fig.add_trace(go.Isosurface(
                x=packageCartesian[0].flatten(),
                y=packageCartesian[1].flatten(),
                z=packageCartesian[2].flatten(),
                value=abs(psi).flatten(),
                #colorscale='RdYlBu_r', #colorscale='Portland',
                isomin=-.75 * abs(psi).min(),
                isomax=.75 * abs(psi).max(),
                surface_count=6,
                opacity=0.5,
                caps=dict(x_show=False, y_show=False, z_show=False),
                **isosurface_kwargs
            ))
        fig.update_layout(
            title = layout_kwargs.get("title", "<b>sp<sup>3</sup>d: Trigonal Bipyramidal</b>"),
            paper_bgcolor=layout_kwargs.get('paper_bgcolor', "black"),
            scene=layout_kwargs.get('scene', dict(
                xaxis=dict(visible=False),
                yaxis=dict(visible=False),
                zaxis=dict(visible=False)
            )),
            font=layout_kwargs.get('font', dict(color="white")),
            **{k: v for k, v in layout_kwargs.items() if k not in ['title', 'paper_bgcolor', 'scene', 'font']}
        )
        fig.show()
    def plot_sp3d2(self,**kwargs):
        PSIP3D2_1,PSIP3D2_2,PSIP3D2_3,PSIP3D2_4,PSIP3D2_5,PSIP3D2_6 = sp3d2()
        packageCartesian=Cartesian_definition()
        layout_kwargs = kwargs.get('layout', {})
        isosurface_kwargs = kwargs.get('isosurface', {})
        fig = go.Figure()
        for psi in [PSIP3D2_1,PSIP3D2_2,PSIP3D2_3,PSIP3D2_4,PSIP3D2_5,PSIP3D2_6]:
            fig.add_trace(go.Isosurface(
                x=packageCartesian[0].flatten(),
                y=packageCartesian[1].flatten(),
                z=packageCartesian[2].flatten(),
                value=abs(psi).flatten(),
                #colorscale='RdYlBu_r', #colorscale='Portland',
                isomin=-.75 * abs(psi).min(),
                isomax=.75 * abs(psi).max(),
                surface_count=6,
                opacity=0.5,
                caps=dict(x_show=False, y_show=False, z_show=False),
                **isosurface_kwargs
            ))
        fig.update_layout(
            title = layout_kwargs.get("title", "<b>sp<sup>3</sup>d<sup>2</sup>: Octahedral</b>"),
            paper_bgcolor=layout_kwargs.get('paper_bgcolor', "black"),
            scene=layout_kwargs.get('scene', dict(
                xaxis=dict(visible=False),
                yaxis=dict(visible=False),
                zaxis=dict(visible=False)
            )),
            font=layout_kwargs.get('font', dict(color="white")),
            **{k: v for k, v in layout_kwargs.items() if k not in ['title', 'paper_bgcolor', 'scene', 'font']}
        )
        fig.show()      