import numpy as np
import  matplotlib.pyplot as plt
import plotly.graph_objects as go
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
    def plot_radial_(self,d,n,l):
        r = np.linspace(0, d, 1000)
        R_function = radial(r, n, l)
        #plt.plot(r, Rnl)
        #plt.plot(r, Rnl**2)
        plt.plot(r, r**2 * R_function**2)
        plt.title(f'$R_{{n, l}}(r)$ distancia = {d}')
        plt.xlabel(r'$r [a_0]$')
        plt.ylabel(r'$R_nl(r) r^2$')
        plt.show()
    def plot_spherical_real(self,package=None,R=None,theta=None,phi=None,fcolors=None,l=None,m=None):
        if package:
            fig = go.Figure(data=[go.Surface(x=package[0]*np.sin(package[1]) * np.cos(package[2]),
                                     y=package[0]*np.sin(package[1]) * np.sin(package[2]),
                                     z=package[0]*np.cos(package[1]),
                                     surfacecolor=package[3],
                                     colorscale='balance')])

    # Show the plot
            fig.update_layout(title=fr'$Y_{package[4], package[5]}$', autosize=False,
                            width=700, height=700,
                            margin=dict(l=65, r=50, b=65, t=90),paper_bgcolor="black",
                            scene=dict(
                            xaxis=dict(visible=False),
                            yaxis=dict(visible=False),
                            zaxis=dict(visible=False)),
                            font=dict(
                            color="white"))
            fig.show()
        else:
            fig = go.Figure(data=[go.Surface(x=R*np.sin(phi) * np.cos(theta),
                                     y=R*np.sin(phi) * np.sin(theta),
                                     z=R*np.cos(phi),
                                     surfacecolor=fcolors,
                                     colorscale='balance')])

    # Show the plot
            fig.update_layout(title=fr'$Y_{l, m}$', autosize=False,
                            width=700, height=700,
                            margin=dict(l=65, r=50, b=65, t=90),paper_bgcolor="black",
                            scene=dict(
                            xaxis=dict(visible=False),
                            yaxis=dict(visible=False),
                            zaxis=dict(visible=False)),
                            font=dict(
                            color="white"))
            fig.show()
    def plot_spherical_imaginary(self,package=None,R=None,theta=None,phi=None,fcolors=None,l=None,m=None):

        if package:
            fig = go.Figure(data=[go.Surface(x=package[0]*np.sin(package[1]) * np.cos(package[2]),
                                     y=package[0]*np.sin(package[1]) * np.sin(package[2]),
                                     z=package[0]*np.cos(package[1]),
                                     surfacecolor=package[3],
                                     colorscale='balance')])

    # Show the plot
            fig.update_layout(title=fr'$Y_{package[4], package[5]}$', autosize=False,
                            width=700, height=700,
                            margin=dict(l=65, r=50, b=65, t=90),paper_bgcolor="black",
                            scene=dict(
                            xaxis=dict(visible=False),
                            yaxis=dict(visible=False),
                            zaxis=dict(visible=False)),
                            font=dict(
                            color="white"))
            fig.show()
        else:
            fig = go.Figure(data=[go.Surface(x=R*np.sin(phi) * np.cos(theta),
                                     y=R*np.sin(phi) * np.sin(theta),
                                     z=R*np.cos(phi),
                                     surfacecolor=fcolors,
                                     colorscale='balance')])

    # Show the plot
            fig.update_layout(title=fr'$Y_{l, m}$', autosize=False,
                            width=700, height=700,
                            margin=dict(l=65, r=50, b=65, t=90),paper_bgcolor="black",
                            scene=dict(
                            xaxis=dict(visible=False),
                            yaxis=dict(visible=False),
                            zaxis=dict(visible=False)),
                            font=dict(
                            color="white"))
            fig.show()
    def plot_wf_2d(self,package=None,psi_sq=None,max_r=None,n=None,l=None,m=None):
        fig, ax = plt.subplots()
        if package:
            ax.contour(package[0], 40, cmap='RdBu', extent=[-package[1], package[1],-package[1], package[1]])
            ax.set_title(r'$|\psi_{{({0}, {1}, {2})}}|^2$'.format(package[2], package[3], package[4]), fontsize=15)
            ax.set_aspect('equal')
            ax.set_xlabel(r'$x$')
            ax.set_ylabel(r'$y$')
            ax.set_xlim(-package[1], package[1])
            ax.set_ylim(-package[1], package[1])
            plt.show()
        else:
            ax.contour(psi_sq, 40, cmap='RdBu', extent=[-max_r, max_r,-max_r, max_r])
            ax.set_title(r'$|\psi_{{({0}, {1}, {2})}}|^2$'.format(n, l, m), fontsize=15)
            ax.set_aspect('equal')
            ax.set_xlabel(r'$x$')
            ax.set_ylabel(r'$y$')
            ax.set_xlim(-max_r, max_r)
            ax.set_ylim(-max_r, max_r)
            plt.show()
    def plot_wf_3d(self,psi):
        x,y,z = Cartesian_definition()
        fig= go.Figure(data=go.Isosurface(
        x=x.flatten(),
        y=y.flatten(),
        z=z.flatten(),
        value=abs(psi).flatten(),
        colorscale='RdBu',
        isomin=-.75*abs(psi).min(),
        isomax=.75*abs(psi).max(),
        surface_count=6,
        opacity=0.5,
        caps=dict(x_show=False,y_show=False,z_show=False)
        ))
        fig.update_layout(paper_bgcolor="black",scene=dict(
                        xaxis=dict(visible=False),
                        yaxis=dict(visible=False),
                        zaxis=dict(visible=False)),
                        font=dict(
                        color="white"))
        fig.show()
    def plot_sp(self,packageCartesian=None,x=None,y=None,z=None):
        if packageCartesian:
            PSI1, PSI2 = sp()
            # Grafico
            fig = go.Figure()
            Iso1 = go.Isosurface(
                x=packageCartesian[0].flatten(),
                y=packageCartesian[1].flatten(),
                z=packageCartesian[2].flatten(),
                value=abs(PSI1).flatten(),
                colorscale='RdYlBu_r',
                isomin=-.75 * abs(PSI1).min(),
                isomax=.75 * abs(PSI1).max(),
                surface_count=6,
                opacity=0.5,
                caps=dict(x_show=False, y_show=False, z_show=False)
            )
            Iso2 = go.Isosurface(
                x=packageCartesian[0].flatten(),
                y=packageCartesian[1].flatten(),
                z=packageCartesian[2].flatten(),
                value=abs(PSI2).flatten(),
                colorscale='RdYlBu_r',
                isomin=-.75 * abs(PSI2).min(),
                isomax=.75 * abs(PSI2).max(),
                surface_count=6,
                opacity=0.5,
                caps=dict(x_show=False, y_show=False, z_show=False)
            )
            fig.update_layout(
            title=fr'$sp:  Lineal$',paper_bgcolor="black",
            scene=dict(
                xaxis=dict(visible=False),
                yaxis=dict(visible=False),
                zaxis=dict(visible=False)
            ),
            font=dict(color="white")
            )
            fig.add_trace(Iso1)
            fig.add_trace(Iso2)
            fig.show()
        else:    
            PSI1, PSI2 = sp()
            # Grafico
            fig = go.Figure()
            Iso1 = go.Isosurface(
                x=x.flatten(),
                y=y.flatten(),
                z=z.flatten(),
                value=abs(PSI1).flatten(),
                colorscale='RdYlBu_r',
                isomin=-.75 * abs(PSI1).min(),
                isomax=.75 * abs(PSI1).max(),
                surface_count=6,
                opacity=0.5,
                caps=dict(x_show=False, y_show=False, z_show=False)
            )
            Iso2 = go.Isosurface(
                x=x.flatten(),
                y=y.flatten(),
                z=z.flatten(),
                value=abs(PSI2).flatten(),
                colorscale='RdYlBu_r',
                isomin=-.75 * abs(PSI2).min(),
                isomax=.75 * abs(PSI2).max(),
                surface_count=6,
                opacity=0.5,
                caps=dict(x_show=False, y_show=False, z_show=False)
            )
            fig.update_layout(
                title=fr'$sp:  Lineal$',paper_bgcolor="black",
                scene=dict(
                    xaxis=dict(visible=False),
                    yaxis=dict(visible=False),
                    zaxis=dict(visible=False)
                ),
                font=dict(color="white")
            )
            fig.add_trace(Iso1)
            fig.add_trace(Iso2)
            fig.show()
    def plot_sp2(self,packageCartesian=None,x=None,y=None,z=None):
        PSI2_1, PSI2_2, PSI2_3 = sp2()
        # Grafico
        if packageCartesian:
            fig = go.Figure()
            Iso1 = go.Isosurface(
                x=packageCartesian[0].flatten(),
                y=packageCartesian[1].flatten(),
                z=packageCartesian[2].flatten(),
                value=abs(PSI2_1).flatten(),
                colorscale='RdYlBu_r',
                isomin=-.75 * abs(PSI2_1).min(),
                isomax=.75 * abs(PSI2_1).max(),
                surface_count=4,
                opacity=0.3,
                caps=dict(x_show=False, y_show=False, z_show=False)
            )
            Iso2 = go.Isosurface(
                x=packageCartesian[0].flatten(),
                y=packageCartesian[1].flatten(),
                z=packageCartesian[2].flatten(),
                value=abs(PSI2_2).flatten(),
                colorscale='RdYlBu_r',
                isomin=-.75 * abs(PSI2_2).min(),
                isomax=.75 * abs(PSI2_2).max(),
                surface_count=4,
                opacity=0.3,
                caps=dict(x_show=False, y_show=False, z_show=False)
            )
            Iso3 = go.Isosurface(
                x=packageCartesian[0].flatten(),
                y=packageCartesian[1].flatten(),
                z=packageCartesian[2].flatten(),
                value=abs(PSI2_3).flatten(),
                colorscale='RdYlBu_r',
                isomin=-.75 * abs(PSI2_3).min(),
                isomax=.75 * abs(PSI2_3).max(),
                surface_count=4,
                opacity=0.3,
                caps=dict(x_show=False, y_show=False, z_show=False)
            )
            fig.update_layout(
                title=fr'$sp^2: Trigonal Planar$',
                template='plotly_dark',
                autosize=False,
                width=700,
                height=700,
                margin=dict(l=65, r=50, b=65, t=90),
                scene=dict(
                    xaxis=dict(showgrid=False, visible=False),
                    yaxis=dict(showgrid=False, visible=False),
                    zaxis=dict(showgrid=False, visible=False)
                )
            )
            fig.add_trace(Iso1)
            fig.add_trace(Iso2)
            fig.add_trace(Iso3)
            fig.show()
        else:
            fig = go.Figure()
            Iso1 = go.Isosurface(
                x=x.flatten(),
                y=y.flatten(),
                z=z.flatten(),
                value=abs(PSI2_1).flatten(),
                colorscale='RdYlBu_r',
                isomin=-.75 * abs(PSI2_1).min(),
                isomax=.75 * abs(PSI2_1).max(),
                surface_count=4,
                opacity=0.3,
                caps=dict(x_show=False, y_show=False, z_show=False)
            )
            Iso2 = go.Isosurface(
                x=x.flatten(),
                y=y.flatten(),
                z=z.flatten(),
                value=abs(PSI2_2).flatten(),
                colorscale='RdYlBu_r',
                isomin=-.75 * abs(PSI2_2).min(),
                isomax=.75 * abs(PSI2_2).max(),
                surface_count=4,
                opacity=0.3,
                caps=dict(x_show=False, y_show=False, z_show=False)
            )
            Iso3 = go.Isosurface(
                x=x.flatten(),
                y=y.flatten(),
                z=z.flatten(),
                value=abs(PSI2_3).flatten(),
                colorscale='RdYlBu_r',
                isomin=-.75 * abs(PSI2_3).min(),
                isomax=.75 * abs(PSI2_3).max(),
                surface_count=4,
                opacity=0.3,
                caps=dict(x_show=False, y_show=False, z_show=False)
            )
            fig.update_layout(
                title=fr'$sp^2: Trigonal Planar$',
                template='plotly_dark',
                autosize=False,
                width=700,
                height=700,
                margin=dict(l=65, r=50, b=65, t=90),
                scene=dict(
                    xaxis=dict(showgrid=False, visible=False),
                    yaxis=dict(showgrid=False, visible=False),
                    zaxis=dict(showgrid=False, visible=False)
                )
            )
            fig.add_trace(Iso1)
            fig.add_trace(Iso2)
            fig.add_trace(Iso3)
            fig.show()
    def plot_sp3(self,packageCartesian=None,x=None,y=None,z=None):
        PSI3_1, PSI3_2, PSI3_3, PSI3_4 = sp3()
        if packageCartesian:
            # Grafico
            fig = go.Figure()
            Iso1 = go.Isosurface(
                x=packageCartesian[0].flatten(),
                y=packageCartesian[1].flatten(),
                z=packageCartesian[2].flatten(),
                value=abs(PSI3_1).flatten(),
                colorscale='RdYlBu_r',
                isomin=-.75 * abs(PSI3_1).min(),
                isomax=.75 * abs(PSI3_1).max(),
                surface_count=6,
                opacity=0.6,
                caps=dict(x_show=False, y_show=False, z_show=False)
            )
            Iso2 = go.Isosurface(
                x=packageCartesian[0].flatten(),
                y=packageCartesian[1].flatten(),
                z=packageCartesian[2].flatten(),
                value=abs(PSI3_2).flatten(),
                colorscale='RdYlBu_r',
                isomin=-.75 * abs(PSI3_2).min(),
                isomax=.75 * abs(PSI3_2).max(),
                surface_count=6,
                opacity=0.3,
                caps=dict(x_show=False, y_show=False, z_show=False)
            )
            Iso3 = go.Isosurface(
                x=packageCartesian[0].flatten(),
                y=packageCartesian[1].flatten(),
                z=packageCartesian[2].flatten(),
                value=abs(PSI3_3).flatten(),
                colorscale='RdYlBu_r',
                isomin=-.75 * abs(PSI3_3).min(),
                isomax=.75 * abs(PSI3_3).max(),
                surface_count=6,
                opacity=0.3,
                caps=dict(x_show=False, y_show=False, z_show=False)
            )
            Iso4 = go.Isosurface(
                x=packageCartesian[0].flatten(),
                y=packageCartesian[1].flatten(),
                z=packageCartesian[2].flatten(),
                value=abs(PSI3_4).flatten(),
                colorscale='RdYlBu_r',
                isomin=-.75 * abs(PSI3_4).min(),
                isomax=.75 * abs(PSI3_4).max(),
                surface_count=6,
                opacity=0.3,
                caps=dict(x_show=False, y_show=False, z_show=False)
            )
            fig.update_layout(
                title=fr'$sp^3: Tetrahedral$',
                template='plotly_dark',
                autosize=False,
                width=700,
                height=700,
                margin=dict(l=65, r=50, b=65, t=90),
                scene=dict(
                    xaxis=dict(showgrid=False, visible=False),
                    yaxis=dict(showgrid=False, visible=False),
                    zaxis=dict(showgrid=False, visible=False)
                )
            )
            fig.add_trace(Iso1)
            fig.add_trace(Iso2)
            fig.add_trace(Iso3)
            fig.add_trace(Iso4)
            fig.show()
    # Grafico
        else:

            fig = go.Figure()
            Iso1 = go.Isosurface(
                x=x.flatten(),
                y=y.flatten(),
                z=z.flatten(),
                value=abs(PSI3_1).flatten(),
                colorscale='RdYlBu_r',
                isomin=-.75 * abs(PSI3_1).min(),
                isomax=.75 * abs(PSI3_1).max(),
                surface_count=6,
                opacity=0.6,
                caps=dict(x_show=False, y_show=False, z_show=False)
            )
            Iso2 = go.Isosurface(
                x=x.flatten(),
                y=y.flatten(),
                z=z.flatten(),
                value=abs(PSI3_2).flatten(),
                colorscale='RdYlBu_r',
                isomin=-.75 * abs(PSI3_2).min(),
                isomax=.75 * abs(PSI3_2).max(),
                surface_count=6,
                opacity=0.3,
                caps=dict(x_show=False, y_show=False, z_show=False)
            )
            Iso3 = go.Isosurface(
                x=x.flatten(),
                y=y.flatten(),
                z=z.flatten(),
                value=abs(PSI3_3).flatten(),
                colorscale='RdYlBu_r',
                isomin=-.75 * abs(PSI3_3).min(),
                isomax=.75 * abs(PSI3_3).max(),
                surface_count=6,
                opacity=0.3,
                caps=dict(x_show=False, y_show=False, z_show=False)
            )
            Iso4 = go.Isosurface(
                x=x.flatten(),
                y=y.flatten(),
                z=z.flatten(),
                value=abs(PSI3_4).flatten(),
                colorscale='RdYlBu_r',
                isomin=-.75 * abs(PSI3_4).min(),
                isomax=.75 * abs(PSI3_4).max(),
                surface_count=6,
                opacity=0.3,
                caps=dict(x_show=False, y_show=False, z_show=False)
            )
            fig.update_layout(
                title=fr'$sp^3: Tetrahedral$',
                template='plotly_dark',
                autosize=False,
                width=700,
                height=700,
                margin=dict(l=65, r=50, b=65, t=90),
                scene=dict(
                    xaxis=dict(showgrid=False, visible=False),
                    yaxis=dict(showgrid=False, visible=False),
                    zaxis=dict(showgrid=False, visible=False)
                )
            )
            fig.add_trace(Iso1)
            fig.add_trace(Iso2)
            fig.add_trace(Iso3)
            fig.add_trace(Iso4)
            fig.show()
    def plot_sp2d(self,packageCartesian=None,x=None,y=None,z=None):
        PSI2D_1,PSI2D_2,PSI2D_3,PSI2D_4 = sp2d()
        if packageCartesian:
            fig = go.Figure()
            Iso1= go.Isosurface(
                x=packageCartesian[0].flatten(),
                y=packageCartesian[1].flatten(),
                z=packageCartesian[2].flatten(),
                value=abs(PSI2D_1).flatten(),
                colorscale='RdYlBu_r',
                isomin=-.75*abs(PSI2D_1).min(),
                isomax=.75*abs(PSI2D_1).max(),
                surface_count=6,
                opacity=0.5,
                caps=dict(x_show=False,y_show=False,z_show=False)
            )
            Iso2= go.Isosurface(
                x=packageCartesian[0].flatten(),
                y=packageCartesian[1].flatten(),
                z=packageCartesian[2].flatten(),
                value=abs(PSI2D_2).flatten(),
                colorscale='RdYlBu_r',
                isomin=-.75*abs(PSI2D_2).min(),
                isomax=.75*abs(PSI2D_2).max(),
                surface_count=6,
                opacity=0.5,
                caps=dict(x_show=False,y_show=False,z_show=False)
            )
            Iso3= go.Isosurface(
                x=packageCartesian[0].flatten(),
                y=packageCartesian[1].flatten(),
                z=packageCartesian[2].flatten(),
                value=abs(PSI2D_3).flatten(),
                colorscale='RdYlBu_r',
                isomin=-.75*abs(PSI2D_3).min(),
                isomax=.75*abs(PSI2D_3).max(),
                surface_count=6,
                opacity=0.5,
                caps=dict(x_show=False,y_show=False,z_show=False)
            )
            Iso4= go.Isosurface(
                x=packageCartesian[0].flatten(),
                y=packageCartesian[1].flatten(),
                z=packageCartesian[2].flatten(),
                value=abs(PSI2D_4).flatten(),
                colorscale='RdYlBu_r',
                isomin=-.75*abs(PSI2D_4).min(),
                isomax=.75*abs(PSI2D_4).max(),
                surface_count=6,
                opacity=0.5,
                caps=dict(x_show=False,y_show=False,z_show=False)
            )
            fig.update_layout(title=fr'$sp^2d: Squere Planar,$',template = 'plotly_dark',autosize=False,
                                    width=700, height=700,
                                    margin=dict(l=65, r=50, b=65, t=90),
                                    scene=dict(xaxis=dict(showgrid=False, visible=False),
                                    yaxis=dict(showgrid=False, visible=False),
                                    zaxis=dict(showgrid=False, visible=False)))
            fig.add_trace(Iso1)
            fig.add_trace(Iso2)
            fig.add_trace(Iso3)
            fig.add_trace(Iso4)
            fig.show()
        else:

            fig = go.Figure()
            Iso1= go.Isosurface(
                x=x.flatten(),
                y=y.flatten(),
                z=z.flatten(),
                value=abs(PSI2D_1).flatten(),
                colorscale='RdYlBu_r',
                isomin=-.75*abs(PSI2D_1).min(),
                isomax=.75*abs(PSI2D_1).max(),
                surface_count=6,
                opacity=0.5,
                caps=dict(x_show=False,y_show=False,z_show=False)
            )
            Iso2= go.Isosurface(
                x=x.flatten(),
                y=y.flatten(),
                z=z.flatten(),
                value=abs(PSI2D_2).flatten(),
                colorscale='RdYlBu_r',
                isomin=-.75*abs(PSI2D_2).min(),
                isomax=.75*abs(PSI2D_2).max(),
                surface_count=6,
                opacity=0.5,
                caps=dict(x_show=False,y_show=False,z_show=False)
            )
            Iso3= go.Isosurface(
                x=x.flatten(),
                y=y.flatten(),
                z=z.flatten(),
                value=abs(PSI2D_3).flatten(),
                colorscale='RdYlBu_r',
                isomin=-.75*abs(PSI2D_3).min(),
                isomax=.75*abs(PSI2D_3).max(),
                surface_count=6,
                opacity=0.5,
                caps=dict(x_show=False,y_show=False,z_show=False)
            )
            Iso4= go.Isosurface(
                x=x.flatten(),
                y=y.flatten(),
                z=z.flatten(),
                value=abs(PSI2D_4).flatten(),
                colorscale='RdYlBu_r',
                isomin=-.75*abs(PSI2D_4).min(),
                isomax=.75*abs(PSI2D_4).max(),
                surface_count=6,
                opacity=0.5,
                caps=dict(x_show=False,y_show=False,z_show=False)
            )
            fig.update_layout(title=fr'$sp^2d: Squere Planar,$',template = 'plotly_dark',autosize=False,
                                    width=700, height=700,
                                    margin=dict(l=65, r=50, b=65, t=90),
                                    scene=dict(xaxis=dict(showgrid=False, visible=False),
                                    yaxis=dict(showgrid=False, visible=False),
                                    zaxis=dict(showgrid=False, visible=False)))
            fig.add_trace(Iso1)
            fig.add_trace(Iso2)
            fig.add_trace(Iso3)
            fig.add_trace(Iso4)
            fig.show()
    def plot_sp3d(self,packageCartesian=None,x=None,y=None,z=None): 
        PSIP3D_1, PSIP3D_2, PSIP3D_3, PSIP3D_4, PSIP3D_5 = sp3d()
        if packageCartesian:
            fig = go.Figure()
            Iso1 = go.Isosurface(
                x=packageCartesian[0].flatten(),
                y=packageCartesian[1].flatten(),
                z=packageCartesian[2].flatten(),
                value=abs(PSIP3D_1).flatten(),
                colorscale='Portland',
                isomin=-.75 * abs(PSIP3D_1).min(),
                isomax=.75 * abs(PSIP3D_1).max(),
                surface_count=6,
                opacity=0.5,
                caps=dict(x_show=False, y_show=False, z_show=False)
            )
            Iso2 = go.Isosurface(
                x=packageCartesian[0].flatten(),
                y=packageCartesian[1].flatten(),
                z=packageCartesian[2].flatten(),
                value=abs(PSIP3D_2).flatten(),
                colorscale='RdYlBu_r',
                isomin=-.75 * abs(PSIP3D_2).min(),
                isomax=.75 * abs(PSIP3D_2).max(),
                surface_count=6,
                opacity=0.5,
                caps=dict(x_show=False, y_show=False, z_show=False)
            )
            Iso3 = go.Isosurface(
                x=packageCartesian[0].flatten(),
                y=packageCartesian[1].flatten(),
                z=packageCartesian[2].flatten(),
                value=abs(PSIP3D_3).flatten(),
                colorscale='RdYlBu_r',
                isomin=-.75 * abs(PSIP3D_3).min(),
                isomax=.75 * abs(PSIP3D_3).max(),
                surface_count=6,
                opacity=0.5,
                caps=dict(x_show=False, y_show=False, z_show=False)
            )
            Iso4 = go.Isosurface(
                x=packageCartesian[0].flatten(),
                y=packageCartesian[1].flatten(),
                z=packageCartesian[2].flatten(),
                value=abs(PSIP3D_4).flatten(),
                colorscale='RdYlBu_r',
                isomin=-.75 * abs(PSIP3D_4).min(),
                isomax=.75 * abs(PSIP3D_4).max(),
                surface_count=6,
                opacity=0.5,
                caps=dict(x_show=False, y_show=False, z_show=False)
            )
            Iso5 = go.Isosurface(
                x=packageCartesian[0].flatten(),
                y=packageCartesian[1].flatten(),
                z=packageCartesian[2].flatten(),
                value=abs(PSIP3D_5).flatten(),
                colorscale='RdYlBu_r',
                isomin=-.75 * abs(PSIP3D_5).min(),
                isomax=.75 * abs(PSIP3D_5).max(),
                surface_count=6,
                opacity=0.5,
                caps=dict(x_show=False, y_show=False, z_show=False)
            )
            fig.update_layout(
                title=fr'$sp^3d$ Trigonal Bipyramidal',
                template='plotly_dark',
                autosize=False,
                width=700,
                height=700,
                margin=dict(l=65, r=50, b=65, t=90),
                scene=dict(
                    xaxis=dict(showgrid=False, visible=False),
                    yaxis=dict(showgrid=False, visible=False),
                    zaxis=dict(showgrid=False, visible=False)
                )
            )
            fig.add_trace(Iso1)
            fig.add_trace(Iso2)
            fig.add_trace(Iso3)
            fig.add_trace(Iso4)
            fig.add_trace(Iso5)
            fig.show()
        else:

            fig = go.Figure()
            Iso1 = go.Isosurface(
                x=x.flatten(),
                y=y.flatten(),
                z=z.flatten(),
                value=abs(PSIP3D_1).flatten(),
                colorscale='Portland',
                isomin=-.75 * abs(PSIP3D_1).min(),
                isomax=.75 * abs(PSIP3D_1).max(),
                surface_count=6,
                opacity=0.5,
                caps=dict(x_show=False, y_show=False, z_show=False)
            )
            Iso2 = go.Isosurface(
                x=x.flatten(),
                y=y.flatten(),
                z=z.flatten(),
                value=abs(PSIP3D_2).flatten(),
                colorscale='RdYlBu_r',
                isomin=-.75 * abs(PSIP3D_2).min(),
                isomax=.75 * abs(PSIP3D_2).max(),
                surface_count=6,
                opacity=0.5,
                caps=dict(x_show=False, y_show=False, z_show=False)
            )
            Iso3 = go.Isosurface(
                x=x.flatten(),
                y=y.flatten(),
                z=z.flatten(),
                value=abs(PSIP3D_3).flatten(),
                colorscale='RdYlBu_r',
                isomin=-.75 * abs(PSIP3D_3).min(),
                isomax=.75 * abs(PSIP3D_3).max(),
                surface_count=6,
                opacity=0.5,
                caps=dict(x_show=False, y_show=False, z_show=False)
            )
            Iso4 = go.Isosurface(
                x=x.flatten(),
                y=y.flatten(),
                z=z.flatten(),
                value=abs(PSIP3D_4).flatten(),
                colorscale='RdYlBu_r',
                isomin=-.75 * abs(PSIP3D_4).min(),
                isomax=.75 * abs(PSIP3D_4).max(),
                surface_count=6,
                opacity=0.5,
                caps=dict(x_show=False, y_show=False, z_show=False)
            )
            Iso5 = go.Isosurface(
                x=x.flatten(),
                y=y.flatten(),
                z=z.flatten(),
                value=abs(PSIP3D_5).flatten(),
                colorscale='RdYlBu_r',
                isomin=-.75 * abs(PSIP3D_5).min(),
                isomax=.75 * abs(PSIP3D_5).max(),
                surface_count=6,
                opacity=0.5,
                caps=dict(x_show=False, y_show=False, z_show=False)
            )
            fig.update_layout(
                title=fr'$sp^3d$ Trigonal Bipyramidal',
                template='plotly_dark',
                autosize=False,
                width=700,
                height=700,
                margin=dict(l=65, r=50, b=65, t=90),
                scene=dict(
                    xaxis=dict(showgrid=False, visible=False),
                    yaxis=dict(showgrid=False, visible=False),
                    zaxis=dict(showgrid=False, visible=False)
                )
            )
            fig.add_trace(Iso1)
            fig.add_trace(Iso2)
            fig.add_trace(Iso3)
            fig.add_trace(Iso4)
            fig.add_trace(Iso5)
            fig.show()
    def plot_sp3d2(self,packageCartesian=None,x=None,y=None,z=None):
      PSIP3D2_1,PSIP3D2_2,PSIP3D2_3,PSIP3D2_4,PSIP3D2_5,PSIP3D2_6 = sp3d2()
      if packageCartesian:
        fig = go.Figure()
        Iso1= go.Isosurface(
            x=packageCartesian[0].flatten(),
            y=packageCartesian[1].flatten(),
            z=packageCartesian[2].flatten(),
            value=abs(PSIP3D2_1).flatten(),
            colorscale='RdYlBu_r',
            isomin=-.75*abs(PSIP3D2_1).min(),
            isomax=.75*abs(PSIP3D2_1).max(),
            surface_count=6,
            opacity=0.5,
            caps=dict(x_show=False,y_show=False,z_show=False)
        )
        Iso2= go.Isosurface(
            x=packageCartesian[0].flatten(),
            y=packageCartesian[1].flatten(),
            z=packageCartesian[2].flatten(),
            value=abs(PSIP3D2_2).flatten(),
            colorscale='RdYlBu_r',
            isomin=-.75*abs(PSIP3D2_2).min(),
            isomax=.75*abs(PSIP3D2_2).max(),
            surface_count=6,
            opacity=0.5,
            caps=dict(x_show=False,y_show=False,z_show=False)
        )
        Iso3= go.Isosurface(
            x=packageCartesian[0].flatten(),
            y=packageCartesian[1].flatten(),
            z=packageCartesian[2].flatten(),
            value=abs(PSIP3D2_3).flatten(),
            colorscale='RdYlBu_r',
            isomin=-.75*abs(PSIP3D2_3).min(),
            isomax=.75*abs(PSIP3D2_3).max(),
            surface_count=6,
            opacity=0.5,
            caps=dict(x_show=False,y_show=False,z_show=False)
        )
        Iso4= go.Isosurface(
            x=packageCartesian[0].flatten(),
            y=packageCartesian[1].flatten(),
            z=packageCartesian[2].flatten(),
            value=abs(PSIP3D2_4).flatten(),
            colorscale='RdYlBu_r',
            isomin=-.75*abs(PSIP3D2_4).min(),
            isomax=.75*abs(PSIP3D2_4).max(),
            surface_count=6,
            opacity=0.5,
            caps=dict(x_show=False,y_show=False,z_show=False)
        )
        Iso5= go.Isosurface(
            x=packageCartesian[0].flatten(),
            y=packageCartesian[1].flatten(),
            z=packageCartesian[2].flatten(),
            value=abs(PSIP3D2_5).flatten(),
            colorscale='RdYlBu_r',
            isomin=-.75*abs(PSIP3D2_5).min(),
            isomax=.75*abs(PSIP3D2_5).max(),
            surface_count=6,
            opacity=0.5,
            caps=dict(x_show=False,y_show=False,z_show=False)
        )
        Iso6= go.Isosurface(
            x=packageCartesian[0].flatten(),
            y=packageCartesian[1].flatten(),
            z=packageCartesian[2].flatten(),
            value=abs(PSIP3D2_6).flatten(),
            colorscale='RdYlBu_r',
            isomin=-.75*abs(PSIP3D2_6).min(),
            isomax=.75*abs(PSIP3D2_6).max(),
            surface_count=6,
            opacity=0.5,
            caps=dict(x_show=False,y_show=False,z_show=False)
        )
        fig.update_layout(title=fr'$sp^3d^2: Octahedral$',template = 'plotly_dark',autosize=False,
                                    width=700, height=700,
                                    margin=dict(l=65, r=50, b=65, t=90),
                                    scene=dict(xaxis=dict(showgrid=False, visible=False),
                                    yaxis=dict(showgrid=False, visible=False),
                                    zaxis=dict(showgrid=False, visible=False)))
        fig.add_trace(Iso1)
        fig.add_trace(Iso2)
        fig.add_trace(Iso3)
        fig.add_trace(Iso4)
        fig.add_trace(Iso5)
        fig.add_trace(Iso6)
        fig.show()      
      else:
            
        fig = go.Figure()
        Iso1= go.Isosurface(
            x=x.flatten(),
            y=y.flatten(),
            z=z.flatten(),
            value=abs(PSIP3D2_1).flatten(),
            colorscale='RdYlBu_r',
            isomin=-.75*abs(PSIP3D2_1).min(),
            isomax=.75*abs(PSIP3D2_1).max(),
            surface_count=6,
            opacity=0.5,
            caps=dict(x_show=False,y_show=False,z_show=False)
        )
        Iso2= go.Isosurface(
            x=x.flatten(),
            y=y.flatten(),
            z=z.flatten(),
            value=abs(PSIP3D2_2).flatten(),
            colorscale='RdYlBu_r',
            isomin=-.75*abs(PSIP3D2_2).min(),
            isomax=.75*abs(PSIP3D2_2).max(),
            surface_count=6,
            opacity=0.5,
            caps=dict(x_show=False,y_show=False,z_show=False)
        )
        Iso3= go.Isosurface(
            x=x.flatten(),
            y=y.flatten(),
            z=z.flatten(),
            value=abs(PSIP3D2_3).flatten(),
            colorscale='RdYlBu_r',
            isomin=-.75*abs(PSIP3D2_3).min(),
            isomax=.75*abs(PSIP3D2_3).max(),
            surface_count=6,
            opacity=0.5,
            caps=dict(x_show=False,y_show=False,z_show=False)
        )
        Iso4= go.Isosurface(
            x=x.flatten(),
            y=y.flatten(),
            z=z.flatten(),
            value=abs(PSIP3D2_4).flatten(),
            colorscale='RdYlBu_r',
            isomin=-.75*abs(PSIP3D2_4).min(),
            isomax=.75*abs(PSIP3D2_4).max(),
            surface_count=6,
            opacity=0.5,
            caps=dict(x_show=False,y_show=False,z_show=False)
        )
        Iso5= go.Isosurface(
            x=x.flatten(),
            y=y.flatten(),
            z=z.flatten(),
            value=abs(PSIP3D2_5).flatten(),
            colorscale='RdYlBu_r',
            isomin=-.75*abs(PSIP3D2_5).min(),
            isomax=.75*abs(PSIP3D2_5).max(),
            surface_count=6,
            opacity=0.5,
            caps=dict(x_show=False,y_show=False,z_show=False)
        )
        Iso6= go.Isosurface(
            x=x.flatten(),
            y=y.flatten(),
            z=z.flatten(),
            value=abs(PSIP3D2_6).flatten(),
            colorscale='RdYlBu_r',
            isomin=-.75*abs(PSIP3D2_6).min(),
            isomax=.75*abs(PSIP3D2_6).max(),
            surface_count=6,
            opacity=0.5,
            caps=dict(x_show=False,y_show=False,z_show=False)
        )
        fig.update_layout(title=fr'$sp^3d^2: Octahedral$',template = 'plotly_dark',autosize=False,
                                    width=700, height=700,
                                    margin=dict(l=65, r=50, b=65, t=90),
                                    scene=dict(xaxis=dict(showgrid=False, visible=False),
                                    yaxis=dict(showgrid=False, visible=False),
                                    zaxis=dict(showgrid=False, visible=False)))
        fig.add_trace(Iso1)
        fig.add_trace(Iso2)
        fig.add_trace(Iso3)
        fig.add_trace(Iso4)
        fig.add_trace(Iso5)
        fig.add_trace(Iso6)
        fig.show()