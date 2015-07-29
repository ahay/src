"""
Renders some contoured and colormapped images of a scalar value field.
 - Left-drag pans the plot.
 - Mousewheel up and down zooms the plot in and out.
 - Pressing "z" brings up the Zoom Box, and you can click-drag a rectangular
   region to zoom.  If you use a sequence of zoom boxes, pressing alt-left-arrow   and alt-right-arrow moves you forwards and backwards through the "zoom
   history".
"""

# Major library imports
from numpy import cos, linspace, log, meshgrid, pi, sin

# Enthought library imports
from enable.api import Component, ComponentEditor
from traits.api import HasTraits, Instance
from traitsui.api import Item, Group, View

# Chaco imports
from chaco.api import ArrayPlotData, ColorBar, gmt_drywet, \
                                 HPlotContainer, LinearMapper, Plot
from chaco.tools.api import PanTool, ZoomTool


#===============================================================================
# # Create the Chaco plot.
#===============================================================================
def _create_plot_component():

    # Create a scalar field to colormap
    x_extents = (-2*pi, 2*pi)
    y_extents = (-1.5*pi, 1.5*pi)
    xs = linspace(-2*pi, 2*pi, 200)
    ys = linspace(-1.5*pi, 1.5*pi, 100)
    x, y = meshgrid(xs,ys)
    zs = sin(log(abs((x+1)**4)+0.05))*cos(y)*1.1*(-y) + \
            sin(((x+1)**2 + y**2)/4)

    # Create a plot data obect and give it this data
    pd = ArrayPlotData()
    pd.set_data("imagedata", zs)

    # Create the left plot, a colormap and simple contours
    lplot = Plot(pd)
    lplot.img_plot("imagedata",
                   name="cm_plot",
                   xbounds=x_extents,
                   ybounds=y_extents,
                   colormap=gmt_drywet)
    lplot.contour_plot("imagedata",
                       type="line",
                       xbounds=x_extents,
                       ybounds=y_extents)

    # Tweak some of the plot properties
    lplot.title = "Colormap and contours"
    lplot.padding = 20
    lplot.bg_color = "white"
    lplot.fill_padding = True

    # Add some tools to the plot
    zoom = ZoomTool(lplot, tool_mode="box", always_on=False)
    lplot.overlays.append(zoom)
    lplot.tools.append(PanTool(lplot, constrain_key="shift"))

    # Right now, some of the tools are a little invasive, and we need the
    # actual CMapImage object to give to them
    cm_plot = lplot.plots["cm_plot"][0]

    # Create the colorbar, handing in the appropriate range and colormap
    colormap = cm_plot.color_mapper
    colorbar = ColorBar(index_mapper=LinearMapper(range=colormap.range),
                        color_mapper=colormap,
                        plot=cm_plot,
                        orientation='v',
                        resizable='v',
                        width=30,
                        padding=20)
    colorbar.padding_top = lplot.padding_top
    colorbar.padding_bottom = lplot.padding_bottom

    # Create the left plot, contours of varying color and width
    rplot = Plot(pd, range2d=lplot.range2d)
    rplot.contour_plot("imagedata",
                       type="line",
                       xbounds=x_extents,
                       ybounds=y_extents,
                       bgcolor="black",
                       levels=15,
                       styles="solid",
                       widths=list(linspace(4.0, 0.1, 15)),
                       colors=gmt_drywet)

    # Add some tools to the plot
    zoom = ZoomTool(rplot, tool_mode="box", always_on=False)
    rplot.overlays.append(zoom)
    rplot.tools.append(PanTool(rplot, constrain_key="shift"))

    # Tweak some of the plot properties
    rplot.title = "Varying contour lines"
    rplot.padding = 20
    rplot.bg_color = "white"
    rplot.fill_padding = True

    # Create a container and add our plots
    container = HPlotContainer(padding=40, fill_padding=True,
                               bgcolor = "white", use_backbuffer=True)
    container.add(colorbar)
    container.add(lplot)
    container.add(rplot)
    return container

#===============================================================================
# Attributes to use for the plot view.
size=(950,650)
title="Some contour plots"
bg_color="lightgray"

#===============================================================================
# # Demo class that is used by the demo.py application.
#===============================================================================
class Demo(HasTraits):
    plot = Instance(Component)

    traits_view = View(
                    Group(
                        Item('plot', editor=ComponentEditor(size=size,
                                                            bgcolor=bg_color),
                             show_label=False),
                        orientation = "vertical"),
                    resizable=True, title=title
                    )
    def _plot_default(self):
         return _create_plot_component()

demo = Demo()

if __name__ == "__main__":
    demo.configure_traits()

# EOF
