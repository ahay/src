#!/usr/bin/env python
"""
Draws an simple RGB image
 - Left-drag pans the plot.
 - Mousewheel up and down zooms the plot in and out.
 - Pressing "z" brings up the Zoom Box, and you can click-drag a rectangular
   region to zoom.  If you use a sequence of zoom boxes, pressing alt-left-arrow
   and alt-right-arrow moves you forwards and backwards through the "zoom
   history".
"""

# Major library imports
from numpy import zeros, uint8

# Enthought library imports
from enable.api import Component, ComponentEditor
from traits.api import HasTraits, Instance
from traitsui.api import Item, Group, View

# Chaco imports
from chaco.api import ArrayPlotData, Plot
from chaco.tools.api import PanTool, ZoomTool
from chaco.tools.image_inspector_tool import ImageInspectorTool, \
     ImageInspectorOverlay

#===============================================================================
# # Create the Chaco plot.
#===============================================================================
def _create_plot_component():

    # Create some RGBA image data
    image = zeros((200,400,4), dtype=uint8)
    image[:,0:40,0] += 255     # Vertical red stripe
    image[0:25,:,1] += 255     # Horizontal green stripe; also yellow square
    image[-80:,-160:,2] += 255 # Blue square
    image[:,:,3] = 255

    # Create a plot data obect and give it this data
    pd = ArrayPlotData()
    pd.set_data("imagedata", image)

    # Create the plot
    plot = Plot(pd, default_origin="top left")
    plot.x_axis.orientation = "top"
    img_plot = plot.img_plot("imagedata")[0]

    # Tweak some of the plot properties
    plot.bgcolor = "white"

    # Attach some tools to the plot
    plot.tools.append(PanTool(plot, constrain_key="shift"))
    plot.overlays.append(ZoomTool(component=plot,
                                    tool_mode="box", always_on=False))

    imgtool = ImageInspectorTool(img_plot)
    img_plot.tools.append(imgtool)
    plot.overlays.append(ImageInspectorOverlay(component=img_plot,
                                               image_inspector=imgtool))
    return plot

#===============================================================================
# Attributes to use for the plot view.
size = (600, 600)
title="Simple image plot"
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

#--EOF---
