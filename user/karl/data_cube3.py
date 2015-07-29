"""
Allows isometric viewing of a 3D data cube.

Click or click-drag in any data window to set the slice to view.
"""

# Outstanding TODOs:
#  - need to add line inspectors to side and bottom plots, and synchronize
#    with center plot
#  - need to set the various image plots to use the same colormap instance,
#    and that colormap's range needs to be set to min/max of the entire cube
#  - refactor create_window() so there is less code duplication
#  - try to eliminate the use of model.xs, ys, zs in favor of bounds tuples
from numpy import amin, amax, zeros, fromfile, transpose, uint8

# Standard library imports
import os, sys, shutil

# Major library imports
from numpy import arange, linspace, nanmin, nanmax, newaxis, pi, sin, cos

# Enthought library imports
from chaco.api import ArrayPlotData, Plot, GridPlotContainer, \
                                 BaseTool, DataRange1D
from chaco.default_colormaps import *
from chaco.tools.api import LineInspector, ZoomTool
from enable.example_support import DemoFrame, demo_main
from enable.api import Window
from traits.api import Any, Array, Bool, Callable, CFloat, CInt, \
        Event, Float, HasTraits, Int, Trait, on_trait_change

import m8r

class Model(HasTraits):
    npts_x = CInt(256)
    npts_y = CInt(256)
    npts_z = CInt(109)

    min_x = CFloat(-2*pi)
    max_x = CFloat(2*pi)
    min_y = CFloat(-2*pi)
    max_y = CFloat(2*pi)
    min_z = CFloat(-pi)
    max_z = CFloat(pi)

    vals = Array

    minval = Float
    maxval = Float

    model_changed = Event

    def __init__(self, *args, **kwargs):
        super(Model, self).__init__(*args, **kwargs)
        self.m8rInput()

    def m8rInput(self):
        #model=m8r.Input('model.rsf')
        # get files names from command line
        filenames=[]
        for parameter in sys.argv[1:]:
            print "processing parameter",parameter
            if parameter.find("=")==-1 :
                print "no = in parameter",parameter,"must be a file name"
                filenames.append(parameter)
        if len(filenames)<1:
            print "just to help me test, if there are no files in the list, "
            print "I will append the file foldplot1.rsf"
            filenames.append('foldplot1.rsf')

        model=m8r.Input(filenames[0])
        self.vals = model[:,:,:]
        # this makes seg fault in nanmin self.vals = model[:,:,::-1]
        npts_x, npts_y, npts_z = self.vals.shape
        print "vals.shape=",self.vals.shape

        self.min_x=model.float("o3")
        self.min_y=model.float("o2")
        self.min_z=model.float("o1")
        print "min x,y,z=",self.min_x,self.min_y,self.min_z
        
        dx=model.float("d3")
        dy=model.float("d2")
        dz=model.float("d1")
        print "dx,dy,dz=",dx,dy,dz

        self.max_x=self.min_x+dx*(npts_x-1)
        self.max_y=self.min_y+dy*(npts_y-1)
        self.max_z=self.min_z+dz*(npts_z-1)

        print "compute linspace"

        print "compute min"
        self.minval = nanmin(self.vals)
        print "compute max"
        self.maxval = nanmax(self.vals)
        print "set model changed"
        self.model_changed = True
        print "leaving m8rInput"
                
class ImageIndexTool(BaseTool):
    """ A tool to set the slice of a cube based on the user's mouse movements
    or clicks.
    """

    # This callback will be called with the index into self.component's
    # index and value:
    #     callback(tool, x_index, y_index)
    # where *tool* is a reference to this tool instance.  The callback
    # can then use tool.token.
    callback = Any()

    # This callback (if it exists) will be called with the integer number
    # of mousewheel clicks
    wheel_cb = Any()

    # This token can be used by the callback to decide how to process
    # the event.
    token  = Any()

    # Whether or not to update the slice info; we enter select mode when
    # the left mouse button is pressed and exit it when the mouse button
    # is released
    # FIXME: This is not used right now.
    select_mode = Bool(False)

    def normal_left_down(self, event):
        self._update_slices(event)

    def normal_right_down(self, event):
        self._update_slices(event)

    def normal_mouse_move(self, event):
        if event.left_down or event.right_down:
            self._update_slices(event)

    def _update_slices(self, event):
            plot = self.component
            ndx = plot.map_index((event.x, event.y),
                                 threshold=5.0, index_only=True)
            if ndx:
                self.callback(self, *ndx)

    def normal_mouse_wheel(self, event):
        if self.wheel_cb is not None:
            self.wheel_cb(self, event.mouse_wheel)


class PlotFrame(DemoFrame):

    # These are the indices into the cube that each of the image plot views
    # will show; the default values are non-zero just to make it a little
    # interesting.
    slice_x = 10
    slice_y = 10
    slice_z = 10

    num_levels = Int(15)
    colormap = Any
    colorcube = Any

    #---------------------------------------------------------------------------
    # Private Traits
    #---------------------------------------------------------------------------

    _cmap = Trait(jet, Callable)

    def _index_callback(self, tool, x_index, y_index):
        plane = tool.token
        if plane == "xy":
            self.slice_x = x_index
            self.slice_y = y_index
        elif plane == "yz":
            # transposed because the plot is oriented vertically
            self.slice_z = x_index
            self.slice_y = y_index
        elif plane == "xz":
            self.slice_x = x_index
            self.slice_z = y_index
        else:
            warnings.warn("Unrecognized plane for _index_callback: %s" % plane)
        self._update_images()
        self.center.invalidate_and_redraw()
        self.right.invalidate_and_redraw()
        self.bottom.invalidate_and_redraw()
        return

    def _wheel_callback(self, tool, wheelamt):
        plane_slice_dict = {"xy": ("slice_z", 2),
                            "yz": ("slice_x", 0),
                            "xz": ("slice_y", 1)}
        attr, shape_ndx = plane_slice_dict[tool.token]
        val = getattr(self, attr)
        max = self.model.vals.shape[shape_ndx]
        if val + wheelamt > max:
            setattr(self, attr, max-1)
        elif val + wheelamt < 0:
            setattr(self, attr, 0)
        else:
            setattr(self, attr, val + wheelamt)

        self._update_images()
        self.center.invalidate_and_redraw()
        self.right.invalidate_and_redraw()
        self.bottom.invalidate_and_redraw()
        return

    def _create_window(self):
        self.model = model = Model()
        self.cmap = jet
        #self._update_model(self.cmap)

        # Create the plot
        self.plotdata = ArrayPlotData()
        self._update_images()

        # Center Plot
        centerplot = Plot(self.plotdata, padding=0)
        imgplot = centerplot.img_plot("xy",
                                xbounds=(model.min_x, model.max_x),
                                ybounds=(model.min_y, model.max_y),
                                colormap=self.cmap)[0]
        self._add_plot_tools(imgplot, "xy")
        self.center = imgplot

        # Right Plot
        rightplot = Plot(self.plotdata, width=150, resizable="v", padding=0)
        rightplot.value_range = centerplot.value_range
        imgplot = rightplot.img_plot("yz",
                                xbounds=(model.min_z, model.max_z),
                                ybounds=(model.min_y, model.max_y),
                                colormap=self.cmap)[0]
        self._add_plot_tools(imgplot, "yz")
        self.right = imgplot

        # Bottom Plot.  Seismic plot axis1 (depth) down into earth 
        #               i.e. z is depth, to altitude.
        bottomplot = Plot(self.plotdata, height=150, resizable="h", 
                          padding=0, origin="top left")
        bottomplot.index_range = centerplot.index_range
        imgplot = bottomplot.img_plot("xz",
                                xbounds=(model.min_x, model.max_x),
                                ybounds=(model.min_z, model.max_z),
                                colormap=self.cmap)[0]
        self._add_plot_tools(imgplot, "xz")
        self.bottom = imgplot

        # Create Container and add all Plots
        container = GridPlotContainer(padding=20, fill_padding=True,
                                      bgcolor="white", use_backbuffer=True,
                                      shape=(2,2), spacing=(20,20))
        container.add(centerplot)
        container.add(rightplot)
        container.add(bottomplot)

        self.container = container
        return Window(self, -1, component=container)

    def _add_plot_tools(self, imgplot, token):
        """ Add LineInspectors, ImageIndexTool, and ZoomTool to the image plots. """

        imgplot.overlays.append(ZoomTool(component=imgplot, tool_mode="box",
                                           enable_wheel=False, always_on=False))
        imgplot.overlays.append(LineInspector(imgplot, axis="index_y", color="white",
            inspect_mode="indexed", write_metadata=True, is_listener=True))
        imgplot.overlays.append(LineInspector(imgplot, axis="index_x", color="white",
            inspect_mode="indexed", write_metadata=True, is_listener=True))
        imgplot.tools.append(ImageIndexTool(imgplot, token=token,
            callback=self._index_callback, wheel_cb=self._wheel_callback))

    def _update_images(self):
        """ Updates the image data in self.plotdata to correspond to the
        slices given.
        """
        range = DataRange1D(low=self.model.minval,
                            high=self.model.maxval)
        self.colormap = self.cmap(range)


        slicexy=self.model.vals[:, :, self.slice_z]
        colorslicexy=(self.colormap.map_screen(slicexy) * 255).astype(uint8)
        # Transposed required because img_plot() expects data in row-major order
        self.plotdata.set_data("xy", transpose(colorslicexy,(1,0,2)))

        
        slicexz=self.model.vals[:, self.slice_y, :]
        colorslicexz=(self.colormap.map_screen(slicexz) * 255).astype(uint8)
        # Transposed required because img_plot() expects data in row-major order
        self.plotdata.set_data("xz", transpose(colorslicexz,(1,0,2)))


        sliceyz=self.model.vals[self.slice_x, :, :]
        colorsliceyz=(self.colormap.map_screen(sliceyz) * 255).astype(uint8)
        #print "colorsliceyz.shape=",colorsliceyz.shape
        #print "type(colorsliceyz)=",type(colorsliceyz)
        #print "type(colorsliceyz[0,0,0])=",type(colorsliceyz[0,0,0])
        self.plotdata.set_data("yz", colorsliceyz)

if __name__ == "__main__":
    demo_main(PlotFrame, size=(800,700), title="Cube analyzer")

