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

# Will hold the path that the user chooses to download to. Will be an empty
# string if the user decides to download to the current directory.
dl_path = ''

# Determines if the script should ask the user if they would like to remove the
# downloaded data.  This defaults to False, because data deletion is
# irreversible, and in the worst case, the user will have to remove it
# manually themselves.
run_cleanup = False

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

    xs = Array
    ys = Array
    vals = Array

    minval = Float
    maxval = Float

    model_changed = Event

    def __init__(self, *args, **kwargs):
        super(Model, self).__init__(*args, **kwargs)
        self.compute_model()

    @on_trait_change("npts_+", "min_+", "max_+")
    def compute_model(self):
        def vfunc(x, y, z):
            return sin(x*z) * cos(y)*sin(z) + sin(0.5*z)

        # Create the axes
        self.xs = linspace(self.min_x, self.max_x, self.npts_x)
        self.ys = linspace(self.min_y, self.max_y, self.npts_y)
        self.zs = linspace(self.min_z, self.max_z, self.npts_z)

        # Generate a cube of values by using newaxis to span new dimensions
        self.vals = vfunc(self.xs[:, newaxis, newaxis],
                          self.ys[newaxis, :, newaxis],
                          self.zs[newaxis, newaxis, :])

        self.minval = nanmin(self.vals)
        self.maxval = nanmax(self.vals)
        self.model_changed = True


class BrainModel(Model):
    def __init__(self, *args, **kwargs):
        download_data()
        super(BrainModel, self).__init__(*args, **kwargs)

    def compute_model(self):
        global dl_path
        mrbrain_path = os.path.join(dl_path, 'voldata', 'MRbrain.')
        nx = 256
        ny = 256
        nz = 109
        full_arr = zeros((nx, ny, nz), dtype='f')
        for i in range(1, 110):
            arr = fromfile(mrbrain_path + str(i), dtype='>u2')
            arr.shape = (256,256)
            full_arr[:,:,i-1] = arr
        self.vals = full_arr

        # Create the axes
        self.xs = arange(nx)
        self.ys = arange(ny)
        self.zs = arange(nz)

        # Generate a cube of values by using newaxis to span new dimensions
        self.minval = nanmin(self.vals)
        self.maxval = nanmax(self.vals)
        self.model_changed = True


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
        # Create the model
        #try:
        #    self.model = model = BrainModel()
        #    cmap = bone
        #except SystemExit:
        #    sys.exit()
        #except:
        #    print "Unable to load BrainModel, using generated data cube."
        self.model = model = Model()
        cmap = jet
        self._update_model(cmap)

        datacube = self.colorcube

        # Create the plot
        self.plotdata = ArrayPlotData()
        self._update_images()

        # Center Plot
        centerplot = Plot(self.plotdata, padding=0)
        imgplot = centerplot.img_plot("xy",
                                xbounds=(model.xs[0], model.xs[-1]),
                                ybounds=(model.ys[0], model.ys[-1]),
                                colormap=cmap)[0]
        self._add_plot_tools(imgplot, "xy")
        self.center = imgplot

        # Right Plot
        rightplot = Plot(self.plotdata, width=150, resizable="v", padding=0)
        rightplot.value_range = centerplot.value_range
        imgplot = rightplot.img_plot("yz",
                                xbounds=(model.zs[0], model.zs[-1]),
                                ybounds=(model.ys[0], model.ys[-1]),
                                colormap=cmap)[0]
        self._add_plot_tools(imgplot, "yz")
        self.right = imgplot

        # Bottom Plot
        bottomplot = Plot(self.plotdata, height=150, resizable="h", padding=0)
        bottomplot.index_range = centerplot.index_range
        imgplot = bottomplot.img_plot("xz",
                                xbounds=(model.xs[0], model.xs[-1]),
                                ybounds=(model.zs[0], model.zs[-1]),
                                colormap=cmap)[0]
        self._add_plot_tools(imgplot, "xz")
        self.bottom = imgplot

        # Create Container and add all Plots
        container = GridPlotContainer(padding=20, fill_padding=True,
                                      bgcolor="white", use_backbuffer=True,
                                      shape=(2,2), spacing=(12,12))
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

    def _update_model(self, cmap):
        range = DataRange1D(low=amin(self.model.vals),
                            high=amax(self.model.vals))
        self.colormap = cmap(range)
        self.colorcube = (self.colormap.map_screen(self.model.vals) * 255).astype(uint8)

    def _update_images(self):
        """ Updates the image data in self.plotdata to correspond to the
        slices given.
        """
        cube = self.colorcube
        pd = self.plotdata
        # These are transposed because img_plot() expects its data to be in
        # row-major order
        pd.set_data("xy", transpose(cube[:, :, self.slice_z], (1,0,2)))
        pd.set_data("xz", transpose(cube[:, self.slice_y, :], (1,0,2)))
        pd.set_data("yz", cube[self.slice_x, :, :])

def download_data():
    global dl_path, run_cleanup

    print 'Please enter the location of the "voldata" subdirectory containing'
    print 'the data files for this demo, or enter a path to download to (7.8MB).'
    print 'Press <ENTER> to download to the current directory.'
    dl_path = raw_input('Path: ').strip().rstrip("/").rstrip("\\")

    if not dl_path.endswith("voldata"):
        voldata_path = os.path.join(dl_path, 'voldata')
    else:
        voldata_path = dl_path
    tar_path = os.path.join(dl_path, 'MRbrain.tar.gz')

    data_good = True
    try:
        for i in range(1,110):
            if not os.path.isfile(os.path.join(voldata_path, "MRbrain.%d" % i)):
                data_good = False
                break
        else:
            data_good = True
    except:
        data_good = False

    if not data_good:
        import urllib
        import tarfile

        if len(dl_path) > 0 and not os.path.exists(dl_path):
            print 'The given path does not exist.'
            run_cleanup = False
            sys.exit()

        if not os.path.isabs(dl_path):
            print 'Downloading to: ' + os.path.join(os.getcwd(), dl_path)
        else:
            print 'Downloading to: ' + dl_path

        try:
            # download and extract the file
            print "Downloading data, Please Wait (7.8MB)"
            opener = urllib.urlopen('http://www-graphics.stanford.edu/data/voldata/MRbrain.tar.gz')
        except:
            print 'Download error. Opening backup data.'
            run_cleanup = False
            raise

        try:
            open(tar_path, 'wb').write(opener.read())
        except:
            print 'Cannot write to the destination directory specified. ' \
                  'Opening backup data.'
            run_cleanup = False
            raise

        tar_file = tarfile.open(tar_path)
        try:
            os.mkdir(voldata_path)
        except:
            pass
        tar_file.extractall(voldata_path)
        tar_file.close()
        os.unlink(tar_path)
    else:
        print 'Previously downloaded data detected.'

def cleanup_data():
    global dl_path

    answer = raw_input('Remove downloaded files? [Y/N]: ')
    if answer.lower() == 'y':
        try:
            shutil.rmtree(os.path.join(dl_path, 'voldata'))
        except:
            pass

if __name__ == "__main__":
    demo_main(PlotFrame, size=(800,700), title="Cube analyzer")
    if run_cleanup:
        cleanup_data()

