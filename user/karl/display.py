from traits.api import HasTraits, Instance, Button, Enum, Int, Float
from traitsui.api import View, Item, Group
from chaco.api import HPlotContainer, Plot, ArrayPlotData, DataRange1D
from enable.api import ColorTrait
from enable.component_editor import ComponentEditor
from numpy import linspace, nanmin, nanmax, transpose, uint8, min, max
from enable.api import BaseTool
import sys
import m8r
from chaco.default_colormaps import *
" derived from /home/karl/learn/python/chaco/button.py "

class SeisData(HasTraits):
    def __init__(self, filename):
        super(SeisData, self).__init__()
        
        self.model=m8r.Input(filename)
        self.vals = self.model[:,:,:]
        print "vals.shape=",self.vals.shape
        for i in range(0,self.vals.shape[1]):
            print "max(vals[:,%d,:])="%i,max(self.vals[:,i,:])

        self.dim=len(self.vals.shape)
        self.axis_start=self.read_axis_float_info("o")
        print "self.axis_start=",self.axis_start
        self.axis_delta=self.read_axis_float_info("d")
        print "self.axis_delta=",self.axis_delta

        self.axis_end=[]
        for i in range(0,self.dim):
            self.axis_end.append(self.axis_start[i]+
                                 (self.vals.shape[i]-1)*self.axis_delta[i])
        print "self.axis_end=",self.axis_end
        
        print "compute min"
        self.minval = nanmin(self.vals)
        print "compute max"
        self.maxval = nanmax(self.vals)
        print "set model changed"
        self.model_changed = True
        print "leaving m8rInput"

    def read_axis_float_info(self,letter):
        list_floats=[]
        for i in range(1,self.dim+1):
            print "get parameter",letter+"%d"%i
            list_floats.append(self.model.float(letter+"%d"%i,1))
        list_floats.reverse()
        return tuple(list_floats)

class ContainerExample(HasTraits):
    plot = Instance(Plot)
    display_button = Button()
    display_button1 = Button()

    traits_view = View(
        Group(Item('display_button', show_label=False),
              Item('display_button1', show_label=False),
              orientation="horizontal"),
        Item('plot',editor=ComponentEditor(), show_label=False),
        width=1000, height=600, resizable=True, title="Chaco Plot",
        buttons=["prev","next","new data"])

    def __init__(self):

        super(ContainerExample, self).__init__()

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

        self.seis_data_0=SeisData(filenames[0])
        self.cmap = jet

        self.displayParameters=DisplayParameters()

        self.slice_y=self.displayParameters.slice_y
        print "self.slice_y=",self.slice_y
        self.arrayPlotData=ArrayPlotData()
        self._update_images()
        bottomplot = Plot(self.arrayPlotData, padding=0, origin="top left")
        imgplot = bottomplot.img_plot("xz",
                                xbounds=(self.seis_data_0.axis_start[0], self.seis_data_0.axis_end[0]),
                                ybounds=(self.seis_data_0.axis_start[2], self.seis_data_0.axis_end[2]),
                                colormap=self.cmap)[0]
        self.bottom = imgplot

        print "x range=",self.seis_data_0.axis_start[0], self.seis_data_0.axis_end[0]
        print "z range=",self.seis_data_0.axis_start[2], self.seis_data_0.axis_end[2]
        #print "colormap.shape=',colormap.shape)
        x = linspace(-14, 14, 100)

        y = x**3
        plotdata = ArrayPlotData(x=x, y=y)

        scatter = Plot(plotdata, origin="top left")
        scatter.plot(("x", "y"), type="scatter", color='blue')
        scatter.title="x**3"
        self.scatter=scatter
        #scatter.origin="top left"
        print "make container"
        container = HPlotContainer(bottomplot)
        print "container.add"
        container.add(scatter)
        #container.spacing = 0

#        scatter.padding_right = 0

#        line.padding_left = 0
#        line.y_axis.orientation = "right"

        self.plot = bottomplot #container  
                               # kls karl !! When I change plot = Instance(Plot) and this line to
                               # self.plot=container I get no img_plotplot to
#        self.plot.tools.append(CustomTool(self.plot))

        self.displayParameters=DisplayParameters()

        imgplot.tools.append(CustomTool(self.plot))

    def _update_images(self):
        range = DataRange1D(low=self.seis_data_0.minval,
                            high=self.seis_data_0.maxval)
        print "self.seis_data_0.min/max=",self.seis_data_0.minval,self.seis_data_0.maxval
        self.colormap = self.cmap(range)

        print "self.slice_y=",self.slice_y
        slicexz=self.seis_data_0.vals[:, self.slice_y, :]
        print "min/max slicexz=",min(slicexz),max(slicexz)
        colorslicexz=(self.colormap.map_screen(slicexz) * 255).astype(uint8)
        # Transposed required because img_plot() expects data in row-major order
#        self.arrayPlotData.set_data("xz", transpose(colorslicexz,(1,0,2)))
        self.arrayPlotData.set_data("xz", colorslicexz)
        print "colorslicexz.shape=",colorslicexz.shape
        print "type(colorslicexz)=",type(colorslicexz)
        print "type(colorslicexz[0,0,0])=",type(colorslicexz[0,0,0])
        print "colorslicexz=",colorslicexz
        print "min(colorslicexz[0,:,:,:]=",min(colorslicexz[:,:,0])
        print "min(colorslicexz[1,:,:,:]=",min(colorslicexz[:,:,1])
        print "min(colorslicexz[2,:,:,:]=",min(colorslicexz[:,:,2])
        print "min(colorslicexz[3,:,:,:]=",min(colorslicexz[:,:,3])

    def _marker_size_changed(self):
        self.scatter.marker_size = self.marker_size
    def _color_changed(self):
        self.scatter.marker_size = self.marker_size
    def _display_button_fired(self):
        print "button pushed"
        self.displayParameters.edit_traits()
        
    def _next_fired(self):
        print "next button pushed"

class DisplayParameters(HasTraits):
    gain = Float(10)
    vaxis = Int(0)
    haxis = Int(1)
    slice_y = Int(63)
    
class CustomTool(BaseTool):
    def normal_mouse_move(self, event):
        print "type event=",type(event),"Screen point:", event.x, event.y
        print "Data:", self.component.map_data((event.x, event.y))
        #print "current_pointer_position", event.current_pointer_position
        #print "scale_xy", event.scale_xy
        #print "event", event


if __name__ == "__main__":
    myContainerExample=ContainerExample()
    myContainerExample.configure_traits()
