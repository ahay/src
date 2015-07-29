#!/usr/bin/env python
from traits.api import HasTraits, Instance, Button, Enum, Int, Float, Range
from traitsui.api import View, Item, Group
from chaco.api import HPlotContainer, Plot, ArrayPlotData, DataRange1D
from chaco.tools.api import PanTool, ZoomTool
from enable.api import ColorTrait
from enable.component_editor import ComponentEditor
from numpy import linspace, nanmin, nanmax, transpose, uint8, min, max, ones
from enable.api import BaseTool
import sys
import m8r
from chaco.default_colormaps import *
" derived from /home/karl/learn/python/chaco/button.py "

# add the unzoom button function karl kls

class SeisData(HasTraits):
    def __init__(self, filename):
        super(SeisData, self).__init__()
        
        self.model=m8r.Input(filename)
        self.vals = self.model[:,:,:]

        print "input file shape=",self.vals.shape
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
        
        print "compute min/max"
        
        max_n_samples=100
        inc=ones(self.dim,dtype=int)
        last=ones(self.dim,dtype=int)
        for i in range(0,self.dim):
            inc[i]=self.vals.shape[i]/100
            #print "self.vals.shape=", self.vals.shape ,"inc=",inc
            if(inc[i]<1):
                inc[i]=1
            last[i]=(self.vals.shape[i]/inc[i]-1)*inc[i]
            #print "self.vals.shape=", self.vals.shape
            #print "inc=",inc,"last=",last
        subsetvals=self.vals[:last[0]:inc[0],:last[1]:inc[1],:last[2]:inc[2]]
        #print "subsetvals.shape=",subsetvals.shape
        mymin=min(subsetvals)
        mymax=max(subsetvals)
        print "min/max",mymin,mymax
        self.maxval = max([abs(mymin),mymax])
        self.minval = -self.maxval
        print "min=",self.minval,"max=",self.maxval

        print "leaving m8rInput"

    def read_axis_float_info(self,letter):
        list_floats=[]
        for i in range(1,self.dim+1):
            # print "get parameter",letter+"%d"%i
            list_floats.append(self.model.float(letter+"%d"%i,1))
        list_floats.reverse()
        return tuple(list_floats)

class ContainerExample(HasTraits):
    plot = Instance(HPlotContainer)
    display_button = Button()
    display_button1 = Button()
    prev = Button()
    next = Button()
    unzoom = Button()

    traits_view = View(
        Group(Item('display_button', show_label=False),
              Item('display_button1', show_label=False),
              Item('prev', show_label=False),
              Item('next', show_label=False),
              Item('unzoom', show_label=False),
              orientation="horizontal"),
        Item('plot',editor=ComponentEditor(), show_label=False),
        # put the info box that lists the mouse location tuple kls karl
        width=1000, height=600, resizable=True, title="Madagascar Comparison",
        # How do I activate these? buttons=["do_nothing","do_nothing_1"]
        )

    def __init__(self):

        super(ContainerExample, self).__init__()

        # make a list of the files names on the command line
        # command line entries without = are assumed to be files
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

        # there should be a try while these files are opened to identify 
        # bad file names and make error messages.

        self.seis_data=[]
        for filename in filenames:
            print "  "
            print "+++++++++++++++++++++++++++++++++++++++++++"
            print "+++++++++++++++++++++++++++++++++++++++++++"
            print "open the file",filename
            self.seis_data.append(SeisData(filename))

        # kls files should be tested to make sure they are all same shape
        print "number files=",len(filenames)
        self.displayParameters=DisplayParameters()
        self.cmap = jet  # kls this should come from display parameters

        self.slice_y=self.displayParameters.slice_y # kls this should be pointer


        print "self.slice_y=",self.slice_y
        self.arrayPlotDatas=[]
        for filename in filenames:
            self.arrayPlotDatas.append(ArrayPlotData())
        self._update_images()

        x_extents=(self.seis_data[0].axis_start[1], 
                   self.seis_data[0].axis_end[1])
        y_extents=(self.seis_data[0].axis_start[2], 
                   self.seis_data[0].axis_end[2])
        bottomplot = Plot(self.arrayPlotDatas[0], origin="top left")
        self.bottomplot=bottomplot
        imgplot = bottomplot.img_plot("xz",
                                xbounds=x_extents,
                                ybounds=y_extents,
                                colormap=self.cmap)[0]
        self.bottom = imgplot

        plotright = Plot(self.arrayPlotDatas[1], origin="top left",
                         range2d=bottomplot.range2d)
        imgplotright = plotright.img_plot("xz",
                                xbounds=x_extents,
                                ybounds=y_extents,
                                colormap=self.cmap)[0]
        self.right = imgplotright

        container = HPlotContainer(fill_padding=True,
                               bgcolor = "white", use_backbuffer=True)
        container.add(bottomplot)
        container.add(plotright)

        self.plot = container

        self.displayParameters=DisplayParameters()
        self.slice_y=self.displayParameters.slice_y

        imgplot.tools.append(CustomTool(imgplot))
        imgplot.tools.append(PanTool(imgplot, constrain_key="shift"))
        imgplot.overlays.append(ZoomTool(component=imgplot,
                                    tool_mode="box", always_on=False))
        imgplotright.tools.append(PanTool(imgplotright, constrain_key="shift"))
        imgplotright.overlays.append(ZoomTool(component=self.right,
                                    tool_mode="box", always_on=False))

    def _update_image(self,number):
        if self.displayParameters.gain != 0:
            rgain=1./self.displayParameters.gain
        else:
            rgain=1
        print "rgain=",rgain
        range = DataRange1D(low=self.seis_data[number].minval*rgain,
                            high=self.seis_data[number].maxval*rgain)
        self.colormap = self.cmap(range)

        slice=transpose(self.seis_data[number].vals[self.slice_y, :, :])
        self.slice=slice
        colorslice=(self.colormap.map_screen(slice) * 255).astype(uint8)
        # Transposed required because img_plot() expects data in row-major order
#        self.arrayPlotData.set_data("xz", colorslicexz)
        self.arrayPlotDatas[number].set_data("xz", colorslice)

    def _update_images(self):
        for imageindex in range(len(self.seis_data)):
            self._update_image(imageindex)
        
    def _marker_size_changed(self):
        self.scatter.marker_size = self.marker_size

    def _color_changed(self):
        self.scatter.marker_size = self.marker_size

    def _display_button_fired(self):
        print "Display button pushed"
        self.displayParameters.edit_traits()
        self._update_images()

    def _prev_fired(self):
        print "prev button pushed"
        slice_y = self.slice_y - self.displayParameters.slice_inc 
        if(slice_y < 0):
            slice_y =  self.seis_data[0].vals.shape[0]-1
        print "after decrement slice_y=",slice_y
        self.slice_y=slice_y
        self._update_images()

    def _next_fired(self):
        print "next button pushed"
        slice_y = self.slice_y + self.displayParameters.slice_inc 
        print "shape=",self.seis_data[0].vals.shape
        if(slice_y >= self.seis_data[0].vals.shape[0]):
            slice_y = 0
        print "after increment slice_y=",slice_y
        self.slice_y=slice_y
        self._update_images()

    def _unzoom_fired(self):
        print "unzoom button pushed"
        print "self.bottomplot.range2d=",self.bottomplot.range2d
        print "xmin/xmax=", \
            self.bottomplot.range2d.x_range.low, \
            self.bottomplot.range2d.x_range.high
        print "ymin/ymax=", \
            self.bottomplot.range2d.y_range.low, \
            self.bottomplot.range2d.y_range.high

        self.bottomplot.range2d.x_range.low=self.seis_data[0].axis_start[1]
        self.bottomplot.range2d.x_range.high=self.seis_data[0].axis_end[1]
        self.bottomplot.range2d.y_range.low=self.seis_data[0].axis_start[2]
        self.bottomplot.range2d.y_range.high=self.seis_data[0].axis_end[2]
        

class DisplayParameters(HasTraits):
    gain = Range(low=.0001,high=10000, value=1.0)
    vaxis = Int(0)
    haxis = Int(1)
    slice_y = Int(63)
    slice_inc = Int(10)
    
class CustomTool(BaseTool):
    def normal_mouse_move(self, event):
        #  print "type event=",type(event),"Screen point:", event.x, event.y
        #print "Data:", self.component.map_data((event.x, event.y))
        x=int(round(self.component.x_mapper.map_data(event.x)))
        y=int(round(self.component.y_mapper.map_data(event.y)))
        #  print "Data:", x, y #, self.slice(x,y)
        #print "current_pointer_position", event.current_pointer_position
        #print "scale_xy", event.scale_xy
        #print "event", event


if __name__ == "__main__":
    myContainerExample=ContainerExample()
    myContainerExample.configure_traits()
