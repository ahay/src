"""
This example creates a simple 1D function examiner, illustrating the use of
ChacoPlotEditors for displaying simple plot relations, as well as TraitsUI
integration. Any 1D numpy/scipy.special function should work in the function
text box.
 - Left-drag pans the plot.
 - Mousewheel up and down zooms the plot in and out.
 - Pressing "z" brings up the Zoom Box, and you can click-drag a rectangular
   region to zoom.  If you use a sequence of zoom boxes, pressing alt-left-arrow
   and alt-right-arrow moves you forwards and backwards through the "zoom
   history".
"""

# Major library imports
from numpy import linspace, pi

# Enthought library imports
from traits.api import Array, Dict, Enum, HasTraits, Str
from traitsui.api import Item, View

# Chaco imports
from chaco.chaco_plot_editor import ChacoPlotEditor, \
                                                ChacoPlotItem


class Foo(HasTraits):

    # Public Traits
    xdata = Array
    plot_type = Enum("scatter", "line")
    eq = Str("sin(x)")

    # Default TraitsUI view
    traits_view = View(
                       ChacoPlotItem("xdata", "_ydata",
                                      type_trait="plot_type",

                                      # Basic axis and label properties
                                      show_label=False,
                                      resizable=True,
                                      orientation="h",
                                      x_label = "Index data",
                                      y_label = "Value data",

                                      # Plot properties
                                      color = "green",
                                      bgcolor = "white",

                                      # Specific to scatter plot
                                      marker = "circle",
                                      marker_size = 2,
                                      outline_color = "none",

                                      # Border, padding properties
                                      border_visible=True,
                                      border_width=1,
                                      padding_bg_color = "lightgray"),
                Item("plot_type"),
                Item("eq"),
                resizable=True,
                width=500, height=500)


    # Private Traits
    _d = Dict
    _ydata = Array

    def __init__(self, **kwtraits):
        super(Foo, self).__init__(**kwtraits)
        self._d = dict(x=self.xdata)
        exec "from scipy import *" in self._d
        exec "from scipy.special import *" in self._d
        self._ydata = eval(self.eq, self._d)

    def _eq_changed(self, old, new):
        try:
            self._ydata = eval(new, self._d)
        except:
            pass

#===============================================================================
# # demo object that is used by the demo.py application.
#===============================================================================
demo = Foo(xdata=linspace(-2*pi, 2*pi ,100), eq="sin(x)")

if __name__ == "__main__":
    demo.edit_traits(kind="modal")
