#!/usr/bin/env python
import os
from traits.api import HasTraits, Enum, Range, on_trait_change
from traitsui.api import View, Item

class Wavelet(HasTraits):
    wtype = Enum(['Bi-orthogonal','Linear','Haar'])
    pclip = Range(0.0, 100.0, value=50.0)

    traits_view = View(
       Item('wtype',label='Wavelet Type'),
       Item('pclip',label='Thresholding'),
       buttons=['OK']
    )

    @on_trait_change('wtype,pclip')
    def scons(self):
        'Get parameters from GUI and pass them to SCons'
        command = 'scons -Q type=%s pclip=%d view' % (self.wtype.lower(),
                                                      self.pclip)
        os.system(command)

if __name__ == "__main__":
    w = Wavelet()
    w.configure_traits()
    
