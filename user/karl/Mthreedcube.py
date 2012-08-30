#!/usr/bin/env python
''' Interactively displays a 3D cube of RSF data using Python + MayaVi2 + VTK.  

REQUIRES the PYTHON API, NUMPY AND SCIPY, MAyaVi2, VTK
'''
# Import RSF API
try:
    import rsf.api as rsf
    import numpy
    import scipy
    from enthought.mayavi import mlab
except Exception, e:
    import sys
    print '\nERROR: NEED PYTHON API, NUMPY, SCIPY, MAYAVI\n', e
    sys.exit(1)

class Header:
    def __init__(self, RSFFile):
        self.dims = {}
        self.shape = []
        for i in range(1,4):
            nkey = 'n%d'%i
            dkey = 'd%d'%i
            okey = 'o%d'%i
            self.dims[nkey] = RSFFile.int(nkey)
            self.dims[dkey] = RSFFile.float(dkey)
            self.dims[okey] = RSFFile.float(okey)
            self.shape.append(RSFFile.int(nkey))
        self.shape.reverse()
# end class def    

par = rsf.Par()

fin = rsf.Input()
header = Header(fin)

data = numpy.zeros(header.shape,'f')

fin.read(data)

mlab.pipeline.image_plane_widget(mlab.pipeline.scalar_field(data),
    plane_orientation='x_axes',
    slice_index=10)
    
mlab.pipeline.image_plane_widget(mlab.pipeline.scalar_field(data),
    plane_orientation='y_axes',
    slice_index=10)

mlab.pipeline.image_plane_widget(mlab.pipeline.scalar_field(data),
    plane_orientation='z_axes',
    slice_index=10)
#mlab.outline()
mlab.show()
