#!/usr/bin/env python
''' Interactively displays a 3D cube of RSF data using Python + MayaVi2 + VTK.  

REQUIRES NUMPY, SCIPY, and MAyaVi2
'''
# Import RSF API
try:
    import rsf.api as rsf
    import numpy
    import scipy
    from mayavi import mlab
except Exception, e:
    import sys
    print '\nERROR: NEED NUMPY, SCIPY, and MAyaVi2\n', e
    sys.exit(1)

par = rsf.Par()
fin = rsf.Input()

data = numpy.zeros(fin.shape(),'f')

fin.read(data)

mlab.figure()

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
