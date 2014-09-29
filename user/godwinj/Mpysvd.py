#!/usr/bin/env python
''' Perform SVD on a matrix using SCIPY.  

REQUIRES the PYTHON API, NUMPY AND SCIPY
'''
import sys

# Import RSF API
try:
    import rsf.api as rsf
    import numpy
    import scipy
except Exception, e:
    print \
'''ERROR: NEED PYTHON API, NUMPY, SCIPY '''
    sys.exit(1)

# Initialize RSF command line parser    
par = rsf.Par()
# Read command line variables

# Declare input and outputs
fin = rsf.Input()    # no argument means stdin
fout = rsf.Output()  # no argument means stdout

# Declare optional inputs/outputs
vectors = par.bool  ("vectors", False     ) # Output singular vectors?
if vectors:
    left    = par.string("left") # File to store left singular vectors
    if not left:
        sys.stderr.write('Need left=\n')
        sys.exit(2)
    right   = par.string("right") # File to store right singular vectors
    if not right:
        sys.stderr.write('Need right=\n')
        sys.exit(2)
    lout = rsf.Output(left)  # left singular vectors
    rout = rsf.Output(right) # right singular vectors

# Get dimensions of input header or output header
n1 = fin.int('n1')
n2 = fin.int('n2')

data = numpy.zeros((n2,n1),'f') # Note, we reverse array dims

# Read our input data
fin.read(data)

# Perform our SVD
u,l,v = numpy.linalg.svd(data)

# Setup output headers
fout.put('n2',1)

# Write output data
fout.write(l)

if vectors:
    lout.put('n1',n2)
    lout.put('n2',n2)
    lout.write(u)
    lout.close()

    rout.put('n1',n1)
    rout.put('n2',n1)
    rout.write(v)
    rout.close()

# Clean up files
fout.close()
fin.close()
