#!/bin/env python
# Print a list of main programs

import glob

mains = []
for main in glob.glob('M*.c'):
    mains.append(main[1:-2])
if mains:
    mains.sort()
    print "C:", ' '.join(mains)

mains = []
for main in glob.glob('M*.cc'):
    mains.append(main[1:-3])
if mains:
    mains.sort()
    print "C++:", ' '.join(mains)

mains = []
for main in glob.glob('M*.cu'):
    mains.append(main[1:-3])
if mains:
    mains.sort()
    print "CUDA:", ' '.join(mains)

mains = []
for main in glob.glob('M*.f90'):
    mains.append(main[1:-4])
if mains:
    mains.sort()
    print "F90:", ' '.join(mains)

mains = []
for main in glob.glob('M*.py'):
    mains.append(main[1:-3])
if mains:
    mains.sort()
    print "Python:", ' '.join(mains)
