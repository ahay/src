#!/bin/env python
# List of unused main programs
import glob

import rsf.doc
import rsf.prog

progs = rsf.doc.progs

for main in glob.glob('M*.c'):
    sfmain = 'sf'+main[1:-2]
    if not sfmain in progs.keys():
        print main

for main in glob.glob('M*.cc'):
    sfmain = 'sf'+main[1:-3]
    if not sfmain in progs.keys():
        print main

for main in glob.glob('M*.cu'):
    sfmain = 'sf'+main[1:-3]
    if not sfmain in progs.keys():
        print main

for main in glob.glob('M*.f90'):
    sfmain = 'sf'+main[1:-4]
    if not sfmain in progs.keys():
        print main
