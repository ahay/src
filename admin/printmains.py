#!/bin/env python
# Print a list of main programs

import glob
mains = []
for main in glob.glob('M*.c'):
    mains.append(main[1:-2])
mains.sort()
print ' '.join(mains)
