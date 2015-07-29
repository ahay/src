#!/bin/env python
# Make a list of unused subroutines in the current directory
import re, glob

cfiles = glob.glob('*.c')

includes = {}
for cfile in cfiles:
    if 'T' == cfile[0]:
        continue
    copen = open(cfile)
    for line in copen.readlines():
        gotit = re.search(r'#include\s+\"([^\.]+)\.h\"',line)
        if gotit:
            include = gotit.group(1)+'.c'
            if include != cfile and not includes.get(include):
                includes[include] = 1
    copen.close()

for cfile in cfiles:
    if 'M' != cfile[0] and 'T' != cfile[0] and not includes.get(cfile):
        print cfile
        
    
