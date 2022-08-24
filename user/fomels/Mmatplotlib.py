#!/usr/bin/env python
'Plotting RSF files with matplotlib'

##   Copyright (C) 2010 University of Texas at Austin
##  
##   This program is free software; you can redistribute it and/or modify
##   it under the terms of the GNU General Public License as published by
##   the Free Software Foundation; either version 2 of the License, or
##   (at your option) any later version.
##  
##   This program is distributed in the hope that it will be useful,
##   but WITHOUT ANY WARRANTY; without even the implied warranty of
##   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
##   GNU General Public License for more details.
##  
##   You should have received a copy of the GNU General Public License
##   along with this program; if not, write to the Free Software
##   Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA

import sys
import m8r
import numpy
import matplotlib.pyplot as plt

# self-documentation
if len(sys.argv) < 2:
    sys.stderr.write('Usage: %s <matplotlib function> <plot options> [format=eps] < inp.rsf [ > out.eps] \n\n' % sys.argv[0])
    sys.exit(1)

# matplotlib command
plot = sys.argv[1]

# default picture format
pformat = 'png'
grid = False
title= ''
xlabel= ''
ylabel= ''

# build parameter dictionary
args = {}
for a in sys.argv[2:]:
    key = a.split('=')[0]
    val =  a.replace(key+'=','')
    if key == 'format':
        pformat = val
    elif key == 'grid':
        if val == '1' or val[0] == 'y' or val[0] == 'Y':
            grid = True
    elif key == 'title':
        title = val
    elif key == 'xlabel':
        xlabel =val
    elif key == 'ylabel':
        ylabel =val
    else:
        args[key] = val

# read data
inp  = m8r.Input()
n1 = inp.int("n1")
d1 = inp.float("d1")
o1 = inp.float("o1")
n2 = inp.size(1)
if n2 > 1:
    d2 = inp.float("d2")
    o2 = inp.float("o2")
else:
    d2 = 1.0
    o2 = 0.0

data = numpy.zeros([n1,n2],'f')
inp.read(data)
inp.close()

x1 = numpy.transpose(numpy.tile(numpy.arange(o1, o1+n1*d1, d1,dtype='f'),(n2,1)))
x2 = numpy.tile(numpy.arange(o2, o2+n2*d2, d2,dtype='f'),(n1,1))

# recognize the plotting type
if plot == 'imshow':
    plt.imshow.__call__(data,**args)
elif plot == 'plot':
    plt.plot.__call__(x1,data,**args)
elif plot == 'semilogx':
    plt.semilogx.__call__(x1,data,**args)
elif plot == 'semilogy':
    plt.semilogy.__call__(x1,data,**args)
elif plot == 'loglog':
    plt.loglog.__call__(x1,data,**args)
elif plot == 'plot_surface':
    fig = plt.figure()
    ax = fig.add_subplot(111, projection='3d')
    ax.plot_surface.__call__(x1,x2,data,**args)
else:
    sys.stderr.write('Unrecognized plotting function "%s" \n\n' % plot)
    sys.exit(2)

if grid:
    plt.grid()
if title:
    plt.title(title)
if xlabel:
    plt.xlabel(xlabel)
if ylabel:
    plt.ylabel(ylabel)

# check if standard output
if sys.stdout.isatty():
    plt.show()
else:
    if sys.version_info[0] >= 3:
        plt.savefig(sys.stdout.buffer,format=pformat)
    else:
        plt.savefig(sys.stdout,format=pformat)
sys.exit(0)
