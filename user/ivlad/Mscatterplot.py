#!/usr/bin/env python
'''Scatter plot a 3-columns dataset with symbols colored according to value
Input should have three columns (i.e. n2=3): value, x, y (i.e. same as output
of sfsparsify). Essentially a m8r wrapper for Matplotlib's pyplot.scatter. See
http://matplotlib.sourceforge.net/api/pyplot_api.html#matplotlib.pyplot.scatter
and
http://matplotlib.sourceforge.net/examples/pylab_examples/scatter_star_poly.html

Quick test to check that it works:
sfspike n1=5 n2=3 nsp=3 k1=1,3,4 k2=1,2,3 mag=1,4,2 |\
sfsparsify nonzero=3 | sfscatterplot
'''

# Copyright (C) 2011 Ioan Vlad
#
# This program is free software; you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation; either version 2 of the License, or
# (at your option) any later version.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with this program; if not, write to the Free Software
# Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA

import numpy
import rsf.api as rsf
import matplotlib.pyplot as pyplot

par = rsf.Par()
input  = rsf.Input()
output = rsf.Output()
assert input.type == 'float'

n1 = input.int("n1")

size = par.int("size", 80)
 
v = numpy.zeros(n1, 'f')
x = numpy.zeros(n1, 'f')
y = numpy.zeros(n1, 'f')

input.read(v)
input.read(x)
input.read(y)

pyplot.scatter(x,y,s=size, c=v, marker=(0,3))
pyplot.show()
