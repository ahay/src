#!/usr/bin/env python
'''Displays a 2-D RSF file with Seismic Unix's ximage
Test with:
sfspike n1=5 n2=3 nsp=3 k1=1,3,4 k2=1,2,3 > junk.rsf
sfximage inp=junk.rsf par="perc=100 cmap=rgb1 legend=1"
You should see a picture with blue background and red blobs.

See also sfimage.
'''
# Copyright (C) 2010 Ioan Vlad
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

try: # Give precedence to local version
    import ivlad
except: # Use distributed version
    import rsfuser.ivlad as ivlad

###############################################################################

def main(par):

    inp = par.string('inp') # Input file
    ipar = par.string('par','') # ximage params that can't be found in RSF headr
    verb = par.bool('verb',False)

    prog = 'ximage' 
    tpar = 'n1=n1 n2=n2 f1=o1 f2=o2 d1=d1 d2=d2 label1=label1 label2=label2'
    vflag = ivlad.switch(par.bool('verb',False), 'y', 'n')

    cmd = 'sfwuab inp=%s prog=%s tpar="%s" ipar="%s" verb=%s' %\
        (inp, prog, tpar, ipar, vflag)

    ivlad.exe(cmd, verb)
    
    return ivlad.unix_success

##############################################

ivlad.run(__name__, main, ['inp'])
