#! /usr/bin/env python
'''Copies header of float file, fills output binary with zeros
Usage: sfzcp file1.rsf file2.rsf'''

# Copyright (C) 2009-2010 Ioan Vlad
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

import os, sys

try: # Give precedence to local version
    import ivlad, m8rex
except: # Use distributed version
    import rsf.user.ivlad as ivlad
    import rsf.user.m8rex as m8rex

###############################################################################

def main(par):

    f1 = sys.argv[1] + ' '
    f2 = sys.argv[2]

    f1_type = ivlad.getout('sfgettype', stdin=f1)
    if f1_type != 'SF_FLOAT':
        raise m8rex.TypeHandlingNotImplemented(f1_type)

    dims_str = ivlad.getout('sffiledims','parform=n',f1)
    ndims = int(dims_str.split(':')[0])

    cmd = 'sfspike mag=0 '
    fields = 'n o d unit label'.split()

    for i in range(ndims):
        for z in fields:            
            cmd += ivlad.getout('sfget',['parform=y',z+str(i+1)],f1) + ' '

    cmd += ' >' + f2

    ivlad.exe(cmd)

    return ivlad.unix_success

###############################################################################

if __name__ == '__main__':
    ivlad.run(main, nminarg=2, nmaxarg=2)
