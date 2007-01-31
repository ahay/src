#! /usr/bin/env python
"""NAME
	sfpclip
DESCRIPTION
	Percentile clip. Shell for sfclip(sfquantile(input)).
SYNOPSIS
	sfpclip inp= out= verb=n pclip=99.0
PARAMETERS
	string inp=		Input file
	string out=		Output file
	bool   verb=n [y/n] 	If y, print system commands, outputs
	float  pclip=99.0	Percentile clip
SOURCE
	user/ivlad/pclip.py
"""
# Copyright (C) 2007 Ioan Vlad
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

from os import system
import sys

try:
    import rsf
except:
    import rsfbak as rsf

def main(argv=sys.argv):

    success = 0
    error   = 1
    par = rsf.Par(argv)

    inp = par.string('inp')
    out = par.string('out')
    if inp == None or out == None:
        sys.stderr.write(__doc__)
        return error

    verb  = par.bool('verb', False)
    pclip = par.float('pclip',99)

    if pclip <0 or pclip>100:
        sys.stderr.write('pclip must be between 0 and 100\n')
        return error

    inp = '< ' + inp
    out = ' > '+ out

    command = inp + ' sfquantile pclip=' + str(pclip)

    if verb:
        print(command)

    clip = getoutput( command )

    if not clip:
        sys.stderr.write('sfquantile did not return anything!\n')

    if verb:
        print(clip)

    command = inp + ' sfclip clip=' + clip + out

    if verb:
        print(command)

    system(command)

    return success

if __name__ == '__main__':
    sys.exit(main()) # Exit with the success or error code returned by main