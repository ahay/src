#! /usr/bin/env python
"""
NAME
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
	user/ivlad/Mpclip.py
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


import sys

try:
    import rsf
except: # Madagascar's Python API not installed
    import rsfbak as rsf

try: # Give precedence to local version
    from ivlad import send_to_os
except: # Use distributed version
    from rsfuser.ivlad import send_to_os

def main(argv=sys.argv):

    # Constants
    success = 0
    error   = 1

    # Parse arguments into a parameter table
    par = rsf.Par(argv)

    inp = par.string('inp')
    out = par.string('out')
    if (inp == None) or (out == None):
        sys.stderr.write(__doc__) # self-doc
        return error

    verb = par.bool('verb', False)
    pclip = par.float('pclip',99)

    if pclip <0 or pclip>100:
        sys.stderr.write('pclip must be between 0 and 100\n')
        return error

    clip = send_to_os('sfquantile', arg='pclip='+str(pclip), stdin=inp, want='stdout', verb=verb)

    if not clip:
        sys.stderr.write('sfquantile did not return anything!\n')
        return error

    send_to_os('sfclip', arg='clip='+clip, stdin=inp, stdout=out, verb=verb)

    return success

if __name__ == '__main__':
    sys.exit(main()) # Exit with the success or error code returned by main