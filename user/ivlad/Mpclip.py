#! /usr/bin/env python
'Percentile clip. Shell for sfclip(sfquantile(input)).'

# Copyright (C) 2007, 2009 Ioan Vlad
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

import sys, rsfprog, os

try:
    import rsf
except: # Madagascar's Python API not installed
    import rsfbak as rsf

try: # Give precedence to local version
    import ivlad
except: # Use distributed version
    import rsfuser.ivlad as ivlad

###############################################################################

# Use exceptions for elegant unit testing

class Error(Exception):
    '''Base class for exceptions in this module.'''
    def __str__(self):
        return 'Pclip exception: ' + self.message

class PclipWrongVal(Error):
    '''exception: pclip out of range'''
    def __init__(self):
        self.message = 'pclip must be between 0 and 100'

class NoReturnFromQuantile(Error):
    '''exception: sfquantile failed'''
    def __init__(self):
        self.message = 'sfquantile did not return anything'

# Find programs in RSF
RSFROOT = os.environ.get('RSFROOT')
bindir = os.path.join(RSFROOT,'bin')
sfquantile = os.path.join(bindir,'sfquantile')
sfclip = os.path.join(bindir,'sfclip')

###############################################################################

def main(argv=sys.argv):

    # Parse arguments into a parameter table
    par = rsf.Par(argv)

    inp = par.string('inp') # input file
    out = par.string('out') # output file
    if None in (inp, out):
        rsfprog.selfdoc()
        return

    verb = par.bool('verb', False) # if y, print system commands, outputs
    pclip = par.float('pclip',99)  # percentile clip

    if pclip <0 or pclip>100:
        raise PclipWrongVal

    clip = send_to_os(sfquantile, arg='pclip='+str(pclip),
                      stdin=inp, want='stdout', verb=verb)

    if not clip:
        raise NoReturnFromQuantile

    send_to_os(sfclip, arg='clip='+clip, stdin=inp, stdout=out, verb=verb)

    return ivlad.unix_success

###############################################################################

if __name__ == '__main__':

    try:
        status = main()
    except Error, e:
        print e.message
        status = ivlad.unix_error

    sys.exit(status)
