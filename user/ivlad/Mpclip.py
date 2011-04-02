#! /usr/bin/env python
'Percentile clip. Shell for sfclip(sfquantile(input)).'

# Copyright (C) 2007-2010 Ioan Vlad
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
    import ivlad, m8rex, sf
except: # Use distributed version
    import rsf.user.ivlad as ivlad
    import rsf.user.m8rex as m8rex
    import rsf.user.sf as sf

from rsf.prog import RSFROOT

###############################################################################

def main(par):

    inp = par.string('inp') # input file
    out = par.string('out') # output file

    verb = par.bool('verb', False) # if y, print system commands, outputs
    pclip = par.float('pclip',99)  # percentile clip

    if pclip < 0 or pclip > 100:
        raise m8rex.ParamOutOfRange('pclip',0,100)

    prog_nm_root = os.path.join(RSFROOT,'bin','sf')
    sfquantile = prog_nm_root + 'quantile'
    sfclip     = prog_nm_root + 'clip'

    clip = ivlad.getout('sfquantile', 'pclip='+str(pclip), inp, verb, True)
    sf.clip(inp, out, clip, verb)

    return ivlad.unix_success

###############################################################################

if __name__ == '__main__':
    ivlad.run(main, ['inp','out'])
