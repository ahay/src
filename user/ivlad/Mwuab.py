#!/usr/bin/env python
'Wrapper for Utilities Acting on Binaries'

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
    import rsf.user.ivlad as ivlad

###############################################################################

def main(par):

    prog = par.string('prog') # Non-madagascar utility
    inp  = par.string('inp') # Input file
    tpar_str = par.string('tpar') # Translated params, i.e.: "ni1=n1 ni2=n2"
    ipar_str = par.string('ipar') # Independent params, i.e. "perc=100 cmap=rgb"
    verb = par.bool('verb',False)

    tpar_out = ''
    if tpar_str != None:
        for tpar in tpar_str.split():
            thispar = tpar.split('=')
            tpar_out += thispar[0] + '=' + \
                ivlad.getout('sfget', ['parform=n', thispar[1]], inp, verb) +' '
           
    data = ivlad.getout('sfin',['info=n',inp])

    ivlad.exe(prog + ' <' + data + ' ' + tpar_out + ' ' + ipar_str, verb)
    
    return ivlad.unix_success

##############################################

if __name__ == '__main__':
    ivlad.run(main, ['prog', 'inp'])
