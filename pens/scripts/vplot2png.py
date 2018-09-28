#!/usr/bin/env python
##   Copyright (C) 2008 University of Texas at Austin
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
import tempfile, os, sys, re
import rsf.vplot2eps

def convert(vpl,png,**kw):
    'Vplot to PNG'
    eps = tempfile.mktemp()
    rsf.vplot2eps.convert(vpl,eps,
                          options=kw.get('options','color=y fat=1 fatmult=1.5'))

    option = kw.get('plotoption','')
    pstoimg = kw.get('pstoimg','pstoimg')
    command = 'PAPERSIZE=ledger %s %s -out %s' \
              + ' -type png -interlaced -antialias -quiet -crop a %s' 
    command = command % (pstoimg,eps,png,option)
    fail = os.system(command)
    os.unlink(eps)
    if fail:
        raise RuntimeError('cannot run "%s" ' % pstoimg)

if __name__ == "__main__":
    # own user interface instead of that provided by RSF's Python API
    # because this script has users that do not have RSF
    argc = len(sys.argv)
    prog = sys.argv.pop(0)
    
    if argc < 2:
        print('''
        Usage: %s [options] file.vpl [file.png]
        Converts vplot to PNG.
        [options] are passed to pspen.
        ''' % prog)
        sys.exit(2)

    png = sys.argv.pop()
    if png[-3:] == 'png':
        vpl = sys.argv.pop()
    else:
        vpl = png
        png = re.sub('\.[^\.]+$','.png',vpl)


    if sys.argv:
        convert(vpl,png,options=string.join(sys.argv,' '))
    else:
        convert(vpl,png)
    sys.exit(0)
  
