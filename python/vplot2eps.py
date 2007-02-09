#!/usr/bin/env python
##   Copyright (C) 2006 University of Texas at Austin
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
import os, sys, string, re, tempfile

top = os.environ.get('RSFROOT')
bindir = os.path.join(top,'bin')
vppen = os.path.join(bindir,'vppen')    
pspen = os.path.join(bindir,'pspen')

def convert(vplot,eps,
            options='color=n fat=1 fatmult=1.5 invras=y',
            psborder=0,
            ppi=72 # points per inch resolution
            ):
    "Convert vplot to EPS"
    space = float(os.environ.get('PSBORDER',psborder))
    opts = os.environ.get(os.path.splitext(os.path.basename(eps))[0]+'.pspen',
                          os.environ.get('PSTEXPENOPTS',options))
    print opts
    # get bounding box
    getbb = vppen + ' big=n stat=l %s < %s | head -1' % (opts,vplot)
    
    out = os.popen(getbb)
    head = string.split(out.read())
    out.close() 

    bb = map(lambda x: int((float(head[x])-space)*ppi),[7, 12, 9, 14])

    out = open(eps,"w")
    out.write("%\!PS-Adobe-2.0 EPSF-2.0\n")
    out.write("%%%%BoundingBox: %d %d %d %d\n" % tuple(bb))
    
    name = tempfile.mktemp()
    command = pspen + ' size=a tex=y %s < %s > %s' % (opts,vplot,name)
    os.system(command)
    
    ps = open(name,'r')        
    out.write(ps.read())
    ps.close()
    os.unlink(name)
    
    out.write("\n")
    out.close()

if __name__ == "__main__":
    # own user interface instead of that provided by RSF's Python API
    # because this script has users that do not have RSF
    argc = len(sys.argv)
    prog = sys.argv.pop(0)
    
    if argc < 2:
        print '''
        Usage: %s [options] file.vpl [file.eps]
        Converts vplot to encapsulated postscript.
        [options] are passed to pspen.
        ''' % prog
        sys.exit(2)

    eps = sys.argv.pop()
    if eps[-2:] == 'ps':
        vpl = sys.argv.pop()
    else:
        vpl = eps
        eps = re.sub('\.[^\.]+$','.eps',vpl)

    try:
        if sys.argv:
            convert(vpl,eps,string.join(sys.argv,' '))
        else:
            convert(vpl,eps)
        sys.exit(0)
    except:
        sys.exit(1)
