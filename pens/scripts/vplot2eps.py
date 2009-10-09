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
            psborder=0.1,
            ppi=72,                     # points per inch resolution
            xbmin=None, xbmax=None,     # options to override default bounding box
            ybmin=None, ybmax=None
            ):
    "Convert vplot to EPS"
    
    space = float(os.environ.get('PSBORDER',psborder))
    opts = os.environ.get('PSTEXPENOPTS',options)

    # Get bounding box info from vplot file
    getbb = vppen + ' big=n stat=l %s < %s | head -1' % (opts,vplot)
    
    out = os.popen(getbb)
    head = string.split(out.read())
    out.close() 
    
    # Use parameters for bounding box if supplied.
    # Otherwise use bounding box from vplot file.
    if xbmin==None: xbbm = head[7]
    else:           xbbm = xbmin
    if xbmax==None: xbbp = head[9]
    else:           xbbp = xbmax

    if ybmin==None: ybbm = head[12]
    else:           ybbm = ybmin
    if ybmax==None: ybbp = head[14]
    else:           ybbp = ybmax

    # Compute bounding box
    bbm = map(lambda x: (float(x)-space)*ppi,[xbbm,ybbm])
    bbp = map(lambda x: (float(x)+space)*ppi,[xbbp,ybbp])

    # Round to integer
    ibbm = map(int,bbm)
    ibbp = map(lambda x: int(x)+1,bbp)
    
    out = open(eps,"w")
    out.write("%!PS-Adobe-2.0 EPSF-2.0\n")
    out.write("%%%%BoundingBox: %d %d %d %d\n" % tuple(ibbm+ibbp))
    out.write("%%%%HiResBoundingBox: %g %g %g %g\n" % tuple(bbm+bbp))
    
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
Usage:
%s [xbmin= xbmax= ybmin= ybmax=] [options] file.vpl [file.eps]

Converts vplot to encapsulated postscript.
[xbmin= xbmax= ybmin= ybmax=] overrides the default bounding box.
[options] are passed to pspen.
        ''' % prog
        sys.exit(2)

    eps = sys.argv.pop()
    if eps[-2:] == 'ps':
        vpl = sys.argv.pop()
    else:
        vpl = eps
        eps = re.sub('\.[^\.]+$','.eps',vpl)

    options = []
    xbmin=None; xbmax=None; ybmin=None; ybmax=None
    for item in sys.argv:
        if   item.split('=')[0] == 'xbmin': xbmin = item.split('=')[1] 
        elif item.split('=')[0] == 'xbmax': xbmax = item.split('=')[1] 
        elif item.split('=')[0] == 'ybmin': ybmin = item.split('=')[1] 
        elif item.split('=')[0] == 'ybmax': ybmax = item.split('=')[1]
        else:                               options.append(item)
            
    try:
        if options:
            convert(vpl,eps,options=string.join(options,' '),
                    xbmin=xbmin, xbmax=xbmax, ybmin=ybmin, ybmax=ybmax)
        else:
            convert(vpl,eps,
                    xbmin=xbmin, xbmax=xbmax, ybmin=ybmin, ybmax=ybmax)
    except:
        print 'Failed to convert %s to %s' % (vpl,eps)
        sys.exit(1)
