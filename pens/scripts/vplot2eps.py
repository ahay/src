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
import rsf.prog

top = rsf.prog.RSFROOT
bindir = os.path.join(top,'bin')
vppen = os.path.join(bindir,'vppen')    
pspen = os.path.join(bindir,'pspen')

def convert(vplot,eps,
            options='color=n fat=1 fatmult=1.5 invras=y',
            psborder=0.1,
            ppi=72,                     # points per inch resolution
            xbmin=None, xbmax=None,     # options to override default bounding box
            ybmin=None, ybmax=None,
            cropshift=False
            ):
    "Convert vplot to EPS"
    
    space = float(os.environ.get('PSBORDER',psborder))
    opts = os.environ.get('PSTEXPENOPTS',options)

    # Get bounding box info from vplot file
    getbb = vppen + ' big=y stat=l %s < %s | head -1' % (opts,vplot)
    
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
    bm = [float(x)-space for x in [xbbm,ybbm]]
    bp = [float(x)+space for x in [xbbp,ybbp]]
    
    bbm = [x*ppi for x in bm]
    bbp = [x*ppi for x in bp]

    # Round to integer
    ibbm = list(map(int,bbm))
    ibbp = [int(x)+1 for x in bbp]
    
    out = open(eps,"w")
    out.write("%!PS-Adobe-2.0 EPSF-2.0\n")
    out.write("%%%%BoundingBox: %d %d %d %d\n" % tuple(ibbm+ibbp))
    out.write("%%%%HiResBoundingBox: %g %g %g %g\n" % tuple(bbm+bbp))

    opts += ' xwmin=%g ywmin=%g xwmax=%g ywmax=%g' % tuple(bm+bp)
    
    name = tempfile.mktemp()
    command = pspen + ' size=a tex=y %s < %s > %s' % (opts,vplot,name)
    os.system(command)

    ps = open(name,'r')        

    if (cropshift):
        # For multipage EPS files, the afterward conversion to PDF will cause
        # pages appear to be shifted (except for the very first one).
        # It is probably due to a bug in GS, which happens when the latter is
        # called with -dEPSCrop option from within epstopdf. We apply a 
        # PS translate here in the beginning of every page to counteract that
        # bug.
        for line in ps.readlines ():
            if line.find('showpage') != -1:
                out.write(line)
                out.write('%g %g translate\n' % (-bbm[0],-bbm[1]))
            else:
                out.write(line)
    else:
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
        print('''
Usage:
%s [xbmin= xbmax= ybmin= ybmax=] [cropshift=n] [options] file.vpl [file.eps]

Converts vplot to encapsulated postscript.
[xbmin= xbmax= ybmin= ybmax=] overrides the default bounding box.
[cropshift=y] applies PS translate for every page (except for the first one)
[options] are passed to pspen.
        ''' % prog)
        sys.exit(2)

    eps = sys.argv.pop()
    if eps[-2:] == 'ps':
        vpl = sys.argv.pop()
    else:
        vpl = eps
        eps = re.sub('\.[^\.]+$','.eps',vpl)

    options = []
    xbmin=None; xbmax=None; ybmin=None; ybmax=None
    crshift=False
    for item in sys.argv:
        if   item.split('=')[0] == 'xbmin': xbmin = item.split('=')[1] 
        elif item.split('=')[0] == 'xbmax': xbmax = item.split('=')[1] 
        elif item.split('=')[0] == 'ybmin': ybmin = item.split('=')[1] 
        elif item.split('=')[0] == 'ybmax': ybmax = item.split('=')[1]
        elif (item.split('=')[0] == 'cropshift' and
              item.split('=')[1] == 'y'):   crshift=True
        else:                               options.append(item)
    try:
        if options:
            convert(vpl,eps,options=string.join(options,' '),
                    xbmin=xbmin, xbmax=xbmax, ybmin=ybmin, ybmax=ybmax,
                    cropshift=crshift)
        else:
            convert(vpl,eps,
                    xbmin=xbmin, xbmax=xbmax, ybmin=ybmin, ybmax=ybmax,
                    cropshift=crshift)
    except:
        print('Failed to convert %s to %s' % (vpl,eps))
        sys.exit(1)
