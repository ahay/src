#!/usr/bin/env python

# Copyright (C) 2009 University of Texas at Austin
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
import os, sys, tempfile
import rsf.vplot2eps
import rsf.vplot2png
import rsf.vplot2gif
import rsf.vplot2avi
import rsf.prog

pens = {
    'vpl': 'vp',
    'eps': 'ps',
    'ps': 'ps',
    'ppm': 'ppm',
    'tiff': 'tiff',
    'tif': 'tiff',
    'jpeg': 'jpeg',
    'jpg': 'jpeg',
    'png': 'png',
    'gif': 'gd',
    'mpeg': 'gd',
    'mpg': 'gd',
    'pdf': 'ps',
    'svg': 'svg',
    'avi': 'ppm'
    }

formats = list(pens.keys())
formats.sort()

def exists(pen):
    '''check if a given pen exists'''
    bindir = os.path.join(rsf.prog.RSFROOT,'bin')
    exe = os.path.join(bindir,pen+'pen')
    if os.path.isfile(exe) and os.path.getsize(exe) > 500:
        return exe
    else:
        return None

def which(prog):
    '''find a program in path'''
    path = os.environ.get('PATH','')
    path = path.split(os.pathsep)

    for d in path:
        if sys.platform == 'cygwin' and d.find('WINDOWS') >= 0:
            continue
        exe = os.path.join(d, prog)
        if os.path.isfile(exe):
            return os.path.normpath(exe)
    return None

def convert(vpl,out,format,pen,args,verb=True):
    global pens, formats

    if not format in formats:
        print('Unknown format "%s" ' % format)
        sys.exit(1)
    if not os.path.isfile(vpl):
        print("\"%s\" is not a file" % vpl)
        sys.exit(1)
    if not pen:
        pen = pens[format]

    exe = exists(pen)
    convert = None

    if not exe:

        # Offer alternatives
        if format == 'png':
            mypens = ('gd','png','ps')
        elif format == 'gif':
            mypens = ('gd','ppm')
        elif format == 'jpeg' or format == 'jpg':
            mypens = ('jpeg','gd')
        elif format == 'pdf':
            mypens = ('ps','pdf')
        else:
            convert = which('convert')
            if convert:
                mypens = ('tiff','ppm','png','jpeg','ps')
            else:
                mypens = ()

        for other in mypens:
            if other != pen:
                exe = exists(other)
                if exe:
                    break
        if not exe:
            print('Unsupported program "%spen".' % pen)
            sys.exit(4)

        if convert:
            format2 = format
            if other == 'ps':
                format = 'eps'
            else:
                format = other
            out2 = out
            out = tempfile.mktemp(suffix='.'+format)

        pen = other

    if pen == 'gd':
        args += ' type=%s' % format

    if format == 'eps':
        crshift = False
        if args.find('cropshift=y') != -1:
           crshift = True;
        rsf.vplot2eps.convert(vpl,out,
                              options='color=n fat=1 fatmult=1.5 ' + args,
                              cropshift=crshift)
    elif format == 'png' and pen == 'ps':
        rsf.vplot2png.convert(vpl,out,
                              options='color=y fat=1 fatmult=1.5 ' + args)
    elif format == 'gif' and pen == 'ppm':
        rsf.vplot2gif.convert(vpl,out,args)
    elif format == 'avi':
        if not which('ffmpeg'):
            print("Conversion failed. Please install ffmpeg.")
            sys.exit(5)
        rsf.vplot2avi.convert(vpl,out)
    elif format == 'pdf' and pen == 'ps':
        eps = tempfile.mktemp(suffix='.eps')
        rsf.vplot2eps.convert(vpl,eps,
                              options='color=y fat=1 fatmult=1.5 ' + args)
        epstopdf = which('epstopdf') or which('a2ping') or which('convert')
        if epstopdf:
            command = 'LD_LIBRARY_PATH=%s GS_OPTIONS="%s" %s %s --outfile=%s' \
                % (os.environ.get('LD_LIBRARY_PATH',''),
                   os.environ.get('GS_OPTIONS',''),epstopdf,eps,out)
            print(command)
            fail = os.system(command)
        else:
            fail = True

        os.unlink(eps)

        if fail:
            raise RuntimeError('Cannot convert eps to pdf.') 
    else:
        # default behavior
        run = '%s %s %s > %s' % (exe,args,vpl,out)
        if verb:
            print(run)
        return os.system(run)

    if convert:
        run = '%s %s %s:%s' % (convert,out,format2,out2)
        if verb:
            print(run)
        return os.system(run)

if __name__ == "__main__":
    # own user interface instead of that provided by RSF's Python API
    # because this script may have users that do not have RSF

    usage = '''
Usage: %s [format=] file.vpl [file2.vpl file3.vpl | file.other]

Supported formats: %s
    ''' % (sys.argv[0],' '.join(formats))

    format = None
    pen = None

    files = []
    args = []

    for arg in sys.argv[1:]:
        if '=' in arg:
            if arg[:7] == 'format=':
                format = arg[7:].lower()
            elif arg[:4] == 'pen=':
                pen = arg[4:]
            else:
                args.append(arg)
        else:
            files.append(arg)
    args = ' '.join(args)

    if format != None:
        if not files:
            print(usage)
            sys.exit(1)
        for infile in files:
            # attach format as suffix
            if format == os.path.splitext(infile)[1][1:]:
                print(usage)
                sys.exit(1)
            outfile = os.path.splitext(infile)[0]+'.'+format
            convert(infile,outfile,format,pen,args)
    else:
        if len(files) !=2:
            print(usage)
            sys.exit(2)
        infile = files[0]
        outfile = files[1]
        # get format from suffix
        format = os.path.splitext(outfile)[1][1:].lower()
        if format == 'vpl':
            print('Trying to convert vpl to vpl will only destroy the input')
            sys.exit(3)
        else:
            convert(infile,outfile,format,pen,args)

    sys.exit(0)
