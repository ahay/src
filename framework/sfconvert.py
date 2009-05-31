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
import os, sys

def convert(infile,outfile,format,args):
    if not os.path.isfile(infile):
        print "\"%s\" is not a file" % infile
        sys.exit(1)
    if format == 'vpl':
        run = 'vppen %s %s > %s' % (args,infile,outfile)
        print run
        os.system(run)
    else:
        print 'unsupported format "%s" ' % format

if __name__ == "__main__":
    # own user interface instead of that provided by RSF's Python API
    # because this script may have users that do not have RSF

    usage = '''
Usage: %s [format=] file.vpl [file2.vpl file3.vpl | file.other]

Supported formats: vpl, ps, eps, pdf, svg, ppm, png, jpeg, tiff
    ''' % sys.argv[0]

    format = None
    files = []
    args = []

    for arg in sys.argv[1:]:
        if '=' in arg:
            if arg[:7] == 'format=':
                format = arg[7:].lower()
            else:
                args.append(arg)
        else:
            files.append(arg)
    args = ' '.join(args)

    if format:
        if not files:
            print usage
            sys.exit(1)
        for infile in files:
            # attach format as suffix
            outfile = os.path.splitext(infile)[0]+'.'+format
            convert(infile,outfile,format,args)
    else:
        if len(files) !=2:
            print usage
            sys.exit(2)
        infile = files[0]
        outfile = files[1]
        # get format from suffix
        format = os.path.splitext(outfile)[1][1:].lower()
        convert(infile,outfile,format,args)

    sys.exit(0)
