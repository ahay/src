#!/usr/bin/env python

# Copyright (C) 2007 University of Texas at Austin
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

import os, sys, time
import rsf.vplot2gif

def convert(infile,outfile):
    gif = 'tmp%s.gif' %  time.time()    
    rsf.vplot2gif.convert(infile,gif)
    run = 'ffmpeg -i %s %s' % (gif,outfile)
    os.system(run)
    os.unlink(gif)

if __name__ == "__main__":
    # own user interface instead of that provided by RSF's Python API
    # because this script has users that do not have RSF
    argc = len(sys.argv)

    if argc < 2:
        print "No input"
        sys.exit(1)

    infile = sys.argv[1]

    if not os.path.isfile(infile):
        print "\"%s\" is not a file" % infile
        sys.exit(1)

    if argc < 3:
        outfile = os.path.splitext(infile)[0]+'.avi'
    else:
        outfile = sys.argv[2]

    convert(infile,outfile);
    sys.exit(0)
