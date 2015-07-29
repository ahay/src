#!/usr/bin/env python

##   Copyright (C) 1987 The Board of Trustees of Stanford University
##   Copyright (C) 2004 University of Texas at Austin
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

# Modified from the original C-shell version by Joe Dellinger

import sys, os, time, string, re, shutil
import rsf.prog

def convert(infile,outfile,args=''):
    spacing = float(os.environ.get('GIFBORDER',0.1))
    ppi = int(os.environ.get('PPI',75))
    delay = int(os.environ.get('GIFDELAY',100))

    bindir = os.path.join(rsf.prog.RSFROOT,'bin')
    vppen  = os.path.join(bindir,'vppen') + ' ' + args
    ppmpen = os.path.join(bindir,'ppmpen')
    
    # Use vppen to find out how big and where on the page the plot is.
    stats = os.popen(vppen + ' size=a stat=l < %s' % infile)
    lines = stats.readlines()
    stats.close()

    stat = string.split(lines[0])

    # find the number of frames
    retot = re.compile('Total\s+(\d+)')
    match = retot.search(lines[-1])
    if match:
        frames = int(match.group(1))
    else:
        frames = 1

    if frames > 1:
        import commands
        which_gifsicle = commands.getoutput('which gifsicle')
        if which_gifsicle[:18] == 'which: no gifsicle':
            sys.stderr.write('Missing program: gifsicle\n')
            sys.exit(1)

    xmin = float(stat[7]) - spacing
    xmax = float(stat[9]) + spacing
    xcen = (xmin+xmax)/2.

    ymin = float(stat[12]) - spacing
    ymax = float(stat[14]) + spacing
    ycen = (ymin+ymax)/2.

    width  = int((xmax-xmin)*ppi+0.9999)
    height = int((ymax-ymin)*ppi+0.9999)

    sys.stderr.write('''
    %s will be %d pixels wide, %d pixels tall,
    at %d pixels per inch, with borders %g inches wide.
    ''' % (outfile,width,height,ppi,spacing))

    random = time.time()
    run = vppen + ' size=a outN=vppen.%%d.%s < %s >/dev/null' % (random,infile)
    os.system(run)

    gifs = []
    for i in range(frames):
        vppen = 'vppen.%d.%s' % (i,random)
        gif = '%s.%d' % (outfile,i)
        gifs.append(gif)

        run = ppmpen + ' break=i n1=%d n2=%d ppi=%d size=a ' \
              'xcenter=%g ycenter=%g %s | ' \
              'ppmquant 256 | ppmtogif > %s' % \
              (width,height,ppi,xcen,ycen,vppen,gif)
        os.system (run)
        os.unlink(vppen)

    if frames > 1:
        if not os.path.isdir(outfile): # if not directory
            # combine frames into an animated gif (requires gifsicle)
            gifsicle = 'gifsicle --merge --loopcount=forever --optimize'
            run = '%s --delay=%d %s > %s' % (gifsicle,int(delay),
                                             string.join(gifs),outfile)
            os.system (run)
            map(os.unlink,gifs)
    else:
        shutil.move(gifs[0],outfile)

if __name__ == "__main__":
    # own user interface instead of that provided by RSF's Python API
    # because this script has users that do not have RSF
    argc = len(sys.argv)

    if argc < 2:
        print '''
        vplot2gif myplot.vpl [myplot.gif]

        Convert Vplot format to GIF format at 75 dots per inch, ideal
	for WWW documents. Default output is input file name with
	suffix changed to ".gif".

	vplot2gif finds the smallest bounding box containing the input
	plot, and makes a gif output just big enough to contain that
	image, plus a border of .25 inches all around. (The position
	of the plot on the vplot virtual page is irrelevant!)

	You can override the 75. dots per inch by setting the
	environment variable PPI.

	You can override the .25 inch border by setting the
	environment variable GIFBORDER. (Note, GIFBORDER will only
	come out in physical units if PPI happens to be accurate for
	your display device.)

        In case of multiple frames, vplot2gif creates an animated
        GIF. The time delay between movie frames is 100 ms (controlled
        by GIFDELAY environmental variable).

        If the output file is a directory, the individual frames are
        not combined into an animation.
        '''

        sys.exit(1)

    infile = sys.argv[1]

    if not os.path.isfile(infile):
        print "\"%s\" is not a file" % infile
        sys.exit(1)

    if argc < 3:
        narg = 2
        outfile = os.path.splitext(infile)[0]+'.gif'
    else:
        narg = 3
        outfile = sys.argv[2]

    convert(infile,outfile,' '.join(sys.argv[narg:]))

    sys.exit(0)
