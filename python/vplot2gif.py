#!/usr/bin/env python

import sys, os, time, string, re, commands
 
def convert(infile,outfile):
    spacing = os.environ.get('GIFBORDER',0.25)
    ppi = os.environ.get('PPI',75)
    scale = os.environ.get('PPMSCALE',1)
    delay = os.environ.get('GIFDELAY',100)
    
    stat = string.split(
        commands.getoutput("vppen vpstyle=n stat=l < %s | head -1" % infile))
    total = commands.getoutput("vppen vpstyle=n stat=l < %s | tail -1" % infile)
    
    retot = re.compile('Total\s+(\d+)')
    match = retot.search(total)
    if match:
        frames = int(match.group(1))
    else:
        frames = 1
    
    xmin = float(stat[7]) - spacing
    xmax = float(stat[9]) + spacing
    xcen = (xmin+xmax)/2.

    ymin = float(stat[12]) - spacing
    ymax = float(stat[14]) + spacing
    ycen = (ymin+ymax)/2.

    width  = int((xmax-xmin)*ppi+0.5)
    height = int((ymax-ymin)*ppi+0.5)
    
    cwidth  = int((xmax-xmin-2*spacing)*ppi+0.5)
    cheight = int((ymax-ymin-2*spacing)*ppi+0.5)
    
    random = time.time()

    run = "vppen vpstyle=n outN=vppen.%%d.%s < %s >/dev/null" % (random,infile)
    sys.stderr.write(run+"\n")
    os.system(run)

    gifs = []
    for i in range(frames):
        vppen = 'vppen.%d.%s' % (i,random)
        gif = '%s.%d' % (outfile,i)
        gifs.append(gif)

        run = "ppmpen vpstyle=n break=i n2=%d n1=%d ppi=%d " \
              "xcenter=%d ycenter=%d %s | pnmscale %g  | " \
              "pnmquant 256 | pnmcrop | ppmtogif -interlace > %s" % \
              (height,width,ppi,xcen,ycen,vppen,scale,gif)
        sys.stderr.write(run+"\n")
        os.system (run)
        os.unlink(vppen)

    if outfile[-1] != '/':
        gifsicle = 'gifsicle --merge --loopcount=forever --optimize'
        run = '%s --delay=%d %s > %s' % (gifsicle,int(delay),
                                         string.join(gifs),outfile)
        print run
        os.system (run)

    map(os.unlink,gifs)

if __name__ == "__main__":
    argc = len(sys.argv)

    if argc < 2:
        print "No input"
        sys.exit(1)

    infile = sys.argv[1]
        
    if not os.path.isfile(infile):
        print "\"%s\" is not a file" % infile
        sys.exit(1)

    if argc < 3:
        outfile = os.path.splitext(infile)[0]+'.gif'
    else:
        outfile = sys.argv[2]

    convert(infile,outfile);
   
    sys.exit(0)

    
