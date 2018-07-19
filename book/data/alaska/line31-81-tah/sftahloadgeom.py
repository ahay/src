#!/usr/bin/env python

import numpy as np
import m8r
import sys

import math
 
class InterpText:
    def __init__(self):
        self.xin=[]
        self.yin=[]
    def append(self,x,y):
#        print "append x=" + str(x) + "y=" + str(y)
        self.xin.append(float(x))
        self.yin.append(float(y))
    def read(self,fileName):
#        print "start InterText line loop"
        for self.line in open(fileName):
#            print "line=" + self.line
            tokens=self.line.split()
            self.append(float(tokens[0]),float(tokens[1]))
    def linearInterp(self,x):
        x=float(x)
#        print "compare x=" + str(x) 
#        print "xin[0]=" + xin[0]
        if(x<=self.xin[0]):
            return self.yin[0]
        #print "len(self.xin)=" + str(len(self.xin));
        for i in range(0,len(self.xin)-1):
            #print "test x<=self.xin[i+1]" + str(x) + "<=" + str(self.xin[i+1])
            if(x<=self.xin[i+1]):
                #print "x=" + str(x) + "i=" + str(i)
                return self.yin[i]+ \
                    (x-self.xin[i])*(self.yin[i+1]-self.yin[i])/ \
                      (self.xin[i+1]-self.xin[i])
        return self.yin[-1]

# create InterpText object to interpolate spnElev.txt
spnElev = InterpText()
spnElev.read('spnElev.txt')
# create InterpText object to interpolate recnoSpn.txt
recnoSpn = InterpText()
recnoSpn.read('recnoSpn.txt')

# set up to read parameters and read and write the seismic data
par = m8r.Par()
inp  = m8r.Input()
output = m8r.Output("out")

# get locations of the header keys used to initialize geometry
#input headers
indx_fldr =inp.get_segy_keyindx('fldr')
indx_tracf=inp.get_segy_keyindx('tracf')
#output headers
indx_ep   = inp.get_segy_keyindx('ep')
indx_sx    =inp.get_segy_keyindx('sx')
indx_sy    =inp.get_segy_keyindx('sy')
indx_gx    =inp.get_segy_keyindx('gx')
indx_gy    =inp.get_segy_keyindx('gy')
indx_cdp   =inp.get_segy_keyindx('cdp')
indx_offset=inp.get_segy_keyindx('offset')
indx_sdepth=inp.get_segy_keyindx('sdepth')
indx_selev =inp.get_segy_keyindx('selev')
indx_gelev =inp.get_segy_keyindx('gelev')
indx_sstat =inp.get_segy_keyindx('sstat')
indx_gstat =inp.get_segy_keyindx('gstat')
indx_tstat =inp.get_segy_keyindx('tstat')
indx_cdpx  =inp.get_segy_keyindx('cdpx')
indx_cdpy  =inp.get_segy_keyindx('cdpy')

#sys.stderr.write('indx_fldr=%s\n'%repr(indx_fldr))
#sys.stderr.write('indx_ep=%s\n'%repr(indx_ep))
#sys.stderr.write('indx_cdp=%s\n'%repr(indx_cdp))

# loop over input traces.  Read trace and header, compute headers, write
while True:
    eof_get_tah, intrace, inheader = inp.get_tah()
    if eof_get_tah:
        break

    fldr =inheader[indx_fldr ]
    tracf=inheader[indx_tracf]

    # set up header values with nulls.
    ep=-999999
    sx    =-999999999
    sy    =-999999999
    gx    =-999999999
    gy    =-999999999
    cdp   =-999999
    outtracf=-999999
    offset=-999999
    shotSpn=-999999999
    gSpn=-999999999
    sdepth=-999
    selev=-999
    gelev=-999
    sstat=-999
    gstat=-999
    tstat=-999
    cdpx=-999999999
    cdpy=-999999999
    # compute the shot number corresponding to this fldr (or recno in observer
    # log terminology
    shot=long(round(recnoSpn.linearInterp(fldr)))
    
    # If shot and tracf are valid, interpolate header values.
    # The recnoSpn.txt file sets shot(fldr) values to -999999.  This
    # causes the headers to be set to dummy values and later they will be
    # removed using suwind/sftahwindow.
    if(shot>0 and tracf<97):
        ep=shot
        # artificial coordinates are computed for all the traces.  This avoids
        # a lot of problems with projection from lat/long to coordinates.  
        # That standard projection for the National Petroleum Reserve Alaska
        # produces xy coordinates that are distorted enough that 
        # sqrt(deltax^2 + deltay^2) is a couple percent different from the
        # distances measured on the spheroid.  (You cannot project from a 
        # spheroid to a flat map without distortion)  These distortions would
        # make problems in the (shot,group) to (cdp, offset) map.  This is 
        # all avoided my artificial coordinates for seismic processing then
        # loading real world coordinates in the processed CDP data.

        # smaller tracf are groups (receivers) with small x coordinates than
        # the shot.
        r_offset=(float(tracf)-((96.0+1.0)/2.0))*110.0
        offset=long(round(r_offset))
        # shot x is shot number * 440 shot interval
        r_sx=float(ep)*440.0
        sx=long(round(r_sx))
        # shot and group y coordinates are all 0
        sy=0
        gx=long(round(r_sx+r_offset))
        gy=0
        # cdp computed from shot/group midpoint 
        cdp=long(round((float(sx+gx)/2.0-55.0/2.0)/55.0-752+101))
        outtracf=long(round((gx/440.0+55.0/440.0)*4.0))
        shotSpn=float(sx)/440.0
        gSpn=float(gx)/440.0

        # elevation at the shot
        selev=long(round(spnElev.linearInterp(shotSpn)))
        # elevation at the group (receiver)
        gelev=long(round(spnElev.linearInterp(gSpn)))

        # observer log says shot holes 82.5 ft before spn 123, then changes to
        # 67.5 ft.
        if(shotSpn<123):
            sdepth=long(round(82.5))
        else:
            sdepth=long(round(67.5))

        # Compute shot and receiver statics to sea level using 10000 ft/s
        # Convert to ms.
        sstat=long(round((-selev+sdepth)*1000.0/10000.0))
        gstat=long(round(-gelev*1000.0/10000.0))
        tstat=long(round((-gelev-selev+sdepth)*1000.0/10000.0))

        cdpx =long(round((sx+gx)/2))
        cdpy =long(round((sy+gy)/2))
    # Code paths for bad and good traces merge again here.
    # Print the header values to standard output.
    #keys='ep,sx,sy,gx,gy,cdp,tracf,offset,' 
    #  'selev,gelev,sstat,gstat,tstat'

    inheader[indx_ep    ]=ep
    inheader[indx_sx    ]=sx
    inheader[indx_sy    ]=sy
    inheader[indx_gx    ]=gx
    inheader[indx_gy    ]=gy
    inheader[indx_cdp   ]=cdp
    inheader[indx_tracf ]=outtracf
    inheader[indx_offset]=offset
    inheader[indx_selev ]=selev
    inheader[indx_gelev ]=gelev
    inheader[indx_sstat ]=sstat
    inheader[indx_gstat ]=gstat
    inheader[indx_tstat ]=tstat
    inheader[indx_cdpx  ]=cdpx
    inheader[indx_cdpy  ]=cdpy

    output.put_tah(intrace,inheader)
