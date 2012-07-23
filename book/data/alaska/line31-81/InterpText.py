#!/usr/bin/env python
"""
custom code used to build headers for alaska line 31-81.  This hardwires things
like the shot, group, and cdp intervals, number groups and much more.  You may
you is as a guide line to write custom code for other data.  I found the
InterpText class reuseable.

usage:
Create input files:
hdrfile.txt output from <line.su sugethw output=geom key=tracl,fldr,tracf
spnElev.txt    Column 1 is the spn.  Column 2 is elevation.  Typed from info in 
               surveyors log
recnoSpn.txt   Column 1 is fldr. Column 2 is spn.  fldr is called recno in the
               observers log.  File was typed using info in the observers log
 
The program is run by:
./InterpText.py > hdrfile1.txt

The hdrfile1.txt file is processed using a2b and sushw to load the headers in
file line.su and output allshots.su 
"""
 
import math
 
class InterpText:
    """

    functions to read (xin,yin) pairs and linearly interpolate y(xout)
    xin is assumed to be monotonically increasing

    """

    def read(self,fileName):
        """

        Read pairs from fileName.  Create xin,yin arrays.  xin must be 
        monotonically increasing.  There is no validation for two floats
        per line or xin monotopically increasing

        """

        self.xin=[]
        self.yin=[]
        for self.line in open(fileName):
            tokens=self.line.split()
            self.xin.append(float(tokens[0]))
            self.yin.append(float(tokens[1]))
#        print self.xin
#        print self.yin
    def linearInterp(self,x):
        """

        Find smallest such that x is in the interval [xin[i],xin[i+1]]
        Linear interpolate y(x).  Copy first or last yin for values outside
        the range of xin array.

        """

        if(x<=self.xin[0]):
            return self.yin[0]
        #print "len(self.xin)=" + str(len(self.xin));
        for i in range(0,len(self.xin)-1):
            #print "test x<=self.xin[i+1]" + str(x) + "<=" + str(self.xin[i+1])
            if(x<=self.xin[i+1]):
                #print "x=" + str(x) + "i=" + str(i)
                return self.yin[i]+(x-self.xin[i])*                           \
                       (self.yin[i+1]-self.yin[i])/(self.xin[i+1]-self.xin[i])
        return self.yin[len(self.yin)-1]

import sys

fspnElev = sys.argv[1]

# create InterpText object to interpolate spnElev.txt
spnElev = InterpText()
spnElev.read(fspnElev)

frecnoSpn = sys.argv[2]

# create InterpText object to interpolate recnoSpn.txt
recnoSpn = InterpText()
recnoSpn.read(frecnoSpn)

hdrfile = sys.argv[3]

for line in open(hdrfile):  # read each line in the hdrfile.txt file
    tokens=line.split()
    # get the 3 values on the input record.  (no error checking!!)
    tracl =long(tokens[0])
    fldr  =long(tokens[1])
    tracf =long(tokens[2]) 
    # set up header values with nulls.
    ep=-999999
    sx    =-999999999
    sy    =-999999999
    gx    =-999999999
    gy    =-999999999
    cdp   =-999999
    sutracf=-999999
    offset=-999999
    shotSpn=-999999999
    gSpn=-999999999
    sdepth=-999
    selev=-999
    gelev=-999
    sstat=-999
    gstat=-999
    tstat=-999
    # compute the shot number corresponding to this fldr (or recno in observer
    # log terminology
    shot=long(round(recnoSpn.linearInterp(fldr)))
    
    # If shot and tracf are valid, interpolate header values.
    # The recnoSpn.txt file sets shot(fldr) values to -999999.  This
    # causes the headers to be set to dummy values and later they will be
    # removed using suwind.
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
        roffset=(float(tracf)-((96.0+1.0)/2.0))*110.0
        offset=long(round(roffset))
        # shot x is shot number * 440 shot interval
        r_sx=float(ep)*440.0
        sx=long(round(r_sx))
        # shot and group y coordinates are all 0
        sy=0
        gx=long(round(r_sx+roffset))
        gy=0
        # cdp computed from shot/group midpoint 
        cdp=long(round((float(sx+gx)/2.0-55.0/2.0)/55.0-752+101))
        sutracf=long(round((gx/440.0+55.0/440.0)*4.0))
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

    # Code paths for bad and good traces merge again here.
    # Print the header values to standard output.
    print str(ep)              \
         + " " + str(sx)      \
         + " " + str(sy)      \
         + " " + str(gx)      \
         + " " + str(gy)      \
         + " " + str(cdp)     \
         + " " + str(sutracf) \
         + " " + str(offset)  \
         + " " + str(selev)   \
         + " " + str(gelev)   \
         + " " + str(sstat)   \
         + " " + str(gstat)   \
         + " " + str(tstat)


 
