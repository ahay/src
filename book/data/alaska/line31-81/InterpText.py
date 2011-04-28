#!/usr/bin/env python

import math
 
#for ep in range(100, 112):
#    os.system("echo ep=" + str(ep) )

class InterpText:
    def read(self,fileName):
        self.xin=[]
        self.yin=[]
        for self.line in open(fileName):
            tokens=self.line.split()
            self.xin.append(float(tokens[0]))
            self.yin.append(float(tokens[1]))
#        print self.xin
#        print self.yin
    def linearInterp(self,x):
        if(x<=self.xin[0]):
            return self.yin[0]
        #print "len(self.xin)=" + str(len(self.xin));
        for i in range(0,len(self.xin)-1):
            #print "test x<=self.xin[i+1]" + str(x) + "<=" + str(self.xin[i+1])
            if(x<=self.xin[i+1]):
                #print "x=" + str(x) + "i=" + str(i)
                return self.yin[i]+(x-self.xin[i])*(self.yin[i+1]-self.yin[i])/(self.xin[i+1]-self.xin[i])
        return self.yin[len(self.yin)-1]

spnElev = InterpText()
spnElev.read('spnElev.txt')

recnoSpn = InterpText()
recnoSpn.read('recnoSpn.txt')

for line in open('hdrfile.txt'):
    tokens=line.split()
    tracl =long(tokens[0])
    fldr  =long(tokens[1])
    tracf =long(tokens[2]) 
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
    shot=long(round(recnoSpn.linearInterp(fldr)))
#    str(shot)
    if(shot>0 and tracf<97):
        str(shot)
        ep=shot
        roffset=(float(tracf)-((96.0+1.0)/2.0))*110.0
        offset=long(round(roffset))
        r_sx=float(ep)*440.0
        sx=long(round(r_sx))
        sy=0
        gx=long(round(r_sx+roffset))
        gy=0
        cdp=long(round((float(sx+gx)/2.0-55.0/2.0)/55.0-752+101))
        sutracf=long(round((gx/440.0+55.0/440.0)*4.0))
        shotSpn=float(sx)/440.0
        gSpn=float(gx)/440.0
#             elevation at the shot
        selev=long(round(spnElev.linearInterp(shotSpn)))
#             elevation at the receiver
        gelev=long(round(spnElev.linearInterp(gSpn)))
        if(shotSpn<123):
            sdepth=long(round(82.5))
        else:
            sdepth=long(round(67.5))
        sstat=long(round((-selev+sdepth)*1000.0/10000.0))
        gstat=long(round(-gelev*1000.0/10000.0))
        tstat=long(round((-gelev-selev+sdepth)*1000.0/10000.0))

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



#         + " " +  str(shotSpn)
#         + " " +  str(gSpn   )
 
