#!/usr/bin/env python
import math
 
class InterpText:
    def __init__(self):
        self.xin=[]
        self.yin=[]
    def append(self,x,y):
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

class InterpText2:
    def __init__(self):
        self.xin=[]
        self.yin=[]        
        self.zin=[]
    def read(self,fileName):
        for self.line in open(fileName):
            tokens=self.line.split()
            self.xin.append(float(tokens[0]))
            self.yin.append(float(tokens[1]))
            self.zin.append(float(tokens[2]))
 #       print "xin="
 #       print self.xin
 #       print "yin="
 #       print self.yin
 #       print "zin="
 #       print self.zin

    def linearInterp(self,x):
        x=float(x)
        if x<=self.xin[0] :
            y=self.yin[0]
            z=self.zin[0]
        if x>=self.xin[-1] :
            y=self.yin[-1]
            z=self.zin[-1]
        if x>self.xin[0] and x<self.xin[-1] :
            #print "len(self.xin)=" + str(len(self.xin));
            for i in range(0,len(self.xin)-1):
                #print "test x<=self.xin[i+1]" + str(x) + "<=" + str(self.xin[i+1])
                if(x<=self.xin[i+1]):
                    #print "x=" + str(x) + "i=" + str(i)
                    y=self.yin[i]+ \
                        (x-self.xin[i])*(self.yin[i+1]-self.yin[i])/ \
                          (self.xin[i+1]-self.xin[i])
                    z=self.zin[i]+ \
                        (x-self.xin[i])*(self.zin[i+1]-self.zin[i])/ \
                          (self.xin[i+1]-self.xin[i])
                    break
        return (y,z)

spnElev = InterpText()
spnElev.read('spnElev.txt')
#for shot in range(99,223):
#    print "shot=" + str(shot) + " spnElev=" + str(spnElev.linearInterp(shot))

recnoSpn = InterpText()
recnoSpn.read('recnoSpn.txt')

spnCable = InterpText2()
spnCable.read('spnCable.txt')

spnOffset = InterpText2()
spnOffset.read('spnOffset.txt')

#for shot in range(99,223):
#    print "shot=" + str(shot) + " Cable=" + str(spnCable.linearInterp(shot))


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
    if(shot>0 and tracf<97):
        ep=shot
        # use spnCable.txt file to figure out group location
        (spngroup1,spngroup96)=spnCable.linearInterp(shot)
        cable=InterpText()
        cable.append(1,spngroup1)
        cable.append(96,spngroup96)
        r_gx=cable.linearInterp(tracf)*440.0
        r_gy=0.0
        # source x include the shift from spnOffset.txt
        (deltaEasting,deltaNorthing)=spnOffset.linearInterp(shot)
        # spns increase in the (0,-1) direction xline dir is (1,0)
        r_sx=-1.0*deltaNorthing+float(ep)*440.0
        r_sy= 1.0*deltaEasting + 0.0
        dx=r_sx-r_gx
        dy=r_sy-r_gy
        roffset=math.sqrt(dx*dx+dy*dy)
        offset=long(round(roffset))
        if dx>0 : offset=-offset
        sx=long(round(r_sx))
        sy=long(round(r_sy))
        gx=long(round(r_gx))
        gy=long(round(r_gy))
        cdp=long(round((float(r_sx+r_gx)/2.0-55.0/2.0)/55.0-752+101))
        sutracf=long(round((r_gx/440.0+55.0/440.0)*4.0))
        shotSpn=float(r_sx)/440.0
        gSpn=float(r_gx)/440.0
#             elevation at the shot
        selev=long(round(spnElev.linearInterp(shotSpn)))
#             elevation at the receiver
        gelev=long(round(spnElev.linearInterp(gSpn)))
        sdepth=long(round(82.5))
        sstat=long(round((-selev+sdepth)*1000.0/10000.0))
        gstat=long(round(-gelev*1000.0/10000.0))
        tstat=long(round((-gelev-selev+sdepth)*1000.0/10000.0))

    print str(ep)             \
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

 
