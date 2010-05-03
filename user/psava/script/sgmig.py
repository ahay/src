from rsf.proj import *

def wflds(wfld,cmps,par):
    Flow(wfld,cmps,
         '''
         fft1 inv=n |
         window squeeze=n n1=%(nw)d min1=%(ow)g |
         transp plane=13 memsize=500 |
         spray axis=2 n=1 o=0 d=1 |
         put label1=mx label2=my label3=hx label4=w
         ''' % par )
    
def datum(wfl1,slow,wfl0,par):
    Flow(wfl1,[wfl0,slow],
         '''
         %(CAM)s mode=d
         slo=${SOURCES[1]}
         ''' % par)

def image(imag,slow,data,par):
    Flow(imag,[data,slow],
         '''
         %(CAM)s mode=m
         slo=${SOURCES[1]}
         ''' % par)

def cimage(imag,slow,data,par):
    Flow(imag,[data,slow],
         '''
         %(CAM)s mode=m
         slo=${SOURCES[1]}
         <   ${SOURCES[0]}
         ''' % par, stdin=0)
