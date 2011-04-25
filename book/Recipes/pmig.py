try:
    from rsf.cluster import *
except:
    from rsf.proj import *
import string, sys
import version

def pmig(data,hdip,ydip,frame1=1,frame2=1,frame3=1):

    if version.old_version():
        return # think how to do it better

    def Grey3(name,title,color='I',extra=''):
        Flow('bar'+name,name,'bar')
        Result(name,[name,'bar'+name],
               '''
               byte | transp plane=23 | 
               grey3 title="%s" bar=${SOURCES[1]}
               frame1=%d frame2=%d frame3=%d color=%s 
               flat=y point1=0.75 point2=0.75 %s
               ''' % (title,frame1,frame2,frame3,color,extra))

    Grey3(data,'Data')
    Grey3(hdip,'Offset Slope',  'j','scalebar=y bartype=v')
    Grey3(ydip,'Midpoint Slope','j','scalebar=y bartype=v')

    nmo = data+'-nmo'
    vnmo = data+'-vnmo'
    Flow([nmo,vnmo],[data,hdip],'pnmo dip=${SOURCES[1]} vel=${TARGETS[1]}')
    Grey3(nmo,'Oriented NMO')

    stk = data+'-stk'
    Flow(stk,nmo,'stack')
    Result(stk,'grey title="Oriented Stack" wheretitle=t wherexlabel=b')

    for case in (data,hdip,ydip):
        Flow(case+'2',case,'transp plane=23 memsize=500')

    mig = data+'-mig'
    Flow(mig,[data+'2',hdip+'2',ydip+'2'],
         '''
         pmig hdip=${SOURCES[1]} xdip=${SOURCES[2]} |
         transp plane=23 memsize=500
         ''')
    Grey3(mig,'Oriented Prestack Time Migration')
    Result(mig+'2',mig,
           '''
           stack | 
           grey title="Oriented Prestack Time Migration" 
           wheretitle=t wherexlabel=b
           ''')
    
    mzo = data+'-mzo'
    Flow(mzo,[data+'2',hdip+'2',ydip+'2'],
         '''
         pmig hdip=${SOURCES[1]} xdip=${SOURCES[2]} mzo=y |
         transp plane=23 memsize=500
         ''')
    Grey3(mzo,'Oriented Migration to Zero Offset')
    Result(mzo+'2',mzo,
           'stack | grey title="Oriented DMO Stack" wheretitle=t wherexlabel=b')
    
