from rsfproj import *
import string, sys

def pmig(data,hdip,ydip,frame1=1,frame2=1,frame3=1):

    version = string.replace(sys.version,"+","")
    version = string.split(string.split(version)[0], ".")
    if map(int, version) < [2, 2, 0]:
        return # think how to do it better

    def grey3(title,color='I'):
        return '''byte bar=bar.rsf | transp plane=23 | grey3 title="%s"
        frame1=%d frame2=%d frame3=%d color=%s flat=y point1=0.75 point2=0.75
        ''' % (title,frame1,frame2,frame3,color)

    Result(data,grey3('Data'))
    Result(hdip,grey3('Offset Slope','j')   + ' scalebar=y bartype=v')
    Result(ydip,grey3('Midpoint Slope','j') + ' scalebar=y bartype=v')

    nmo = data+'-nmo'
    vnmo = data+'-vnmo'
    Flow([nmo,vnmo],[data,hdip],'pnmo dip=${SOURCES[1]} vel=${TARGETS[1]}')
    Result(nmo,grey3('Oriented NMO'))

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
    Result(mig,grey3('Oriented Prestack Time Migration'))
    Result(mig+'2',mig,
           'stack | grey title="Oriented Prestack Time Migration" wheretitle=t wherexlabel=b')
    
    mzo = data+'-mzo'
    Flow(mzo,[data+'2',hdip+'2',ydip+'2'],
         '''
         pmig hdip=${SOURCES[1]} xdip=${SOURCES[2]} mzo=y |
         transp plane=23 memsize=500
         ''')
    Result(mzo,grey3('Oriented Migration to Zero Offset'))
    Result(mzo+'2',mzo,
           'stack | grey title="Oriented DMO Stack" wheretitle=t wherexlabel=b')
    
