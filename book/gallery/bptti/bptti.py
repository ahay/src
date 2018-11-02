from rsf.proj import *
import sys

# Fetch will not work unless you do the following:
# 1. Download ModelParams.tar.gz from http://www.freeusp.org/2007_BP_Ani_Vel_Benchmark/
# 2. Put in under $DATAPATH/BP
tgz = 'ModelParams.tar.gz'

Fetch(tgz,'BP',top=os.environ.get('DATAPATH'),server='local')

pars = Split('epsilon delta vp theta')

sgy = {}
for par in pars:
    sgy[par] = os.path.join('ModelParams',par.capitalize() + '_Model.sgy')

zcat = WhereIs('gzcat') or WhereIs('zcat')

Flow(list(sgy.values()),tgz,
     zcat + ' $SOURCE | tar -xvf -',stdin=0,stdout=-1)

def getmod(par):
    if par in pars:
        Flow([par,par+'.asc',par+'.bin'],sgy[par],
             '''
             segyread hfile=${TARGETS[1]} bfile=${TARGETS[2]} read=d |
             put
             o2=0 d2=0.00625 label2=Distance unit2=km
             o1=0 d1=0.00625 label1=Depth unit1=km %s |
             transp plane=12 
             ''' % ('','| scale dscale=0.001')[par=='vp'])
    elif par=='vx':
        Flow(par,'vp epsilon',
             '''
             math e=${SOURCES[1]} output="input*sqrt(1+2*e)"
             ''')
    elif par=='eta':
        Flow(par,'epsilon delta',
             '''
             math e=${SOURCES[0]} d=${SOURCES[1]} output="(e-d)/(1+2*d)"
             ''')
    else:
        print('Unknown parameter', par)
        sys.exit(0)

# Download from http://www.freeusp.org/2007_BP_Ani_Vel_Benchmark/Anisotropic_FD_Model_Shots_part1.sgy.gz

shots =[]
tshots = []
for s in range(4):
    segy = 'Anisotropic_FD_Model_Shots_part%d.sgy' % (s+1)
    segyz = segy+'.gz'

    Fetch(segyz,'BP',top=os.environ.get('DATAPATH'),server='local')
    Flow(segy,segyz,zcat + ' $SOURCE',stdin=0)

    shot = 'shot%d' % s
    Flow([shot,'t'+shot,shot+'.asc'],segy,
         'segyread tfile=${TARGETS[1]} hfile=${TARGETS[2]}')

    shots.append(shot)
    tshots.append('t'+shot)

def getshots(myshots,mytshots):
    Flow(myshots,shots,'cat axis=2 ${SOURCES[1:4]}')
    Flow(mytshots,tshots,'cat axis=2 ${SOURCES[1:4]}')
