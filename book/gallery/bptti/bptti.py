from rsf.proj import *
import sys

# Download from http://www.freeusp.org/2007_BP_Ani_Vel_Benchmark/
tgz = 'ModelParams.tar.gz'

Fetch(tgz,'BP',top=os.environ.get('DATAPATH'),server='local')

pars = Split('epsilon delta vp theta')

sgy = {}
for par in pars:
    sgy[par] = os.path.join('ModelParams',par.capitalize() + '_Model.sgy')

zcat = WhereIs('gzcat') or WhereIs('zcat')

Flow(sgy.values(),tgz,
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
        print 'Unknown parameter', par
        sys.exit(0)

# Download from http://www.freeusp.org/2007_BP_Ani_Vel_Benchmark/Anisotropic_FD_Model_Shots_part1.sgy.gz
