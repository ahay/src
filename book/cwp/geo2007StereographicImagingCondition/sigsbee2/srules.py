from rsf.proj import *
import fdmod,spmig,zomig

# ------------------------------------------------------------
def param():
    par = {
        'kt':25,    # wavelet delay
        'nx':800, 'ox':10.025, 'dx':0.025, 'lx':'position','ux':'kft',
        'nz':400,  'oz':6,     'dz':0.025, 'lz':'depth','uz':'kft',
        'nt':2000, 'ot':0,     'dt':0.008, 'lt':'time','ut':'s',
        'nw':200,  'ow':1,'jw':1,'verb':'y','nrmax':3,'dtmax':0.0001,
        'pmx':100,'tmx':100,
        #
        'oanga':0,   'nanga':60,  'danga':3.,    # SS angle 1
        'oangb':0,   'nangb':1,   'dangb':5.,    # SS angle 2
                  'nl':100, 'dl':1,     # SS line length
        'sig':1
        }

    par['xsou1']=16
    par['xsou2']=24

    par['ow']=2
    par['kw']=par['nt']/2+1    
    par['dw']=1/(par['nt']*par['dt'])
    par['fw']=int(par['ow']/par['dw'])

    par['tcut']=6
    return par

# ------------------------------------------------------------
def data():
    
    # velocity
    vstr = 'sigsbee2a_stratigraphy.sgy'
    Fetch(vstr,'sigsbee')
    Flow('zvstr tzvstr ./shead ./bshead',vstr,
         '''
         segyread
         tape=$SOURCE
         tfile=${TARGETS[1]}
         hfile=${TARGETS[2]}
         bfile=${TARGETS[3]}
         ''',stdin=0)

    # reflectivity
    refl = 'sigsbee2a_reflection_coefficients.sgy'
    Fetch(refl,'sigsbee')
    Flow('zrefl tzrefl',refl,
         '''
         segyread
         tape=$SOURCE
         tfile=${TARGETS[1]}
         hfile=/dev/null
         bfile=/dev/null
         ''',stdin=0)

# ------------------------------------------------------------
def run(par):
    # ------------------------------------------------------------
    fdmod.point('ss1',par['xsou1'],par['oz'],par)
    fdmod.point('ss2',par['xsou2'],par['oz'],par)
    Plot('ss1','window      |' + fdmod.ssplot('plotcol=5',par))
    Plot('ss2','window      |' + fdmod.ssplot('plotcol=5',par))

    # ------------------------------------------------------------
    # velocity
    Plot(  'vel',fdmod.cgrey('bias=4.8 allpos=y pclip=99 color=j',par))
    Result('vel','vel ss1 ss2','Overlay')

    # slowness
    Flow('slo','vel',
         '''
         math output=1/input |
         transp |
         transp plane=23 |
         put o2=0 d2=1 label2=y
         ''')
    Result('slo','window | transp |'
           + fdmod.cgrey('allpos=y pclip=95 bias=0.125',par))

    # reflectivity
    Flow('ref','del',
         '''
         transp |
         transp plane=23 |
         put o2=0 d2=1 label2=y
         ''')
    Result('ref','window | transp |' + fdmod.cgrey('pclip=99',par))
    
    # ------------------------------------------------------------
    # wavelet
    fdmod.wavelet('wav',8,par)
    Result('wav','window n1=500 |' + fdmod.waveplot('',par))

    # ------------------------------------------------------------
    for i in ('1','2'):
        Flow('spk'+i,'wav',
             '''
             pad beg2=%d n2out=%d |
             put o2=%g d2=%g
             ''' % ( (par['xsou'+i]-par['ox'])/par['dx'],par['nx'],par['ox'],par['dx']) )
        Result('spk'+i,fdmod.dgrey('pclip=100',par))
        
    # ------------------------------------------------------------
    for i in ('1','2'):
        # source wavefield (from time to frequency)
        zomig.wflds('dds'+i,'spk'+i,par)
        
        # wavefield extrapolation MODELING
        spmig.modelPW3('ddr'+i,'slo','dds'+i,'ref',par)
        
    # ------------------------------------------------------------
    # shots 1 and 2
    Flow('dds0','dds1 dds2','add ${SOURCES[1]}')
    Flow('ddr0','ddr1 ddr2','add ${SOURCES[1]}') # both shots

    for i in ('1','2','0'):
        Result('dds'+i,'window | real | smooth rect1=5 | sfgrey pclip=100')
        Result('ddr'+i,'window | real | smooth rect1=5 | sfgrey pclip=100')

    # ------------------------------------------------------------
    # available dds[0,1,2], ddr[0,1,2], slo
    # ------------------------------------------------------------
    for i in ('1','2','0'):
        # recorded data (from frequency to time) 
        Flow('ttr'+i,'ddr'+i,
             '''
             window |
             transp |
             pad beg1=%(fw)d n1out=%(kw)d |
             fft1 inv=y opt=n
             ''' % par)
        Result('ttr'+i,fdmod.dgrey('pclip=100 min1=1 max1=6 screenratio=0.5 screenht=7',par))
    
        # wavefield extrapolation MIGRATION
        spmig.imagePW3('ii'+i,'cc'+i,'slo','dds'+i,'ddr'+i,par)
        Plot(  'ii'+i,'window min1=%(lox)d max1=%(hix)d | transp |' %par
               + fdmod.cgrey('pclip=99',par))
        Result('ii'+i,['ii'+i,'ss1','ss2'],'Overlay')

    Flow('ii','ii1 ii2','add ${SOURCES[1]}')
    Plot('ii','window min1=%(lox)d max1=%(hix)d | transp |' %par
         + fdmod.cgrey('pclip=100',par))
    Result('ii',['ii','ss1','ss2'],'Overlay')
    
    # ------------------------------------------------------------
    # datuming
    zomig.Cwfone3('wfs','dds0','slo',par) # source wavefield for one shot
    zomig.Awfone3('wfr','ddr0','slo',par) # receiver wavefield for two shots

    # ------------------------------------------------------------
    # data in the time domain
    for k in ('s','r'):    
        Flow('q'+k,'wf'+k,
             '''
             window min1=%(lox)d max1=%(hix)d j1=2 j3=2 |
             transp plane=23 memsize=500 |
             transp plane=12 memsize=500 |
             pad beg1=%(fw)d n1out=%(kw)d |
             fft1 inv=y opt=n |
             window max1=%(tcut)g
             ''' % par)
        
        Plot('q'+k,'window j3=10 |' +
             fdmod.dgrey('gainpanel=a pclip=99',par),view=1)

    Flow(  'qi','qs qr','add mode=p ${SOURCES[1]}')
    Plot('qi','window j3=10 |'
           + fdmod.dgrey('gainpanel=a pclip=100',par),view=1)
        
    # SIC
    Flow('kk',['qs','qr'],
         '''
         sic3d ur=${SOURCES[1]} nbuf=500 verb=y stack=n
         oanga=%(oanga)g nanga=%(nanga)d danga=%(danga)g
         oangb=%(oangb)g nangb=%(nangb)d dangb=%(dangb)g
         nl=%(nl)d dl=%(dl)g
         sig=%(sig)g
         ''' % par)
    
    Plot('kk','transp plane=23 | stack | transp |'
           + fdmod.cgrey('pclip=100',par))
    Result('kk',['kk','ss1','ss2'],'Overlay')
