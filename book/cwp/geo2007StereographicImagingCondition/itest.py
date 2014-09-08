from rsf.proj import *
import fdmod
from math import *

# ------------------------------------------------------------
def param():
    par = {
        'nt':1000, 'dt':0.0005,'ot':0, 'lt':'time', 'ut':'s',     
        'nx':1000, 'dx':2.0,   'ox':0, 'lx':'position', 'ux':'m',
        'nz':200,  'dz':2.0,   'oz':0, 'lz':'depth', 'uz':'m',
        'kt':50,    # wavelet delay
        'frq':50,
        'jsnap':4,
        'vprt':0,
        'test':0,
        'nb':100,
        'oanga':0,   'nanga':60, 'danga':3.,  # SS angle 1
        'oangb':0,   'nangb':2,  'dangb':1,   # SS angle 2
                     'nl':100,   'dl':1.0,    # SS line length
        'sig':1                               # SS gaussian decay
        }

    par['kx']=par['nx']/2-75
    par['mx']=par['nx']/2-50
    par['kz']=par['nz']/2-25
    par['mz']=par['nz']/2+25
        
    # add F-D modeling parameters
    fdmod.param(par)
    return par

# ------------------------------------------------------------
# sources
def onesou(ss,par):
    fdmod.point(ss+'1',par['ox']+5*par['nx']/10*par['dx'],0,par)
    Flow(ss,ss+'1','window')

def twosou(ss,par):
    fdmod.point(ss+'0',par['ox']+3.5*par['nx']/10*par['dx'],0,par)
    fdmod.point(ss+'1',par['ox']+5.0*par['nx']/10*par['dx'],0,par)    
    Flow(ss,[ss+'0',ss+'1'],'cat axis=2 space=n ${SOURCES[0:2]}', stdin=0)

def thrsou(ss,par):
    fdmod.point(ss+'0',par['ox']+3*par['nx']/10*par['dx'],0,par)
    fdmod.point(ss+'1',par['ox']+5*par['nx']/10*par['dx'],0,par)
    fdmod.point(ss+'2',par['ox']+6*par['nx']/10*par['dx'],0,par)
    Flow(ss,[ss+'0',ss+'1',ss+'2'],'cat axis=2 space=n ${SOURCES[0:3]}', stdin=0)

# ------------------------------------------------------------
# reflectors
def flat(refl,par):
    Flow(refl,None,
         '''
         spike nsp=1 mag=1
         n1=%(nz)d o1=%(oz)g d1=%(dz)g k1=151
         n2=%(nx)d o2=%(ox)g d2=%(dx)g
         ''' % par)

# segment
def segment(refl,par):
    Flow(refl,None,
         '''
         spike nsp=1 mag=1
         n1=%(nz)d o1=%(oz)g d1=%(dz)g k1=151
         n2=%(nx)d o2=%(ox)g d2=%(dx)g k2=475 l2=525 |
         smooth rect2=20
         ''' % par)

def rand(refl,par):
    Flow('rand',None,
         '''
         spike nsp=1 mag=1 n1=%d k1=1 |
         noise |
         window n2=1 |
         pad n2out=%d |         
         clip clip=0.25 |
         transp |
         put n1=%d n2=1
         ''' %(par['nz']/20,20,par['nz']))
    Flow('refl','rand',
         '''
         window f1=100 |
         pad beg1=100 |
         spray axis=2 n=%(nx)d o=%(ox)g  d=%(dx)g |
         put                   o1=%(oz)d d1=%(dz)g |
         scale axis=123
         ''' %par)


#    Flow(refl,None,
#         '''
#         spike nsp=1 mag=0
#         n1=%(nz)d o1=%(oz)g d1=%(dz)g |
#         noise |
#         scale axis=123 |
#         clip clip=0.25 |
#         scale axis=123 |
#         window f1=100 |
#         pad beg1=100 |
#         spray axis=2 n=%(nx)d o=%(ox)g d=%(dx)g
#         ''' % par)

def wave(refl,par):
    k=0.2
    l=par['xmax']
    d=300
    layers = (
        ((0.000*l,d),
         (0.250*l,d*(1-k)),
         (0.500*l,d*(1+k)),
         (0.750*l,d*(1-k)),
         (1.000*l,d)),
        ((0,1.5),(l,1.5))
        )
    
    vi=0
    vt=1.0
    vels = "%s,%s,%s" % (vi,vt,vi)
    drvs = "%s,%s" % (0,0)
    
    dim1 = 'd1=%(dx)g o1=%(ox)g n1=%(nx)d' % par
    dim2 = 'd2=%(dz)g o2=%(oz)g n2=%(nz)d' % par
    
    for i in range(len(layers)):
        
        inp = 'inp%d' % (i+1)
        Flow('./'+inp+'.rsf',None,
             'echo %s in=$TARGET data_format=ascii_float n1=2 n2=%d' % \
             (string.join(map(lambda x: string.join(map(str,x)),layers[i]),
                          ' '),
              len(layers[i])))
        
    Flow(refl+'lay1','inp1','dd form=native | spline %s fp=%s' % (dim1,drvs))
    Flow(refl+'lay2',refl+'lay1','add add=%(dz)f' % par)
    Flow(refl+'layers',[refl+'lay1',refl+'lay2'],'cat axis=2 ${SOURCES[1:2]}')
    Flow(refl,refl+'layers',
         '''
         unif2 v00=%s n1=%d d1=%g o1=%g
         ''' % (vels,par['nz'],par['dz'],par['oz']) )

def dips(refl,ang,par):

    l=par['xmax']
    d = 300
    f = ang;
    k = (l/(2*d)) * sin(pi*f/180.);
    
    layers = (
        ((0.000*l,d*(1+k)),
         (1.000*l,d*(1-k))),
        ((0,1.5),(l,1.5))
        )

    vi=0
    vt=1.0
    vels = "%s,%s,%s" % (vi,vt,vi)

    t=-2*d*k/l;
    drvs = "%s,%s" % (t,t)
    
    dim1 = 'd1=%(dx)g o1=%(ox)g n1=%(nx)d' % par
    dim2 = 'd2=%(dz)g o2=%(oz)g n2=%(nz)d' % par
    
    for i in range(len(layers)):
        
        inp = 'inp%d' % (i+1)
        Flow('./'+inp+'.rsf',None,
             'echo %s in=$TARGET data_format=ascii_float n1=2 n2=%d' % \
             (string.join(map(lambda x: string.join(map(str,x)),layers[i]),
                          ' '),
              len(layers[i])))
        
    Flow(refl+'lay1','inp1','dd form=native | spline %s fp=%s' % (dim1,drvs))
    Flow(refl+'lay2',refl+'lay1','add add=%(dz)f' % par)
    Flow(refl+'layers',[refl+'lay1',refl+'lay2'],'cat axis=2 ${SOURCES[1:2]}')
    Flow(refl,refl+'layers',
         '''
         unif2 v00=%s n1=%d d1=%g o1=%g
         ''' % (vels,par['nz'],par['dz'],par['oz']) )

# ------------------------------------------------------------
def run(par):
    # experiments
    fdmod.horizontal('rr',0,par)

    Plot('rr','window j2=10|' + fdmod.rrplot('',par))
    Plot('ss','window      |' + fdmod.ssplot('',par))
    Plot('sx','window      |' + fdmod.ssplot('plotcol=5',par))
    
    # wavelet
    fdmod.wavelet('wav_',par['frq'],par)
    Flow('wav','wav_','transp')
    Result('wav','window n2=200 |' + fdmod.waveplot('',par))
    
    # velocity
    Flow('vbck',None,
         '''
         math n1=%(nz)d o1=%(oz)g d1=%(dz)g output="2000" |
         spray axis=2 n=%(nx)d o=%(ox)g d=%(dx)g
         ''' % par)
    Flow('vprt',None,
         '''
         spike nsp=1 mag=1
         n1=%(nz)d o1=%(oz)g d1=%(dz)g k1=%(kz)d l1=%(mz)d
         n2=%(nx)d o2=%(ox)g d2=%(dx)g k2=%(kx)d l2=%(mx)d |
         smooth rect1=25 rect2=25 repeat=3 |
         scale axis=123 |
         scale rscale=%(vprt)g
         ''' % par)
    Flow(  'velo','vbck vprt','add ${SOURCES[1]}')
    Plot(  'velo',fdmod.cgrey('allpos=y bias=1200 pclip=100 color=g',par))
    Result('velo',['velo','ss','sx'],'Overlay')
    
    # density
    Flow('dens','velo','math output=1')
    
    # reflector
    Plot('refl','refl velo',
         '''
         depth2time velocity=${SOURCES[1]} dt=%(dt)g nt=%(nt)d |
         scale rscale=-1 |
         ricker1 frequency=%(frq)g |
         time2depth velocity=${SOURCES[1]} dz=%(dz)g nz=%(nz)d |
         ''' % par + fdmod.cgrey('pclip=100',par))
    Result('refl',['refl','ss','sx'],'Overlay')

    # mask
    Flow('mask',None,
         '''
         spike nsp=1 mag=1
         n1=%(nx)d o1=%(ox)g d1=%(dx)g k1=101 l1=900
         n2=%(nt)d o2=%(ot)g d2=%(dt)g |
         smooth rect1=100 |
         scale axis=123
         ''' % par)
    Result('mask','transp |' + fdmod.dgrey('allpos=y pclip=100',par))
        
    # F-D modeling (born)
    fdmod.lwefd1(
        'do','wo',
        'dd','wd',
        'wav','velo','dens','refl','ss','rr','jsnap=100',par)
    Result('do','transp | window min1=0.25 |' + fdmod.dgrey('min1=0.25 pclip=100',par))
    Result('dd','transp | window min1=0.25 |' + fdmod.dgrey('min1=0.25 pclip=100',par))
    Plot('wo',fdmod.wgrey('',par),view=1)
    Plot('wd',fdmod.wgrey('',par),view=1)

    # source data and wavefield
    fdmod.awefd1(
        'ds','ws',
        'wav','velo','dens','sx','rr','',par)
    Plot('ws','window j3=20 |' + fdmod.wgrey('',par),view=1)

    # receiver wavefield
    Flow('du','dd mask',
         'add mode=p ${SOURCES[1]} | reverse which=2 opt=i verb=y')
    fdmod.awefd(
        'dx','wx',
        'du','velo','dens','rr','rr','',par)
    Flow('dr','dx','reverse which=2 opt=i verb=y')
    Flow('wr','wx','reverse which=4 opt=i verb=y')
    Plot('wr','window j3=20 |' + fdmod.wgrey('',par),view=1)

    for i in range(0,par['nt']/100,1):
        fdmod.wframe('wo'+'-'+str(i),'wo',i,'pclip=99.9',par)
        fdmod.wframe('wd'+'-'+str(i),'wd',i,'pclip=100',par)

    for i in range(0,par['nt']/100,1):
        fdmod.wframe('wx'+'-'+str(i),'wx',i*25,'pclip=99.9',par)
    
    # ------------------------------------------------------------
    minx=500
    maxx=1500
    minz=par['oz']+par['nz']*par['dz']/2
    numz=par['nz']/2

    mint=0.1
    numt=150
    maxt=mint+numt*par['dt']*par['jsnap']

    # wavefield
    for i in ('s','r'):
        Flow('u'+i,'w'+i,
             '''
             window min1=%(zmin)g max1=%(zmax)g min2=%(xmin)g max2=%(xmax)g |
             scale axis=123
             ''' % par)
        Plot('u'+i,'window j3=10 |' + fdmod.wgrey('pclip=99',par),view=1)

        for k in range(0,par['nt']/par['jsnap'],25):
            fdmod.wframe('u'+i+'-'+str(k/25),'u'+i,k,'pclip=99',par)
        
        # windowed wavefields
        Flow('p'+i,'u'+i,
             '''
             window min1=%g n1=%g min2=%g max2=%g min3=%g n3=%g
             ''' % (minz,numz,
                    minx,maxx,
                    mint,numt))
        Flow('q'+i,'p'+i,'transp plane=13 memsize=500')
        Flow('o'+i,'q'+i,'transp plane=23 memsize=500')
        
    Flow(  'qi','qs qr','add mode=p ${SOURCES[1]}')
    Flow(  'oi','os or','add mode=p ${SOURCES[1]}')
    
    for i in ('s','r','i'):
        Plot('q'+i,'window j3=10 |'
             + fdmod.dgrey('gainpanel=a pclip=100',par),view=1)    
        Plot('o'+i,'window j3=10 | transp |'
             + fdmod.egrey('gainpanel=a pclip=100',par),view=1)
        
        Flow(['q'+i+'plt','q'+i+'bar'],'q'+i,
             'byte bar=${TARGETS[1]} gainpanel=a pclip=100')
        for k in range(10):
            Result('q'+i+'plt'+str(k),['q'+i+'plt','q'+i+'bar'],
                   'window n3=1 f3=%d |' % (10*k) +
                   fdmod.dgrey(
                '''
                bar=${SOURCES[1]} min1=%g max1=%g min2=%g max2=%g
                labelsz=8 labelfat=3 screenratio=1.5
                ''' %(mint,maxt,minx,maxx),par))
        
    # cut along the reflectors
    Flow('cut','refl',
         '''
         window min1=%g n1=%g min2=%g max2=%g |
         spray axis=3 n=%d o=%g d=%g |
         transp plane=13 memsize=500
         ''' % (minz,numz,
                minx,maxx,
                numt,0.1,par['dt']*4) )
    
    for i in ('s','r'):
        Flow('c'+i,['q'+i,'cut'],
             '''
             add mode=p ${SOURCES[1]} |
             transp plane=23 |
             stack
             ''')
        Result('c'+i,fdmod.dgrey('',par))

        Flow(  'f'+i,'q'+i,'window n3=1 min3=300')
        Result('f'+i,fdmod.dgrey('',par))

    # ------------------------------------------------------------
    # conventional IC
    Flow(  'ii',['ps','pr'],'ic ur=${SOURCES[1]} version=0 nbuf=500 verb=y')
    Plot(  'ii',fdmod.cgrey('pclip=99.9',par))
    Result('ii',['ii','ss','sx'],'Overlay')

# ------------------------------------------------------------
def stereo2d(par):
    # 2D stereographic IC
    Flow('jj',['qs','qr'],
         '''
         sic ur=${SOURCES[1]} nbuf=500 verb=y
         oa=%(oa)g na=%(na)d da=%(da)g
         nl=%(nl)d dl=%(dl)g
         sig=%(sig)g
         ''' % par)
    Plot(  'jj','transp | ' + fdmod.cgrey('pclip=99.9',par))
    Result('jj',['jj','jj','ss','sx'],'Overlay')
    
# ------------------------------------------------------------
def stereo3d(par):
    # 3D stereographic IC
    Flow('kk',['qs','qr'],
         '''
         sic3d ur=${SOURCES[1]} nbuf=500 verb=y stack=n
         oanga=%(oanga)g nanga=%(nanga)d danga=%(danga)g
         oangb=%(oangb)g nangb=%(nangb)d dangb=%(dangb)g
         nl=%(nl)d dl=%(dl)g
         sig=%(sig)g
         ''' % par )    
    Plot(  'kk','transp plane=23 | stack | transp | '
           + fdmod.cgrey('pclip=99.9',par))
    Result('kk',['kk','ss','sx'],'Overlay')
    
# ------------------------------------------------------------
def test(par):

    for i in ('s','r'):
        Flow(  'v'+i,'u'+i,'window n3=1 f3=400')
        Result('v'+i,fdmod.cgrey('',par))
        
        Flow(  'k'+i,'f'+i,
               '''
               lstk verb=y
               oa=%(oa)d na=%(na)d da=%(da)g
               nl=%(nl)d dl=%(dl)g sig=%g
               ''' % par )
        Result('k'+i,fdmod.dgrey('gainpanel=a pclip=100',par))
        
        jx = 20
        Result('k'+i+'all','k'+i,
               '''
               window j2=%d |
               transp plane=23 |
               put n3=1 n2=%d |
               grey pclip=100
               ''' % (jx,par['nx']/2/jx*na) )
        
        Flow('t'+i,'k'+i,'byte gainpanel=a pclip=100 | put label1=x label2=t label3=a')
        
        Result('ff','fs fr','cat axis=3 space=n ${SOURCES[1]} | transp | grey pclip=100 gainpanel=e')
        Result('tt','ts tr',
               '''
               cat axis=3 space=n ${SOURCES[1]} |
               window j1=10 |
               transp plane=12 |
               transp plane=23 |
               grey title=""
               ''')
        
#        for i in ('s','r'):
#            Flow(['t'+i+'plt','t'+i+'bar'],'k'+i,'byte bar=${TARGETS[1]} gainpanel=a pclip=99.9')
#            for k in range(20):
#                Result('t'+i+'plt'+str(k),['t'+i+'plt','t'+i+'bar'],
#                       'window n1=1 f1=%d | transp |' % (20*k) +
#                       'grey title="" screenratio=0.2 screenht=3 label2=x unit2=m label1=p unit1="#" ')
                
        Flow('fi','fs fr','add mode=p ${SOURCES[1]}')
        Flow('ki','ks kr',
             '''
             add mode=p ${SOURCES[1]} |
             transp plane=23 memsize=1000 |
             stack
             ''')
        
        Result('fi',fdmod.dgrey('',par))
        Result('ki',fdmod.dgrey('',par))
        
# ------------------------------------------------------------
def debug(par):
    proj = Project()
    proj.Program('SICt.x',Split('SICt.c intutil.c'))

    Flow('k1',['qs','qr'],
         '''
         sic ur=${SOURCES[1]} nbuf=500 verb=y
         oa=%(oa)g na=%(na)d da=%(da)g
         nl=%(nl)d dl=%(dl)g
         sig=%(sig)g
         ''' % par)
    
    Flow('k2',['ps','pr'],
         '''
         sic3d
         ur=${SOURCES[1]} verb=y
         oanga=%(oanga)g nanga=%(nanga)d danga=%(danga)g
         oangb=%(oangb)g nangb=%(nangb)d dangb=%(dangb)g
         nl=%(nl)d dl=%(dl)g
         sig=%(sig)g
         ''' % par )
    Flow('k3',['qs','qr','SICt.x'],
         '''
         ./SICt.x
         ur=${SOURCES[1]} verb=y
         oanga=%(oanga)g nanga=%(nanga)d danga=%(danga)g
         oangb=%(oangb)g nangb=%(nangb)d dangb=%(dangb)g
         nl=%(nl)d dl=%(dl)g
         sig=%(sig)g
         ''' % par )
    Flow('k4',['qs','qr','SICt.x'],
         '''
         ./SICt.x
         ur=${SOURCES[1]} verb=y
         oa=%(oa)g na=%(na)d da=%(da)g
         nl=%(nl)d dl=%(dl)g
         sig=%(sig)g
         ''' % par )
    
    Plot(  'k1','transp |' + fdmod.cgrey('pclip=100',par))
    Plot(  'k2',             fdmod.cgrey('pclip=100',par))
    Plot(  'k3','transp |' + fdmod.cgrey('pclip=100',par))
    Plot(  'k4','transp |' + fdmod.cgrey('pclip=100',par))
    
    for id in ('4','3'):
        Result('k'+id,['k'+id,'ss','sx'],'Overlay')

def seislet(uu,us,ur,n1,n2,par):

    rect1=10
    rect2=10
    p0=0
    pmin=-100
    eps=0.1

    for i in (us,ur):
        Flow(i+'dat',i,'pad n2out=%d ' % n2)
        Flow(i+'dip',i+'dat',
             'dip rect1=%d rect2=%d p0=%g pmin=%g' % (rect1,rect2,p0,pmin))
        Result(i+'dip','grey color=j scalebar=y')
        
        Flow(i+'sei',[i+'dat',i+'dip'],
             'seislet dip=${SOURCES[1]} eps=%g adj=y inv=y' % eps)
        Flow(i+'inv',[i+'sei',i+'dip'],
             'seislet dip=${SOURCES[1]} eps=%g' % eps)
        Result(i+'inv','grey')

    Flow(uu+'sei',[us+'sei',ur+'sei'],'add mode=p ${SOURCES[1]}')
    Flow(uu+'dat',[uu+'sei',us+'dip'],'seislet dip=${SOURCES[1]} eps=%g' % eps)

    for i in (us,ur,uu):
        Result(i+'sei','put o2=0 d2=1 | grey unit2=')
        Result(i+'dat','grey')

def dwt(uu,us,ur,par):

    for i in (us,ur):
        Flow(i+'dwt',i,'transp | dwt | transp')
#    Flow(uu+'dwt',[us+'dwt',ur+'dwt'],'add mode=p ${SOURCES[1]}')

    Flow(us+'msk',us+'dwt','envelope | scale axis=123')
    Flow(uu+'dwt',[us+'msk',ur+'dwt'],'add mode=p ${SOURCES[1]}')


    for i in (us,ur,uu):
        Result(i+'dwt',
               'put o2=0 d2=1 | grey  title="DWT" label2=Scale unit2=')

        Flow(i+'inv',i+'dwt','transp | dwt adj=y inv=y | transp')
        Result(i+'inv','grey title="INV"')

    Flow(uu,uu+'inv','window')


def new(par):

    fdmod.point3('ssold0',par['ox']+3*par['nx']/10*par['dx'],par['oz'],1,par)
    fdmod.point3('ssold1',par['ox']+5*par['nx']/10*par['dx'],par['oz'],1,par)
    fdmod.point3('ssold2',par['ox']+6*par['nx']/10*par['dx'],par['oz'],1,par)
    Flow('ssold',['ssold0','ssold1','ssold2'],'cat axis=2 space=n ${SOURCES[0:3]}', stdin=0)    
    Flow('wavold','wav','window')
    fdmod.lmodel(
        'doold','woold',
        'ddold','wdold',
        'wavold','velo','refl',
        'ssold','rr','jsnap=100 nbz=100 nbx=100 ',par)
    
    Result('woold',fdmod.wgrey('',par))
    Result('wdold',fdmod.wgrey('',par))

    fdmod.point('ssnew0',par['ox']+3*par['nx']/10*par['dx'],par['oz'],par)
    fdmod.point('ssnew1',par['ox']+5*par['nx']/10*par['dx'],par['oz'],par)
    fdmod.point('ssnew2',par['ox']+6*par['nx']/10*par['dx'],par['oz'],par)
    Flow('ssnew',['ssnew0','ssnew1','ssnew2'],'cat axis=2 space=n ${SOURCES[0:3]}', stdin=0)    
    Flow('wavnew','wav','window squeeze=n')
    fdmod.lwefd1(
        'donew','wonew',
        'ddnew','wdnew',
        'wavnew','velo','dens','refl',
        'ssnew','rr','jsnap=100 nb=100 ',par)
    
    Result('wonew',fdmod.wgrey('',par))
    Result('wdnew',fdmod.wgrey('',par))

