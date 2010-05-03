from rsf.proj import *
import zomig,spmig,fdmod

def param(par):
    par['verb']='y'
    par['eps']=0.1
    par['nrmax']=3
    par['dtmax']=0.00005
    par['tmx']=16

    par['lx']='Position'
    par['ly']=''
    par['lz']='Depth'
    par['lt']='Time'
    
    # ------------------------------------------------------------
    par['xpad']=par['nx']/2.
    par['xsou']=par['ox'] + par['xpad'] * par['dx']
    par['xoff']=          - par['xpad'] * par['dx']
    par['zsou']=par['oz']
    
    par['dw']= 1/(par['nt']*par['dt'])
    par['begw']=par['ow']/par['dw']
    par['padw']=par['nt']/2+1

    par['kt']=50
    
def test(par):
    # ------------------------------------------------------------
    # source coordinate
    fdmod.point('ss',par['xsou'],par['zsou'],par)
    Plot('ss',fdmod.ssplot('',par))
    Result('vel',['vel','ss'],'Overlay')
    Result('ref',['ref','ss'],'Overlay')
    Result('img',['img','ss'],'Overlay')

    # slowness
    Flow('slo','vel',
         '''
         transp |
         math "output=1/input" |
         spray axis=2 n=1 |
         put label2=y
         ''' % par )
    Result('slo','window | transp |' + fdmod.cgrey('allpos=y bias=0.065',par))

    # migration wavelet
    Flow('wvl',None,
         '''
         spike nsp=1 mag=1
         n1=%(nt)d d1=%(dt)g o1=0        k1=1
         n2=1      d2=%(dx)g o2=%(xsou)g |
         scale axis=123 |
         put label1=t label2=x label3=y
         ''' % par)

    # modeling wavelet (time domain)
    Flow('wav',None,
         '''
         spike nsp=1 mag=1
         n1=%(nt)d d1=%(dt)g o1=0   k1=%(kt)d
         n2=1      d2=%(dx)g o2=%(xsou)g |
         ricker1 frequency=%(frq)g |
         scale axis=123 |
         put label1=t label2=x label3=y
         ''' % par)
    Result('wav','window n1=200 |' + fdmod.waveplot('',par))

    # modeling wavelet (frequency domain)
    Flow('sou','wav',
         '''
         fft1 |
         window squeeze=n n1=%(nw)d min1=%(ow)g |
         pad beg2=%(xpad)d n2out=%(nx)d |
         put label1=w label2=x label3=y |
         transp memsize=250 plane=12 |
         transp memsize=250 plane=23 
         ''' % par)

    # global slowness perturnation
    Flow('ds',None,
         '''
         spike nsp=1 mag=0.00005
         n1=%(nz)d d1=%(dz)g o1=%(oz)g
         n2=%(nx)d d2=%(dx)g o2=%(ox)g |
         put label1=z label2=x label3=y |
         smooth rect1=1 rect2=1 |
         transp plane=12 |
         transp plane=23 |
         rtoc
         ''' % par)
    Plot('ds','window | real | transp |'+ fdmod.cgrey('pclip=99.9',par))

    # ------------------------------------------------------------
    # zero-offset
    # ------------------------------------------------------------
    # data
    zomig.model3('Zmod','slo','img',par)
    Flow('Zdat','Zmod',
         '''
         transp plane=23 |
         transp plane=12 |
         pad beg1=%(begw)d n1out=%(padw)d |
         fft1 inv=y |
         put o1=0 label1=t label2=x
         ''' % par)
    Result('Zdat',fdmod.dgrey('screenratio=0.5 screenht=7',par))
    
    # wavefield
    zomig.Awftwo3('woz','Zmod','slo',par)

    # migration
    zomig.image3('Zimg','slo','Zmod',par)
    Plot(  'Zimg','window | transp |' + fdmod.cgrey('',par))
    Result('Zimg','Zimg','Overlay')

    # WEMVA zero-offset
    zomig.s2i(  'ds' ,'ZFds','woz','slo',par) # forward   F(ds)
    zomig.i2s('ZFds','ZAFds','woz','slo',par) # adjoint A(F(ds))
    Result( 'ZFds','window | real | transp |'+ fdmod.cgrey('',par))
    Result('ZAFds','window | real | transp |'+ fdmod.cgrey('pclip=98',par))

    # ------------------------------------------------------------
    # shot-record 
    # ------------------------------------------------------------
    # data
    spmig.modelPW3('Smod','slo','sou','ref',par)
    Flow('Sdat','Smod',
         '''
         transp plane=23 |
         transp plane=12 |
         pad beg1=%(begw)d n1out=%(padw)d |
         fft1 inv=y |
         window f1=%(kt)d |
         pad n1out=%(nt)d |
         put o1=0 o2=%(xoff)g o3=%(xsou)g
         ''' % par)
    Result('Sdat',
           fdmod.dgrey('min2=%g max2=%g label2="Offset" screenratio=0.5 screenht=7'
                       %(par['xoff'],-par['xoff']),par))
    
    # wavefields
    spmig.wflds  ('dos','dor','wvl','Sdat',par)
    zomig.Cwfone3('wos','dos','slo',par) #   source
    zomig.Awfone3('wor','dor','slo',par) # receiver
    
    # migration
    spmig.imagePW3('Simg','cig','slo','dos','dor',par)
    Plot(  'Simg','window | transp |'+ fdmod.cgrey('',par))
    Result('Simg','Simg ss','Overlay')

    # WEMVA shot-record
    spmig.s2i(  'ds', 'SFds','wos','wor','slo',par) # forward   F(ds)
    spmig.i2s('SFds','SAFds','wos','wor','slo',par) # adjoint A(F(ds))
    Result( 'SFds','window | real | transp |'+ fdmod.cgrey('',par))
    Result('SAFds','window | real | transp |'+ fdmod.cgrey('pclip=98',par))

    # ------------------------------------------------------------
    for ispk in range(par['nspk']):
        i = par['ospk'] + ispk
        tag = str(i)

        par['xx']=par['nx']/2+i*par['jspk']
        par['xs']=par['xx']    # x start
        par['xe']=par['xx']    # x end
        par['zs']=par['nz']/2  # z start
        par['ze']=par['nz']/2  # z end

        # slowness perturbation
        Flow('ds'+tag,None,
             '''
             spike nsp=1 mag=1
             n1=%(nz)d d1=%(dz)g o1=%(oz)g k1=%(zs)d l1=%(ze)d
             n2=%(nx)d d2=%(dx)g o2=%(ox)g k2=%(xs)d l2=%(xe)d |
             put label1=z label2=x label3=y |
             smooth rect1=6 rect2=2 repeat=3 |
             scale axis=123 |
             scale rscale=0.00005 |
             transp plane=12 |
             transp plane=23 |
             rtoc
             ''' % par)
        Result('ds'+tag,'window | real | transp |'+ fdmod.cgrey('',par))

        # image perturbation
        Flow('di'+tag,'ZFds msk',
             '''
             window squeeze=n n1=1 f1=%(xx)d |
             pad beg1=%(xx)d n1out=%(nx)d |
             put o1=%(ox)g |
             math m=${SOURCES[1]} output="input*m"
             ''' %par)
        Result('di'+tag,'window | real | transp | smooth rect2=3 repeat=3 |'+ fdmod.cgrey('',par))

        # WEMVA zero-offset
        zomig.s2i(  'ds'+tag, 'ZFds'+tag,'woz','slo',par) # forward   F(ds)
        zomig.i2s('ZFds'+tag,'ZAFds'+tag,'woz','slo',par) # adjoint A(F(ds))
        zomig.i2s(  'di'+tag, 'ZAdi'+tag,'woz','slo',par) # adjoint   A(di)
        zomig.s2i('ZAdi'+tag,'ZFAdi'+tag,'woz','slo',par) # forward F(A(di))

        Result(   'ZFds'+tag,'window | real | transp |'+ fdmod.cgrey('',par))
        Result(  'ZAFds'+tag,'window | real | transp |'+ fdmod.cgrey('pclip=99.9',par))
        Result(   'ZAdi'+tag,'window | real | transp |'+ fdmod.cgrey('pclip=99.9',par))
        Result(  'ZFAdi'+tag,'window | real | transp |'+ fdmod.cgrey('',par))

        # WEMVA shot-record
        spmig.s2i(  'ds'+tag, 'SFds'+tag,'wos','wor','slo',par) # forward   F(ds)
        spmig.i2s('SFds'+tag,'SAFds'+tag,'wos','wor','slo',par) # adjoint A(F(ds))
        spmig.i2s(  'di'+tag, 'SAdi'+tag,'wos','wor','slo',par) # adjoint   A(di)
        spmig.s2i('SAdi'+tag,'SFAdi'+tag,'wos','wor','slo',par) # forward F(A(di))

        Result(   'SFds'+tag,'window | real | transp |'+ fdmod.cgrey('',par))
        Result(  'SAFds'+tag,'window | real | transp |'+ fdmod.cgrey('pclip=99.9',par))
        Result(   'SAdi'+tag,'window | real | transp |'+ fdmod.cgrey('pclip=99.9',par))
        Result(  'SFAdi'+tag,'window | real | transp |'+ fdmod.cgrey('',par))
        
