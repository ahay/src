from rsfproj import *
import zomig,spmig,pplot

def igrey(custom,par):
    return '''
    grey  labelrot=n wantaxis=y wanttitle=y title="" pclip=100
    label1="z" label2="x"
    %s %s
    ''' % (par['greypar'],custom)

def test(par):
    par['zmin']=par['oz']
    par['zmax']=par['oz'] + (par['nz']-1) * par['dz']
    par['xmin']=par['ox']
    par['xmax']=par['ox'] + (par['nx']-1) * par['dx']
    
    par['xpad']=par['nx']/2.
    par['xsou']=par['ox'] + par['xpad'] * par['dx']
    par['xoff']=          - par['xpad'] * par['dx']
    
    par['dw']= 1/(par['nt']*par['dt'])

    par['begw']=par['ow']/par['dw']
    par['padw']=par['nt']/2+1

    par['kt']=125

    # slowness
    Flow('slo','vel',
         '''
         transp |
         math "output=1/input" |
         spray axis=2 n=1 |
         put label2=y
         ''' % par )
    Result('slo','slo','window | transp |' + igrey('color=j allpos=y',par))

    # migration wavelet
    Flow('wvl',None,
         '''
         spike nsp=1 mag=1
         n1=%(nt)d d1=%(dt)g o1=0        k1=1
         n2=1      d2=%(dx)g o2=%(xsou)g |
         put label1=t label2=x label3=y
         ''' % par)

    # modeling wavelet (time domain)
    Flow('wav',None,
         '''
         spike nsp=1 mag=1
         n1=%(nt)d d1=%(dt)g o1=0        k1=%(kt)d
         n2=1      d2=%(dx)g o2=%(xsou)g |
         ricker1 frequency=%(frq)g |
         put label1=t label2=x label3=y
         ''' % par)
    Result('wav','wav','graph title=" "')

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

    # ------------------------------------------------------------
    # DATA/IMAGE (by exploding reflector modeling)
    # ------------------------------------------------------------
    Flow('vvv','ref',
         '''
         window | transp |
         math output=%(vel)g
         ''' % par )

    Flow('img',['ref','vvv'],
         '''
         window | transp |
         depth2time velocity=${SOURCES[1]} dt=%(dt)g nt=%(nt)d |
         ricker1 frequency=%(frq)g |
         time2depth velocity=${SOURCES[1]} dz=%(dz)g nz=%(nz)d |
         transp |
         spray axis=2 n=1 o=0 d=1 |
         put label1=x label2=y label3=z 
         ''' % par )    
    Result('img','img','window | transp |' + igrey('',par))
    
    zomig.model('modZ','slo','img',par)
    zomig.image('iii','slo','modZ',par)
    
    Flow('zoff','modZ',
         '''
         transp plane=23 |
         transp plane=12 |
         pad beg1=%(begw)d n1out=%(padw)d |
         fft1 inv=y |
         put o1=0 label1=t label2=x
         ''' % par)
    Result('zoff','zoff','grey pclip=100')
    
    # wavefields
    zomig.wflds('doz','zoff',par)
    zomig.Awftwo('woz','doz','slo',par)

    # migration
    zomig.image('imgZ','slo','doz',par)
    Plot('imgZ','imgZ','window | transp |' + igrey('title="image" pclip=99',par))
    # ------------------------------------------------------------

    # ------------------------------------------------------------
    # DATA/IMAGE (by shot-record modeling)
    # ------------------------------------------------------------
    spmig.model('modS','slo','sou','ref',par)
    Flow('sone','modS',
         '''
         transp plane=23 |
         transp plane=12 |
         pad beg1=%(begw)d n1out=%(padw)d |
         fft1 inv=y |
         window f1=%(kt)d |
         pad n1out=%(nt)d |
         put o1=0 o2=%(xoff)g o3=%(xsou)g
         ''' % par)
    Result('sone','sone','grey pclip=100')
    
    # wavefields
    spmig.wflds('dos','dor','wvl','sone',par)
    zomig.Cwfone('wos','dos','slo',par) # source
    zomig.Awfone('wor','dor','slo',par) # receiver
    
    # migration
    spmig.imagePW('imgS','cig','slo','dos','dor',par)
    Plot('imgS','imgS','window | transp |'+ igrey('title="img" pclip=99',par))
    
    # ------------------------------------------------------------
    # slowness anomalies at different locations
    for i in range(3):
        par['xs']=par['nx']/2+i*20 - 1 # x start
        par['xe']=par['nx']/2+i*20 + 1 # x end

        par['zs']=par['nz']/2 - 3
        par['ze']=par['nz']/2 + 3
        
        dsl = 'dsl' + str(i)
        Flow(dsl,None,
             '''
             spike nsp=1 mag=0.00005
             n1=%(nz)d d1=%(dz)g o1=%(oz)g k1=%(zs)d l1=%(ze)d
             n2=%(nx)d d2=%(dx)g o2=%(ox)g k2=%(xs)d l2=%(xe)d |
             put label1=z label2=x label3=y |
             smooth rect1=1 rect2=1 |
             transp plane=12 |
             transp plane=23 |
             rtoc
             ''' % par)
        Plot(dsl,dsl,'window | real | transp |'+ igrey('title="dsl" color=F',par))

        # WEMVA zero-offset
        scatZ= 'scatZ'+ str(i)
        dimZ = 'dimZ' + str(i)
        bslZ = 'bslZ' + str(i)
        zomig.s2i(dsl ,dimZ,'woz','slo',par) # forward
        zomig.i2s(dimZ,bslZ,'woz','slo',par) # adjoint
        
        Plot(dimZ,dimZ,'window | real | transp |'+ igrey('title="dim" pclip=99',par))
        Plot(bslZ,bslZ,'window | real | transp |'+ igrey('title="bsl" color=F pclip=99.9',par))
        pplot.p2x2(scatZ,'imgZ',dsl,dimZ,bslZ,0.5,0.55,-10,-12)

        # WEMVA shot-record
        scatS= 'scatS'+ str(i)
        dimS = 'dimS' + str(i)
        bslS = 'bslS' + str(i)
        spmig.s2i(dsl, dimS,'wos','wor','slo',par) # forward
        spmig.i2s(dimS,bslS,'wos','wor','slo',par) # adjoint
        
        Plot(dimS,dimS,'window | real | transp |'+ igrey('title="dim" pclip=99',par))
        Plot(bslS,bslS,'window | real | transp |'+ igrey('title="bsl" color=F pclip=99.9',par))
        pplot.p2x2(scatS,'imgS',dsl,dimS,bslS,0.5,0.55,-10,-12)
        
    # ------------------------------------------------------------
    
    par['mx']=par['nx']/2
    Flow('dim','dimZ0',
         '''
         window n1=1 f1=%(mx)d |
         spray axis=1 n=%(nx)d o=%(ox)g d=%(dx)g |
         spray axis=2 n=1 o=0 d=1
         ''' % par)
    Plot('dim','dim','window | real | transp |'+ igrey('pclip=100',par))
    
    Flow('don','dimZ0',
         '''
         window | transp | stack |
         transp plane=12 |
         transp plane=23 |
         put o1=%(ox)g d1=%(dx)g
         ''' % par)

    zomig.i2s('dim','bakZ','woz','slo',par)
    Plot('bakZ','bakZ','window | real | transp |'+ igrey('color=F pclip=99',par))
    Result('kerZ',['imgZ','dim','bakZ'],'OverUnderAniso')

    spmig.i2s('dim','bakS','wos','wor','slo',par)
    Plot('bakS','bakS','window | real | transp |'+ igrey('color=F pclip=99',par))
    Result('kerS',['imgS','dim','bakS'],'OverUnderAniso')

    # ------------------------------------------------------------
    # sensitivity kernels
    for j in range(3):
        par['xx']=par['nx']/2+j*20
        
        dim = 'dim' + str(j)
        Flow(dim,'don','pad beg1=%(xx)d n1out=%(nx)d | put o1=%(ox)d' % par)
        Plot(dim,dim,'window | real | transp |'+ igrey('pclip=100',par))

        bakZ = 'bakZ' + str(j)
        zomig.i2s(dim,bakZ,'woz','slo',par)
        Plot(bakZ,bakZ,'window | real | transp |'+ igrey('color=F pclip=99.9',par))
        kerZ = 'kerZ' + str(j)
        Result(kerZ,['imgZ',dim,bakZ],'OverUnderAniso')

        bakS = 'bakS' + str(j)
        spmig.i2s(dim,bakS,'wos','wor','slo',par)
        Plot(bakS,bakS,'window | real | transp |'+ igrey('color=F pclip=99.9',par))
        kerS = 'kerS' + str(j)
        Result(kerS,['imgS',dim,bakS],'OverUnderAniso')
