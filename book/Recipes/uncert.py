from rsf.proj import *
import velcon

def uncert(data,        # data name
           nv,          # continuation steps
           v0,          # initial velocity
           dv,          # velocity step
           nx,          # lateral dimension
           nh,          # number of offsets
           padt,        # time padding
           padt2,       # extra time padding
           padx=None,   # lateral padding
           v1=None,     # other velocity
           n1=None,     # time extent
           dt=0.004,    # time sampling
           dx=0.01,     # midpoint sampling
           units='km',  # lateral units
           vslope=None, # semblance muting
           vx0=0,       # semblance muting
           x0=0,        # lateral origin
           rect1=10,    # vertical smoothing
           rect2=10):   # lateral  smoothing

    velcon.velcon(data,nv,v0,dv,nx,nh,padt,padt2,padx,v1,n1,dt,dx,units,vslope,vx0,x0,rect1,rect2)
    
    vlf=data+'-vlf'
    vlf2=data+'-vlf1'
    npk = data+'-npk'

    # To estimate uncertainty: measure dt/dv, measure dv, multiply
    ref = data+'-ref'
    Flow(ref,[vlf,npk],
         '''
         transp | refer ref=${SOURCES[1]} | window n1=9 f1=%d | transp
         ''' % (nv-5))

    if vslope:
        refer = '''
        mutter x0=%g v0=%g half=n |
        transp | refer ref=${SOURCES[1]} | transp
        ''' % (vx0,vslope)
    else:
        refer = '''
        transp | refer ref=${SOURCES[1]} | transp
        '''
    
    Flow(ref+'2',[vlf2,npk],refer)

    dtdv = data+'-dtdv'
    Flow(dtdv,ref,
         '''
         dip n4=0 rect1=%d rect2=5 rect3=%d |
         window n2=1 min2=0 |
         scale dscale=%g
         ''' % (rect1,rect2,dt/dv))
    Flow(dtdv+'0',[dtdv,npk],
         'math vel=${SOURCES[1]} output="vel*x1*input*input" ')
    Result(dtdv,
           '''
           grey title="Structural Sensitivity in T"
           label1=Time unit1=s label2="Lateral Position" unit2=%s
           scalebar=y color=j allpos=y barlabel="dt/dv" barunit="s\^\s75 2\s100 \_/%s"
           ''' % (units,units))

    dxdv = data+'-dxdv'
    Flow(dxdv,ref,
         '''
         transp plane=13 |
         dip n4=0 rect1=%d rect2=5 rect3=%d |
         window n2=1 min2=0 |
         scale dscale=%g | transp
         ''' % (rect1,rect2,dx/dv))
    Result(dxdv,
           '''
           grey title="Structural Sensitivity in X"
           label1=Time unit1=s label2="Lateral Position" unit2=%s
           scalebar=y color=j barlabel="dx/dv" barunit=s
           ''' % units)

    ddv = data+'-ddv'
    Flow(ref+'3',ref+'2','stack norm=n')
    Flow(ddv,[ref+'2',ref+'3'],
         '''
         math output="x2*x2*input" |
         stack norm=n |
         add mode=d ${SOURCES[1]} |
         math output="sqrt(input)"
         ''')   
    Result(ddv,
           '''
           grey title="Velocity Uncertainty" allpos=y
           label1=Time unit1=s label2="Lateral Position" unit2=%s
           scalebar=y color=j barlabel=Velocity barunit="%s/s"
           ''' % (units,units))

    unc = data+'-unc'
    Flow(unc,[dtdv,ddv],'math dv=${SOURCES[1]} output="abs(0.5*input*dv)" ')
    Result(unc,
           '''
           grey title="Structural Uncertainty"
           scalebar=y color=j allpos=y
           label1=Time unit1=s label2="Lateral Position" unit2=%s
           barlabel="Vertical Uncertainty (s)"
           ''' % units)
    Flow(unc+'2',[dxdv,ddv],
         'math dv=${SOURCES[1]} output="abs(0.5*input*dv)" ')
    Result(unc+'2',
           '''
           grey title="Structural Uncertainty"
           scalebar=y color=j allpos=y
           label1=Time unit1=s label2="Lateral Position" unit2=%s
           barlabel="Lateral Uncertainty (%s)"
           ''' % (units,units))
    
    ddip=data+'-ddip'
    Flow(ddip,data,
         '''
         window |
         dip rect1=%d rect2=%d rect3=1
         ''' % (rect1,rect2))
    dpwd=data+'-dpwd'
    Flow(dpwd,[data,ddip],'window | pwd dip=${SOURCES[1]}')

    data2 = data+'2'
    Flow(data2,dpwd,'window n4=1 | spray axis=2 n=1')

    mig2=data2+'-mig'
    Flow(mig2,data2,'preconstkirch vel=%g' % v0)

    cip2=data2+'-cip'
    Flow(cip2,mig2,velcon.mig2cip)

    if padx:
        pad2=data2+'-pad'
        Flow(pad2,cip2,'pad n3=%d' % padx)
    else:
        pad2=cip2

    ckx2=data2+'-ckx'
    vlf2=data2+'-vlf'
    Flow(ckx2,pad2,'cosft sign3=1')
    Flow(vlf2,ckx2,
         '''
         fourvc nv=%d dv=%g v0=%g pad=%d pad2=%d |
         cosft sign3=-1 | window n3=%d
         ''' % (nv,dv,v0,padt,padt2,nx))

    
    vlf22=data2+'-vlf2'
    Flow(vlf22,pad2,
         '''
         transp plane=23 memsize=500 |
         fourvc2 nv=%d dv=%g v0=%g pad=%d pad2=%d |
         window n2=%d | transp plane=23 memsize=500
         ''' % (nv,dv,v0,padt,padt2,nx))

    
    foc=data2+'-foc'
    Flow(foc,[vlf2,vlf22],
         '''
         focus rect1=%d rect3=%d |
         math vlf=${SOURCES[1]} output=vlf/input
         ''' % (2*rect1,2*rect2))
