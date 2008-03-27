from rsfproj import *
import dix

mig2cip = None

def velcon(data,        # data name
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
    '''Velocity continuation'''
    global mig2cip

    vm = v0+0.5*nv*dv
    
    mig=data+'-mig'
    Flow(mig,data,'preconstkirch vel=%g' % v0,split=[0],axis=[4,nh])

    if n1:
        mig2cip = '''
        transp plane=24 memsize=500 | halfint inv=1 adj=1 | window n1=%d
        ''' % n1
    else:
        mig2cip = '''
        transp plane=24 memsize=500 | halfint inv=1 adj=1
        '''
        n1=100

    cip=data+'-cip'
    Flow(cip,mig,mig2cip)

    if padx:
        pad=data+'-pad'
        Flow(pad,cip,'pad n3=%d' % padx)
    else:
        pad=cip

    ckx=data+'-ckx'
    vlf=data+'-vlf'
    vlf2=data+'-vlf2'
    Flow(ckx,pad,'cosft sign3=1')
    Flow(vlf,ckx,
         '''
         fourvc nv=%d dv=%g v0=%g pad=%d pad2=%d |
         cosft sign3=-1 | window n3=%d
         ''' % (nv,dv,v0,padt,padt2,nx))

    Flow(vlf2,pad,
         '''
         transp plane=23 memsize=500 |
         fourvc2 nv=%d dv=%g v0=%g pad=%d pad2=%d |
         window n2=%d | transp plane=23 memsize=500
         ''' % (nv,dv,v0,padt,padt2,nx))
    
    if v1:
        Flow(mig+'1',data,'preconstkirch vel=%g' % v1,split=[0],axis=[4,nh])
        Flow(cip+'1',mig+'1',mig2cip)

        migr = data+'-migr'

        Flow(migr,cip,'stack norm=y')
        Plot(migr,'grey title=Migration0')

        Flow(migr+'1',cip+'1','stack norm=y')
        Plot(migr+'1','grey title=Migration1')

        vlfr = data+'-vlfr'

        Flow(vlfr,vlf,'window n2=1 min2=%g' % v1)
        Plot(vlfr,'grey title="Velocity Continuation 0 -> 1" ')

        Result(migr,[migr,migr+'1',vlfr],'SideBySideAniso')


    if vslope:
        pick = '''
        mutter x0=%g v0=%g half=n |
        pick rect1=%d rect2=%d | transp plane=23 memsize=500
        ''' % (vx0,vslope,rect1,rect2)
    else:
        pick = '''
        pick rect1=%d rect2=%d | transp plane=23 memsize=500
        ''' % (rect1,rect2)

    npk = data+'-npk'
    Flow(npk,vlf2,pick)
    Plot(npk, 	 
         ''' 	 
         grey pclip=100 color=j bias=%g 	 
         scalebar=y title="Picked Migration Velocity" 	 
         label1=Time unit1=s label2="Lateral Position" unit2=%s 	 
         barlabel=Velocity barunit="%s/s" barreverse=y 	 
         ''' % (vm,units,units))
        
    fmg = data+'-fmg'
    Flow(fmg,[vlf,npk],'slice pick=${SOURCES[1]} | transp plane=23')

    Result(fmg,
           'grey title=Slice label1=Time unit1=s label2="Lateral Position" unit2=%s ' % units)

    agc = data+'-agc'
    Flow(agc,fmg,'agc rect1=200')
    Plot(fmg,agc,
         '''
         grey title="Time-Migrated Image" label1="Time" unit1=s label2="Lateral Position" unit2=%s
         ''' % units)
    Result(agc,
           '''
           grey title=Picked pclip=98 label1="Time" unit1=s label2="Lateral Position" unit2=%s
           ''' % units)

    Result(fmg+'2',[fmg,npk],'SideBySideAniso',vppen='txscale=1.2')

    Flow(agc+'2',[vlf2,npk],'slice pick=${SOURCES[1]} | window')
    
