from rsfproj import *

def velcon(data,        # data name
           nv,          # continuation steps
           v0,          # initial velocity
           dv,          # velocity step
           nx,          # lateral dimension
           padt,        # time padding
           padt2,       # extra time padding
           padx=None,   # lateral padding
           v1=None,     # other velocity
           n1=None,     # time extent
           rect1=10,    # vertical smoothing
           rect2=10):   # lateral  smoothing
    '''Velocity continuation'''
    
    mig=data+'-mig'
    Flow(mig,data,'preconstkirch vel=%g' % v0)

    if n1:
        mig2cip = '''
        transp plane=24 memsize=500 | halfint inv=1 adj=1 | window n1=%d
        ''' % n1
    else:
        mig2cip = '''
        transp plane=24 memsize=500 | halfint inv=1 adj=1
        '''

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
        Flow(mig+'1',data,'preconstkirch vel=%g' % v1)
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

    npk = data+'-npk'
    Flow(npk,vlf,
         '''
         pick rect1=%d rect2=%d | transp plane=23 memsize=500
         ''' % (rect1,rect2))
    Plot(npk,
         '''
         grey pclip=100 color=j bias=%g 
         scalebar=y title="Picked RMS Velocity"
         label1="Time (s)" label2="Lateral (km)"
         barlabel="Velocity (km/s)"
         ''' % (v0+0.5*nv*dv))
    for line in (0,1):
        Plot(npk+str(line),npk,
             '''
             contour nc=40 wanttitle=n plotcol=%d plotfat=%d
             scalebar=y wantaxis=n barlabel="Velocity (km/s)"
             ''' % ((0,10),(7,1))[line])
    Result(npk,[npk,npk+'0',npk+'1'],'Overlay')

    fmg = data+'-fmg'
    Flow(fmg,[vlf,npk],'slice pick=${SOURCES[1]} | transp plane=23')
    Result(fmg,'grey title=Slice label1="Time (s)" label2="Lateral (km)" ')

    agc = data+'-agc'
    Flow(agc,fmg,'agc rect1=200')
    Result(agc,
           'grey title=Picked pclip=98 label1="Time (s)" label2="Lateral (km)" ')
    
    Plot(fmg,agc,
         'grey title="Time-Migrated Image" label1="Time (s)" label2="Lateral (km)" ')
    Result(fmg+'2',[fmg,npk],'SideBySideIso',vppen='txscale=1.2')

    
