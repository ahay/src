from rsfproj import *

def stack(name,
          v0,
          nv,
          dv,
          rect1=10,
          rect2=10,
          vslope=None,
          f1=1,
          nout=2048,
          vx0=None):

    scn = name+'-scn'
    vel = name+'-vel'
    
    Flow(scn,name,'mutter v0=%g | vscan semblance=y v0=%g nv=%d dv=%g' % (v0,v0,nv,dv))

    if vslope:
        pick = 'mutter x0=%g v0=%g half=n | pick rect1=%d rect2=%d | window' % (vx0,vslope,rect1,rect2)
    else:
        pick = 'pick rect1=%d rect2=%d | window' % (rect1,rect2)
        
    Flow(vel,scn,pick)

    Result(vel,'grey color=j scalebar=y title="RMS velocity" bias=%g' % (v0+0.5*nv*dv))

    nmo = name+'-nmo'
    stk = name+'-stk'

    Flow(nmo,[name,vel],'mutter v0=%g | nmo velocity=${SOURCES[1]}' % v0)
    Flow(stk,nmo,'stack')
    Result(stk,'grey title="NMO Stack" ')

    Flow(stk+'2',nmo,
         '''
         window f1=%d | logstretch nout=%d |
         fft1 | transp plane=13 |
         finstack |
         transp |
         fft1 inv=y | window n1=%d |
         logstretch inv=y | pad beg1=%d
         ''' % (f1,nout,nout,f1))
    Result(stk+'2','grey title="DMO Stack" ')

    dip = name+'-dip'
    Flow(dip,stk+'2','dip rect1=%d rect2=%d' % (rect1,rect2))
    Result(dip,'grey title=Slope color=j scalebar=y')

    pwd=name+'-pwd'
    Flow(pwd,[stk+'2',dip],'pwd dip=${SOURCES[1]}')
    Result(pwd,'grey title=Diffractions')

    vlf=name+'-vlf'
    Flow(vlf,pwd,
         '''
         pad n2=521 | cosft sign2=1 | spray axis=2 n=1 o=0 d=1 |
         stolt vel=1.5 nf=4 |
         fourvc nv=80 dv=0.0125 v0=1.3 |
         cosft sign3=-1 |
         window n3=250
         ''')
