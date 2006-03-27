from rsfproj import *
import string, sys


def stack(name,
          v0,
          nv,
          dv,
          nx,
          padx,
          nt,
          tmin=0,
          tmax=10,
          rect1=10,
          rect2=10,
          srect1=1,
          srect2=3,
          vslope=None,
          units='km',
          f1=1,
          j3=1,
          dx=1,
          x0=0,
          nout=2048,
          vx0=None):


    version = string.replace(sys.version,"+","")
    version = string.split(string.split(version)[0], ".")
    if map(int, version) < [2, 2, 0]:
        return # think how to do it better

    scn = name+'-scn'
    vel = name+'-vel'
    
    Flow(scn,name,
         'mutter v0=%g | vscan semblance=y v0=%g nv=%d dv=%g' % (v0,v0,nv,dv))

    if vslope:
        pick = 'mutter x0=%g v0=%g half=n | ' % (vx0,vslope)
    else:
        pick = ''

    pick = pick + 'pick rect1=%d rect2=%d | window' % (rect1,rect2)
        
    def grey(title):
        return '''
        window n1=%d |
        grey title="%s" 
        label1="Time (s)" label2="Lateral (%s)"
        ''' % (nt,title,units)

    def velgrey(title):
        return grey(title) + '''
        color=j scalebar=y bias=%g barlabel="Velocity (%s/s)"
        barreverse=y
        ''' % (v0+0.5*nv*dv,units)

    Flow(vel,scn,pick)   
    Result(vel,velgrey('RMS Velocity'))
 
    nmo = name+'-nmo'
    stk = name+'-stk'

    Flow(nmo,[name,vel],'mutter v0=%g | nmo velocity=${SOURCES[1]}' % v0)
    Flow(stk,nmo,'stack')
    Result(stk,grey('NMO Stack'))

    Flow(stk+'2',nmo,
         '''
         window f1=%d | logstretch nout=%d |
         fft1 | transp plane=13 memsize=500 |
         finstack |
         transp memsize=500 |
         fft1 inv=y | window n1=%d |
         logstretch inv=y | pad beg1=%d
         ''' % (f1,nout,nout,f1))
    Result(stk+'2',grey('DMO Stack'))

    dip = name+'-dip'
    Flow(dip,stk+'2','dip rect1=%d rect2=%d' % (rect1,rect2))
    Result(dip,grey('Dominant Slope') + \
           ' color=j scalebar=y barlabel="Slope (samples)" ')

    pwd=name+'-pwd'
    Flow(pwd,[stk+'2',dip],'pwd dip=${SOURCES[1]}')
    Result(pwd,grey('Separated Diffractions'))

    shp=name+'-shp'
    Flow(shp,[stk+'2',dip],
         '''
         pwdsmooth2 dip=${SOURCES[1]} rect1=%d rect2=%d |
         add ${SOURCES[0]} scale=-1,1
         ''' % (srect1,srect2))
    Result(shp,grey('Diffractions'))

#    dips = name+'-dips'
#    Flow(dips,nmo,'dip rect1=%d rect2=2 rect3=%d' % (rect1,rect2))

#    pwds = name+'-pwds'
#    Flow(pwds,[nmo,dips],'pwd dip=${SOURCES[1]}')


#    difs = name+'-difs'
#    Flow(difs,pwds,
#         '''
#         window n4=1 f1=%d | logstretch nout=%d |
#         fft1 | transp plane=13 |
#         finstack |
#         transp |
#         fft1 inv=y | window n1=%d |
#         logstretch inv=y | pad beg1=%d
#         ''' % (f1,nout,nout,f1))
#    Result(difs,grey('DMO Stack of Diffractions'))

#    dif = name+'-dif'
#    Flow(dif,[pwds,vel],'window f4=1 | inmo velocity=${SOURCES[1]}')

    velcon = '''
    pad n2=%d | cosft sign2=1 | 
    stolt vel=%g nf=4 | spray axis=2 n=1 o=0 d=1 |
    fourvc nv=%d dv=%g v0=%g |
    cosft sign3=-1 | 
    window n3=%d |
    put o3=%g
    ''' % (padx,v0,nv,dv,v0,nx,x0)

    vlf=name+'-vlf'
    Flow(vlf,shp,velcon)

    if j3 > 1:
        focus = 'window j3=%d | ' % j3
    else:
        focus = ''
        
    focus = focus + '''    
    focus rect1=%d rect3=%d |
    math output="1/abs(input)" |
    cut max1=%g | cut min1=%g 
    ''' % (2*rect1,2*rect2,tmin,tmax)
 
    foc=name+'-foc'
    Flow(foc,vlf,focus)

    pik=name+'-pik'

    if j3 > 1:
        pick2 = pick + ''' |
        transp |
        spline n1=%d d1=%g o1=%g |
        transp
        ''' % (nx,dx,x0)
    else:
        pick2 = pick
        
    Flow(pik,foc,pick2)
    Result(pik,velgrey('Migration Velocity'))

    slc=name+'-slc'
    Flow(slc,[vlf,pik],'slice pick=${SOURCES[1]}')

    Result(slc,grey('Migrated Diffractions'))

    Flow(vlf+'2',stk+'2',velcon)
    Flow(slc+'2',[vlf+'2',pik],'slice pick=${SOURCES[1]}')

    Result(slc+'2',grey('Migrated Stack'))

    Flow(slc+'1',slc+'2','agc rect1=200')

    Result(slc+'1',grey('Migrated Stack'))
