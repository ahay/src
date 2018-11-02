from __future__ import absolute_import
from rsf.proj import *
import fdd,pot,pplot

s_coefm = '-0.000348,+0.004785,-0.039876,+0.598144'
s_coefp = '-0.598144,+0.039876,-0.004785,+0.000348'

c_coefm ='-0.0027,+0.0285,-0.1495,+0.5981'
c_coefp ='-0.5981,+0.1495,-0.0285,+0.0027'

#c_coefm = '0,0,0,+1'
#c_coefp = '-1,0,0,0'

# ------------------------------------------------------------
def opplot(title,custom,par):
    return '''
    grey title="%s"
    screenratio=1 pclip=100
    wantscalebar=n wantaxis=n
    %s
    ''' % (title,par['labelattr'] +' '+custom)

# ------------------------------------------------------------
def filter(coef,par):
    return'''
    spike nsp=9 mag=%s
    n1=9 o1=-4 d1=1 k1=1,2,3,4,5,6,7,8,9
    ''' % coef
# ------------------------------------------------------------
def fdx(par):
    return list(filter('%s' %('0,'+s_coefm+','  +s_coefp     ),par)) + '| transp'
def cdx(par):
    return list(filter('%s' %(     c_coefm+',0,'+c_coefp     ),par)) + '| transp'
def bdx(par):
    return list(filter('%s' %(     s_coefm+','  +s_coefp+',0'),par)) + '| transp'
# ------------------------------------------------------------
def fdz(par):
    return list(filter('%s' %('0,'+s_coefm+','  +s_coefp     ),par))
def cdz(par):
    return list(filter('%s' %(     c_coefm+',0,'+c_coefp     ),par))
def bdz(par):
    return list(filter('%s' %(     s_coefm+','  +s_coefp+',0'),par))
# ------------------------------------------------------------
def ddapply(out,inp,flt,stat,par):
    Flow(out,[inp,flt],
         'convolve2 verb=y stat=%s flt=${SOURCES[1]}' % stat)
# ------------------------------------------------------------

def ddapply3d(out,inp,flt,stat,par):
    Flow(out,[inp,flt],
         'boxfilter verb=y stat=%s flt=${SOURCES[1]}' % stat)
# ------------------------------------------------------------
def derivatives(par):
    Flow('fdx',None,fdd.fdx(par))
    Flow('fdz',None,fdd.fdz(par))

    Flow('bdx',None,fdd.bdx(par))
    Flow('bdz',None,fdd.bdz(par))

    Flow('dxC',None,fdd.cdx(par))
    Flow('dzC',None,fdd.cdz(par))

# ------------------------------------------------------------
def oneirST(inp,zdel,xdel,m,n,custom,par):

    par['nzspk']=2*m+1
    par['nxspk']=2*n+1
    par['fzspk']=(par['nzspk']-1)/2+1
    par['fxspk']=(par['nxspk']-1)/2+1

    Flow(inp,None,
         '''
         spike nsp=1 mag=1
         n1=%(nzspk)d o1=%(oz)g d1=%(dz)g k1=%(fzspk)d
         n2=%(nxspk)d o2=%(ox)g d2=%(dx)g k2=%(fxspk)d
         ''' % par)

    for j in (zdel,xdel):
        ddapply(inp+'-'+j,inp,j,'y',par)


    Flow(inp+'-all',[inp+'-'+zdel,inp+'-'+xdel],
         '''
         cat axis=3 space=n ${SOURCES[1]} |
         scale axis=123 |
         byte gainpanel=a pclip=98
         %s
         ''' % custom)
    Plot(inp+'-'+zdel,inp+'-all','window n3=1 f3=0 |'
         + opplot("\F5 L\_z",custom,par))
    Plot(inp+'-'+xdel,inp+'-all','window n3=1 f3=1 |'
         + opplot("\F5 L\_x",custom,par))
#    Plot(inp+'-'+zdel,inp+'-all','window n3=1 f3=0 |' + opplot('',''))
#    Plot(inp+'-'+xdel,inp+'-all','window n3=1 f3=1 |' + opplot('',''))

    if('ys' not in par): par['ys']=0.75
    if('xs' not in par): par['xs']=0.75
    if('xc' not in par): par['xc']=-8.25
    pplot.p1x2(inp,inp+'-'+zdel,inp+'-'+xdel,par['ys'],par['xs'],par['xc'])

#    pplot.p2x1(inp,inp+'-'+zdel,inp+'-'+xdel,0.5,0.5,-9.00)

# ------------------------------------------------------------

def oneirST_box(inp,zdel,xdel,m,n,box,custom,par):

    par['nzspk']=2*m+1
    par['nxspk']=2*n+1
    par['fzspk']=(par['nzspk']-1)/2+1
    par['fxspk']=(par['nxspk']-1)/2+1

    Flow(inp,None,
         '''
         spike nsp=1 mag=1
         n1=%(nzspk)d o1=%(oz)g d1=%(dz)g k1=%(fzspk)d
         n2=%(nxspk)d o2=%(ox)g d2=%(dx)g k2=%(fxspk)d
         ''' % par)

    for j in (zdel,xdel):
        ddapply(inp+'-'+j,inp,j,'y',par)

    Flow(inp+'-all',[inp+'-'+zdel,inp+'-'+xdel],
         '''
         cat axis=3 space=n ${SOURCES[1]} |
         scale axis=123 |
         byte gainpanel=a pclip=99.9|
         put o1=0 d1=1 o2=0 d2=1 label1=Sample# label2=Sample# unit1= unit2=
         %s
         ''' % custom)
    Plot(inp+'-'+zdel+'-',inp+'-all','window n3=1 f3=0 |'
         + opplot("\F5 L\_z ",' wantaxis=y   o2num=0 d2num=10 n2tic=6 o1num=0 d1num=10 n1tic=6 '+custom,par))
    Plot(inp+'-'+xdel+'-',inp+'-all','window n3=1 f3=1 |'
         + opplot("\F5 L\_x ",' wantaxis=y   o2num=0 d2num=10 n2tic=6 o1num=0 d1num=10 n1tic=6 unit1= label1= wantaxis1= '+custom,par))

#    Plot(inp+'-'+zdel,inp+'-all','window n3=1 f3=0 |' + opplot('',''))
#    Plot(inp+'-'+xdel,inp+'-all','window n3=1 f3=1 |' + opplot('',''))

    Plot(inp+'-'+zdel,[inp+'-'+zdel+'-',box],'Overlay')
    Plot(inp+'-'+xdel,[inp+'-'+xdel+'-',box],'Overlay')

    if('ys' not in par): par['ys']=0.75
    if('xs' not in par): par['xs']=0.75
    if('xc' not in par): par['xc']=-8.25
    pplot.p1x2(inp,inp+'-'+zdel,inp+'-'+xdel,par['ys'],par['xs'],par['xc'])

#    pplot.p2x1(inp,inp+'-'+zdel,inp+'-'+xdel,0.5,0.5,-9.00)

# ------------------------------------------------------------
def oneirNS(inp,zdel,xdel,m,n,custom,par):

    for     jx in range(3):
        fx = (jx+1) * par['nx']/4 + 1
        for jz in range(3):
            fz = (jz+1) * par['nz']/4 + 1
            tag = str(jx)+str(jz)

            Flow(zdel+tag,zdel,'window n3=1 f3=%d n4=1 f4=%d' %(fz,fx) )
            Flow(xdel+tag,xdel,'window n3=1 f3=%d n4=1 f4=%d' %(fz,fx) )

            oneirST(inp+tag,zdel+tag,xdel+tag,m,n,custom,par)

# ------------------------------------------------------------

def oneirNS_rot(inp,zdel,xdel,nu,m,n,custom,par):

    Flow('sinnu','nu','math output="sin(input*3.1415926/180)" ')
    nz=par['nzspk']
    nx=par['nxspk']
    print(nz,nx)
    for     jx in range(3):
        fx = (jx+1) * par['nx']/4 + 1
        for jz in range(3):
            fz = (jz+1) * par['nz']/4 + 1
            tag = str(jx)+str(jz)

            Flow(zdel+tag+'-t',zdel,'window n3=1 f3=%d n4=1 f4=%d' %(fz,fx) )
            Flow(xdel+tag+'-t',xdel,'window n3=1 f3=%d n4=1 f4=%d' %(fz,fx) )

            Flow(zdel+tag,['sinnu',zdel+tag+'-t',xdel+tag+'-t'],
                 '''
                 window n1=1 f1=%d n2=1 f2=%d |
                 spray axis=1 n=%d  |
                 spray axis=2 n=%d  |
                 math output="z*sqrt(1-input*input)-x*input"
                      z=${SOURCES[1]} x=${SOURCES[2]}
                 '''%(fz,fx,nz,nx) )
            Flow(xdel+tag,['sinnu',zdel+tag+'-t',xdel+tag+'-t'],
                 '''
                 window n1=1 f1=%d n2=1 f2=%d |
                 spray axis=1 n=%d  |
                 spray axis=2 n=%d  |
                 math output="z*input+x*sqrt(1-input*input)"
                      z=${SOURCES[1]} x=${SOURCES[2]}
                 '''%(fz,fx,nz,nx) )

            oneirST(inp+tag,zdel+tag,xdel+tag,m,n,custom,par)


# ------------------------------------------------------------
# construct derivatives using different sigma in gaussian taper
def separatorD(zdel,xdel,inp,ccc,stat,domain,tapertype,sigma,order,n,m,par):

    Flow([zdel+'-tmp',xdel+'-tmp'],
         [inp,ccc],
         ''' sfederiv2d
         verb=y stat=%s domain=%s ompnth=8
         tapertype=%s sig=%f order=%d
         ccc=${SOURCES[1]}
         zdel=${TARGETS[0]}
         xdel=${TARGETS[1]}
         ''' % (stat,domain,tapertype,sigma,order))

    Flow(zdel,[zdel+'-tmp'],
         '''
         window f1=%d n1=%d f2=%d n2=%d |
         put o1=%g d1=1 o2=%g d2=1
         ''' %(par['spikez']-n-1,2*n+1,
               par['spikex']-m-1,2*m+1,
               -n,-m))

    Flow(xdel,[xdel+'-tmp'],
         '''
         window f1=%d n1=%d f2=%d n2=%d |
         put o1=%g d1=1 o2=%g d2=1
         ''' %(par['spikez']-m-1,2*m+1,
               par['spikex']-n-1,2*n+1,
               -n,-m))



def fftshift(shift,original,par):
    n1=par['nz']/2
    n2=par['nx']/2
    Flow(original+'q1',original,'window n1=%d n2=%d f1=0 f2=0'%(n1,n2))
    Flow(original+'q2',original,'window n1=%d n2=%d f1=0 f2=%d'%(n1,n2,n2))
    Flow(original+'q3',original,'window n1=%d n2=%d f1=%d f2=0'%(n1,n2,n1))
    Flow(original+'q4',original,'window n1=%d n2=%d f1=%d f2=%d'%(n1,n2,n1,n2))

    Flow(original+'q24',[original+'q4',original+'q2'],'cat ${SOURCES[1]} axis=1 space=n')
    Flow(original+'q13',[original+'q3',original+'q1'],'cat ${SOURCES[1]} axis=1 space=n')

    Flow(shift,[original+'q24',original+'q13'],'cat ${SOURCES[1]} axis=2 space=n')
#Flow(shift+'-tmp',[original+'q24',original+'q13'],'cat ${SOURCES[1]} axis=2 space=n')
#    a1=2*m+1
#    a2=2*n+1
#    b1=n1-m
#    b2=n2-n
#    Flow(shift,shift+'-tmp','window n1=%d n2=%d f1=%d f2=%d'%(a1,a2,b1,b2))












def KtoX(zderX,xderX,zderK,xderK,par):

    Flow(zderX+'-t',zderK,
         '''rtoc |
        math output="input*sqrt(-1)*(-1)" |
        fft3 axis=1 opt=n inv=y |
        fft3 axis=2 opt=n inv=y |
        real
              ''')
    Result(zderX,'grey color=E pclip=99 screenratio=1')

    fftshift(zderX,zderX+'-t',par)
    Result(zderX,'grey color=E pclip=99 screenratio=1')


    Flow(xderX+'-t',xderK,
         '''rtoc |
        math output="input*sqrt(-1)*(-1)" |
        fft3 axis=1 opt=n inv=y |
        fft3 axis=2 opt=n inv=y |
        real
              ''')
    Result(xderX,'grey color=E pclip=99 screenratio=1')

    fftshift(xderX,xderX+'-t',par)
    Result(xderX,'grey color=E pclip=99 screenratio=1')




def SepK(up,us,uz,ux,kz,kx,par):

#    Flow(uz+'k',uz,'rtoc | fft3 axis=1 opt=n pad=1 | fft3 axis=2 opt=n pad=1')
#    Flow(ux+'k',ux,'rtoc | fft3 axis=1 opt=n pad=1 | fft3 axis=2 opt=n pad=1')
#wavefield*separator in k domain
    Flow(up+'tmp1',[kz, uz+'k'],
         '''
     rtoc |
     math uz=${SOURCES[1]} output="input*sqrt(-1)*uz*(-1)" |
     fft3 axis=1 inv=y opt=n |
     fft3 axis=2 inv=y opt=n
     ''')
    Flow(up+'tmp2',[kx, ux+'k'],
         '''
     rtoc |
     math ux=${SOURCES[1]} output="input*sqrt(-1)*ux*(-1)" |
     fft3 axis=1 inv=y opt=n |
     fft3 axis=2 inv=y opt=n
     ''')
#sum all dimensions
    Flow(up,[up+'tmp1',up+'tmp2'],
         '''
     math x=${SOURCES[0]} y=${SOURCES[1]} output="x+y"|
     real |
     put d1=%(dz)f d2=%(dx)f o1=%(oz)f o2=%(ox)f
     '''%par)

#    Result(up,'grey color=E screenratio=1')




    Flow(us+'tmp1',[kz, ux+'k'],
         '''
     rtoc |
     math uz=${SOURCES[1]} output="input*sqrt(-1)*uz*(-1)" |
     fft3 axis=1 inv=y opt=n |
     fft3 axis=2 inv=y opt=n
     ''')
    Flow(us+'tmp2',[kx, uz+'k'],
         '''
     rtoc |
     math ux=${SOURCES[1]} output="input*sqrt(-1)*ux*(-1)" |
     fft3 axis=1 inv=y opt=n |
     fft3 axis=2 inv=y opt=n
     ''')

    Flow(us,[us+'tmp1',us+'tmp2'],
         '''
     math x=${SOURCES[0]} y=${SOURCES[1]} output="x-y" |
     real |
     put d1=%(dz)f d2=%(dx)f  o1=%(oz)f o2=%(ox)f
     '''%par)

#    Result(us,'  grey color=E screenratio=1')




def oneirST_S(inp,zdel,xdel,m,n,custom,par):

    par['nzspk']=2*m+1
    par['nxspk']=2*n+1
    par['fzspk']=(par['nzspk']-1)/2+1
    par['fxspk']=(par['nxspk']-1)/2+1

    Flow(inp,None,
         '''
         spike nsp=1 mag=1
         n1=%(nzspk)d o1=%(oz)g d1=%(dz)g k1=%(fzspk)d
         n2=%(nxspk)d o2=%(ox)g d2=%(dx)g k2=%(fxspk)d
         ''' % par)

    for j in (zdel,xdel):
##        ddapply(inp+'-'+j,inp,j,'y',par)
        Flow(inp+'-'+j,inp,'math output="input"')

    Flow(inp+'-all',[inp+'-'+zdel,inp+'-'+xdel],
         '''
         cat axis=3 space=n ${SOURCES[1]} |
         scale axis=123 |
         byte gainpanel=a pclip=98
         %s
         ''' % custom)
    Plot(inp+'-'+zdel,inp+'-all','window n3=1 f3=0 |'
         + opplot("\F15 \\v148 \F3 /\F15 \\v148 \F3 z",custom,par))
    Plot(inp+'-'+xdel,inp+'-all','window n3=1 f3=1 |'
         + opplot("\F15 \\v148 \F3 /\F15 \\v148 \F3 x",custom,par))
#    Plot(inp+'-'+zdel,inp+'-all','window n3=1 f3=0 |' + opplot('',''))
#    Plot(inp+'-'+xdel,inp+'-all','window n3=1 f3=1 |' + opplot('',''))

    if('ys' not in par): par['ys']=0.75
    if('xs' not in par): par['xs']=0.75
    if('xc' not in par): par['xc']=-8.25
    pplot.p1x2(inp,inp+'-'+zdel,inp+'-'+xdel,par['ys'],par['xs'],par['xc'])

#    pplot.p2x1(inp,inp+'-'+zdel,inp+'-'+xdel,0.5,0.5,-9.00)
def KtoX2(zderX,xderX,zderK,xderK,par):

    Flow(zderX+'-t',zderK,
         '''rtoc |
        math output="input*sqrt(-1)*(-1)" |
        fft3d axis=1 cnt=y inv=y |
        fft3d axis=2 cnt=y inv=y |
        real
              ''')
    Result(zderX+'-t','grey color=E pclip=99 screenratio=1')

    fftshift(zderX,zderX+'-t',par)
    Result(zderX,'grey color=E pclip=99 screenratio=1')


    Flow(xderX+'-t',xderK,
         '''rtoc |
        math output="input*sqrt(-1)*(-1)" |
        fft3d axis=1 cnt=y inv=y |
        fft3d axis=2 cnt=y inv=y |
        real
              ''')
    Result(xderX+'-t','grey color=E pclip=99 screenratio=1')

    fftshift(xderX,xderX+'-t',par)
    Result(xderX,'grey color=E pclip=99 screenratio=1')



# construct derivatives
#def separatorP(zdel,xdel,inp,ccc,nu,stat,domain,n,m,par):
#
#    Flow([zdel+'-tmp',xdel+'-tmp'],
#         [inp,ccc,nu,
#          'Code/EDERIV2-k.x'],
#         '''
#         ${SOURCES[3]} verb=y stat=%s domain=%s ompnth=1
#         ccc=${SOURCES[1]}
#          nu=${SOURCES[2]}
#         zdel=${TARGETS[0]}
#         xdel=${TARGETS[1]}
#         ''' % (stat,domain))
#
#    Flow(zdel,[zdel+'-tmp'],
#         '''
#         window f1=%d n1=%d f2=%d n2=%d |
#         put o1=%g d1=1 o2=%g d2=1
#         ''' %(par['spikez']-n-1,2*n+1,
#               par['spikex']-m-1,2*m+1,
#               -n,-m))
#
#    Flow(xdel,[xdel+'-tmp'],
#         '''
#         window f1=%d n1=%d f2=%d n2=%d |
#         put o1=%g d1=1 o2=%g d2=1
#         ''' %(par['spikez']-m-1,2*m+1,
#               par['spikex']-n-1,2*n+1,
#               -n,-m))
## ------------------------------------------------------------








##############################################################
# 3D


#def separator3(zdel,xdel,ydel,inp,ccc,nu,stat,domain,n,m,l,par):
#
#    Flow([zdel,xdel,ydel],
#         [inp,ccc,nu,
#          'Code/EDERIV3.x'],
#         '''
#         ${SOURCES[3]} verb=y stat=%s domain=%s ompnth=8
#         ccc=${SOURCES[1]}
#          nu=${SOURCES[2]}
#         zdel=${TARGETS[0]}
#         xdel=${TARGETS[1]}
#         ydel=${TARGETS[2]}
#         ''' % (stat,domain))
#


def separator3PSS(pzdel,pxdel,pydel,
                  vzdel,vxdel,vydel,
                  hzdel,hxdel,hydel,
                  inp,ccc,nu,stat,domain,n,m,l,par):

    Flow([pzdel,pxdel,pydel,
          vzdel,vxdel,vydel,
          hzdel,hxdel,hydel],
         [inp,ccc,nu],
          '''
         ederiv3dfilters verb=y stat=%s domain=%s ompnth=8
         ccc=${SOURCES[1]}
          nu=${SOURCES[2]}
         pzdel=${TARGETS[0]}  pxdel=${TARGETS[1]}  pydel=${TARGETS[2]}
         vzdel=${TARGETS[3]}  vxdel=${TARGETS[4]}  vydel=${TARGETS[5]}
         hzdel=${TARGETS[6]}  hxdel=${TARGETS[7]}  hydel=${TARGETS[8]}
         ''' % (stat,domain))


def separator3TTI(zdel,xdel,ydel,inp,ccc,stat,domain,n,m,l,par):

    Flow([zdel,xdel,ydel],
         [inp,ccc],
         '''
         ederiv3d verb=y stat=%s domain=%s ompnth=8
         ccc=${SOURCES[1]}
         zdel=${TARGETS[0]}
         xdel=${TARGETS[1]}
         ydel=${TARGETS[2]}
         ''' % (stat,domain))


def SepK3(up,us,uz,ux,kz,kx,par):

    Flow(uz+'k',uz,'rtoc | fft3 axis=1 opt=n pad=1 | fft3 axis=2 opt=n pad=1| fft3 axis=3 opt=n pad=1')
    Flow(ux+'k',ux,'rtoc | fft3 axis=1 opt=n pad=1 | fft3 axis=2 opt=n pad=1| fft3 axis=3 opt=n pad=1')
#wavefield*separator in k domain
    Flow(up+'tmp1',[kz, uz+'k'],
         '''
     rtoc |
     math uz=${SOURCES[1]} output="input*sqrt(-1)*uz*(-1)" |
     fft3 axis=1 inv=y opt=n |
     fft3 axis=2 inv=y opt=n |
     fft3 axis=3 inv=y opt=n
     ''')
    Flow(up+'tmp2',[kx, ux+'k'],
         '''
     rtoc |
     math ux=${SOURCES[1]} output="input*sqrt(-1)*ux*(-1)" |
     fft3 axis=1 inv=y opt=n |
     fft3 axis=2 inv=y opt=n |
     fft3 axis=3 inv=y opt=n
     ''')
#sum all dimensions
    Flow(up,[up+'tmp1',up+'tmp2'],
         '''
     math x=${SOURCES[0]} y=${SOURCES[1]} output="x+y"|
     real |
     put d1=%(dz)f d2=%(dx)f d3=%(dy)f
     '''%par)





    Flow(us+'tmp1',[kz, ux+'k'],
         '''
     rtoc |
     math uz=${SOURCES[1]} output="input*sqrt(-1)*uz*(-1)" |
     fft3 axis=1 inv=y opt=n |
     fft3 axis=2 inv=y opt=n |
     fft3 axis=3 inv=y opt=n
     ''')
    Flow(us+'tmp2',[kx, uz+'k'],
         '''
     rtoc |
     math ux=${SOURCES[1]} output="input*sqrt(-1)*ux*(-1)" |
     fft3 axis=1 inv=y opt=n |
     fft3 axis=2 inv=y opt=n |
     fft3 axis=3 inv=y opt=n
     ''')

    Flow(us,[us+'tmp1',us+'tmp2'],
         '''
     math x=${SOURCES[0]} y=${SOURCES[1]} output="x-y" |
     real |
     put d1=%(dz)f d2=%(dx)f d3=%(dy)f
     '''%par)


# ------------------------------------------------------------
# 3d vti separation
def SepKP3(up,usv,ush,uz,ux,uy,kz,kx,ky,par):

    Flow(uz+'k',uz,'rtoc | fft3 axis=1 opt=n pad=1 | fft3 axis=2 opt=n pad=1| fft3 axis=3 opt=n pad=1')
    Flow(ux+'k',ux,'rtoc | fft3 axis=1 opt=n pad=1 | fft3 axis=2 opt=n pad=1| fft3 axis=3 opt=n pad=1')
    Flow(uy+'k',uy,'rtoc | fft3 axis=1 opt=n pad=1 | fft3 axis=2 opt=n pad=1| fft3 axis=3 opt=n pad=1')

    Flow(kx+'k',kx,'rtoc | math output="input*sqrt(-1)*(-1)"   ')
    Flow(ky+'k',ky,'rtoc | math output="input*sqrt(-1)*(-1)"  ')
    Flow(kz+'k',kz,'rtoc | math output="input*sqrt(-1)*(-1)"  ')
#wavefield*separator in k domain
    Flow(up,[ux+'k', uy+'k',uz+'k',kx+'k',ky+'k',kz+'k'],
         '''
     math  ux=${SOURCES[0]} uy=${SOURCES[1]} uz=${SOURCES[2]}
           kx=${SOURCES[3]} ky=${SOURCES[4]} kz=${SOURCES[5]}
           output="(kx*ux+ky*uy+kz*uz)" |
     fft3 axis=1 inv=y opt=n |
     fft3 axis=2 inv=y opt=n |
     fft3 axis=3 inv=y opt=n |
     real |  put d1=%(dz)f d2=%(dx)f d3=%(dy)f
     '''%par)

# ------------------------------------------------------------
#   us-z
    Flow(ush,[ux+'k',uy+'k',kx+'k',ky+'k'],
         '''
     math ux=${SOURCES[0]} uy=${SOURCES[1]} dx=${SOURCES[2]} dy=${SOURCES[3]}
          output="(uy*dx-ux*dy)  "  |
     fft3 axis=1 inv=y opt=n |
     fft3 axis=2 inv=y opt=n |
     fft3 axis=3 inv=y opt=n |
     real | put d1=%(dz)f d2=%(dx)f d3=%(dy)f
     '''%par)
    Flow(usv+'-x',[uy+'k',uz+'k',ky+'k',kz+'k',kx+'k'],
         '''
     math uy=${SOURCES[0]} uz=${SOURCES[1]} dy=${SOURCES[2]} dz=${SOURCES[3]} dx=${SOURCES[4]}
          output="(uz*dy-uy*dz)*dy  "  |
     fft3 axis=1 inv=y opt=n |
     fft3 axis=2 inv=y opt=n |
     fft3 axis=3 inv=y opt=n |
     real | put d1=%(dz)f d2=%(dx)f d3=%(dy)f
     '''%par)
    Flow(usv+'-y',[ux+'k',uz+'k',kx+'k',kz+'k',ky+'k'],
         '''
     math ux=${SOURCES[0]} uz=${SOURCES[1]} dx=${SOURCES[2]} dz=${SOURCES[3]} dy=${SOURCES[4]}
          output="(ux*dz-uz*dx)*dx "  |
     fft3 axis=1 inv=y opt=n |
     fft3 axis=2 inv=y opt=n |
     fft3 axis=3 inv=y opt=n |
     real | put d1=%(dz)f d2=%(dx)f d3=%(dy)f
     '''%par)


    Flow(usv,[usv+'-x',usv+'-y',kx,ky],
          '''math a=${SOURCES[0]} b=${SOURCES[1]}
              output="a-b" ''')

# ------------------------------------------------------------
# ------------------------------------------------------------
def SepKP3TTI(up,usv,ush,uz,ux,uy,kz,kx,ky,nu,alpha,par):

    kp=up+'k'
    kh=ush+'k'
    kv=usv+'k'

# ------------------------------------------------------------
##
 #  wavefield
 ##

#    Flow(ux+'k',ux,'rtoc | fft3 axis=1 opt=n pad=1 | fft3 axis=2 opt=n pad=1| fft3 axis=3 opt=n pad=1')
#    Flow(uy+'k',uy,'rtoc | fft3 axis=1 opt=n pad=1 | fft3 axis=2 opt=n pad=1| fft3 axis=3 opt=n pad=1')
#    Flow(uz+'k',uz,'rtoc | fft3 axis=1 opt=n pad=1 | fft3 axis=2 opt=n pad=1| fft3 axis=3 opt=n pad=1')

##
 #  P wave polarization
 ##


    Flow(kp+'xk',kx,'rtoc | math output="input*sqrt(-1)*(-1)"   ')
    Flow(kp+'yk',ky,'rtoc | math output="input*sqrt(-1)*(-1)"  ')
    Flow(kp+'zk',kz,'rtoc | math output="input*sqrt(-1)*(-1)"  ')


# ------------------------------------------------------------
##
 # SH polarization
 ##

    Flow(kh+'xk',[kx,ky,kz,nu,alpha],
         '''
     math  kx=${SOURCES[0]} ky=${SOURCES[1]} kz=${SOURCES[2]}
           nu=${SOURCES[3]} alpha=${SOURCES[4]}
           output="-cos(nu)*ky-sin(nu)*sin(alpha)*kz" |
           rtoc | math output="input*sqrt(-1)*(-1) "
        ''')
    Flow(kh+'yk',[kx,ky,kz,nu,alpha],
         '''
     math  kx=${SOURCES[0]} ky=${SOURCES[1]} kz=${SOURCES[2]}
           nu=${SOURCES[3]} alpha=${SOURCES[4]}
           output=" cos(nu)*kx+sin(nu)*cos(alpha)*kz" |
           rtoc | math output="input*sqrt(-1)*(-1) "
        ''')
    Flow(kh+'zk',[kx,ky,kz,nu,alpha],
         '''
     math  kx=${SOURCES[0]} ky=${SOURCES[1]} kz=${SOURCES[2]}
           nu=${SOURCES[3]} alpha=${SOURCES[4]}
           output=" -sin(nu)*cos(alpha)*ky+sin(nu)*sin(alpha)*kx" |
           rtoc | math output="input*sqrt(-1)*(-1) "
        ''')

##
 # SV polarization
 ## do not multiply with -i, it is squared and makes it -1
    Flow(kv+'xk',[kx,ky,kz,nu,alpha],
         '''
     math  kx=${SOURCES[0]} ky=${SOURCES[1]} kz=${SOURCES[2]}
           nu=${SOURCES[3]} alpha=${SOURCES[4]}
           output="-cos(nu)*kx*kz-sin(nu)*cos(alpha)*(ky*ky+kz*kz)+sin(nu)*sin(alpha)*kx*ky " |
           rtoc
        ''')
    Flow(kv+'yk',[kx,ky,kz,nu,alpha],
         '''
     math  kx=${SOURCES[0]} ky=${SOURCES[1]} kz=${SOURCES[2]}
           nu=${SOURCES[3]} alpha=${SOURCES[4]}
           output="-cos(nu)*ky*kz-sin(nu)*sin(alpha)*(kx*kx+kz*kz)+sin(nu)*cos(alpha)*kx*ky"|
           rtoc
        ''')
    Flow(kv+'zk',[kx,ky,kz,nu,alpha],
         '''
     math  kx=${SOURCES[0]} ky=${SOURCES[1]} kz=${SOURCES[2]}
           nu=${SOURCES[3]} alpha=${SOURCES[4]}
           output="  cos(nu)*(ky*ky+kx*kx)+sin(nu)*cos(alpha)*kx*kz+sin(nu)*sin(alpha)*ky*kz " |
           rtoc
        ''')

#wavefield*separator in k domain
# ------------------------------------------------------------

    Flow(up,[ux+'k', uy+'k',uz+'k',kp+'xk',kp+'yk',kp+'zk'],
         '''
     math  ux=${SOURCES[0]} uy=${SOURCES[1]} uz=${SOURCES[2]}
           kx=${SOURCES[3]} ky=${SOURCES[4]} kz=${SOURCES[5]}
           output="kx*ux+ky*uy+kz*uz" |
     fft3 axis=1 inv=y opt=n |
     fft3 axis=2 inv=y opt=n |
     fft3 axis=3 inv=y opt=n |
     real |  put d1=%(dz)f d2=%(dx)f d3=%(dy)f
     '''%par)

    Flow(usv,[ux+'k', uy+'k',uz+'k',kv+'xk',kv+'yk',kv+'zk'],
         '''
     math  ux=${SOURCES[0]} uy=${SOURCES[1]} uz=${SOURCES[2]}
           kx=${SOURCES[3]} ky=${SOURCES[4]} kz=${SOURCES[5]}
           output="kx*ux+ky*uy+kz*uz" |
     fft3 axis=1 inv=y opt=n |
     fft3 axis=2 inv=y opt=n |
     fft3 axis=3 inv=y opt=n |
     real |  put d1=%(dz)f d2=%(dx)f d3=%(dy)f
     '''%par)

    Flow(ush,[ux+'k', uy+'k',uz+'k',kh+'xk',kh+'yk',kh+'zk'],
         '''
     math  ux=${SOURCES[0]} uy=${SOURCES[1]} uz=${SOURCES[2]}
           kx=${SOURCES[3]} ky=${SOURCES[4]} kz=${SOURCES[5]}
           output="kx*ux+ky*uy+kz*uz" |
     fft3 axis=1 inv=y opt=n |
     fft3 axis=2 inv=y opt=n |
     fft3 axis=3 inv=y opt=n |
     real |  put d1=%(dz)f d2=%(dx)f d3=%(dy)f
     '''%par)

# ------------------------------------------------------------

def RevLeft(out,inp,par):
    par['halfnx']=par['nx']/2

    Flow(inp+'l',inp,'window n2=%(halfnx)d f2=0'%par)
    Flow(inp+'r',inp,'window n2=%(halfnx)d f2=%(halfnx)d'%par)
    Flow(out,[inp+'l', inp+'r'],'scale rscale=-1 | cat ${SOURCES[1]} axis=2 space=n')


def RevLeftC(out,inp,par):

    par['halfnx']=par['nx']/2
    par['halfny']=par['ny']/2

    par['r0']=par['nx']/2+1
    par['re']=par['nx']/2-1

    Flow(inp+'l',inp,'window n2=%(halfnx)d f2=0'%par)

    Flow(inp+'c-tmp',inp,'window n2=1 f2=%(halfnx)d'%par)
    Flow(inp+'c-tmp-b',inp+'c-tmp','window n2=%(halfny)d f2=0'%par)
    Flow(inp+'c-tmp-t',inp+'c-tmp','window n2=%(halfny)d f2=%(halfny)d'%par)
    Flow(inp+'c',[inp+'c-tmp-b', inp+'c-tmp-t'],'scale rscale=-1 | cat ${SOURCES[1]} axis=2 space=n| transp plane=23')

    Flow(inp+'r',inp,'window n2=%(re)d f2=%(r0)d '%par)

    Flow(out,[inp+'l',inp+'c',inp+'r'],'scale rscale=-1 | cat ${SOURCES[1:3]} axis=2 space=n')


# ------------------------------------------------------------
# 3d vti separation
def SepKS3(usx,usy,usz,uz,ux,uy,kz,kx,ky,par):

    Flow(uz+'k',uz,'rtoc | fft3 axis=1 opt=n pad=1 | fft3 axis=2 opt=n pad=1| fft3 axis=3 opt=n pad=1')
    Flow(ux+'k',ux,'rtoc | fft3 axis=1 opt=n pad=1 | fft3 axis=2 opt=n pad=1| fft3 axis=3 opt=n pad=1')
    Flow(uy+'k',uy,'rtoc | fft3 axis=1 opt=n pad=1 | fft3 axis=2 opt=n pad=1| fft3 axis=3 opt=n pad=1')

    Flow(kx+'k',kx,'rtoc | math output="input*sqrt(-1)*(-1)"   ')
    Flow(ky+'k',ky,'rtoc | math output="input*sqrt(-1)*(-1)"  ')
    Flow(kz+'k',kz,'rtoc | math output="input*sqrt(-1)*(-1)"  ')

# ------------------------------------------------------------
#   3 components of curl(wavefield)
    Flow(usz,[ux+'k',uy+'k',kx+'k',ky+'k'],
         '''
     math ux=${SOURCES[0]} uy=${SOURCES[1]} dx=${SOURCES[2]} dy=${SOURCES[3]}
          output="uy*dx-ux*dy  "  |
     fft3 axis=1 inv=y opt=n |
     fft3 axis=2 inv=y opt=n |
     fft3 axis=3 inv=y opt=n |
     real | put d1=%(dz)f d2=%(dx)f d3=%(dy)f
     '''%par)
    Flow(usx,[uy+'k',uz+'k',ky+'k',kz+'k'],
         '''
     math uy=${SOURCES[0]} uz=${SOURCES[1]} dy=${SOURCES[2]} dz=${SOURCES[3]}
          output="(uz*dy-uy*dz)  "  |
     fft3 axis=1 inv=y opt=n |
     fft3 axis=2 inv=y opt=n |
     fft3 axis=3 inv=y opt=n |
     real | put d1=%(dz)f d2=%(dx)f d3=%(dy)f
     '''%par)
    Flow(usy,[uz+'k',ux+'k',kz+'k',kx+'k'],
         '''
     math uz=${SOURCES[0]} ux=${SOURCES[1]} dz=${SOURCES[2]} dx=${SOURCES[3]}
          output="(ux*dz-uz*dx) "  |
     fft3 axis=1 inv=y opt=n |
     fft3 axis=2 inv=y opt=n |
     fft3 axis=3 inv=y opt=n |
     real | put d1=%(dz)f d2=%(dx)f d3=%(dy)f
     '''%par)





# 3d vti separation
def SepKSproj3(usx,usy,usz,uz,ux,uy,kz,kx,ky,par):

    Flow(uz+'k',uz,'rtoc | fft3 axis=1 opt=n pad=1 | fft3 axis=2 opt=n pad=1| fft3 axis=3 opt=n pad=1')
    Flow(ux+'k',ux,'rtoc | fft3 axis=1 opt=n pad=1 | fft3 axis=2 opt=n pad=1| fft3 axis=3 opt=n pad=1')
    Flow(uy+'k',uy,'rtoc | fft3 axis=1 opt=n pad=1 | fft3 axis=2 opt=n pad=1| fft3 axis=3 opt=n pad=1')

    Flow(kx+'k',kx,'rtoc | math output="input*sqrt(-1)*(-1)"   ')
    Flow(ky+'k',ky,'rtoc | math output="input*sqrt(-1)*(-1)"  ')
    Flow(kz+'k',kz,'rtoc | math output="input*sqrt(-1)*(-1)"  ')

# ------------------------------------------------------------
#   3 components of curl(wavefield)
    Flow(usx,[ux+'k',uy+'k',uz+'k',kx+'k',ky+'k',kz+'k'],
         '''
     math ux=${SOURCES[0]} uy=${SOURCES[1]} uz=${SOURCES[2]} dx=${SOURCES[3]} dy=${SOURCES[4]} dz=${SOURCES[5]}
          output="(dy*dy+dz*dz)*ux+(-dx*dy)*uy+(-dx*dz)*uz "  |
     fft3 axis=1 inv=y opt=n |
     fft3 axis=2 inv=y opt=n |
     fft3 axis=3 inv=y opt=n |
     real | put d1=%(dz)f d2=%(dx)f d3=%(dy)f
     '''%par)
    Flow(usy,[ux+'k',uy+'k',uz+'k',kx+'k',ky+'k',kz+'k'],
         '''
     math ux=${SOURCES[0]} uy=${SOURCES[1]} uz=${SOURCES[2]} dx=${SOURCES[3]} dy=${SOURCES[4]} dz=${SOURCES[5]}
          output="(-dx*dy)*ux+(dx*dx+dz*dz)*uy+(-dy*dz)*uz  "  |
     fft3 axis=1 inv=y opt=n |
     fft3 axis=2 inv=y opt=n |
     fft3 axis=3 inv=y opt=n |
     real | put d1=%(dz)f d2=%(dx)f d3=%(dy)f
     '''%par)
    Flow(usz,[ux+'k',uy+'k',uz+'k',kx+'k',ky+'k',kz+'k'],
         '''
     math ux=${SOURCES[0]} uy=${SOURCES[1]} uz=${SOURCES[2]} dx=${SOURCES[3]} dy=${SOURCES[4]} dz=${SOURCES[5]}
          output="(-dx*dz)*ux+(-dy*dz)*uy+(dx*dx+dy*dy)*uz  "  |
     fft3 axis=1 inv=y opt=n |
     fft3 axis=2 inv=y opt=n |
     fft3 axis=3 inv=y opt=n |
     real | put d1=%(dz)f d2=%(dx)f d3=%(dy)f
     '''%par)



##############################################################
#make s wave amplitudes more accurate
def filters(f1,f2,dz,par):

    Flow([f1+'tmp',f2+'tmp'],[dz, 'Code/weights.x'],
         '''
          Code/weights.x  theta0=30
              f1=${TARGETS[0]} f2=${TARGETS[1]}
         ''')
    Flow(f1,f1+'tmp','rtoc')
    Flow(f2,f2+'tmp','rtoc')
#f1 for sh
#f2 for s

def makeSops(px,py,pz,shx,shy,shz,svx,svy,svz,dx,dy,dz,f1,f2):

    p=1

    Flow(px,[dx,dy,dz],'math output="dx/(sqrt(dx*dx+dy*dy+dz*dz)+.01)" dx=${SOURCES[0]} dy=${SOURCES[1]} dz=${SOURCES[2]}|rtoc |math output="I*input"  ')
    Flow(py,[dx,dy,dz],'math output="dy/(sqrt(dx*dx+dy*dy+dz*dz)+.01)" dx=${SOURCES[0]} dy=${SOURCES[1]} dz=${SOURCES[2]}|rtoc |math output="I*input"  ')
    Flow(pz,[dx,dy,dz],'math output="dz/(sqrt(dx*dx+dy*dy+dz*dz)+.01)" dx=${SOURCES[0]} dy=${SOURCES[1]} dz=${SOURCES[2]}|rtoc |math output="I*input"  ')
# ------------------------------------------------------------

    Flow('tmpa1',[dx,dy,dz],'''math output="(dy*dy+dz*dz)/((dx*dx+dy*dy+dz*dz)+.001)"
                      dx=${SOURCES[0]} dy=${SOURCES[1]} dz=${SOURCES[2]} |
                      rtoc''')
    Flow('tmpb1',[dx, dy],'math dx=${SOURCES[0]}  dy=${SOURCES[1]} output="dy/(sqrt(dx*dx+dy*dy)+.01)" |rtoc|math output="I*input" ')
#    Flow(shx,'tmpa1 tmpb1','add ${SOURCES[1]}  scale=1,1 ')
    Flow(shx,'tmpa1 tmpb1 f1 f2','math s=${SOURCES[0]} sh=${SOURCES[1]} f1=${SOURCES[2]} f2=${SOURCES[3]} output="s*f2*%f+sh*f1*(2-%f)" '%(p,p))


    Flow('tmpa2',[dx,dy,dz],'''math output="-dx*dy/((dx*dx+dy*dy+dz*dz)+.001)"
                      dx=${SOURCES[0]} dy=${SOURCES[1]} dz=${SOURCES[2]} |
                      rtoc''')
    Flow('tmpb2',[dx, dy],'math dx=${SOURCES[0]}  dy=${SOURCES[1]} output="dx/(sqrt(dx*dx+dy*dy)+.01)" |rtoc|math output="-I*input" ')
#    Flow(shy,'tmpa2 tmpb2','add ${SOURCES[1]}  scale=1,1 ')
    Flow(shy,'tmpa2 tmpb2 f1 f2','math s=${SOURCES[0]} sh=${SOURCES[1]} f1=${SOURCES[2]} f2=${SOURCES[3]} output="s*f2*%f+sh*f1*(2-%f)" '%(p,p))

    Flow('tmpa3',[dx,dy,dz],'''math output="-dx*dz/(sqrt(dx*dx+dy*dy+dz*dz)+.001)"
                      dx=${SOURCES[0]}  dy=${SOURCES[1]} dz=${SOURCES[2]} |
                      rtoc''')
    Flow('tmpb3',dx,'rtoc | math output="0*sqrt(-1)" ')
#    Flow(shz,'tmpa3 tmpb3','add ${SOURCES[1]}  scale=1,1 ')
    Flow(shz,'tmpa3 tmpb3 f1 f2','math s=${SOURCES[0]} sh=${SOURCES[1]} f1=${SOURCES[2]} f2=${SOURCES[3]} output="s*f2*%f+sh*f1*(2-%f)" '%(p,p))

# ------------------------------------------------------------
# SV separators
    Flow('tmpa4','tmpa2','math output="input" ')
    Flow('tmpb4',[dx,dy,dz],'''math dx=${SOURCES[0]}  dy=${SOURCES[1]} dz=${SOURCES[2]}
              output="-dx*dz/(sqrt(dx*dx+dy*dy)+.01)/(sqrt(dx*dx+dy*dy+dz*dz)+.001)" |rtoc ''')
#    Flow(svx,'tmpa4 tmpb4','add ${SOURCES[1]}  scale=1,1 ')
    Flow(svx,'tmpa4 tmpb4 f1 f2','math s=${SOURCES[0]} sh=${SOURCES[1]} f1=${SOURCES[2]} f2=${SOURCES[3]} output="s*f2*%f+sh*f1*(2-%f)" '%(p,p))



    Flow('tmpa5',[dx,dy,dz],'''math output="(dx*dx+dz*dz)/((dx*dx+dy*dy+dz*dz)+.001)"
                      dx=${SOURCES[0]} dy=${SOURCES[1]} dz=${SOURCES[2]} |
                      rtoc''')
    Flow('tmpb5',[dx,dy,dz],'''math dx=${SOURCES[0]}  dy=${SOURCES[1]} dz=${SOURCES[2]}
                    output="-dy*dz/(sqrt(dx*dx+dy*dy)+.01)/(sqrt(dx*dx+dy*dy+dz*dz)+.001)" | rtoc ''')
#    Flow(svy,'tmpa5 tmpb5','add ${SOURCES[1]}  scale=1,1 ')
    Flow(svy,'tmpa5 tmpb5 f1 f2','math s=${SOURCES[0]} sh=${SOURCES[1]} f1=${SOURCES[2]} f2=${SOURCES[3]} output="s*f2*%f+sh*f1*(2-%f)" '%(p,p))




    Flow('tmpa6',[dx,dy,dz],'''math output="-dy*dz/((dx*dx+dy*dy+dz*dz)+.001)"
                      dx=${SOURCES[0]} dy=${SOURCES[1]} dz=${SOURCES[2]} |
                      rtoc''')
    Flow('tmpb6',[dx,dy,dz],'''math dx=${SOURCES[0]}  dy=${SOURCES[1]} dz=${SOURCES[2]}
                    output="sqrt(dx*dx+dy*dy)/(sqrt(dx*dx+dy*dy+dz*dz)+.001)" | rtoc ''')
#    Flow(svz,'tmpa6 tmpb6','add ${SOURCES[1]}  scale=1,1 ')
    Flow(svz,'tmpa6 tmpb6 f1 f2','math s=${SOURCES[0]} sh=${SOURCES[1]} f1=${SOURCES[2]} f2=${SOURCES[3]} output="s*f2*%f+sh*f1*(2-%f)" '%(p,p))

# ------------------------------------------------------------

def SepKSS(p,sh,sv,uz,ux,uy,pkx,pky,pkz,shkx,shky,shkz,svkx,svky,svkz,par):

    Flow(uz+'k',uz,'rtoc | fft3 axis=1 opt=n pad=1 | fft3 axis=2 opt=n pad=1| fft3 axis=3 opt=n pad=1')
    Flow(ux+'k',ux,'rtoc | fft3 axis=1 opt=n pad=1 | fft3 axis=2 opt=n pad=1| fft3 axis=3 opt=n pad=1')
    Flow(uy+'k',uy,'rtoc | fft3 axis=1 opt=n pad=1 | fft3 axis=2 opt=n pad=1| fft3 axis=3 opt=n pad=1')

    Flow(sh,[ux+'k',uy+'k',uz+'k',shkx,shky,shkz],
         '''
     math ux=${SOURCES[0]} uy=${SOURCES[1]} uz=${SOURCES[2]}
          dx=${SOURCES[3]} dy=${SOURCES[4]} dz=${SOURCES[5]}
          output="dx*ux+dy*uy+dz*uz "  |
     fft3 axis=1 inv=y opt=n |
     fft3 axis=2 inv=y opt=n |
     fft3 axis=3 inv=y opt=n |
     real | put d1=%(dz)f d2=%(dx)f d3=%(dy)f
     '''%par)

    Flow(sv,[ux+'k',uy+'k',uz+'k',svkx,svky,svkz],
         '''
     math ux=${SOURCES[0]} uy=${SOURCES[1]} uz=${SOURCES[2]}
          dx=${SOURCES[3]} dy=${SOURCES[4]} dz=${SOURCES[5]}
          output="dx*ux+dy*uy+dz*uz "  |
     fft3 axis=1 inv=y opt=n |
     fft3 axis=2 inv=y opt=n |
     fft3 axis=3 inv=y opt=n |
     real | put d1=%(dz)f d2=%(dx)f d3=%(dy)f
     '''%par)

    Flow(p,[ux+'k',uy+'k',uz+'k',pkx,pky,pkz],
         '''
     math ux=${SOURCES[0]} uy=${SOURCES[1]} uz=${SOURCES[2]}
          dx=${SOURCES[3]} dy=${SOURCES[4]} dz=${SOURCES[5]}
          output="dx*ux+dy*uy+dz*uz "  |
     fft3 axis=1 inv=y opt=n |
     fft3 axis=2 inv=y opt=n |
     fft3 axis=3 inv=y opt=n |
     real | put d1=%(dz)f d2=%(dx)f d3=%(dy)f
     '''%par)



def SeperatorK3(kz,kx,ky,nu,alpha,index,index2,par):

    kp='kp'+index+'_nu'+index2
    kh='kh'+index+'_nu'+index2
    kv='kv'+index+'_nu'+index2

    print(kp)
##
 #  P wave polarization
 ##


    Flow(kp+'xk',kx,'rtoc | math output="input*sqrt(-1)*(-1)"   ')
    Flow(kp+'yk',ky,'rtoc | math output="input*sqrt(-1)*(-1)"  ')
    Flow(kp+'zk',kz,'rtoc | math output="input*sqrt(-1)*(-1)"  ')


# ------------------------------------------------------------
##
 # SH polarization
 ##
# SH polarization is not influenced by anisotropy. It is decouple and polarized
# orthogonal to symmetry planes, therefore, can just nxk to find SH polarization

    Flow(kh+'xk',[kx,ky,kz,nu,alpha],
         '''
     math  kx=${SOURCES[0]} ky=${SOURCES[1]} kz=${SOURCES[2]}
           nu=${SOURCES[3]} alpha=${SOURCES[4]}
           output="-cos(nu)*ky-sin(nu)*sin(alpha)*kz" |
           rtoc | math output="input*sqrt(-1)*(-1) "
        ''')
    Flow(kh+'yk',[kx,ky,kz,nu,alpha],
         '''
     math  kx=${SOURCES[0]} ky=${SOURCES[1]} kz=${SOURCES[2]}
           nu=${SOURCES[3]} alpha=${SOURCES[4]}
           output=" cos(nu)*kx+sin(nu)*cos(alpha)*kz" |
           rtoc | math output="input*sqrt(-1)*(-1) "
        ''')
    Flow(kh+'zk',[kx,ky,kz,nu,alpha],
         '''
     math  kx=${SOURCES[0]} ky=${SOURCES[1]} kz=${SOURCES[2]}
           nu=${SOURCES[3]} alpha=${SOURCES[4]}
           output=" -sin(nu)*cos(alpha)*ky+sin(nu)*sin(alpha)*kx" |
           rtoc | math output="input*sqrt(-1)*(-1) "
        ''')
#
#    Flow(kh+'xk',[nu,alpha],
#         '''
#     math  nu=${SOURCES[0]} alpha=${SOURCES[1]}
#           output="-cos(nu)*(x3-120)-sin(nu)*sin(alpha)*(x1-120)" |
#           rtoc | math output="input*sqrt(-1)*(-1) "
#        ''')
#    Flow(kh+'yk',[nu,alpha],
#         '''
#           math  nu=${SOURCES[0]} alpha=${SOURCES[1]}
#           output=" cos(nu)*(x2-120)+sin(nu)*cos(alpha)*(x1-120)" |
#           rtoc | math output="input*sqrt(-1)*(-1) "
#        ''')
#    Flow(kh+'zk',[nu,alpha],
#         '''
#           math  nu=${SOURCES[0]} alpha=${SOURCES[1]}
#           output=" -sin(nu)*cos(alpha)*(x3-120)+sin(nu)*sin(alpha)*(x2-120)" |
#           rtoc | math output="input*sqrt(-1)*(-1) "
#        ''')



##
 # SV polarization
 ## do not multiply with -i, it is squared and makes it -1
    Flow(kv+'xk',[kx,ky,kz,nu,alpha],
         '''
     math  kx=${SOURCES[0]} ky=${SOURCES[1]} kz=${SOURCES[2]}
           nu=${SOURCES[3]} alpha=${SOURCES[4]}
           output="-cos(nu)*kx*kz-sin(nu)*cos(alpha)*(ky*ky+kz*kz)+sin(nu)*sin(alpha)*kx*ky " |
           rtoc
        ''')
    Flow(kv+'yk',[kx,ky,kz,nu,alpha],
         '''
     math  kx=${SOURCES[0]} ky=${SOURCES[1]} kz=${SOURCES[2]}
           nu=${SOURCES[3]} alpha=${SOURCES[4]}
           output="-cos(nu)*ky*kz-sin(nu)*sin(alpha)*(kx*kx+kz*kz)+sin(nu)*cos(alpha)*kx*ky"|
           rtoc
        ''')
    Flow(kv+'zk',[kx,ky,kz,nu,alpha],
         '''
     math  kx=${SOURCES[0]} ky=${SOURCES[1]} kz=${SOURCES[2]}
           nu=${SOURCES[3]} alpha=${SOURCES[4]}
           output="  cos(nu)*(ky*ky+kx*kx)+sin(nu)*cos(alpha)*kx*kz+sin(nu)*sin(alpha)*ky*kz " |
           rtoc
        ''')


def SeperateInK3(up,usv,ush,uz,ux,uy,kz,kx,ky,nu,alpha,index,index2,par):

#wavefield*separator in k domain
# ------------------------------------------------------------
    kp='kp'+index+'_nu'+index2
    kh='kh'+index+'_nu'+index2
    kv='kv'+index+'_nu'+index2
    Flow(up,[ux+'k', uy+'k',uz+'k',kp+'xk',kp+'yk',kp+'zk'],
         '''
     math  ux=${SOURCES[0]} uy=${SOURCES[1]} uz=${SOURCES[2]}
           kx=${SOURCES[3]} ky=${SOURCES[4]} kz=${SOURCES[5]}
           output="kx*ux+ky*uy+kz*uz" |
     fft3 axis=1 inv=y opt=n |
     fft3 axis=2 inv=y opt=n |
     fft3 axis=3 inv=y opt=n |
     real |  put d1=%(dz)f d2=%(dx)f d3=%(dy)f
     '''%par)

    Flow(usv,[ux+'k', uy+'k',uz+'k',kv+'xk',kv+'yk',kv+'zk'],
         '''
     math  ux=${SOURCES[0]} uy=${SOURCES[1]} uz=${SOURCES[2]}
           kx=${SOURCES[3]} ky=${SOURCES[4]} kz=${SOURCES[5]}
           output="kx*ux+ky*uy+kz*uz" |
     fft3 axis=1 inv=y opt=n |
     fft3 axis=2 inv=y opt=n |
     fft3 axis=3 inv=y opt=n |
     real |  put d1=%(dz)f d2=%(dx)f d3=%(dy)f
     '''%par)

    Flow(ush,[ux+'k', uy+'k',uz+'k',kh+'xk',kh+'yk',kh+'zk'],
         '''
     math  ux=${SOURCES[0]} uy=${SOURCES[1]} uz=${SOURCES[2]}
           kx=${SOURCES[3]} ky=${SOURCES[4]} kz=${SOURCES[5]}
           output="kx*ux+ky*uy+kz*uz" |
     fft3 axis=1 inv=y opt=n |
     fft3 axis=2 inv=y opt=n |
     fft3 axis=3 inv=y opt=n |
     real |  put d1=%(dz)f d2=%(dx)f d3=%(dy)f
     '''%par)

# ------------------------------------------------------------
