from rsf.proj import *

# ------------------------------------------------------------
def igrey(args,par):
    return '''
    grey labelrot=n %s
    label1="z" label2="x" title=" "
    min2=%g max2=%g min1=%g max1=%g
    ''' % (args,par['xmin'],par['xmax'],par['zmin'],par['zmax'])
def pgraph(args,par):
    return '''
    graph %s
    yreverse=y symbolsz=16 wantaxis=n title=" "
    min1=%g max1=%g min2=%g max2=%g
    ''' % (args,par['xmin'],par['xmax'],par['zmin'],par['zmax'])
def dgrey(args,par):
    return '''
    grey labelrot=n %s
    label1="t" label2="x" title=" "
    min2=%g max2=%g
    ''' % (args,par['xmin'],par['xmax'])

# ------------------------------------------------------------

def velo(velo,par):
    Plot(velo,velo,igrey('pclip=100 bias=900 allpos=y',par))

def coor(s_i,ss,rr,par):
    r_o = par['r_o']
    r_d = par['r_d']
    r_n = par['r_n']

    s_o = par['s_o']
    s_d = par['s_d']
    s_n = par['s_n']
    s_x = s_o + s_i * s_d

    _zs = '_zs' + str(s_i)
    _xs = '_xs' + str(s_i)
    _rs = '_rs' + str(s_i)
    Flow(_zs,None,'math n1=1 d1=0 o1=0 output=0')
    Flow(_xs,None,'math n1=1 d1=0 o1=0 output=%s' % s_x)
    Flow(_rs,None,'math n1=1 d1=0 o1=0 output=1')
    Flow(ss,[_xs,_zs,_rs],
         '''
         cat axis=2 space=n ${SOURCES[0]} ${SOURCES[1]} ${SOURCES[2]} |
         transp
         ''', stdin=0)

    _zr = '_zr' + str(s_i)
    _xr = '_xr' + str(s_i)    
    Flow(_zr,None,'math n1=%d d1=%g o1=%g output="00"' % (r_n,r_d,s_x) )
    Flow(_xr,None,'math n1=%d d1=%g o1=%g output="x1"' % (r_n,r_d,s_x) )
    Flow(rr,[_xr,_zr]     ,
         '''
         cat axis=2 space=n ${SOURCES[0]} ${SOURCES[1]}               |
         transp
         ''', stdin=0)

    Plot(ss,ss,'window n1=2 | dd type=complex | window | ' + pgraph('symbol=* plotcol=2',par))
    Plot(rr,rr,'window n1=2 | dd type=complex | window | ' + pgraph('symbol=. plotcol=1 symbolsz=.01',par))
    
def shot(s_i,data,wfld,wave,velo,ss,rr,par):

    s_o = par['s_o']
    s_d = par['s_d']
    s_n = par['s_n']
    s_x = s_o + s_i * s_d

    par['s_x']=s_x

    _dd = '_dd' + str(s_i)
    Flow([_dd,wfld],[wave,velo,ss,rr],
         '''
         afdm2d
         verb=%(verb)s
         abc=%(abc)s
         free=%(free)s
             snap=%(snap)s jsnap=%(jsnap)d
             nbz=%(nbz)d tz=%(tz)g
             nbx=%(nbx)d tx=%(tx)g
         vel=${SOURCES[1]}
         sou=${SOURCES[2]}
         rec=${SOURCES[3]}
         wfl=${TARGETS[1]}
         ''' % par)

    Flow(data,_dd,
         '''
         window min2=%(ft)g |
         pad n2out=%(nt)d  |
         put o1=0 |
         transp |
         mutter half=n t0=%(ft)g v0=1500 |
         put o3=%(s_x)g d3=%(s_d)g
         ''' % par )

    Plot(data,data,'put o2=%s | ' % s_x + dgrey('pclip=100',par))

