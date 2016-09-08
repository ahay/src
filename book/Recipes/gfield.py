try:    from rsf.cluster import * 
except: from rsf.proj    import *
import math

def execute(rr,seed,par):
    lpar = par.copy()

    if(not lpar.has_key('ff')): lpar['ff']=0
    if(not lpar.has_key('ru')): lpar['ru']=5
    if(not lpar.has_key('rv')): lpar['rv']=1
    if(not lpar.has_key('aa')): lpar['aa']=1
    
    # directions 
    lpar['ux'] = math.cos(math.pi*par['ff']/180.)
    lpar['uz'] = math.sin(math.pi*par['ff']/180.)
    # cannot use math.radians in old Python
    
    lpar['vx'] =  lpar['uz']
    lpar['vz'] = -lpar['ux']

    # double and center the grid
    lpar['nx']=2*(par['nx'] + abs(par['ox'])/par['dx'])
    lpar['nz']=2*(par['nz'] + abs(par['oz'])/par['dz'])
    lpar['ox']=-lpar['nx']*par['dx']/2
    lpar['oz']=-lpar['nz']*par['dz']/2

    lpar['seed']=seed
    
    # IID noise
    Flow(rr+'-n',None,
         '''
         math output="0"
         n1=%(nz)d d1=%(dz)g o1=%(oz)g label1="z" 
         n2=%(nx)d d2=%(dx)g o2=%(ox)g label2="x" |
         noise seed=%(seed)d type=y |
         scale axis=123
         ''' % lpar)
    
    # ------------------------------------------------------------    
    # distance
    Flow(rr+'-lu',rr+'-n','math output="(x2*(%g) + x1*(%g))/(%g)"'% (lpar['ux'],lpar['uz'],lpar['ru']) )
    Flow(rr+'-lv',rr+'-n','math output="(x2*(%g) + x1*(%g))/(%g)"'% (lpar['vx'],lpar['vz'],lpar['rv']) )
    
    Flow(rr+'-l',[rr+'-lu',rr+'-lv'],
         'math lu=${SOURCES[0]} lv=${SOURCES[1]} output="sqrt(lu*lu+lv*lv)"', stdin=0)
    
    # ------------------------------------------------------------
    # covariance
    Flow(  rr+'-c',rr+'-l','math output="exp(-( input^(%(aa)g)) )"' % lpar)
    Flow(  rr+'-d',rr+'-c','rotate rot1=%g rot2=%g' % (par['nz']/2,par['nx']/2))
    
    # ------------------------------------------------------------
    # FFT
    Flow(  rr+'-nf',rr+'-n','rtoc | fft3 opt=y axis=1 | fft3 opt=y axis=2')    
    Flow(  rr+'-cf',rr+'-d','rtoc | fft3 opt=y axis=1 | fft3 opt=y axis=2')    
    Flow(  rr+'-af',rr+'-cf','real | clip2 lower=0 | math output="sqrt(input)" | rtoc')
    
    # ------------------------------------------------------------
    # random field
    Flow(rr+'-k',[rr+'-af',rr+'-nf'],
         '''
         add mode=m ${SOURCES[1]} |
         fft3 opt=y axis=2 inv=y |
         fft3 opt=y axis=1 inv=y |
         real
         ''')

    # ------------------------------------------------------------
    Flow(rr,rr+'-k','window min1=%(oz)g min2=%(ox)g | scale axis=123' % par)

# ------------------------------------------------------------
def run2d(R,par,seed=1.0,aa=0.0,ra=1.0,rb=1.0,alpha=1.0):

    # IID noise
    Flow(R+'-n',None,
         '''
         math output="0"
         n1=%(nz)d d1=%(dz)g o1=%(oz)g label1=%(lz)s unit1=%(uz)s
         n2=%(nx)d d2=%(dx)g o2=%(ox)g label2=%(lx)s unit2=%(ux)s |
         '''%par +
         '''
         noise seed=%g type=y | scale axis=123
         '''%seed)

    # ------------------------------------------------------------    
    # distance
    ax = math.cos(math.pi*(   aa)/180.)
    az = math.sin(math.pi*(   aa)/180.)
    
    bx = math.cos(math.pi*(-90+aa)/180.)
    bz = math.sin(math.pi*(-90+aa)/180.)

#    print ax,az
#    print bx,bz
    
    Flow(R+'-la',R+'-n','math output="( (x2*(%g) + x1*(%g) )/(%g) )^2"'% (ax,az,ra) )
    Flow(R+'-lb',R+'-n','math output="( (x2*(%g) + x1*(%g) )/(%g) )^2"'% (bx,bz,rb) )
    Flow(R+'-l',[R+'-la',R+'-lb'],'add ${SOURCES[1]} | math output="sqrt(input)"')

    # ------------------------------------------------------------
    # covariance
    Flow(  R+'-c',R+'-l','math output="exp(-( input^(%g)) )"'%alpha)

    # ------------------------------------------------------------
    # FFT
    Flow(  R+'-nf',R+'-n','rtoc | fft3 axis=1 | fft3 axis=2')    
    Flow(  R+'-cf',R+'-c','rtoc | fft3 axis=1 | fft3 axis=2')    
    Flow(  R+'-kf',R+'-cf','real | clip2 lower=0 | math output="sqrt(input)" | rtoc')
    
    # ------------------------------------------------------------
    # random field
    Flow(R,[R+'-kf',R+'-nf'],
         '''
         add mode=m ${SOURCES[1]} |
         fft3 axis=2 inv=y |
         fft3 axis=1 inv=y |
         real
         ''')

# ------------------------------------------------------------
def run3d(R,par,seed=1.0,aa=0.0,bb=0.0,ra=1.0,rb=1.0,rc=1.0,alpha=1.0):

    # IID noise
    Flow(R+'-n',None,
         '''
         math output="0"
         n1=%(nz)d d1=%(dz)g o1=%(oz)g label1=%(lz)s unit1=%(uz)s
         n2=%(nx)d d2=%(dx)g o2=%(ox)g label2=%(lx)s unit2=%(ux)s
         n3=%(ny)d d3=%(dy)g o3=%(oy)g label3=%(ly)s unit3=%(uy)s |
         '''%par +
         '''
         noise seed=%g type=y | scale axis=123
         '''%seed)

    # ------------------------------------------------------------    
    # distance
    ax = math.cos(math.pi*(   aa)/180.) * math.cos(math.pi*(   bb)/180.)
    ay = math.cos(math.pi*(   aa)/180.) * math.sin(math.pi*(   bb)/180.)
    az = math.sin(math.pi*(   aa)/180.)

    bx = math.cos(math.pi*(-90+aa)/180.) * math.cos(math.pi*(   bb)/180.)
    by = math.cos(math.pi*(-90+aa)/180.) * math.sin(math.pi*(   bb)/180.)
    bz = math.sin(math.pi*(-90+aa)/180.)

    cx = ay*bz - by*az;
    cy = az*bx - bz*ax;
    cz = ax*by - bx*ay;

    Flow(R+'-la',R+'-n','math output="( (x2*(%g) + x3*(%g) + x1*(%g) )/(%g) )^2"'% (ax,ay,az,rb) )    
    Flow(R+'-lb',R+'-n','math output="( (x2*(%g) + x3*(%g) + x1*(%g) )/(%g) )^2"'% (bx,by,bz,ra) )
    Flow(R+'-lc',R+'-n','math output="( (x2*(%g) + x3*(%g) + x1*(%g) )/(%g) )^2"'% (cx,cy,cz,rc) )
    Flow(R+'-l',[R+'-la',R+'-lb',R+'-lc'],'add ${SOURCES[1]} ${SOURCES[2]} | math output="sqrt(input)"')

    # ------------------------------------------------------------
    # covariance
    Flow(  R+'-c',R+'-l','math output="exp(-( input^(%g)) )"'%alpha)

    # ------------------------------------------------------------
    # FFT
    Flow(  R+'-nf',R+'-n','rtoc | fft3 opt=y axis=1 | fft3 opt=y axis=2 | fft3 opt=y axis=3')    
    Flow(  R+'-cf',R+'-c','rtoc | fft3 opt=y axis=1 | fft3 opt=y axis=2 | fft3 opt=y axis=3')    
    Flow(  R+'-kf',R+'-cf','real | clip2 lower=0 | math output="sqrt(input)" | rtoc')
    
    # ------------------------------------------------------------
    # random field
    Flow(R,[R+'-kf',R+'-nf'],
         '''
         add mode=m ${SOURCES[1]} |
         fft3 opt=y axis=3 inv=y |
         fft3 opt=y axis=2 inv=y |
         fft3 opt=y axis=1 inv=y |
         real
         ''')

