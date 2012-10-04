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

