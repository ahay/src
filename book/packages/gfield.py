from rsfproj import *

def execute(rr,par):
    lpar = par.copy()
    lpar['nx']=2*par['nx']
    lpar['nz']=2*par['nz']
    lpar['ox']=-par['nx']*par['dx']
    lpar['oz']=-par['nz']*par['dz']
    
    # IID noise
    Flow(rr+'-n',None,
         '''
         math output="0"
         n1=%(nz)d d1=%(dz)g o1=%(oz)g label1="z" 
         n2=%(nx)d d2=%(dx)g o2=%(ox)g label2="x" |
         noise type=y |
         scale axis=123
         ''' % lpar)

    # ------------------------------------------------------------    
    # distance
    Flow(rr+'-lu',rr+'-n','math output="(x2*%g + x1*%g)/%g"'% (lpar['ux'],lpar['uz'],lpar['ru']) )
    Flow(rr+'-lv',rr+'-n','math output="(x2*%g + x1*%g)/%g"'% (lpar['vx'],lpar['vz'],lpar['rv']) )
    
    Flow(rr+'-l',[rr+'-lu',rr+'-lv'],
         'math lu=${SOURCES[0]} lv=${SOURCES[1]} output="sqrt(lu*lu+lv*lv)"', stdin=0)
    
    # ------------------------------------------------------------
    # covariance
    Flow(  rr+'-c',rr+'-l','math output="exp(-(input^%(a)g))"' % lpar)
    
    # ------------------------------------------------------------
    # FFT
    Flow(  rr+'-nf',rr+'-n','rtoc | fft3 opt=n axis=1 | fft3 opt=n axis=2')    
    Flow(  rr+'-cf',rr+'-c','rtoc | fft3 opt=n axis=1 | fft3 opt=n axis=2')    
    Flow(  rr+'-af',rr+'-cf','real | clip2 lower=0 | math output="sqrt(input)" | rtoc')
    
    # ------------------------------------------------------------
    # random field
    Flow(rr+'-k',[rr+'-af',rr+'-nf'],
         '''
         add mode=m ${SOURCES[1]} |
         fft3 opt=n axis=1 inv=y |
         fft3 opt=n axis=2 inv=y |
         real
         ''')

    # ------------------------------------------------------------
    Flow(rr,rr+'-k','window min1=0 min2=0 | scale axis=123')


