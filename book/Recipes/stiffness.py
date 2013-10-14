try:    from rsf.cluster import *
except: from rsf.proj    import *
from math import *

# ------------------------------------------------------------
# isotropic stiffness tensor
def iso2d(cc,vp,vs,ro,par):
    
    # Lame parameters
    # lambda = ro * (vp^2 - 2 vs^2)
    Flow(cc+'lambda',[ro,vp,vs],
         '''
         math output="ro*(vp*vp-2*vs*vs)"
         ro=${SOURCES[0]}
         vp=${SOURCES[1]}
         vs=${SOURCES[2]}     
         ''')
    
    #     mu = ro *           vs^2
    Flow(cc+'mu',[ro,vs],
         '''
         math output="ro*vs*vs"
         ro=${SOURCES[0]}
         vs=${SOURCES[1]}
         ''')
    
    # c11 = lambda + 2 mu
    # c33 = lambda + 2 mu
    # c55 = mu
    # c13 = lambda
    Flow(cc+'-11',[cc+'lambda',cc+'mu'],
         'math l=${SOURCES[0]} m=${SOURCES[1]} output="l+2*m"')
    Flow(cc+'-33',[cc+'lambda',cc+'mu'],
         'math l=${SOURCES[0]} m=${SOURCES[1]} output="l+2*m"')
    Flow(cc+'-55',[cc+'lambda',cc+'mu'],
         'math l=${SOURCES[0]} m=${SOURCES[1]} output="m"')
    Flow(cc+'-13',[cc+'lambda',cc+'mu'],
         'math l=${SOURCES[0]} m=${SOURCES[1]} output="l"')

    Flow(cc,[cc+'-11',
             cc+'-33',
             cc+'-55',
             cc+'-13'],
         'cat axis=3 space=n ${SOURCES[1:4]}')

# ------------------------------------------------------------
def iso3d(cc,vp,vs,ro,par):
    
    # Lame parameters
    # lambda = ro * (vp^2 - 2 vs^2)
    Flow(cc+'lambda',[ro,vp,vs],
         '''
         math output="ro*(vp*vp-2*vs*vs)"
         ro=${SOURCES[0]}
         vp=${SOURCES[1]}
         vs=${SOURCES[2]}     
         ''')
    
    #     mu = ro *           vs^2
    Flow(cc+'mu',[ro,vs],
         '''
         math output="ro*vs*vs"
         ro=${SOURCES[0]}
         vs=${SOURCES[1]}
         ''')
    
    # c11 = lambda + 2 mu
    # c22 = lambda + 2 mu
    # c33 = lambda + 2 mu
    Flow(cc+'-11',[cc+'lambda',cc+'mu'],
         'math l=${SOURCES[0]} m=${SOURCES[1]} output="l+2*m"')
    Flow(cc+'-22',cc+'-11','window')
    Flow(cc+'-33',cc+'-11','window')

    # c12 = lambda
    # c13 = lambda
    # c23 = lambda
    Flow(cc+'-12',[cc+'lambda',cc+'mu'],
         'math l=${SOURCES[0]} m=${SOURCES[1]} output="l"')
    Flow(cc+'-13',cc+'-12','window')
    Flow(cc+'-23',cc+'-12','window')

    # c44 = mu
    # c55 = mu
    # c66 = mu
    Flow(cc+'-44',[cc+'lambda',cc+'mu'],
         'math l=${SOURCES[0]} m=${SOURCES[1]} output="m"')
    Flow(cc+'-55',cc+'-44','window')
    Flow(cc+'-66',cc+'-44','window')

    Flow(cc,[cc+'-11',
             cc+'-22',
             cc+'-33',
             cc+'-44',
             cc+'-55',
             cc+'-66',
             cc+'-12',
             cc+'-13',
             cc+'-23'],
         'cat axis=4 space=n ${SOURCES[1:9]}')

# ------------------------------------------------------------
# VTI medium stiffness tensor
def vti2d(cc,vp,vs,ro,epsilon,delta,par):
    
    Flow(cc+'-33',[vp,ro],
         '''
         math output="ro*vp*vp"
         vp=${SOURCES[0]}
         ro=${SOURCES[1]}
         ''')    
    Flow(cc+'-55',[vs,ro],
         '''
         math output="ro*vs*vs"
         vs=${SOURCES[0]}
         ro=${SOURCES[1]}
         ''')
    Flow(cc+'-11',[cc+'-33',epsilon],
         '''
         math output="input*(1+2*eps2)"
         eps2=${SOURCES[1]}
         ''')
    Flow(cc+'-13',[cc+'-33',cc+'-55',delta],
         '''
         math output="sqrt(2*c33*(c33-c55)*del2+(c33-c55)*(c33-c55))-c55"
         c33=${SOURCES[0]}
         c55=${SOURCES[1]}
         del2=${SOURCES[2]}
         ''')
    Flow(cc,[cc+'-11',
             cc+'-33',
             cc+'-55',
             cc+'-13'],
         'cat axis=3 space=n ${SOURCES[1:4]}')

# ------------------------------------------------------------
def vti3d(cc,vp,vs,ro,epsilon,delta,gamma,par):

    Flow(cc+'-33',[ro,vp],
         '''
         math output="ro*(vp*vp)"
         ro=${SOURCES[0]}
         vp=${SOURCES[1]}
         ''')
    
    Flow(cc+'-55',[ro,vs],
         '''
         math output="ro*(vs*vs)"
         ro=${SOURCES[0]}
         vs=${SOURCES[1]}
         ''')

    Flow(cc+'-11',[cc+'-33',epsilon],
         '''
         math output="input*(1+2*eps2)"
         eps2=${SOURCES[1]}
         ''')
    
    Flow(cc+'-22',[cc+'-33',epsilon],
         '''
         math output="input*(1+2*eps1)"
         eps1=${SOURCES[1]}
         ''')

    Flow(cc+'-66',[cc+'-55',gamma],
         '''
         math output="input*(1+2*gam1)"
         gam1=${SOURCES[1]}
         ''')
    
    Flow(cc+'-44',[cc+'-66',gamma],
         '''
         math output="input/(1+2*gam2)"
         gam2=${SOURCES[1]}
         ''')

    Flow(cc+'-13',[cc+'-33',cc+'-55',delta],
         '''
         math output="sqrt(2*c33*(c33-c55)*del2+(c33-c55)*(c33-c55))-c55"
         c33=${SOURCES[0]}
         c55=${SOURCES[1]}
         del2=${SOURCES[2]}
         ''')
    Flow(cc+'-23',[cc+'-33',cc+'-44',delta],
         '''
         math output="sqrt(2*c33*(c33-c44)*del1+(c33-c44)*(c33-c44))-c44"
         c33=${SOURCES[0]}
         c44=${SOURCES[1]}
         del1=${SOURCES[2]}
         ''')    
    Flow(cc+'-12',[cc+'-11',cc+'-66'],
         '''
         math output="c11-2*c66"
         c11=${SOURCES[0]}
         c66=${SOURCES[1]}
         ''')
    
    Flow(cc,[cc+'-11',
             cc+'-22',
             cc+'-33',
             cc+'-44',
             cc+'-55',
             cc+'-66',
             cc+'-12',
             cc+'-13',
             cc+'-23'],
         'cat axis=4 space=n ${SOURCES[1:9]}')

# ------------------------------------------------------------
# TTI stiffness tensor
def tti2d(mm,vp,vs,ro,epsilon,delta,nu,par):
    
    Flow(nu+'-rad',nu,'math output="3.1415*input/180."')
    
    vti2d(mm+'-cc',vp,vs,ro,epsilon,delta,par)
    Flow(mm+'-cc11',mm+'-cc','window n3=1 f3=0')
    Flow(mm+'-cc33',mm+'-cc','window n3=1 f3=1')
    Flow(mm+'-cc55',mm+'-cc','window n3=1 f3=2')
    Flow(mm+'-cc13',mm+'-cc','window n3=1 f3=3')
    
    Flow(mm+'11',[mm+'-cc11',mm+'-cc13',mm+'-cc33',mm+'-cc55',nu+'-rad'],
         '''
         math output="c11*cos(nu)^4+2*(c13+2*c55)*cos(nu)^2*sin(nu)^2+c33*sin(nu)^4"
         c11=${SOURCES[0]}
         c13=${SOURCES[1]}
         c33=${SOURCES[2]}
         c55=${SOURCES[3]}
         nu=${SOURCES[4]}
         ''')
    Flow(mm+'33',[mm+'-cc11',mm+'-cc13',mm+'-cc33',mm+'-cc55',nu+'-rad'],
         '''
         math output="c33*cos(nu)^4+2*(c13+2*c55)*cos(nu)^2*sin(nu)^2+c11*sin(nu)^4"
         c11=${SOURCES[0]}
         c13=${SOURCES[1]}
         c33=${SOURCES[2]}
         c55=${SOURCES[3]}
         nu=${SOURCES[4]}
         ''')
    Flow(mm+'55',[mm+'-cc11',mm+'-cc13',mm+'-cc33',mm+'-cc55',nu+'-rad'],
         '''
         math output="(c11-2*c13+c33+4*c55-(c11-2*c13+c33-4*c55)*cos(4*nu))/8"
         c11=${SOURCES[0]}
         c13=${SOURCES[1]}
         c33=${SOURCES[2]}
         c55=${SOURCES[3]}
         nu=${SOURCES[4]}
         ''')
    Flow(mm+'13',[mm+'-cc11',mm+'-cc13',mm+'-cc33',mm+'-cc55',nu+'-rad'],
         '''
         math output="(c11+6*c13+c33-4*c55-(c11-2*c13+c33-4*c55)*cos(4*nu))/8"
         c11=${SOURCES[0]}
         c13=${SOURCES[1]}
         c33=${SOURCES[2]}
         c55=${SOURCES[3]}
         nu=${SOURCES[4]}
         ''')
    
    Flow(mm,[mm+'11',mm+'33',mm+'55',mm+'13'],
         'cat axis=3 space=n ${SOURCES[1:4]}')
    
def tti2dobsolete(mm,vp,vs,ro,epsilon,delta,nu,par):
    # Is this really 2D TTI?  Jia's code suggests that 6 stiffness coefs are 
    # required here...  Hence, this will be re-routed to her code, which is
    # compliant with the newer elastic wave equation requirements.
    useOLD = False
    if useOLD:
        Flow(nu+'-rad',nu,'math output="3.1415*input/180."')
        
        vti2d(mm+'-cc',vp,vs,ro,epsilon,delta,par)
        Flow(mm+'-cc11',mm+'-cc','window n3=1 f3=0')
        Flow(mm+'-cc33',mm+'-cc','window n3=1 f3=1')
        Flow(mm+'-cc55',mm+'-cc','window n3=1 f3=2')
        Flow(mm+'-cc13',mm+'-cc','window n3=1 f3=3')
            
        Flow(mm+'11',[mm+'-cc11',mm+'-cc13',mm+'-cc33',mm+'-cc55',nu+'-rad'],
             '''
             math output="c11*cos(nu)^4+2*(c13+2*c55)*cos(nu)^2*sin(nu)^2+c33*sin(nu)^4"
             c11=${SOURCES[0]}
             c13=${SOURCES[1]}
             c33=${SOURCES[2]}
             c55=${SOURCES[3]}
             nu=${SOURCES[4]}
             ''')
        Flow(mm+'33',[mm+'-cc11',mm+'-cc13',mm+'-cc33',mm+'-cc55',nu+'-rad'],
             '''
             math output="c33*cos(nu)^4+2*(c13+2*c55)*cos(nu)^2*sin(nu)^2+c11*sin(nu)^4"
             c11=${SOURCES[0]}
             c13=${SOURCES[1]}
             c33=${SOURCES[2]}
             c55=${SOURCES[3]}
             nu=${SOURCES[4]}
             ''')
        Flow(mm+'55',[mm+'-cc11',mm+'-cc13',mm+'-cc33',mm+'-cc55',nu+'-rad'],
             '''
             math output="(c11-2*c13+c33+4*c55-(c11-2*c13+c33-4*c55)*cos(4*nu))/8"
             c11=${SOURCES[0]}
             c13=${SOURCES[1]}
             c33=${SOURCES[2]}
             c55=${SOURCES[3]}
             nu=${SOURCES[4]}
             ''')
        Flow(mm+'13',[mm+'-cc11',mm+'-cc13',mm+'-cc33',mm+'-cc55',nu+'-rad'],
             '''
             math output="(c11+6*c13+c33-4*c55-(c11-2*c13+c33-4*c55)*cos(4*nu))/8"
             c11=${SOURCES[0]}
             c13=${SOURCES[1]}
             c33=${SOURCES[2]}
             c55=${SOURCES[3]}
             nu=${SOURCES[4]}
             ''')
        
        Flow(mm,[mm+'11',mm+'33',mm+'55',mm+'13'],
             'cat axis=3 space=n ${SOURCES[1:4]}')
    else:
        Flow(mm,[vp,vs,ro,epsilon,delta,nu],
            '''
            stiff dim=2
            vp=${SOURCES[0]}
            vs=${SOURCES[1]}
            ro=${SOURCES[2]}
            eps=${SOURCES[3]}
            del=${SOURCES[4]}
            nu=${SOURCES[5]}
            gam=${SOURCES[0]}
            alp=${SOURCES[0]}
            ''')
         
def tti3d(mm,vp,vs,ro,epsilon,delta,gamma,nu,alpha,par):
    # Updated per the new elastic wave equation specifications.
    # Requires sfstiffness.
    Flow(mm,[vp,vs,ro,epsilon,delta,gamma,nu,alpha],
        '''
        stiff dim=3
        vp=${SOURCES[0]}
        vs=${SOURCES[1]}
        ro=${SOURCES[2]}
        eps=${SOURCES[3]}
        del=${SOURCES[4]}
        gam=${SOURCES[5]}
        nu=${SOURCES[6]}
        alp=${SOURCES[7]}
        ''',stdin=0)
    
# ------------------------------------------------------------
def ort2d(cc,vp,vs,ro,eps2,del2,par):

    Flow(cc+'-33',[ro,vp],
         '''
         math output="ro*(vp*vp)"
         ro=${SOURCES[0]}
         vp=${SOURCES[1]}
         ''')
    
    Flow(cc+'-55',[ro,vs],
         '''
         math output="ro*(vs*vs)"
         ro=${SOURCES[0]}
         vs=${SOURCES[1]}
         ''')

    Flow(cc+'-11',[cc+'-33',eps2],
         '''
         math output="input*(1+2*eps2)"
         eps2=${SOURCES[1]}
         ''')

    Flow(cc+'-13',[cc+'-33',cc+'-55',del2],
         '''
         math output="sqrt(2*c33*(c33-c55)*del2+(c33-c55)*(c33-c55))-c55"
         c33=${SOURCES[0]}
         c55=${SOURCES[1]}
         del2=${SOURCES[2]}
         ''')

    Flow(cc,[cc+'-11',
             cc+'-33',
             cc+'-55',
             cc+'-13'],
         'cat axis=3 space=n ${SOURCES[1:4]}')

# ------------------------------------------------------------
def ort3d(cc,vp,vs,ro,eps1,eps2,del1,del2,del3,gam1,gam2,par):

    Flow(cc+'-33',[ro,vp],
         '''
         math output="ro*(vp*vp)"
         ro=${SOURCES[0]}
         vp=${SOURCES[1]}
         ''')
    
    Flow(cc+'-55',[ro,vs],
         '''
         math output="ro*(vs*vs)"
         ro=${SOURCES[0]}
         vs=${SOURCES[1]}
         ''')

    Flow(cc+'-11',[cc+'-33',eps2],
         '''
         math output="input*(1+2*eps2)"
         eps2=${SOURCES[1]}
         ''')
    
    Flow(cc+'-22',[cc+'-33',eps1],
         '''
         math output="input*(1+2*eps1)"
         eps1=${SOURCES[1]}
         ''')

    Flow(cc+'-66',[cc+'-55',gam1],
         '''
         math output="input*(1+2*gam1)"
         gam1=${SOURCES[1]}
         ''')
    
    Flow(cc+'-44',[cc+'-66',gam2],
         '''
         math output="input/(1+2*gam2)"
         gam2=${SOURCES[1]}
         ''')

    Flow(cc+'-13',[cc+'-33',cc+'-55',del2],
         '''
         math output="sqrt(2*c33*(c33-c55)*del2+(c33-c55)*(c33-c55))-c55"
         c33=${SOURCES[0]}
         c55=${SOURCES[1]}
         del2=${SOURCES[2]}
         ''')
    Flow(cc+'-23',[cc+'-33',cc+'-44',del1],
         '''
         math output="sqrt(2*c33*(c33-c44)*del1+(c33-c44)*(c33-c44))-c44"
         c33=${SOURCES[0]}
         c44=${SOURCES[1]}
         del1=${SOURCES[2]}
         ''')    
    Flow(cc+'-12',[cc+'-11',cc+'-66',del3],
         '''
         math output="sqrt(2*c11*(c11-c66)*del3+(c11-c66)*(c11-c66))-c66"
         c11=${SOURCES[0]}
         c66=${SOURCES[1]}
         del3=${SOURCES[2]}
         ''')
    
    Flow(cc,[cc+'-11',
             cc+'-22',
             cc+'-33',
             cc+'-44',
             cc+'-55',
             cc+'-66',
             cc+'-12',
             cc+'-13',
             cc+'-23'],
         'cat axis=4 space=n ${SOURCES[1:9]}')

# ------------------------------------------------------------
def cplot2d(cc,i1,i2,par):

    Flow(cc+'-nul',cc,'window n1=1 n2=1 n3=1 | math output=0')
    Flow(cc+'-c11',cc,'window n1=1 f1=%d n2=1 f2=%d n3=1 f3=0' %(i1,i2))
    Flow(cc+'-c33',cc,'window n1=1 f1=%d n2=1 f2=%d n3=1 f3=1' %(i1,i2))
    Flow(cc+'-c55',cc,'window n1=1 f1=%d n2=1 f2=%d n3=1 f3=2' %(i1,i2))
    Flow(cc+'-c13',cc,'window n1=1 f1=%d n2=1 f2=%d n3=1 f3=3' %(i1,i2))

    Flow(cc+'-row1',[cc+'-c11',cc+'-c13',cc+'-nul'],'cat axis=2 space=n ${SOURCES[1:3]}')
    Flow(cc+'-row2',[cc+'-c13',cc+'-c33',cc+'-nul'],'cat axis=2 space=n ${SOURCES[1:3]}')
    Flow(cc+'-row3',[cc+'-nul',cc+'-nul',cc+'-c55'],'cat axis=2 space=n ${SOURCES[1:3]}')

    Flow(cc+'-all',[cc+'-row1',cc+'-row2',cc+'-row3'],'cat axis=1 space=n ${SOURCES[1:3]}')

    Result(cc,cc+'-all','grey pclip=100 title="" wantaxis=n screenratio=1 allpos=y')

# ------------------------------------------------------------
def cplot3d(cc,i1,i2,i3,par):

    Flow(cc+'-nul',cc,'window n1=1 n2=1 n3=1 n4=1 | math output=0')

    Flow(cc+'-c11',cc,'window n1=1 f1=%d n2=1 f2=%d n3=1 f3=%d n4=1 f4=0' %(i1,i2,i3))
    Flow(cc+'-c22',cc,'window n1=1 f1=%d n2=1 f2=%d n3=1 f3=%d n4=1 f4=1' %(i1,i2,i3))
    Flow(cc+'-c33',cc,'window n1=1 f1=%d n2=1 f2=%d n3=1 f3=%d n4=1 f4=2' %(i1,i2,i3))
    Flow(cc+'-c44',cc,'window n1=1 f1=%d n2=1 f2=%d n3=1 f3=%d n4=1 f4=3' %(i1,i2,i3))
    Flow(cc+'-c55',cc,'window n1=1 f1=%d n2=1 f2=%d n3=1 f3=%d n4=1 f4=4' %(i1,i2,i3))
    Flow(cc+'-c66',cc,'window n1=1 f1=%d n2=1 f2=%d n3=1 f3=%d n4=1 f4=5' %(i1,i2,i3))
    Flow(cc+'-c12',cc,'window n1=1 f1=%d n2=1 f2=%d n3=1 f3=%d n4=1 f4=6' %(i1,i2,i3))
    Flow(cc+'-c13',cc,'window n1=1 f1=%d n2=1 f2=%d n3=1 f3=%d n4=1 f4=7' %(i1,i2,i3))
    Flow(cc+'-c23',cc,'window n1=1 f1=%d n2=1 f2=%d n3=1 f3=%d n4=1 f4=8' %(i1,i2,i3))

    Flow(cc+'-row1',[cc+'-c11',cc+'-c12',cc+'-c13',cc+'-nul',cc+'-nul',cc+'-nul'],'cat axis=2 space=n ${SOURCES[1:6]}')
    Flow(cc+'-row2',[cc+'-c12',cc+'-c22',cc+'-c23',cc+'-nul',cc+'-nul',cc+'-nul'],'cat axis=2 space=n ${SOURCES[1:6]}')
    Flow(cc+'-row3',[cc+'-c13',cc+'-c23',cc+'-c33',cc+'-nul',cc+'-nul',cc+'-nul'],'cat axis=2 space=n ${SOURCES[1:6]}')
    Flow(cc+'-row4',[cc+'-nul',cc+'-nul',cc+'-nul',cc+'-c44',cc+'-nul',cc+'-nul'],'cat axis=2 space=n ${SOURCES[1:6]}')
    Flow(cc+'-row5',[cc+'-nul',cc+'-nul',cc+'-nul',cc+'-nul',cc+'-c55',cc+'-nul'],'cat axis=2 space=n ${SOURCES[1:6]}')
    Flow(cc+'-row6',[cc+'-nul',cc+'-nul',cc+'-nul',cc+'-nul',cc+'-nul',cc+'-c66'],'cat axis=2 space=n ${SOURCES[1:6]}')

    Flow(cc+'-all',[cc+'-row1',cc+'-row2',cc+'-row3',cc+'-row4',cc+'-row5',cc+'-row6'],'cat axis=1 space=n ${SOURCES[1:6]}')

    Result(cc,cc+'-all','grey pclip=100 title="" wantaxis=n screenratio=1 allpos=y color=j')
    
def fcplot3d(cc,i1,i2,i3,par):
    ''' Full 21 coefficient plot for stiffness tensor '''

    Flow(cc+'-nul',cc,'window n1=1 n2=1 n3=1 n4=1 | math output=0')

    Flow(cc+'-c11',cc,'window n1=1 f1=%d n2=1 f2=%d n3=1 f3=%d n4=1 f4=0' %(i1,i2,i3))
    Flow(cc+'-c12',cc,'window n1=1 f1=%d n2=1 f2=%d n3=1 f3=%d n4=1 f4=1' %(i1,i2,i3))
    Flow(cc+'-c13',cc,'window n1=1 f1=%d n2=1 f2=%d n3=1 f3=%d n4=1 f4=2' %(i1,i2,i3))
    Flow(cc+'-c14',cc,'window n1=1 f1=%d n2=1 f2=%d n3=1 f3=%d n4=1 f4=3' %(i1,i2,i3))
    Flow(cc+'-c15',cc,'window n1=1 f1=%d n2=1 f2=%d n3=1 f3=%d n4=1 f4=4' %(i1,i2,i3))
    Flow(cc+'-c16',cc,'window n1=1 f1=%d n2=1 f2=%d n3=1 f3=%d n4=1 f4=5' %(i1,i2,i3))
    Flow(cc+'-c22',cc,'window n1=1 f1=%d n2=1 f2=%d n3=1 f3=%d n4=1 f4=6' %(i1,i2,i3))
    Flow(cc+'-c23',cc,'window n1=1 f1=%d n2=1 f2=%d n3=1 f3=%d n4=1 f4=7' %(i1,i2,i3))
    Flow(cc+'-c24',cc,'window n1=1 f1=%d n2=1 f2=%d n3=1 f3=%d n4=1 f4=8' %(i1,i2,i3))
    Flow(cc+'-c25',cc,'window n1=1 f1=%d n2=1 f2=%d n3=1 f3=%d n4=1 f4=9' %(i1,i2,i3))
    Flow(cc+'-c26',cc,'window n1=1 f1=%d n2=1 f2=%d n3=1 f3=%d n4=1 f4=10' %(i1,i2,i3))
    Flow(cc+'-c33',cc,'window n1=1 f1=%d n2=1 f2=%d n3=1 f3=%d n4=1 f4=11' %(i1,i2,i3))
    Flow(cc+'-c34',cc,'window n1=1 f1=%d n2=1 f2=%d n3=1 f3=%d n4=1 f4=12' %(i1,i2,i3))
    Flow(cc+'-c35',cc,'window n1=1 f1=%d n2=1 f2=%d n3=1 f3=%d n4=1 f4=13' %(i1,i2,i3))
    Flow(cc+'-c36',cc,'window n1=1 f1=%d n2=1 f2=%d n3=1 f3=%d n4=1 f4=14' %(i1,i2,i3))
    Flow(cc+'-c44',cc,'window n1=1 f1=%d n2=1 f2=%d n3=1 f3=%d n4=1 f4=15' %(i1,i2,i3))
    Flow(cc+'-c45',cc,'window n1=1 f1=%d n2=1 f2=%d n3=1 f3=%d n4=1 f4=16' %(i1,i2,i3))
    Flow(cc+'-c46',cc,'window n1=1 f1=%d n2=1 f2=%d n3=1 f3=%d n4=1 f4=17' %(i1,i2,i3))
    Flow(cc+'-c55',cc,'window n1=1 f1=%d n2=1 f2=%d n3=1 f3=%d n4=1 f4=18' %(i1,i2,i3))
    Flow(cc+'-c56',cc,'window n1=1 f1=%d n2=1 f2=%d n3=1 f3=%d n4=1 f4=19' %(i1,i2,i3))
    Flow(cc+'-c66',cc,'window n1=1 f1=%d n2=1 f2=%d n3=1 f3=%d n4=1 f4=20' %(i1,i2,i3))
    Flow(cc+'-row1',[cc+'-c11',cc+'-c12',cc+'-c13',cc+'-c14',cc+'-c15',cc+'-c16'],'cat axis=2 space=n ${SOURCES[1:6]}')
    Flow(cc+'-row2',[cc+'-c12',cc+'-c22',cc+'-c23',cc+'-c24',cc+'-c25',cc+'-c26'],'cat axis=2 space=n ${SOURCES[1:6]}')
    Flow(cc+'-row3',[cc+'-c13',cc+'-c23',cc+'-c33',cc+'-c34',cc+'-c35',cc+'-c36'],'cat axis=2 space=n ${SOURCES[1:6]}')
    Flow(cc+'-row4',[cc+'-c14',cc+'-c24',cc+'-c34',cc+'-c44',cc+'-c45',cc+'-c46'],'cat axis=2 space=n ${SOURCES[1:6]}')
    Flow(cc+'-row5',[cc+'-c15',cc+'-c25',cc+'-c35',cc+'-c45',cc+'-c55',cc+'-c56'],'cat axis=2 space=n ${SOURCES[1:6]}')
    Flow(cc+'-row6',[cc+'-c16',cc+'-c26',cc+'-c36',cc+'-c46',cc+'-c56',cc+'-c66'],'cat axis=2 space=n ${SOURCES[1:6]}')

    Flow(cc+'-all',[cc+'-row1',cc+'-row2',cc+'-row3',cc+'-row4',cc+'-row5',cc+'-row6'],'cat axis=1 space=n ${SOURCES[1:6]}')

    Result(cc,cc+'-all','grey pclip=100 title="" wantaxis=n screenratio=1 allpos=y color=j')    
