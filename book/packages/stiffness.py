from rsfproj import *
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
    # c13 = lambda
    # c44 = mu
    Flow(cc+'-11',[cc+'lambda',cc+'mu'],
         'math l=${SOURCES[0]} m=${SOURCES[1]} output="l+2*m"')
    Flow(cc+'-13',[cc+'lambda',cc+'mu'],
         'math l=${SOURCES[0]} m=${SOURCES[1]} output="l"')
    Flow(cc+'-33',[cc+'lambda',cc+'mu'],
         'math l=${SOURCES[0]} m=${SOURCES[1]} output="l+2*m"')
    Flow(cc+'-55',[cc+'lambda',cc+'mu'],
         'math l=${SOURCES[0]} m=${SOURCES[1]} output="m"')

    Flow(cc+'-15',cc+'-11','math output=0')
    Flow(cc+'-35',cc+'-11','math output=0')

    Flow(cc,[cc+'-11',cc+'-13',cc+'-15',cc+'-33',cc+'-35',cc+'-55'],
         'cat axis=3 space=n ${SOURCES[1:6]}')
    
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
         math output="2*epsilon*c33+c33"
         c33=${SOURCES[0]}
         epsilon=${SOURCES[1]}
         ''')
    Flow(cc+'-13',[cc+'-33',cc+'-55',delta],
         '''
         math output="sqrt(2*c33*(c33-c55)*delta+(c33-c55)*(c33-c55))-c55"
         c33=${SOURCES[0]}
         c55=${SOURCES[1]}
         delta=${SOURCES[2]}
         ''')

    Flow(cc+'-15',cc+'-11','math output=0')
    Flow(cc+'-35',cc+'-11','math output=0')
    
    Flow(cc,[cc+'-11',cc+'-13',cc+'-15',cc+'-33',cc+'-35',cc+'-55'],
         'cat axis=3 space=n ${SOURCES[1:6]}')

# ------------------------------------------------------------
# TTI stiffness tensor
def tti2d(mm,vp,vs,ro,epsilon,delta,nu,par):

    Flow(nu+'-rad',nu,'math output="3.1415*input/180."')
    
    vti2d(mm+'-cc',vp,vs,ro,epsilon,delta,par)
    Flow(mm+'-cc11',mm+'-cc','window n3=1 f3=0')
    Flow(mm+'-cc13',mm+'-cc','window n3=1 f3=1')
    Flow(mm+'-cc33',mm+'-cc','window n3=1 f3=3')
    Flow(mm+'-cc55',mm+'-cc','window n3=1 f3=5')
        
    Flow(mm+'11',[mm+'-cc11',mm+'-cc13',mm+'-cc33',mm+'-cc55',nu+'-rad'],
         '''
         math output="c11*cos(nu)^4+2*(c13+2*c55)*cos(nu)^2*sin(nu)^2+c33*sin(nu)^4"
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
    Flow(mm+'15',[mm+'-cc11',mm+'-cc13',mm+'-cc33',mm+'-cc55',nu+'-rad'],
         '''
         math output="(c11-c33+(c11-2*c13+c33-4*c55)*cos(2*nu))*sin(2*nu)/4"
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
    Flow(mm+'35',[mm+'-cc11',mm+'-cc13',mm+'-cc33',mm+'-cc55',nu+'-rad'],
         '''
         math output="-(-c11+c33+(c11-2*c13+c33-4*c55)*cos(2*nu))*sin(2*nu)/4"
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
    
    Flow(mm,[mm+'11',mm+'13',mm+'15',mm+'33',mm+'35',mm+'55'],
         'cat axis=3 space=n ${SOURCES[1:6]}')
