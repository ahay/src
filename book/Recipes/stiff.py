try:
    from rsf.cluster import *
except:
    from rsf.proj import *
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

# ------------------------------------------------------------
def orth3d(cc,vp,vs,ro,epsilon1,epsilon2,delta1,delta2,delta3,gamma1,gamma2,par):
      
    c33=ro*vp*vp
    c55=ro*vs*vs

    c66=(1+2*gamma1)*c55
    c44=c66/(1+2*gamma2)

    c11=(1+2*epsilon2)*c33
    c22=(1+2*epsilon1)*c33

    c12=sqrt(2*c11*(c11-c66)*delta3+(c11-c66)*(c11-c66))-c66    
    c13=sqrt(2*c33*(c33-c55)*delta2+(c33-c55)*(c33-c55))-c55
    c23=sqrt(2*c33*(c33-c44)*delta1+(c33-c44)*(c33-c44))-c44
    
#    Flow('zero',None,
#         '''
#          spike nsp=1 mag=0.0
#          n1=%(nx)d o1=%(ox)g d1=%(dx)g 
#          n2=%(ny)d o2=%(oy)g d2=%(dy)g 
#          n3=%(nz)d o3=%(oz)g d3=%(dz)g |
#          put label1=%(lx)s label2=%(ly)s label3=%(lz)s 
#              unit1=%(uz)s unit2=%(ux)s  unit3=%(uz)s
#         ''' % par)

    Flow(cc+'-33','zero3d','add add=%f'%c33)
    Flow(cc+'-55','zero3d','add add=%f'%c55)
    Flow(cc+'-66','zero3d','add add=%f'%c66)
    Flow(cc+'-44','zero3d','add add=%f'%c44)
    Flow(cc+'-11','zero3d','add add=%f'%c11)
    Flow(cc+'-13','zero3d','add add=%f'%c13)
    Flow(cc+'-22','zero3d','add add=%f'%c22)
    Flow(cc+'-23','zero3d','add add=%f'%c23)
    Flow(cc+'-12','zero3d','add add=%f'%c12)

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

    Flow(cc+'-66',[cc+'-55',gam1],
         '''
         math output="input*(1+2*gam1)"
         gam1=${SOURCES[1]}
         ''')

    Flow(cc+'-44',[cc+'-66',gam2],
         '''
         math output="input*(1+2*gam2)"
         gam2=${SOURCES[1]}
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

    Flow(cc+'-12',[cc+'-11',cc+'-66',del3],
         '''
         math output="sqrt(2*c11*(c11-c66)*del3+(c11-c66)*(c11-c66))-c66"
         c11=${SOURCES[0]}
         c66=${SOURCES[1]}
         del3=${SOURCES[2]}
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
# ------------------------------------------------------------
# VTI medium stiffness tensor
def vti2d_point(vp,vs,ro,epsilon,delta):
   
    c33=ro*vp*vp
    c55=ro*vs*vs
    c11=2*epsilon*c33+c33
    c13=sqrt(2*c33*(c33-c55)*delta+(c33-c55)*(c33-c55))-c55
    c15=0
    c35=0
    cc=[c11,c13,c15,c33,c35,c55]   
    return cc
# ------------------------------------------------------------
# ------------------------------------------------------------
# TTI stiffness tensor
def tti2d_point(vp,vs,ro,epsilon,delta,nu):

    nu=3.1415*nu/180.
    cc=vti2d_point(vp,vs,ro,epsilon,delta)
    c11=cc[0]
    c13=cc[1]
    c15=cc[2]
    c33=cc[3]
    c35=cc[4]
    c55=cc[5]
        

    
    m11=c11*pow(cos(nu),4) + 2*(c13+2*c55)*pow(cos(nu),2)*pow(sin(nu),2) + c33*pow(sin(nu),4)
    m13=(c11+6*c13+c33-4*c55-(c11-2*c13+c33-4*c55)*cos(4*nu))/8
    m15=(c11-c33+(c11-2*c13+c33-4*c55)*cos(2*nu))*sin(2*nu)/4
    m33=c33*pow(cos(nu),4) + 2*(c13+2*c55)*pow(cos(nu),2)*pow(sin(nu),2) + c11*pow(sin(nu),4)
    m35=-(-c11+c33+(c11-2*c13+c33-4*c55)*cos(2*nu))*sin(2*nu)/4
    m55=(c11-2*c13+c33+4*c55-(c11-2*c13+c33-4*c55)*cos(4*nu))/8
    mm=[m11,m13,m15,m33,m35,m55]
    
    return mm

## 
 # 
# ------------------------------------------------------------
def tti3d_point(vp,vs,ro,epsilon,delta,gamma,nu,alpha):
    
   
#    2d vti coefficients	    
    c33=vp *vp *ro ;
    c55=vs *vs *ro ;
    c11=2*epsilon *c33+c33;
    c13=sqrt(2*c33*(c33-c55)*delta +(c33-c55)*(c33-c55))-c55;
    c66=2*gamma *c55+c55;
    c12=c11-2*c66;

    
    def COSNU4(n): return pow(cos(n),4)
    def COSNU3(n): return pow(cos(n),3)
    def COSNU2(n): return pow(cos(n),2)
    def SINNU4(n): return pow(sin(n),4)
    def SINNU3(n): return pow(sin(n),3)
    def SINNU2(n): return pow(sin(n),2)
    def COS4NU(n): return cos(4*n)
    def COS2NU(n): return cos(2*n)
    def COSNU(n):  return cos(n)
    def SIN4NU(n): return sin(4*n)
    def SIN2NU(n): return sin(2*n)
    def SINNU(n):  return sin(n)    
    nu=3.1415*nu/180.
    
#    2d tti coefficients	    
    m11=c11*COSNU4(nu)+2*(c13+2*c55)*COSNU2(nu)*SINNU2(nu)+c33*SINNU4(nu);
    m12=c12*COSNU2(nu)+c13*SINNU2(nu);
    m13=(c11+6*c13+c33-4*c55-(c11-2*c13+c33-4*c55)*COS4NU(nu))/8;
    m15=(c11-c33+(c11-2*c13+c33-4*c55)*COS2NU(nu))*SIN2NU(nu)/4;
    m22=c11;
    m23=c13*COSNU2(nu)+c12*SINNU2(nu);
    m25=(c12-c13)*COSNU(nu)*SINNU(nu);
    m33=c33*COSNU4(nu)+2*(c13+2*c55)*COSNU2(nu)*SINNU2(nu)+c11*SINNU4(nu);
    m35=-(-c11+c33+(c11-2*c13+c33-4*c55)*COS2NU(nu))*SIN2NU(nu)/4;
    m44=c55*COSNU2(nu)+c66*SINNU2(nu);
    m46=(c66-c55)*COSNU(nu)*SINNU(nu);
    m55=(c11-2*c13+c33+4*c55-(c11-2*c13+c33-4*c55)*COS4NU(nu))/8;
    m66=c66*COSNU2(nu)+c55*SINNU2(nu);
# ------------------------------------------------------------
    nu=3.1415*alpha/180.

    c11=m11*COSNU4(nu)+2*(m12+2*m66)*COSNU2(nu)*SINNU2(nu)+m22*SINNU4(nu);
    c12=(m11+6*m12+m22-4*m66-(m11-2*m12+m22-4*m66)*COSNU4(nu))/8;
    c13=m13*COSNU2(nu)+m23*SINNU2(nu);
    c14=(m15-2*m46)*COSNU2(nu)*SINNU(nu)+m25*SINNU3(nu);
    c15=m15*COSNU3(nu)+(m25+2*m46)*COSNU(nu)*SINNU2(nu);
    c16=(m11-m22-(m11-2*m12+m22-4*m66)*COS2NU(nu))*SIN2NU(nu)/4;
    c22=m22*COSNU4(nu)+2*(m12+2*m66)*COSNU2(nu)*SINNU2(nu)+m11*SINNU4(nu);
    c23=m23*COSNU2(nu)+m13*SINNU2(nu);
    c24=(m25+2*m46)*COSNU2(nu)*SINNU(nu)+m15*SINNU3(nu);
    c25=m25*COSNU3(nu)+(m15-2*m46)*COSNU(nu)*SINNU2(nu);
    c26=(m11-m22-(m11-2*m12+m22-4*m66)*COS2NU(nu))*SIN2NU(nu)/4;
    c33=m33;
    c34=m35*SINNU(nu);
    c35=m35*COSNU(nu);
    c36=(m13-m23)*COSNU(nu)*SINNU(nu);
    c44=m44*COSNU2(nu)+m55*SINNU2(nu);
    c45=(-m44+m55)*COSNU(nu)*SINNU(nu);
    c46=m46*COSNU3(nu)+(m15-m25-m46)*COSNU(nu)*SINNU2(nu);
    c55=m55*COSNU2(nu)+m44*SINNU2(nu);
    c56=m46*SINNU3(nu)+(m15-m25-m46)*COSNU2(nu)*SINNU(nu);
    c66=(m11-2*m12+m22+4*m66-(m11-2*m12+m22-4*m66)*COS4NU(nu))/8;


    cc=[c11, c12, c13, c14, c15, c16, c22, c23, c24, c25, c26, c33, c34, c35, c36, c44, c45, c46, c55, c56, c66]
    
    return cc
