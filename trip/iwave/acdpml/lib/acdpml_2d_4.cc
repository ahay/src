#include "cstd.h"

void acdpml_2d_4(float ** uc,    // current field
	   float ** up,        // previous field
	   float ** csq,       // csq
	   float ** phi1,      // phi1
	   float ** phi0,      // phi0
           float * dp1,        // damping profile zeta_x
           float * dp0,        // damping profile zeta_x
           float *di,
           float dt,
           int * s,            // start index
           int * e,            // end index
           float c0,
           float * c1,
           float * c2,
           int * lbc,
           int * rbc
	   ) {
    int i1, i0;
    //fprintf(stderr, "dp0[%d] = %f\n",e[0], dp0[e[0]]);
    // PML

    //fprintf(stderr, " after computing Du_x Du_z acdpml_2d_2!!!\n");

    // \phi separate loops along boundary !!!!!!!!!!!!!!!! csq and damp profile
    // compute interior of the domain
    // update wavefield up
    for (i1=s[1]; i1<=e[1]; i1++) {
        for (i0=s[0]; i0<=e[0]; i0++) {
            float lap = (c0* uc[i1][i0] +
                         c1[1]*(uc[i1-1][i0  ] + uc[i1+1][i0  ]) +
                         c2[1]*(uc[i1-2][i0  ] + uc[i1+2][i0  ]) +
                         c1[0]*(uc[i1  ][i0-1] + uc[i1  ][i0+1]) +
                         c2[0]*(uc[i1  ][i0-2] + uc[i1  ][i0+2]));
            
            float cff = 1.0/(1.0+(dp1[i1]+dp0[i0])*dt/2.0);
            float cffuc = (2.0-dp1[i1]*dp0[i0]*dt*dt) * cff;
            float cffup = ((dp1[i1]+dp0[i0])/2.0*dt-1.0) * cff;
            float cff1 = dt*dt/2.0/di[1]*cff;
            float cff0 = dt*dt/2.0/di[0]*cff;
            
            up[i1][i0] = cffuc * uc[i1][i0] + cffup * up[i1][i0] + cff * csq[i1][i0] * lap
            +cff1*(phi1[i1][i0-1]+phi1[i1][i0]-phi1[i1-1][i0-1]-phi1[i1-1][i0])
            +cff0*(phi0[i1-1][i0]+phi0[i1][i0]-phi0[i1-1][i0-1]-phi0[i1][i0-1]);
        }
    }
    
    for (i1=s[1]; i1<e[1]; i1++) {
        for (i0=s[0]; i0<e[0]; i0++) {
            float cff1 = (2.0-dt*dp1[i1])/(2.0+dt*dp1[i1]);
            float cff0 = (2.0-dt*dp0[i0])/(2.0+dt*dp0[i0]);
            float tmp = (csq[i1][i0]+csq[i1+1][i0]+csq[i1][i0+1]+csq[i1+1][i0+1])/4.0*2.0*dt;
            float tmpux= (uc[i1+1][i0]+uc[i1+1][i0+1]-
                          uc[i1  ][i0]-uc[i1  ][i0+1]+
                          up[i1+1][i0]+up[i1+1][i0+1]-
                          up[i1  ][i0]-up[i1  ][i0+1])/4.0/di[1]*(dp0[i0]-dp1[i1])/(2.0+dt*dp1[i1]);
            float tmpuz= (uc[i1][i0+1]+uc[i1+1][i0+1]-
                          uc[i1][i0  ]-uc[i1+1][i0  ]+
                          up[i1][i0+1]+up[i1+1][i0+1]-
                          up[i1][i0  ]-up[i1+1][i0  ])/4.0/di[0]*(dp1[i1]-dp0[i0])/(2.0+dt*dp0[i0]);
            
            phi1[i1][i0]= phi1[i1][i0]*cff1 + tmpux*tmp;
            phi0[i1][i0]= phi0[i1][i0]*cff0 + tmpuz*tmp;
        }
    }
    
    i1=s[1]-1;
    i0=s[0]-1;
    float cff1 = (2.0-dt*dp1[i1])/(2.0+dt*dp1[i1]);
    float cff0 = (2.0-dt*dp0[i0])/(2.0+dt*dp0[i0]);
    float tmpux= (uc[i1+1][i0]+uc[i1+1][i0+1]-
                  uc[i1  ][i0]-uc[i1  ][i0+1]+
                  up[i1+1][i0]+up[i1+1][i0+1]-
                  up[i1  ][i0]-up[i1  ][i0+1])/4.0/di[1]*(dp0[i0]-dp1[i1])/(2.0+dt*dp1[i1]);
    float tmpuz= (uc[i1][i0+1]+uc[i1+1][i0+1]-
                  uc[i1][i0  ]-uc[i1+1][i0  ]+
                  up[i1][i0+1]+up[i1+1][i0+1]-
                  up[i1][i0  ]-up[i1+1][i0  ])/4.0/di[0]*(dp1[i1]-dp0[i0])/(2.0+dt*dp0[i0]);
    
    phi1[i1][i0]= phi1[i1][i0]*cff1 + tmpux*2.0*dt*(csq[i1+1][i0+1]);
	phi0[i1][i0]= phi0[i1][i0]*cff0 + tmpuz*2.0*dt*(csq[i1+1][i0+1]);
    
    // compute i1=s[1]-1
    i1=s[1]-1;
    for (i0=s[0]; i0<=e[0]; i0++){
        cff1 = (2.0-dt*dp1[i1])/(2.0+dt*dp1[i1]);
        cff0 = (2.0-dt*dp0[i0])/(2.0+dt*dp0[i0]);
        tmpux= (uc[i1+1][i0]+uc[i1+1][i0+1]-
                uc[i1  ][i0]-uc[i1  ][i0+1]+
                up[i1+1][i0]+up[i1+1][i0+1]-
                up[i1  ][i0]-up[i1  ][i0+1])/4.0/di[1]*(dp0[i0]-dp1[i1])/(2.0+dt*dp1[i1]);
        tmpuz= (uc[i1][i0+1]+uc[i1+1][i0+1]-
                uc[i1][i0  ]-uc[i1+1][i0  ]+
                up[i1][i0+1]+up[i1+1][i0+1]-
                up[i1][i0  ]-up[i1+1][i0  ])/4.0/di[0]*(dp1[i1]-dp0[i0])/(2.0+dt*dp0[i0]);
        phi1[i1][i0]= phi1[i1][i0]*cff1 + tmpux*2.0*dt*(csq[i1+1][i0]);
        
        phi0[i1][i0]= phi0[i1][i0]*cff0 + tmpuz*2.0*dt*(csq[i1+1][i0]);
    }
    // compute i0=s[0]-1
    i0=s[0]-1;
    for (i1=s[1]; i1<=e[1]; i1++){
        cff1 = (2.0-dt*dp1[i1])/(2.0+dt*dp1[i1]);
        cff0 = (2.0-dt*dp0[i0])/(2.0+dt*dp0[i0]);
        tmpux= (uc[i1+1][i0]+uc[i1+1][i0+1]-
                uc[i1  ][i0]-uc[i1  ][i0+1]+
                up[i1+1][i0]+up[i1+1][i0+1]-
                up[i1  ][i0]-up[i1  ][i0+1])/4.0/di[1]*(dp0[i0]-dp1[i1])/(2.0+dt*dp1[i1]);
        tmpuz= (uc[i1][i0+1]+uc[i1+1][i0+1]-
                uc[i1][i0  ]-uc[i1+1][i0  ]+
                up[i1][i0+1]+up[i1+1][i0+1]-
                up[i1][i0  ]-up[i1+1][i0  ])/4.0/di[0]*(dp1[i1]-dp0[i0])/(2.0+dt*dp0[i0]);
        
        phi1[i1][i0] = phi1[i1][i0]*cff1 + tmpux*2.0*dt*(csq[i1][i0+1]);
        
        phi0[i1][i0] = phi0[i1][i0]*cff0 + tmpuz*2.0*dt*(csq[i1][i0+1]);
    }
    // compute i1=e[1]
    i1=e[1];
    for (i0=s[0]; i0<=e[0]; i0++){
        cff1 = (2.0-dt*dp1[i1])/(2.0+dt*dp1[i1]);
        cff0 = (2.0-dt*dp0[i0])/(2.0+dt*dp0[i0]);
        tmpux= (uc[i1+1][i0]+uc[i1+1][i0+1]-
                uc[i1  ][i0]-uc[i1  ][i0+1]+
                up[i1+1][i0]+up[i1+1][i0+1]-
                up[i1  ][i0]-up[i1  ][i0+1])/4.0/di[1]*(dp0[i0]-dp1[i1])/(2.0+dt*dp1[i1]);
        tmpuz= (uc[i1][i0+1]+uc[i1+1][i0+1]-
                uc[i1][i0  ]-uc[i1+1][i0  ]+
                up[i1][i0+1]+up[i1+1][i0+1]-
                up[i1][i0  ]-up[i1+1][i0  ])/4.0/di[0]*(dp1[i1]-dp0[i0])/(2.0+dt*dp0[i0]);
        
        phi1[i1][i0] = phi1[i1][i0]*cff1 + tmpux*2.0*dt*(csq[i1][i0]);
        
        phi0[i1][i0] = phi0[i1][i0]*cff0 + tmpuz*2.0*dt*(csq[i1][i0]);
    }
    // compute i0=e[0]
    i0=e[0];
    for (i1=s[1]; i1<e[1]; i1++) {
        cff1 = (2.0-dt*dp1[i1])/(2.0+dt*dp1[i1]);
        cff0 = (2.0-dt*dp0[i0])/(2.0+dt*dp0[i0]);
        tmpux= (uc[i1+1][i0]+uc[i1+1][i0+1]-
                uc[i1  ][i0]-uc[i1  ][i0+1]+
                up[i1+1][i0]+up[i1+1][i0+1]-
                up[i1  ][i0]-up[i1  ][i0+1])/4.0/di[1]*(dp0[i0]-dp1[i1])/(2.0+dt*dp1[i1]);
        tmpuz= (uc[i1][i0+1]+uc[i1+1][i0+1]-
                uc[i1][i0  ]-uc[i1+1][i0  ]+
                up[i1][i0+1]+up[i1+1][i0+1]-
                up[i1][i0  ]-up[i1+1][i0  ])/4.0/di[0]*(dp1[i1]-dp0[i0])/(2.0+dt*dp0[i0]);
        
        phi1[i1][i0] = phi1[i1][i0]*cff1 + tmpux*2.0*dt*(csq[i1][i0]);
        
        phi0[i1][i0] = phi0[i1][i0]*cff0 + tmpuz*2.0*dt*(csq[i1][i0]);
    }



    // Homogeneous Dirichlet boundary conditions
  if (lbc[1]) {
#pragma ivdep
    for (i0=s[0];i0<=e[0];i0++) {
      up[s[1]-2][i0]=-up[s[1]][i0];
      up[s[1]-1][i0]=0;
    }
  }

  if (rbc[1]) {
#pragma ivdep
    for (i0=s[0];i0<=e[0];i0++) {
      up[e[1]+2][i0]=-up[e[1]][i0];
      up[e[1]+1][i0]=0;
    }
  }

  if (lbc[0]) {
#pragma ivdep
    for (i1=s[1];i1<=e[1];i1++) {
      up[i1][s[0]-2]=-up[i1][s[0]];
      up[i1][s[0]-1]=0;
    }
  }
  if (rbc[0]) {
#pragma ivdep
    for (i1=s[1];i1<=e[1];i1++) {
      up[i1][e[0]+2]=-up[i1][e[0]];
      up[i1][e[0]+1]=0;
    }
  }

}

