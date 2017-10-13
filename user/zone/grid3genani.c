/* 3-D velocity grid for ray tracing in general anisotropic media. */
/*
  Copyright (C) 2004 University of Texas at Austin
  
  This program is free software; you can redistribute it and/or modify
  it under the terms of the GNU General Public License as published by
  the Free Software Foundation; either version 2 of the License, or
  (at your option) any later version.
  
  This program is distributed in the hope that it will be useful,
  but WITHOUT ANY WARRANTY; without even the implied warranty of
  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
  GNU General Public License for more details.
  
  You should have received a copy of the GNU General Public License
  along with this program; if not, write to the Free Software
  Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA
*/

#include <math.h>
#include <rsf.h>

#include "eno3.h"
#include "bond.h"

#ifndef _grid3genani_h

typedef struct grid3genani* grid3genani;
/* abstract data type */
/*^*/

#endif

struct grid3genani {
    eno3 p11, p12, p13, p14, p15, p16, p22, p23, p24, p25, p26, p33, p34, p35, p36, p44, p45, p46, p55, p56, p66, pthe, pphi;
    int n1, n2, n3;
    float o1, d1, o2, d2, o3, d3;
};
/* concrete data type */

grid3genani grid3genani_init (int n1, float o1, float d1 /* first axis */, 
		    int n2, float o2, float d2 /* second axis */,
		    int n3, float o3, float d3 /* third axis */,
		    float *c11                 /* [n1*n2*n3], gets corrupted */, 
		    float *c12                 /* [n1*n2*n3], gets corrupted */, 
		    float *c13                 /* [n1*n2*n3], gets corrupted */, 
		    float *c14                 /* [n1*n2*n3], 0 if NULL */, 
		    float *c15                 /* [n1*n2*n3], 0 if NULL */, 
		    float *c16                 /* [n1*n2*n3], 0 if NULL */, 
		    float *c22                 /* [n1*n2*n3], gets corrupted */, 
		    float *c23                 /* [n1*n2*n3], gets corrupted */, 
		    float *c24                 /* [n1*n2*n3], 0 if NULL */, 
		    float *c25                 /* [n1*n2*n3], 0 if NULL */, 
		    float *c26                 /* [n1*n2*n3], 0 if NULL */, 
		    float *c33                 /* [n1*n2*n3], gets corrupted */, 
		    float *c34                 /* [n1*n2*n3], 0 if NULL */, 
		    float *c35                 /* [n1*n2*n3], 0 if NULL */, 
		    float *c36                 /* [n1*n2*n3], 0 if NULL */, 
		    float *c44                 /* [n1*n2*n3], gets corrupted */, 
		    float *c45                 /* [n1*n2*n3], 0 if NULL */, 
		    float *c46                 /* [n1*n2*n3], 0 if NULL */, 
		    float *c55                 /* [n1*n2*n3], gets corrupted */, 
		    float *c56                 /* [n1*n2*n3], 0 if NULL */, 
		    float *c66                 /* [n1*n2*n3], gets corrupted */,
		    float *theta               /* [n1*n2*n3], 0 if NULL */, 
		    float *phi                 /* [n1*n2*n3], 0 if NULL */, 
		    int order                  /* interpolation order */)
/*< Initialize grid object >*/
{
    grid3genani grd;
    
    grd = (grid3genani) sf_alloc(1,sizeof(*grd));

    grd->n1 = n1; grd->o1 = o1; grd->d1 = d1; // z
    grd->n2 = n2; grd->o2 = o2; grd->d2 = d2; // y
    grd->n3 = n3; grd->o3 = o3; grd->d3 = d3; // x

    grd->p11 = eno3_init (order, n1, n2, n3);
    grd->p12 = eno3_init (order, n1, n2, n3);
    grd->p13 = eno3_init (order, n1, n2, n3);
    grd->p14 = eno3_init (order, n1, n2, n3);
    grd->p15 = eno3_init (order, n1, n2, n3);
    grd->p16 = eno3_init (order, n1, n2, n3);
    grd->p22 = eno3_init (order, n1, n2, n3);
    grd->p23 = eno3_init (order, n1, n2, n3);
    grd->p24 = eno3_init (order, n1, n2, n3);
    grd->p25 = eno3_init (order, n1, n2, n3);
    grd->p26 = eno3_init (order, n1, n2, n3);
    grd->p33 = eno3_init (order, n1, n2, n3);
    grd->p34 = eno3_init (order, n1, n2, n3);
    grd->p35 = eno3_init (order, n1, n2, n3);
    grd->p36 = eno3_init (order, n1, n2, n3);
    grd->p44 = eno3_init (order, n1, n2, n3);
    grd->p45 = eno3_init (order, n1, n2, n3);
    grd->p46 = eno3_init (order, n1, n2, n3);
    grd->p55 = eno3_init (order, n1, n2, n3);
    grd->p56 = eno3_init (order, n1, n2, n3);
    grd->p66 = eno3_init (order, n1, n2, n3);
    grd->pthe = eno3_init (order, n1, n2, n3);
    grd->pphi = eno3_init (order, n1, n2, n3);

/* If NULL -> orthorhombic */
    eno3_set1 (grd->p11, c11);
    eno3_set1 (grd->p12, c12);
    eno3_set1 (grd->p13, c13);
    if (NULL!=c14) eno3_set1 (grd->p14, c14);
    if (NULL!=c15) eno3_set1 (grd->p15, c15);
    if (NULL!=c16) eno3_set1 (grd->p16, c16);
    eno3_set1 (grd->p22, c22);
    eno3_set1 (grd->p23, c23);
    if (NULL!=c24) eno3_set1 (grd->p24, c24);
    if (NULL!=c25) eno3_set1 (grd->p25, c25);
    if (NULL!=c26) eno3_set1 (grd->p26, c26);
    eno3_set1 (grd->p33, c33);
    if (NULL!=c34) eno3_set1 (grd->p34, c34);
    if (NULL!=c35) eno3_set1 (grd->p35, c35);
    if (NULL!=c36) eno3_set1 (grd->p36, c36);
    eno3_set1 (grd->p44, c44);
    if (NULL!=c45) eno3_set1 (grd->p45, c45);
    if (NULL!=c46) eno3_set1 (grd->p46, c46);
    eno3_set1 (grd->p55, c55);
    if (NULL!=c56) eno3_set1 (grd->p56, c56);
    eno3_set1 (grd->p66, c66);
    if (NULL!=theta) eno3_set1 (grd->pthe,theta);
    if (NULL!=phi) eno3_set1 (grd->pphi,phi);
    return grd;
}

void grid3genani_p_rhs(void* par /* grid */, 
		float* xy /* coordinate [6] (z,y,x) */,
		float* rhs  /* right-hand side [6] */,
		int* rungecount /* count of runge-kutta solve */,
		float** oldg /* Polarization from storage*/)
/*< right-hand side for the qP ray tracing system >*/
{
    grid3genani grd;
    float x, y, z, n1, n2, n3, g[3];
    float c11,c12,c13,c14,c15,c16,c22,c23,c24,c25,c26,c33,c34,c35,c36,c44,c45,c46,c55,c56,c66,the,phi;
    float dc11[3],dc12[3],dc13[3],dc14[3],dc15[3],dc16[3],dc22[3],dc23[3],dc24[3],dc25[3],dc26[3],dc33[3],dc34[3],dc35[3],dc36[3],dc44[3],dc45[3],dc46[3],dc55[3],dc56[3],dc66[3],dthe[3],dphi[3];
    int i, j, k;
    
    grd = (grid3genani) par;
    
    x = (xy[2]-grd->o3)/grd->d3; i = x; x -= i;
    y = (xy[1]-grd->o2)/grd->d2; j = y; y -= j;
    z = (xy[0]-grd->o1)/grd->d1; k = z; z -= k;

    eno3_apply(grd->p11, k, j, i, z, y, x, &c11, dc11, BOTH);
    eno3_apply(grd->p12, k, j, i, z, y, x, &c12, dc12, BOTH);
    eno3_apply(grd->p13, k, j, i, z, y, x, &c13, dc13, BOTH);
    eno3_apply(grd->p14, k, j, i, z, y, x, &c14, dc14, BOTH);
    eno3_apply(grd->p15, k, j, i, z, y, x, &c15, dc15, BOTH);
    eno3_apply(grd->p16, k, j, i, z, y, x, &c16, dc16, BOTH);
    eno3_apply(grd->p22, k, j, i, z, y, x, &c22, dc22, BOTH);
    eno3_apply(grd->p23, k, j, i, z, y, x, &c23, dc23, BOTH);
    eno3_apply(grd->p24, k, j, i, z, y, x, &c24, dc24, BOTH);
    eno3_apply(grd->p25, k, j, i, z, y, x, &c25, dc25, BOTH);
    eno3_apply(grd->p26, k, j, i, z, y, x, &c26, dc26, BOTH);
    eno3_apply(grd->p33, k, j, i, z, y, x, &c33, dc33, BOTH);
    eno3_apply(grd->p34, k, j, i, z, y, x, &c34, dc34, BOTH);
    eno3_apply(grd->p35, k, j, i, z, y, x, &c35, dc35, BOTH);
    eno3_apply(grd->p36, k, j, i, z, y, x, &c36, dc36, BOTH);
    eno3_apply(grd->p44, k, j, i, z, y, x, &c44, dc44, BOTH);
    eno3_apply(grd->p45, k, j, i, z, y, x, &c45, dc45, BOTH);
    eno3_apply(grd->p46, k, j, i, z, y, x, &c46, dc46, BOTH);
    eno3_apply(grd->p55, k, j, i, z, y, x, &c55, dc55, BOTH);
    eno3_apply(grd->p56, k, j, i, z, y, x, &c56, dc56, BOTH);
    eno3_apply(grd->p66, k, j, i, z, y, x, &c66, dc66, BOTH);
    eno3_apply(grd->pthe, k, j, i, z, y, x, &the, dthe, BOTH);
    eno3_apply(grd->pphi, k, j, i, z, y, x, &phi, dphi, BOTH);
    
    /* LA PACK to solve Christoffel */
    /* definition for LAPACK SVD ROUTINEs */
    char    jobz='V';  // for SVD 
    char    uplo='U';  // for SVD 
    int     M=3;       // for SVD 
    int     LDA=M;     // for SVD 
    int     LWORK=8*M; // for SVD 
    int     INFO;      // for SVD 
    double  Chr[9], ww[3], work[24];  // Lapack SVD array 
    
    
    float one = sqrt(xy[3]*xy[3] + xy[4]*xy[4] + xy[5]*xy[5]);
    n1 = xy[5]/one; n2 = xy[4]/one; n3 = xy[3]/one;
    
    /*Bond transformation*/
    bond(&phi,&the,&c11,&c12,&c13,&c14,&c15,&c16,&c22,&c23,&c24,&c25,&c26,&c33,&c34,&c35,&c36,&c44,&c45,&c46,&c55,&c56,&c66);
    
    /* Christoffel matrix */
    Chr[0] = c11*n1*n1 + c66*n2*n2 + c55*n3*n3 + 2*c16*n1*n2 + 2*c15*n1*n3 + 2*c56*n2*n3;
    Chr[4] = c66*n1*n1 + c22*n2*n2 + c44*n3*n3 + 2*c26*n1*n2 + 2*c46*n1*n3 + 2*c24*n2*n3;
    Chr[8] = c55*n1*n1 + c44*n2*n2 + c33*n3*n3 + 2*c45*n1*n2 + 2*c35*n1*n3 + 2*c34*n2*n3;
    Chr[1] = Chr[3] = c16*n1*n1 + c26*n2*n2 + c45*n3*n3 + (c12+c66)*n1*n2 + (c14+c56)*n1*n3 + (c25+c46)*n2*n3; 
    Chr[2] = Chr[6] = c15*n1*n1 + c46*n2*n2 + c35*n3*n3 + (c14+c56)*n1*n2 + (c13+c55)*n1*n3 + (c36+c45)*n2*n3; 
    Chr[5] = Chr[7] = c56*n1*n1 + c24*n2*n2 + c34*n3*n3 + (c25+c46)*n1*n2 + (c36+c45)*n1*n3 + (c23+c44)*n2*n3; 

    /* LAPACK's ssyev routine (slow but accurate) */
    dsyev_(&jobz, &uplo, &M, Chr, &LDA, ww, work, &LWORK, &INFO);
    
        
    int quot = (int) floorf(rungecount[0]/3);
    int rem = (int) rungecount[0]-quot*3;
    
    // qP
    g[0] = Chr[6]; g[1] = Chr[7]; g[2] = Chr[8]; // Polarization (g) of x,y,z
    xy[5] = n1/sqrt(ww[2]); xy[4] = n2/sqrt(ww[2]); xy[3] = n3/sqrt(ww[2]); 
    
    if(rem == 0) {
    	oldg[quot][0] = g[0]; oldg[quot][1] = g[1]; oldg[quot][2] = g[2]; // Store for every last solve of RK3 at each timme step
    }
    rungecount[0] += 1;
    
    
// z
rhs[0] = g[2]*(c13*g[0]*xy[5] + c36*g[1]*xy[5] + c35*g[2]*xy[5] + c36*g[0]*xy[4] + c23*g[1]*xy[4] + c34*g[2]*xy[4] + c35*g[0]*xy[3] + c34*g[1]*xy[3] + c33*g[2]*xy[3]) + 
   g[1]*(c14*g[0]*xy[5] + c46*g[1]*xy[5] + c45*g[2]*xy[5] + c46*g[0]*xy[4] + c24*g[1]*xy[4] + c44*g[2]*xy[4] + c45*g[0]*xy[3] + c44*g[1]*xy[3] + c34*g[2]*xy[3]) + 
   g[0]*(c15*g[0]*xy[5] + c56*g[1]*xy[5] + c55*g[2]*xy[5] + c56*g[0]*xy[4] + c25*g[1]*xy[4] + c45*g[2]*xy[4] + c55*g[0]*xy[3] + c45*g[1]*xy[3] + c35*g[2]*xy[3]);

// y
rhs[1] = g[1]*(c12*g[0]*xy[5] + c26*g[1]*xy[5] + c25*g[2]*xy[5] + c26*g[0]*xy[4] + c22*g[1]*xy[4] + c24*g[2]*xy[4] + c25*g[0]*xy[3] + c24*g[1]*xy[3] + c23*g[2]*xy[3]) + 
   g[2]*(c14*g[0]*xy[5] + c46*g[1]*xy[5] + c45*g[2]*xy[5] + c46*g[0]*xy[4] + c24*g[1]*xy[4] + c44*g[2]*xy[4] + c45*g[0]*xy[3] + c44*g[1]*xy[3] + c34*g[2]*xy[3]) + 
   g[0]*(c16*g[0]*xy[5] + c66*g[1]*xy[5] + c56*g[2]*xy[5] + c66*g[0]*xy[4] + c26*g[1]*xy[4] + c46*g[2]*xy[4] + c56*g[0]*xy[3] + c46*g[1]*xy[3] + c36*g[2]*xy[3]);

// x
rhs[2] = g[0]*(c11*g[0]*xy[5] + c16*g[1]*xy[5] + c15*g[2]*xy[5] + c16*g[0]*xy[4] + c12*g[1]*xy[4] + c14*g[2]*xy[4] + c15*g[0]*xy[3] + c14*g[1]*xy[3] + c13*g[2]*xy[3]) + 
   g[2]*(c15*g[0]*xy[5] + c56*g[1]*xy[5] + c55*g[2]*xy[5] + c56*g[0]*xy[4] + c25*g[1]*xy[4] + c45*g[2]*xy[4] + c55*g[0]*xy[3] + c45*g[1]*xy[3] + c35*g[2]*xy[3]) + 
   g[1]*(c16*g[0]*xy[5] + c66*g[1]*xy[5] + c56*g[2]*xy[5] + c66*g[0]*xy[4] + c26*g[1]*xy[4] + c46*g[2]*xy[4] + c56*g[0]*xy[3] + c46*g[1]*xy[3] + c36*g[2]*xy[3]);

// pz
rhs[3] =  (-(dc11[0]*pow(g[0],2)*pow(xy[5],2)) - dc66[0]*pow(g[1],2)*pow(xy[5],2) - 2*dc15[0]*g[0]*g[2]*pow(xy[5],2) - 2*dc56[0]*g[1]*g[2]*pow(xy[5],2) - dc55[0]*pow(g[2],2)*pow(xy[5],2) - 
        2*dc12[0]*g[0]*g[1]*xy[5]*xy[4] - 2*dc66[0]*g[0]*g[1]*xy[5]*xy[4] - 2*dc26[0]*pow(g[1],2)*xy[5]*xy[4] - 2*dc14[0]*g[0]*g[2]*xy[5]*xy[4] - 2*dc56[0]*g[0]*g[2]*xy[5]*xy[4] - 
        2*dc25[0]*g[1]*g[2]*xy[5]*xy[4] - 2*dc46[0]*g[1]*g[2]*xy[5]*xy[4] - 2*dc45[0]*pow(g[2],2)*xy[5]*xy[4] - dc66[0]*pow(g[0],2)*pow(xy[4],2) - 2*dc26[0]*g[0]*g[1]*pow(xy[4],2) - 
        dc22[0]*pow(g[1],2)*pow(xy[4],2) - 2*dc46[0]*g[0]*g[2]*pow(xy[4],2) - 2*dc24[0]*g[1]*g[2]*pow(xy[4],2) - dc44[0]*pow(g[2],2)*pow(xy[4],2) - 2*dc16[0]*g[0]*xy[5]*(g[1]*xy[5] + g[0]*xy[4]) - 
        2*((dc15[0]*pow(g[0],2) + g[1]*((dc14[0] + dc56[0])*g[0] + dc46[0]*g[1]) + ((dc13[0] + dc55[0])*g[0] + (dc36[0] + dc45[0])*g[1])*g[2] + dc35[0]*pow(g[2],2))*xy[5] + 
        (dc56[0]*pow(g[0],2) + g[1]*((dc25[0] + dc46[0])*g[0] + dc24[0]*g[1]) + ((dc36[0] + dc45[0])*g[0] + (dc23[0] + dc44[0])*g[1])*g[2] + dc34[0]*pow(g[2],2))*xy[4])*xy[3] - 
        (dc55[0]*pow(g[0],2) + 2*dc45[0]*g[0]*g[1] + dc44[0]*pow(g[1],2) + 2*dc35[0]*g[0]*g[2] + 2*dc34[0]*g[1]*g[2] + dc33[0]*pow(g[2],2))*pow(xy[3],2))/2;

// py
rhs[4] =  (-(dc11[1]*pow(g[0],2)*pow(xy[5],2)) - dc66[1]*pow(g[1],2)*pow(xy[5],2) - 2*dc15[1]*g[0]*g[2]*pow(xy[5],2) - 2*dc56[1]*g[1]*g[2]*pow(xy[5],2) - dc55[1]*pow(g[2],2)*pow(xy[5],2) - 
        2*dc12[1]*g[0]*g[1]*xy[5]*xy[4] - 2*dc66[1]*g[0]*g[1]*xy[5]*xy[4] - 2*dc26[1]*pow(g[1],2)*xy[5]*xy[4] - 2*dc14[1]*g[0]*g[2]*xy[5]*xy[4] - 2*dc56[1]*g[0]*g[2]*xy[5]*xy[4] - 
        2*dc25[1]*g[1]*g[2]*xy[5]*xy[4] - 2*dc46[1]*g[1]*g[2]*xy[5]*xy[4] - 2*dc45[1]*pow(g[2],2)*xy[5]*xy[4] - dc66[1]*pow(g[0],2)*pow(xy[4],2) - 2*dc26[1]*g[0]*g[1]*pow(xy[4],2) - 
        dc22[1]*pow(g[1],2)*pow(xy[4],2) - 2*dc46[1]*g[0]*g[2]*pow(xy[4],2) - 2*dc24[1]*g[1]*g[2]*pow(xy[4],2) - dc44[1]*pow(g[2],2)*pow(xy[4],2) - 2*dc16[1]*g[0]*xy[5]*(g[1]*xy[5] + g[0]*xy[4]) - 
        2*((dc15[1]*pow(g[0],2) + g[1]*((dc14[1] + dc56[1])*g[0] + dc46[1]*g[1]) + ((dc13[1] + dc55[1])*g[0] + (dc36[1] + dc45[1])*g[1])*g[2] + dc35[1]*pow(g[2],2))*xy[5] + 
        (dc56[1]*pow(g[0],2) + g[1]*((dc25[1] + dc46[1])*g[0] + dc24[1]*g[1]) + ((dc36[1] + dc45[1])*g[0] + (dc23[1] + dc44[1])*g[1])*g[2] + dc34[1]*pow(g[2],2))*xy[4])*xy[3] - 
        (dc55[1]*pow(g[0],2) + 2*dc45[1]*g[0]*g[1] + dc44[1]*pow(g[1],2) + 2*dc35[1]*g[0]*g[2] + 2*dc34[1]*g[1]*g[2] + dc33[1]*pow(g[2],2))*pow(xy[3],2))/2;

// px
rhs[5] =  (-(dc11[2]*pow(g[0],2)*pow(xy[5],2)) - dc66[2]*pow(g[1],2)*pow(xy[5],2) - 2*dc15[2]*g[0]*g[2]*pow(xy[5],2) - 2*dc56[2]*g[1]*g[2]*pow(xy[5],2) - dc55[2]*pow(g[2],2)*pow(xy[5],2) - 
        2*dc12[2]*g[0]*g[1]*xy[5]*xy[4] - 2*dc66[2]*g[0]*g[1]*xy[5]*xy[4] - 2*dc26[2]*pow(g[1],2)*xy[5]*xy[4] - 2*dc14[2]*g[0]*g[2]*xy[5]*xy[4] - 2*dc56[2]*g[0]*g[2]*xy[5]*xy[4] - 
        2*dc25[2]*g[1]*g[2]*xy[5]*xy[4] - 2*dc46[2]*g[1]*g[2]*xy[5]*xy[4] - 2*dc45[2]*pow(g[2],2)*xy[5]*xy[4] - dc66[2]*pow(g[0],2)*pow(xy[4],2) - 2*dc26[2]*g[0]*g[1]*pow(xy[4],2) - 
        dc22[2]*pow(g[1],2)*pow(xy[4],2) - 2*dc46[2]*g[0]*g[2]*pow(xy[4],2) - 2*dc24[2]*g[1]*g[2]*pow(xy[4],2) - dc44[2]*pow(g[2],2)*pow(xy[4],2) - 2*dc16[2]*g[0]*xy[5]*(g[1]*xy[5] + g[0]*xy[4]) - 
        2*((dc15[2]*pow(g[0],2) + g[1]*((dc14[2] + dc56[2])*g[0] + dc46[2]*g[1]) + ((dc13[2] + dc55[2])*g[0] + (dc36[2] + dc45[2])*g[1])*g[2] + dc35[2]*pow(g[2],2))*xy[5] + 
        (dc56[2]*pow(g[0],2) + g[1]*((dc25[2] + dc46[2])*g[0] + dc24[2]*g[1]) + ((dc36[2] + dc45[2])*g[0] + (dc23[2] + dc44[2])*g[1])*g[2] + dc34[2]*pow(g[2],2))*xy[4])*xy[3] - 
        (dc55[2]*pow(g[0],2) + 2*dc45[2]*g[0]*g[1] + dc44[2]*pow(g[1],2) + 2*dc35[2]*g[0]*g[2] + 2*dc34[2]*g[1]*g[2] + dc33[2]*pow(g[2],2))*pow(xy[3],2))/2;


        rhs[3] /= grd->d3; rhs[4] /= grd->d2; rhs[5] /= grd->d1; // Correct normalization from eno

}

void grid3genani_s1_rhs(void* par /* grid */, 
		float* xy /* coordinate [6] */,
		float* rhs  /* right-hand side [6] */,
		int* rungecount /* count of runge-kutta solve */,
		float** oldg /* Polarization from storage*/)
/*< right-hand side for the qS1 ray tracing system >*/
{
    grid3genani grd;
    float x, y, z, n1, n2, n3, g[3];
    float c11,c12,c13,c14,c15,c16,c22,c23,c24,c25,c26,c33,c34,c35,c36,c44,c45,c46,c55,c56,c66,the,phi;
    float dc11[3],dc12[3],dc13[3],dc14[3],dc15[3],dc16[3],dc22[3],dc23[3],dc24[3],dc25[3],dc26[3],dc33[3],dc34[3],dc35[3],dc36[3],dc44[3],dc45[3],dc46[3],dc55[3],dc56[3],dc66[3],dthe[3],dphi[3];
    int i, j, k;
    
    grd = (grid3genani) par;
    
    x = (xy[0]-grd->o1)/grd->d1; i = x; x -= i;
    y = (xy[1]-grd->o2)/grd->d2; j = y; y -= j;
    z = (xy[2]-grd->o3)/grd->d3; k = z; z -= k;

    eno3_apply(grd->p11, i, j, k, x, y, z, &c11, dc11, BOTH);
    eno3_apply(grd->p12, i, j, k, x, y, z, &c12, dc12, BOTH);
    eno3_apply(grd->p13, i, j, k, x, y, z, &c13, dc13, BOTH);
    eno3_apply(grd->p14, i, j, k, x, y, z, &c14, dc14, BOTH);
    eno3_apply(grd->p15, i, j, k, x, y, z, &c15, dc15, BOTH);
    eno3_apply(grd->p16, i, j, k, x, y, z, &c16, dc16, BOTH);
    eno3_apply(grd->p22, i, j, k, x, y, z, &c22, dc22, BOTH);
    eno3_apply(grd->p23, i, j, k, x, y, z, &c23, dc23, BOTH);
    eno3_apply(grd->p24, i, j, k, x, y, z, &c24, dc24, BOTH);
    eno3_apply(grd->p25, i, j, k, x, y, z, &c25, dc25, BOTH);
    eno3_apply(grd->p26, i, j, k, x, y, z, &c26, dc26, BOTH);
    eno3_apply(grd->p33, i, j, k, x, y, z, &c33, dc33, BOTH);
    eno3_apply(grd->p34, i, j, k, x, y, z, &c34, dc34, BOTH);
    eno3_apply(grd->p35, i, j, k, x, y, z, &c35, dc35, BOTH);
    eno3_apply(grd->p36, i, j, k, x, y, z, &c36, dc36, BOTH);
    eno3_apply(grd->p44, i, j, k, x, y, z, &c44, dc44, BOTH);
    eno3_apply(grd->p45, i, j, k, x, y, z, &c45, dc45, BOTH);
    eno3_apply(grd->p46, i, j, k, x, y, z, &c46, dc46, BOTH);
    eno3_apply(grd->p55, i, j, k, x, y, z, &c55, dc55, BOTH);
    eno3_apply(grd->p56, i, j, k, x, y, z, &c56, dc56, BOTH);
    eno3_apply(grd->p66, i, j, k, x, y, z, &c66, dc66, BOTH);
    eno3_apply(grd->pthe, k, j, i, z, y, x, &the, dthe, BOTH);
    eno3_apply(grd->pphi, k, j, i, z, y, x, &phi, dphi, BOTH);
    
    /* LA PACK to solve Christoffel */
    /* definition for LAPACK SVD ROUTINEs */
    char    jobz='V';  // for SVD 
    char    uplo='U';  // for SVD 
    int     M=3;       // for SVD 
    int     LDA=M;     // for SVD 
    int     LWORK=8*M; // for SVD 
    int     INFO;      // for SVD 
    double  Chr[9], ww[3], work[24];  // Lapack SVD array 
    
    float one = sqrt(xy[3]*xy[3] + xy[4]*xy[4] + xy[5]*xy[5]);
    n1 = xy[5]/one; n2 = xy[4]/one; n3 = xy[3]/one;

    /*Bond transformation*/
    bond(&phi,&the,&c11,&c12,&c13,&c14,&c15,&c16,&c22,&c23,&c24,&c25,&c26,&c33,&c34,&c35,&c36,&c44,&c45,&c46,&c55,&c56,&c66);
    
    /* Christoffel matrix */
    Chr[0] = c11*n1*n1 + c66*n2*n2 + c55*n3*n3 + 2*c16*n1*n2 + 2*c15*n1*n3 + 2*c56*n2*n3;
    Chr[4] = c66*n1*n1 + c22*n2*n2 + c44*n3*n3 + 2*c26*n1*n2 + 2*c46*n1*n3 + 2*c24*n2*n3;
    Chr[8] = c55*n1*n1 + c44*n2*n2 + c33*n3*n3 + 2*c45*n1*n2 + 2*c35*n1*n3 + 2*c34*n2*n3;
    Chr[1] = Chr[3] = c16*n1*n1 + c26*n2*n2 + c45*n3*n3 + (c12+c66)*n1*n2 + (c14+c56)*n1*n3 + (c25+c46)*n2*n3; 
    Chr[2] = Chr[6] = c15*n1*n1 + c46*n2*n2 + c35*n3*n3 + (c14+c56)*n1*n2 + (c13+c55)*n1*n3 + (c36+c45)*n2*n3; 
    Chr[5] = Chr[7] = c56*n1*n1 + c24*n2*n2 + c34*n3*n3 + (c25+c46)*n1*n2 + (c36+c45)*n1*n3 + (c23+c44)*n2*n3; 

    /* LAPACK's ssyev routine (slow but accurate) */
    dsyev_(&jobz, &uplo, &M, Chr, &LDA, ww, work, &LWORK, &INFO);
    

    int quot = (int) floorf(rungecount[0]/3);
    int rem = (int) rungecount[0]-quot*3;
    
    // qS1
    g[0] = Chr[3]; g[1] = Chr[4]; g[2] = Chr[5];
    xy[5] = n1/sqrt(ww[1]); xy[4] = n2/sqrt(ww[1]); xy[3] = n3/sqrt(ww[1]); 
 
    
    if(rem == 0) {
    	oldg[quot][0] = g[0]; oldg[quot][1] = g[1]; oldg[quot][2] = g[2]; // Store for every last solve of RK3 at each timme step
    }
    rungecount[0] += 1;


// z
rhs[0] = g[2]*(c13*g[0]*xy[5] + c36*g[1]*xy[5] + c35*g[2]*xy[5] + c36*g[0]*xy[4] + c23*g[1]*xy[4] + c34*g[2]*xy[4] + c35*g[0]*xy[3] + c34*g[1]*xy[3] + c33*g[2]*xy[3]) + 
   g[1]*(c14*g[0]*xy[5] + c46*g[1]*xy[5] + c45*g[2]*xy[5] + c46*g[0]*xy[4] + c24*g[1]*xy[4] + c44*g[2]*xy[4] + c45*g[0]*xy[3] + c44*g[1]*xy[3] + c34*g[2]*xy[3]) + 
   g[0]*(c15*g[0]*xy[5] + c56*g[1]*xy[5] + c55*g[2]*xy[5] + c56*g[0]*xy[4] + c25*g[1]*xy[4] + c45*g[2]*xy[4] + c55*g[0]*xy[3] + c45*g[1]*xy[3] + c35*g[2]*xy[3]);

// y
rhs[1] = g[1]*(c12*g[0]*xy[5] + c26*g[1]*xy[5] + c25*g[2]*xy[5] + c26*g[0]*xy[4] + c22*g[1]*xy[4] + c24*g[2]*xy[4] + c25*g[0]*xy[3] + c24*g[1]*xy[3] + c23*g[2]*xy[3]) + 
   g[2]*(c14*g[0]*xy[5] + c46*g[1]*xy[5] + c45*g[2]*xy[5] + c46*g[0]*xy[4] + c24*g[1]*xy[4] + c44*g[2]*xy[4] + c45*g[0]*xy[3] + c44*g[1]*xy[3] + c34*g[2]*xy[3]) + 
   g[0]*(c16*g[0]*xy[5] + c66*g[1]*xy[5] + c56*g[2]*xy[5] + c66*g[0]*xy[4] + c26*g[1]*xy[4] + c46*g[2]*xy[4] + c56*g[0]*xy[3] + c46*g[1]*xy[3] + c36*g[2]*xy[3]);

// x
rhs[2] = g[0]*(c11*g[0]*xy[5] + c16*g[1]*xy[5] + c15*g[2]*xy[5] + c16*g[0]*xy[4] + c12*g[1]*xy[4] + c14*g[2]*xy[4] + c15*g[0]*xy[3] + c14*g[1]*xy[3] + c13*g[2]*xy[3]) + 
   g[2]*(c15*g[0]*xy[5] + c56*g[1]*xy[5] + c55*g[2]*xy[5] + c56*g[0]*xy[4] + c25*g[1]*xy[4] + c45*g[2]*xy[4] + c55*g[0]*xy[3] + c45*g[1]*xy[3] + c35*g[2]*xy[3]) + 
   g[1]*(c16*g[0]*xy[5] + c66*g[1]*xy[5] + c56*g[2]*xy[5] + c66*g[0]*xy[4] + c26*g[1]*xy[4] + c46*g[2]*xy[4] + c56*g[0]*xy[3] + c46*g[1]*xy[3] + c36*g[2]*xy[3]);

// pz
rhs[3] =  (-(dc11[0]*pow(g[0],2)*pow(xy[5],2)) - dc66[0]*pow(g[1],2)*pow(xy[5],2) - 2*dc15[0]*g[0]*g[2]*pow(xy[5],2) - 2*dc56[0]*g[1]*g[2]*pow(xy[5],2) - dc55[0]*pow(g[2],2)*pow(xy[5],2) - 
        2*dc12[0]*g[0]*g[1]*xy[5]*xy[4] - 2*dc66[0]*g[0]*g[1]*xy[5]*xy[4] - 2*dc26[0]*pow(g[1],2)*xy[5]*xy[4] - 2*dc14[0]*g[0]*g[2]*xy[5]*xy[4] - 2*dc56[0]*g[0]*g[2]*xy[5]*xy[4] - 
        2*dc25[0]*g[1]*g[2]*xy[5]*xy[4] - 2*dc46[0]*g[1]*g[2]*xy[5]*xy[4] - 2*dc45[0]*pow(g[2],2)*xy[5]*xy[4] - dc66[0]*pow(g[0],2)*pow(xy[4],2) - 2*dc26[0]*g[0]*g[1]*pow(xy[4],2) - 
        dc22[0]*pow(g[1],2)*pow(xy[4],2) - 2*dc46[0]*g[0]*g[2]*pow(xy[4],2) - 2*dc24[0]*g[1]*g[2]*pow(xy[4],2) - dc44[0]*pow(g[2],2)*pow(xy[4],2) - 2*dc16[0]*g[0]*xy[5]*(g[1]*xy[5] + g[0]*xy[4]) - 
        2*((dc15[0]*pow(g[0],2) + g[1]*((dc14[0] + dc56[0])*g[0] + dc46[0]*g[1]) + ((dc13[0] + dc55[0])*g[0] + (dc36[0] + dc45[0])*g[1])*g[2] + dc35[0]*pow(g[2],2))*xy[5] + 
        (dc56[0]*pow(g[0],2) + g[1]*((dc25[0] + dc46[0])*g[0] + dc24[0]*g[1]) + ((dc36[0] + dc45[0])*g[0] + (dc23[0] + dc44[0])*g[1])*g[2] + dc34[0]*pow(g[2],2))*xy[4])*xy[3] - 
        (dc55[0]*pow(g[0],2) + 2*dc45[0]*g[0]*g[1] + dc44[0]*pow(g[1],2) + 2*dc35[0]*g[0]*g[2] + 2*dc34[0]*g[1]*g[2] + dc33[0]*pow(g[2],2))*pow(xy[3],2))/2;

// py
rhs[4] =  (-(dc11[1]*pow(g[0],2)*pow(xy[5],2)) - dc66[1]*pow(g[1],2)*pow(xy[5],2) - 2*dc15[1]*g[0]*g[2]*pow(xy[5],2) - 2*dc56[1]*g[1]*g[2]*pow(xy[5],2) - dc55[1]*pow(g[2],2)*pow(xy[5],2) - 
        2*dc12[1]*g[0]*g[1]*xy[5]*xy[4] - 2*dc66[1]*g[0]*g[1]*xy[5]*xy[4] - 2*dc26[1]*pow(g[1],2)*xy[5]*xy[4] - 2*dc14[1]*g[0]*g[2]*xy[5]*xy[4] - 2*dc56[1]*g[0]*g[2]*xy[5]*xy[4] - 
        2*dc25[1]*g[1]*g[2]*xy[5]*xy[4] - 2*dc46[1]*g[1]*g[2]*xy[5]*xy[4] - 2*dc45[1]*pow(g[2],2)*xy[5]*xy[4] - dc66[1]*pow(g[0],2)*pow(xy[4],2) - 2*dc26[1]*g[0]*g[1]*pow(xy[4],2) - 
        dc22[1]*pow(g[1],2)*pow(xy[4],2) - 2*dc46[1]*g[0]*g[2]*pow(xy[4],2) - 2*dc24[1]*g[1]*g[2]*pow(xy[4],2) - dc44[1]*pow(g[2],2)*pow(xy[4],2) - 2*dc16[1]*g[0]*xy[5]*(g[1]*xy[5] + g[0]*xy[4]) - 
        2*((dc15[1]*pow(g[0],2) + g[1]*((dc14[1] + dc56[1])*g[0] + dc46[1]*g[1]) + ((dc13[1] + dc55[1])*g[0] + (dc36[1] + dc45[1])*g[1])*g[2] + dc35[1]*pow(g[2],2))*xy[5] + 
        (dc56[1]*pow(g[0],2) + g[1]*((dc25[1] + dc46[1])*g[0] + dc24[1]*g[1]) + ((dc36[1] + dc45[1])*g[0] + (dc23[1] + dc44[1])*g[1])*g[2] + dc34[1]*pow(g[2],2))*xy[4])*xy[3] - 
        (dc55[1]*pow(g[0],2) + 2*dc45[1]*g[0]*g[1] + dc44[1]*pow(g[1],2) + 2*dc35[1]*g[0]*g[2] + 2*dc34[1]*g[1]*g[2] + dc33[1]*pow(g[2],2))*pow(xy[3],2))/2;

// px
rhs[5] =  (-(dc11[2]*pow(g[0],2)*pow(xy[5],2)) - dc66[2]*pow(g[1],2)*pow(xy[5],2) - 2*dc15[2]*g[0]*g[2]*pow(xy[5],2) - 2*dc56[2]*g[1]*g[2]*pow(xy[5],2) - dc55[2]*pow(g[2],2)*pow(xy[5],2) - 
        2*dc12[2]*g[0]*g[1]*xy[5]*xy[4] - 2*dc66[2]*g[0]*g[1]*xy[5]*xy[4] - 2*dc26[2]*pow(g[1],2)*xy[5]*xy[4] - 2*dc14[2]*g[0]*g[2]*xy[5]*xy[4] - 2*dc56[2]*g[0]*g[2]*xy[5]*xy[4] - 
        2*dc25[2]*g[1]*g[2]*xy[5]*xy[4] - 2*dc46[2]*g[1]*g[2]*xy[5]*xy[4] - 2*dc45[2]*pow(g[2],2)*xy[5]*xy[4] - dc66[2]*pow(g[0],2)*pow(xy[4],2) - 2*dc26[2]*g[0]*g[1]*pow(xy[4],2) - 
        dc22[2]*pow(g[1],2)*pow(xy[4],2) - 2*dc46[2]*g[0]*g[2]*pow(xy[4],2) - 2*dc24[2]*g[1]*g[2]*pow(xy[4],2) - dc44[2]*pow(g[2],2)*pow(xy[4],2) - 2*dc16[2]*g[0]*xy[5]*(g[1]*xy[5] + g[0]*xy[4]) - 
        2*((dc15[2]*pow(g[0],2) + g[1]*((dc14[2] + dc56[2])*g[0] + dc46[2]*g[1]) + ((dc13[2] + dc55[2])*g[0] + (dc36[2] + dc45[2])*g[1])*g[2] + dc35[2]*pow(g[2],2))*xy[5] + 
        (dc56[2]*pow(g[0],2) + g[1]*((dc25[2] + dc46[2])*g[0] + dc24[2]*g[1]) + ((dc36[2] + dc45[2])*g[0] + (dc23[2] + dc44[2])*g[1])*g[2] + dc34[2]*pow(g[2],2))*xy[4])*xy[3] - 
        (dc55[2]*pow(g[0],2) + 2*dc45[2]*g[0]*g[1] + dc44[2]*pow(g[1],2) + 2*dc35[2]*g[0]*g[2] + 2*dc34[2]*g[1]*g[2] + dc33[2]*pow(g[2],2))*pow(xy[3],2))/2;


        rhs[3] /= grd->d3; rhs[4] /= grd->d2; rhs[5] /= grd->d1; // Correct normalization from eno
}

void grid3genani_s2_rhs(void* par /* grid */, 
		float* xy /* coordinate [6] */,
		float* rhs  /* right-hand side [6] */,
		int* rungecount /* count of runge-kutta solve */,
		float** oldg /* Polarization from storage*/)
/*< right-hand side for the qS2 ray tracing system >*/
{
    grid3genani grd;
    float x, y, z, n1, n2, n3, g[3];
    float c11,c12,c13,c14,c15,c16,c22,c23,c24,c25,c26,c33,c34,c35,c36,c44,c45,c46,c55,c56,c66,the,phi;
    float dc11[3],dc12[3],dc13[3],dc14[3],dc15[3],dc16[3],dc22[3],dc23[3],dc24[3],dc25[3],dc26[3],dc33[3],dc34[3],dc35[3],dc36[3],dc44[3],dc45[3],dc46[3],dc55[3],dc56[3],dc66[3],dthe[3],dphi[3];
    int i, j, k;
    
    grd = (grid3genani) par;
    
    x = (xy[0]-grd->o1)/grd->d1; i = x; x -= i;
    y = (xy[1]-grd->o2)/grd->d2; j = y; y -= j;
    z = (xy[2]-grd->o3)/grd->d3; k = z; z -= k;

    eno3_apply(grd->p11, i, j, k, x, y, z, &c11, dc11, BOTH);
    eno3_apply(grd->p12, i, j, k, x, y, z, &c12, dc12, BOTH);
    eno3_apply(grd->p13, i, j, k, x, y, z, &c13, dc13, BOTH);
    eno3_apply(grd->p14, i, j, k, x, y, z, &c14, dc14, BOTH);
    eno3_apply(grd->p15, i, j, k, x, y, z, &c15, dc15, BOTH);
    eno3_apply(grd->p16, i, j, k, x, y, z, &c16, dc16, BOTH);
    eno3_apply(grd->p22, i, j, k, x, y, z, &c22, dc22, BOTH);
    eno3_apply(grd->p23, i, j, k, x, y, z, &c23, dc23, BOTH);
    eno3_apply(grd->p24, i, j, k, x, y, z, &c24, dc24, BOTH);
    eno3_apply(grd->p25, i, j, k, x, y, z, &c25, dc25, BOTH);
    eno3_apply(grd->p26, i, j, k, x, y, z, &c26, dc26, BOTH);
    eno3_apply(grd->p33, i, j, k, x, y, z, &c33, dc33, BOTH);
    eno3_apply(grd->p34, i, j, k, x, y, z, &c34, dc34, BOTH);
    eno3_apply(grd->p35, i, j, k, x, y, z, &c35, dc35, BOTH);
    eno3_apply(grd->p36, i, j, k, x, y, z, &c36, dc36, BOTH);
    eno3_apply(grd->p44, i, j, k, x, y, z, &c44, dc44, BOTH);
    eno3_apply(grd->p45, i, j, k, x, y, z, &c45, dc45, BOTH);
    eno3_apply(grd->p46, i, j, k, x, y, z, &c46, dc46, BOTH);
    eno3_apply(grd->p55, i, j, k, x, y, z, &c55, dc55, BOTH);
    eno3_apply(grd->p56, i, j, k, x, y, z, &c56, dc56, BOTH);
    eno3_apply(grd->p66, i, j, k, x, y, z, &c66, dc66, BOTH);
    eno3_apply(grd->pthe, k, j, i, z, y, x, &the, dthe, BOTH);
    eno3_apply(grd->pphi, k, j, i, z, y, x, &phi, dphi, BOTH);
    
    /* LA PACK to solve Christoffel */
    /* definition for LAPACK SVD ROUTINEs */
    char    jobz='V';  // for SVD 
    char    uplo='U';  // for SVD 
    int     M=3;       // for SVD 
    int     LDA=M;     // for SVD 
    int     LWORK=8*M; // for SVD 
    int     INFO;      // for SVD 
    double  Chr[9], ww[3], work[24];  // Lapack SVD array 
    
    float one = sqrt(xy[3]*xy[3] + xy[4]*xy[4] + xy[5]*xy[5]);
    n1 = xy[5]/one; n2 = xy[4]/one; n3 = xy[3]/one;
    
    /*Bond transformation*/
    bond(&phi,&the,&c11,&c12,&c13,&c14,&c15,&c16,&c22,&c23,&c24,&c25,&c26,&c33,&c34,&c35,&c36,&c44,&c45,&c46,&c55,&c56,&c66);
    
    /* Christoffel matrix */
    Chr[0] = c11*n1*n1 + c66*n2*n2 + c55*n3*n3 + 2*c16*n1*n2 + 2*c15*n1*n3 + 2*c56*n2*n3;
    Chr[4] = c66*n1*n1 + c22*n2*n2 + c44*n3*n3 + 2*c26*n1*n2 + 2*c46*n1*n3 + 2*c24*n2*n3;
    Chr[8] = c55*n1*n1 + c44*n2*n2 + c33*n3*n3 + 2*c45*n1*n2 + 2*c35*n1*n3 + 2*c34*n2*n3;
    Chr[1] = Chr[3] = c16*n1*n1 + c26*n2*n2 + c45*n3*n3 + (c12+c66)*n1*n2 + (c14+c56)*n1*n3 + (c25+c46)*n2*n3; 
    Chr[2] = Chr[6] = c15*n1*n1 + c46*n2*n2 + c35*n3*n3 + (c14+c56)*n1*n2 + (c13+c55)*n1*n3 + (c36+c45)*n2*n3; 
    Chr[5] = Chr[7] = c56*n1*n1 + c24*n2*n2 + c34*n3*n3 + (c25+c46)*n1*n2 + (c36+c45)*n1*n3 + (c23+c44)*n2*n3; 

    /* LAPACK's ssyev routine (slow but accurate) */
    dsyev_(&jobz, &uplo, &M, Chr, &LDA, ww, work, &LWORK, &INFO);
        
    int quot = (int) floorf(rungecount[0]/3);
    int rem = (int) rungecount[0]-quot*3;
    
    // qS2
    g[0] = Chr[0]; g[1] = Chr[1]; g[2] = Chr[2];
    xy[5] = n1/sqrt(ww[0]); xy[4] = n2/sqrt(ww[0]); xy[3] = n3/sqrt(ww[0]); 
    
    if(rem == 0) {
    	oldg[quot][0] = g[0]; oldg[quot][1] = g[1]; oldg[quot][2] = g[2]; // Store for every last solve of RK3 at each timme step
    }
    rungecount[0] += 1;
    
// z
rhs[0] = g[2]*(c13*g[0]*xy[5] + c36*g[1]*xy[5] + c35*g[2]*xy[5] + c36*g[0]*xy[4] + c23*g[1]*xy[4] + c34*g[2]*xy[4] + c35*g[0]*xy[3] + c34*g[1]*xy[3] + c33*g[2]*xy[3]) + 
   g[1]*(c14*g[0]*xy[5] + c46*g[1]*xy[5] + c45*g[2]*xy[5] + c46*g[0]*xy[4] + c24*g[1]*xy[4] + c44*g[2]*xy[4] + c45*g[0]*xy[3] + c44*g[1]*xy[3] + c34*g[2]*xy[3]) + 
   g[0]*(c15*g[0]*xy[5] + c56*g[1]*xy[5] + c55*g[2]*xy[5] + c56*g[0]*xy[4] + c25*g[1]*xy[4] + c45*g[2]*xy[4] + c55*g[0]*xy[3] + c45*g[1]*xy[3] + c35*g[2]*xy[3]);

// y
rhs[1] = g[1]*(c12*g[0]*xy[5] + c26*g[1]*xy[5] + c25*g[2]*xy[5] + c26*g[0]*xy[4] + c22*g[1]*xy[4] + c24*g[2]*xy[4] + c25*g[0]*xy[3] + c24*g[1]*xy[3] + c23*g[2]*xy[3]) + 
   g[2]*(c14*g[0]*xy[5] + c46*g[1]*xy[5] + c45*g[2]*xy[5] + c46*g[0]*xy[4] + c24*g[1]*xy[4] + c44*g[2]*xy[4] + c45*g[0]*xy[3] + c44*g[1]*xy[3] + c34*g[2]*xy[3]) + 
   g[0]*(c16*g[0]*xy[5] + c66*g[1]*xy[5] + c56*g[2]*xy[5] + c66*g[0]*xy[4] + c26*g[1]*xy[4] + c46*g[2]*xy[4] + c56*g[0]*xy[3] + c46*g[1]*xy[3] + c36*g[2]*xy[3]);

// x
rhs[2] = g[0]*(c11*g[0]*xy[5] + c16*g[1]*xy[5] + c15*g[2]*xy[5] + c16*g[0]*xy[4] + c12*g[1]*xy[4] + c14*g[2]*xy[4] + c15*g[0]*xy[3] + c14*g[1]*xy[3] + c13*g[2]*xy[3]) + 
   g[2]*(c15*g[0]*xy[5] + c56*g[1]*xy[5] + c55*g[2]*xy[5] + c56*g[0]*xy[4] + c25*g[1]*xy[4] + c45*g[2]*xy[4] + c55*g[0]*xy[3] + c45*g[1]*xy[3] + c35*g[2]*xy[3]) + 
   g[1]*(c16*g[0]*xy[5] + c66*g[1]*xy[5] + c56*g[2]*xy[5] + c66*g[0]*xy[4] + c26*g[1]*xy[4] + c46*g[2]*xy[4] + c56*g[0]*xy[3] + c46*g[1]*xy[3] + c36*g[2]*xy[3]);

// pz
rhs[3] =  (-(dc11[0]*pow(g[0],2)*pow(xy[5],2)) - dc66[0]*pow(g[1],2)*pow(xy[5],2) - 2*dc15[0]*g[0]*g[2]*pow(xy[5],2) - 2*dc56[0]*g[1]*g[2]*pow(xy[5],2) - dc55[0]*pow(g[2],2)*pow(xy[5],2) - 
        2*dc12[0]*g[0]*g[1]*xy[5]*xy[4] - 2*dc66[0]*g[0]*g[1]*xy[5]*xy[4] - 2*dc26[0]*pow(g[1],2)*xy[5]*xy[4] - 2*dc14[0]*g[0]*g[2]*xy[5]*xy[4] - 2*dc56[0]*g[0]*g[2]*xy[5]*xy[4] - 
        2*dc25[0]*g[1]*g[2]*xy[5]*xy[4] - 2*dc46[0]*g[1]*g[2]*xy[5]*xy[4] - 2*dc45[0]*pow(g[2],2)*xy[5]*xy[4] - dc66[0]*pow(g[0],2)*pow(xy[4],2) - 2*dc26[0]*g[0]*g[1]*pow(xy[4],2) - 
        dc22[0]*pow(g[1],2)*pow(xy[4],2) - 2*dc46[0]*g[0]*g[2]*pow(xy[4],2) - 2*dc24[0]*g[1]*g[2]*pow(xy[4],2) - dc44[0]*pow(g[2],2)*pow(xy[4],2) - 2*dc16[0]*g[0]*xy[5]*(g[1]*xy[5] + g[0]*xy[4]) - 
        2*((dc15[0]*pow(g[0],2) + g[1]*((dc14[0] + dc56[0])*g[0] + dc46[0]*g[1]) + ((dc13[0] + dc55[0])*g[0] + (dc36[0] + dc45[0])*g[1])*g[2] + dc35[0]*pow(g[2],2))*xy[5] + 
        (dc56[0]*pow(g[0],2) + g[1]*((dc25[0] + dc46[0])*g[0] + dc24[0]*g[1]) + ((dc36[0] + dc45[0])*g[0] + (dc23[0] + dc44[0])*g[1])*g[2] + dc34[0]*pow(g[2],2))*xy[4])*xy[3] - 
        (dc55[0]*pow(g[0],2) + 2*dc45[0]*g[0]*g[1] + dc44[0]*pow(g[1],2) + 2*dc35[0]*g[0]*g[2] + 2*dc34[0]*g[1]*g[2] + dc33[0]*pow(g[2],2))*pow(xy[3],2))/2;

// py
rhs[4] =  (-(dc11[1]*pow(g[0],2)*pow(xy[5],2)) - dc66[1]*pow(g[1],2)*pow(xy[5],2) - 2*dc15[1]*g[0]*g[2]*pow(xy[5],2) - 2*dc56[1]*g[1]*g[2]*pow(xy[5],2) - dc55[1]*pow(g[2],2)*pow(xy[5],2) - 
        2*dc12[1]*g[0]*g[1]*xy[5]*xy[4] - 2*dc66[1]*g[0]*g[1]*xy[5]*xy[4] - 2*dc26[1]*pow(g[1],2)*xy[5]*xy[4] - 2*dc14[1]*g[0]*g[2]*xy[5]*xy[4] - 2*dc56[1]*g[0]*g[2]*xy[5]*xy[4] - 
        2*dc25[1]*g[1]*g[2]*xy[5]*xy[4] - 2*dc46[1]*g[1]*g[2]*xy[5]*xy[4] - 2*dc45[1]*pow(g[2],2)*xy[5]*xy[4] - dc66[1]*pow(g[0],2)*pow(xy[4],2) - 2*dc26[1]*g[0]*g[1]*pow(xy[4],2) - 
        dc22[1]*pow(g[1],2)*pow(xy[4],2) - 2*dc46[1]*g[0]*g[2]*pow(xy[4],2) - 2*dc24[1]*g[1]*g[2]*pow(xy[4],2) - dc44[1]*pow(g[2],2)*pow(xy[4],2) - 2*dc16[1]*g[0]*xy[5]*(g[1]*xy[5] + g[0]*xy[4]) - 
        2*((dc15[1]*pow(g[0],2) + g[1]*((dc14[1] + dc56[1])*g[0] + dc46[1]*g[1]) + ((dc13[1] + dc55[1])*g[0] + (dc36[1] + dc45[1])*g[1])*g[2] + dc35[1]*pow(g[2],2))*xy[5] + 
        (dc56[1]*pow(g[0],2) + g[1]*((dc25[1] + dc46[1])*g[0] + dc24[1]*g[1]) + ((dc36[1] + dc45[1])*g[0] + (dc23[1] + dc44[1])*g[1])*g[2] + dc34[1]*pow(g[2],2))*xy[4])*xy[3] - 
        (dc55[1]*pow(g[0],2) + 2*dc45[1]*g[0]*g[1] + dc44[1]*pow(g[1],2) + 2*dc35[1]*g[0]*g[2] + 2*dc34[1]*g[1]*g[2] + dc33[1]*pow(g[2],2))*pow(xy[3],2))/2;

// px
rhs[5] =  (-(dc11[2]*pow(g[0],2)*pow(xy[5],2)) - dc66[2]*pow(g[1],2)*pow(xy[5],2) - 2*dc15[2]*g[0]*g[2]*pow(xy[5],2) - 2*dc56[2]*g[1]*g[2]*pow(xy[5],2) - dc55[2]*pow(g[2],2)*pow(xy[5],2) - 
        2*dc12[2]*g[0]*g[1]*xy[5]*xy[4] - 2*dc66[2]*g[0]*g[1]*xy[5]*xy[4] - 2*dc26[2]*pow(g[1],2)*xy[5]*xy[4] - 2*dc14[2]*g[0]*g[2]*xy[5]*xy[4] - 2*dc56[2]*g[0]*g[2]*xy[5]*xy[4] - 
        2*dc25[2]*g[1]*g[2]*xy[5]*xy[4] - 2*dc46[2]*g[1]*g[2]*xy[5]*xy[4] - 2*dc45[2]*pow(g[2],2)*xy[5]*xy[4] - dc66[2]*pow(g[0],2)*pow(xy[4],2) - 2*dc26[2]*g[0]*g[1]*pow(xy[4],2) - 
        dc22[2]*pow(g[1],2)*pow(xy[4],2) - 2*dc46[2]*g[0]*g[2]*pow(xy[4],2) - 2*dc24[2]*g[1]*g[2]*pow(xy[4],2) - dc44[2]*pow(g[2],2)*pow(xy[4],2) - 2*dc16[2]*g[0]*xy[5]*(g[1]*xy[5] + g[0]*xy[4]) - 
        2*((dc15[2]*pow(g[0],2) + g[1]*((dc14[2] + dc56[2])*g[0] + dc46[2]*g[1]) + ((dc13[2] + dc55[2])*g[0] + (dc36[2] + dc45[2])*g[1])*g[2] + dc35[2]*pow(g[2],2))*xy[5] + 
        (dc56[2]*pow(g[0],2) + g[1]*((dc25[2] + dc46[2])*g[0] + dc24[2]*g[1]) + ((dc36[2] + dc45[2])*g[0] + (dc23[2] + dc44[2])*g[1])*g[2] + dc34[2]*pow(g[2],2))*xy[4])*xy[3] - 
        (dc55[2]*pow(g[0],2) + 2*dc45[2]*g[0]*g[1] + dc44[2]*pow(g[1],2) + 2*dc35[2]*g[0]*g[2] + 2*dc34[2]*g[1]*g[2] + dc33[2]*pow(g[2],2))*pow(xy[3],2))/2;


        rhs[3] /= grd->d3; rhs[4] /= grd->d2; rhs[5] /= grd->d1; // Correct normalization from eno
}



void grid3genani_s1c_rhs(void* par /* grid */, 
		float* xy /* coordinate [6] */,
		float* rhs  /* right-hand side [6] */,
		int* rungecount /* count of runge-kutta solve */,
		float** oldg /* Polarization from storage*/)
/*< right-hand side for the most coupled qS1 ray tracing system >*/
{
    grid3genani grd;
    float x, y, z, n1, n2, n3, g[3];
    float c11,c12,c13,c14,c15,c16,c22,c23,c24,c25,c26,c33,c34,c35,c36,c44,c45,c46,c55,c56,c66,the,phi;
    float dc11[3],dc12[3],dc13[3],dc14[3],dc15[3],dc16[3],dc22[3],dc23[3],dc24[3],dc25[3],dc26[3],dc33[3],dc34[3],dc35[3],dc36[3],dc44[3],dc45[3],dc46[3],dc55[3],dc56[3],dc66[3],dthe[3],dphi[3];
    int i, j, k;
    
    grd = (grid3genani) par;
    
    x = (xy[0]-grd->o1)/grd->d1; i = x; x -= i;
    y = (xy[1]-grd->o2)/grd->d2; j = y; y -= j;
    z = (xy[2]-grd->o3)/grd->d3; k = z; z -= k;

    eno3_apply(grd->p11, i, j, k, x, y, z, &c11, dc11, BOTH);
    eno3_apply(grd->p12, i, j, k, x, y, z, &c12, dc12, BOTH);
    eno3_apply(grd->p13, i, j, k, x, y, z, &c13, dc13, BOTH);
    eno3_apply(grd->p14, i, j, k, x, y, z, &c14, dc14, BOTH);
    eno3_apply(grd->p15, i, j, k, x, y, z, &c15, dc15, BOTH);
    eno3_apply(grd->p16, i, j, k, x, y, z, &c16, dc16, BOTH);
    eno3_apply(grd->p22, i, j, k, x, y, z, &c22, dc22, BOTH);
    eno3_apply(grd->p23, i, j, k, x, y, z, &c23, dc23, BOTH);
    eno3_apply(grd->p24, i, j, k, x, y, z, &c24, dc24, BOTH);
    eno3_apply(grd->p25, i, j, k, x, y, z, &c25, dc25, BOTH);
    eno3_apply(grd->p26, i, j, k, x, y, z, &c26, dc26, BOTH);
    eno3_apply(grd->p33, i, j, k, x, y, z, &c33, dc33, BOTH);
    eno3_apply(grd->p34, i, j, k, x, y, z, &c34, dc34, BOTH);
    eno3_apply(grd->p35, i, j, k, x, y, z, &c35, dc35, BOTH);
    eno3_apply(grd->p36, i, j, k, x, y, z, &c36, dc36, BOTH);
    eno3_apply(grd->p44, i, j, k, x, y, z, &c44, dc44, BOTH);
    eno3_apply(grd->p45, i, j, k, x, y, z, &c45, dc45, BOTH);
    eno3_apply(grd->p46, i, j, k, x, y, z, &c46, dc46, BOTH);
    eno3_apply(grd->p55, i, j, k, x, y, z, &c55, dc55, BOTH);
    eno3_apply(grd->p56, i, j, k, x, y, z, &c56, dc56, BOTH);
    eno3_apply(grd->p66, i, j, k, x, y, z, &c66, dc66, BOTH);
    eno3_apply(grd->pthe, k, j, i, z, y, x, &the, dthe, BOTH);
    eno3_apply(grd->pphi, k, j, i, z, y, x, &phi, dphi, BOTH);
    
    /* LA PACK to solve Christoffel */
    /* definition for LAPACK SVD ROUTINEs */
    char    jobz='V';  // for SVD 
    char    uplo='U';  // for SVD 
    int     M=3;       // for SVD 
    int     LDA=M;     // for SVD 
    int     LWORK=8*M; // for SVD 
    int     INFO;      // for SVD 
    double  Chr[9], ww[3], work[24];  // Lapack SVD array 
    
    float one = sqrt(xy[3]*xy[3] + xy[4]*xy[4] + xy[5]*xy[5]);
    n1 = xy[5]/one; n2 = xy[4]/one; n3 = xy[3]/one;
    
    /*Bond transformation*/
    bond(&phi,&the,&c11,&c12,&c13,&c14,&c15,&c16,&c22,&c23,&c24,&c25,&c26,&c33,&c34,&c35,&c36,&c44,&c45,&c46,&c55,&c56,&c66);
    
    /* Christoffel matrix */
    Chr[0] = c11*n1*n1 + c66*n2*n2 + c55*n3*n3 + 2*c16*n1*n2 + 2*c15*n1*n3 + 2*c56*n2*n3;
    Chr[4] = c66*n1*n1 + c22*n2*n2 + c44*n3*n3 + 2*c26*n1*n2 + 2*c46*n1*n3 + 2*c24*n2*n3;
    Chr[8] = c55*n1*n1 + c44*n2*n2 + c33*n3*n3 + 2*c45*n1*n2 + 2*c35*n1*n3 + 2*c34*n2*n3;
    Chr[1] = Chr[3] = c16*n1*n1 + c26*n2*n2 + c45*n3*n3 + (c12+c66)*n1*n2 + (c14+c56)*n1*n3 + (c25+c46)*n2*n3; 
    Chr[2] = Chr[6] = c15*n1*n1 + c46*n2*n2 + c35*n3*n3 + (c14+c56)*n1*n2 + (c13+c55)*n1*n3 + (c36+c45)*n2*n3; 
    Chr[5] = Chr[7] = c56*n1*n1 + c24*n2*n2 + c34*n3*n3 + (c25+c46)*n1*n2 + (c36+c45)*n1*n3 + (c23+c44)*n2*n3; 

    /* LAPACK's ssyev routine (slow but accurate) */
    dsyev_(&jobz, &uplo, &M, Chr, &LDA, ww, work, &LWORK, &INFO);
    
    
    int quot = (int) floorf(rungecount[0]/3);
    int rem = (int) rungecount[0]-quot*3;
    
/*    sf_warning("count %d quo %d rem %d", rungecount[0], quot,rem);*/
    
    // qS1
    g[0] = Chr[3]; g[1] = Chr[4]; g[2] = Chr[5];
    xy[5] = n1/sqrt(ww[1]); xy[4] = n2/sqrt(ww[1]); xy[3] = n3/sqrt(ww[1]); 
/*    sf_warning("pickg %g %g %g ",g[0],g[1],g[2]);*/
    
    // Switch to qS2 if not coupled
    if (quot > 1 && fabsf(g[0]*oldg[quot-1][0]+g[1]*oldg[quot-1][1]+g[2]*oldg[quot-1][2]) < fabsf(Chr[0]*oldg[quot-1][0]+Chr[1]*oldg[quot-1][1]+Chr[2]*oldg[quot-1][2]) ) {
        g[0] = Chr[0]; g[1] = Chr[1]; g[2] = Chr[2];
        xy[5] = n1/sqrt(ww[0]); xy[4] = n2/sqrt(ww[0]); xy[3] = n3/sqrt(ww[0]); 
/*        sf_warning("swap oldg %g %g %g to %g %g %g \n",oldg[quot-1][0],oldg[quot-1][1],oldg[quot-1][2],g[0],g[1],g[2]);*/
    }

    if(rem == 0) {
    	oldg[quot][0] = g[0]; oldg[quot][1] = g[1]; oldg[quot][2] = g[2]; // Store for every last solve of RK3 at each timme step
    }
    rungecount[0] += 1;

// z
rhs[0] = g[2]*(c13*g[0]*xy[5] + c36*g[1]*xy[5] + c35*g[2]*xy[5] + c36*g[0]*xy[4] + c23*g[1]*xy[4] + c34*g[2]*xy[4] + c35*g[0]*xy[3] + c34*g[1]*xy[3] + c33*g[2]*xy[3]) + 
   g[1]*(c14*g[0]*xy[5] + c46*g[1]*xy[5] + c45*g[2]*xy[5] + c46*g[0]*xy[4] + c24*g[1]*xy[4] + c44*g[2]*xy[4] + c45*g[0]*xy[3] + c44*g[1]*xy[3] + c34*g[2]*xy[3]) + 
   g[0]*(c15*g[0]*xy[5] + c56*g[1]*xy[5] + c55*g[2]*xy[5] + c56*g[0]*xy[4] + c25*g[1]*xy[4] + c45*g[2]*xy[4] + c55*g[0]*xy[3] + c45*g[1]*xy[3] + c35*g[2]*xy[3]);

// y
rhs[1] = g[1]*(c12*g[0]*xy[5] + c26*g[1]*xy[5] + c25*g[2]*xy[5] + c26*g[0]*xy[4] + c22*g[1]*xy[4] + c24*g[2]*xy[4] + c25*g[0]*xy[3] + c24*g[1]*xy[3] + c23*g[2]*xy[3]) + 
   g[2]*(c14*g[0]*xy[5] + c46*g[1]*xy[5] + c45*g[2]*xy[5] + c46*g[0]*xy[4] + c24*g[1]*xy[4] + c44*g[2]*xy[4] + c45*g[0]*xy[3] + c44*g[1]*xy[3] + c34*g[2]*xy[3]) + 
   g[0]*(c16*g[0]*xy[5] + c66*g[1]*xy[5] + c56*g[2]*xy[5] + c66*g[0]*xy[4] + c26*g[1]*xy[4] + c46*g[2]*xy[4] + c56*g[0]*xy[3] + c46*g[1]*xy[3] + c36*g[2]*xy[3]);

// x
rhs[2] = g[0]*(c11*g[0]*xy[5] + c16*g[1]*xy[5] + c15*g[2]*xy[5] + c16*g[0]*xy[4] + c12*g[1]*xy[4] + c14*g[2]*xy[4] + c15*g[0]*xy[3] + c14*g[1]*xy[3] + c13*g[2]*xy[3]) + 
   g[2]*(c15*g[0]*xy[5] + c56*g[1]*xy[5] + c55*g[2]*xy[5] + c56*g[0]*xy[4] + c25*g[1]*xy[4] + c45*g[2]*xy[4] + c55*g[0]*xy[3] + c45*g[1]*xy[3] + c35*g[2]*xy[3]) + 
   g[1]*(c16*g[0]*xy[5] + c66*g[1]*xy[5] + c56*g[2]*xy[5] + c66*g[0]*xy[4] + c26*g[1]*xy[4] + c46*g[2]*xy[4] + c56*g[0]*xy[3] + c46*g[1]*xy[3] + c36*g[2]*xy[3]);

// pz
rhs[3] =  (-(dc11[0]*pow(g[0],2)*pow(xy[5],2)) - dc66[0]*pow(g[1],2)*pow(xy[5],2) - 2*dc15[0]*g[0]*g[2]*pow(xy[5],2) - 2*dc56[0]*g[1]*g[2]*pow(xy[5],2) - dc55[0]*pow(g[2],2)*pow(xy[5],2) - 
        2*dc12[0]*g[0]*g[1]*xy[5]*xy[4] - 2*dc66[0]*g[0]*g[1]*xy[5]*xy[4] - 2*dc26[0]*pow(g[1],2)*xy[5]*xy[4] - 2*dc14[0]*g[0]*g[2]*xy[5]*xy[4] - 2*dc56[0]*g[0]*g[2]*xy[5]*xy[4] - 
        2*dc25[0]*g[1]*g[2]*xy[5]*xy[4] - 2*dc46[0]*g[1]*g[2]*xy[5]*xy[4] - 2*dc45[0]*pow(g[2],2)*xy[5]*xy[4] - dc66[0]*pow(g[0],2)*pow(xy[4],2) - 2*dc26[0]*g[0]*g[1]*pow(xy[4],2) - 
        dc22[0]*pow(g[1],2)*pow(xy[4],2) - 2*dc46[0]*g[0]*g[2]*pow(xy[4],2) - 2*dc24[0]*g[1]*g[2]*pow(xy[4],2) - dc44[0]*pow(g[2],2)*pow(xy[4],2) - 2*dc16[0]*g[0]*xy[5]*(g[1]*xy[5] + g[0]*xy[4]) - 
        2*((dc15[0]*pow(g[0],2) + g[1]*((dc14[0] + dc56[0])*g[0] + dc46[0]*g[1]) + ((dc13[0] + dc55[0])*g[0] + (dc36[0] + dc45[0])*g[1])*g[2] + dc35[0]*pow(g[2],2))*xy[5] + 
        (dc56[0]*pow(g[0],2) + g[1]*((dc25[0] + dc46[0])*g[0] + dc24[0]*g[1]) + ((dc36[0] + dc45[0])*g[0] + (dc23[0] + dc44[0])*g[1])*g[2] + dc34[0]*pow(g[2],2))*xy[4])*xy[3] - 
        (dc55[0]*pow(g[0],2) + 2*dc45[0]*g[0]*g[1] + dc44[0]*pow(g[1],2) + 2*dc35[0]*g[0]*g[2] + 2*dc34[0]*g[1]*g[2] + dc33[0]*pow(g[2],2))*pow(xy[3],2))/2;

// py
rhs[4] =  (-(dc11[1]*pow(g[0],2)*pow(xy[5],2)) - dc66[1]*pow(g[1],2)*pow(xy[5],2) - 2*dc15[1]*g[0]*g[2]*pow(xy[5],2) - 2*dc56[1]*g[1]*g[2]*pow(xy[5],2) - dc55[1]*pow(g[2],2)*pow(xy[5],2) - 
        2*dc12[1]*g[0]*g[1]*xy[5]*xy[4] - 2*dc66[1]*g[0]*g[1]*xy[5]*xy[4] - 2*dc26[1]*pow(g[1],2)*xy[5]*xy[4] - 2*dc14[1]*g[0]*g[2]*xy[5]*xy[4] - 2*dc56[1]*g[0]*g[2]*xy[5]*xy[4] - 
        2*dc25[1]*g[1]*g[2]*xy[5]*xy[4] - 2*dc46[1]*g[1]*g[2]*xy[5]*xy[4] - 2*dc45[1]*pow(g[2],2)*xy[5]*xy[4] - dc66[1]*pow(g[0],2)*pow(xy[4],2) - 2*dc26[1]*g[0]*g[1]*pow(xy[4],2) - 
        dc22[1]*pow(g[1],2)*pow(xy[4],2) - 2*dc46[1]*g[0]*g[2]*pow(xy[4],2) - 2*dc24[1]*g[1]*g[2]*pow(xy[4],2) - dc44[1]*pow(g[2],2)*pow(xy[4],2) - 2*dc16[1]*g[0]*xy[5]*(g[1]*xy[5] + g[0]*xy[4]) - 
        2*((dc15[1]*pow(g[0],2) + g[1]*((dc14[1] + dc56[1])*g[0] + dc46[1]*g[1]) + ((dc13[1] + dc55[1])*g[0] + (dc36[1] + dc45[1])*g[1])*g[2] + dc35[1]*pow(g[2],2))*xy[5] + 
        (dc56[1]*pow(g[0],2) + g[1]*((dc25[1] + dc46[1])*g[0] + dc24[1]*g[1]) + ((dc36[1] + dc45[1])*g[0] + (dc23[1] + dc44[1])*g[1])*g[2] + dc34[1]*pow(g[2],2))*xy[4])*xy[3] - 
        (dc55[1]*pow(g[0],2) + 2*dc45[1]*g[0]*g[1] + dc44[1]*pow(g[1],2) + 2*dc35[1]*g[0]*g[2] + 2*dc34[1]*g[1]*g[2] + dc33[1]*pow(g[2],2))*pow(xy[3],2))/2;

// px
rhs[5] =  (-(dc11[2]*pow(g[0],2)*pow(xy[5],2)) - dc66[2]*pow(g[1],2)*pow(xy[5],2) - 2*dc15[2]*g[0]*g[2]*pow(xy[5],2) - 2*dc56[2]*g[1]*g[2]*pow(xy[5],2) - dc55[2]*pow(g[2],2)*pow(xy[5],2) - 
        2*dc12[2]*g[0]*g[1]*xy[5]*xy[4] - 2*dc66[2]*g[0]*g[1]*xy[5]*xy[4] - 2*dc26[2]*pow(g[1],2)*xy[5]*xy[4] - 2*dc14[2]*g[0]*g[2]*xy[5]*xy[4] - 2*dc56[2]*g[0]*g[2]*xy[5]*xy[4] - 
        2*dc25[2]*g[1]*g[2]*xy[5]*xy[4] - 2*dc46[2]*g[1]*g[2]*xy[5]*xy[4] - 2*dc45[2]*pow(g[2],2)*xy[5]*xy[4] - dc66[2]*pow(g[0],2)*pow(xy[4],2) - 2*dc26[2]*g[0]*g[1]*pow(xy[4],2) - 
        dc22[2]*pow(g[1],2)*pow(xy[4],2) - 2*dc46[2]*g[0]*g[2]*pow(xy[4],2) - 2*dc24[2]*g[1]*g[2]*pow(xy[4],2) - dc44[2]*pow(g[2],2)*pow(xy[4],2) - 2*dc16[2]*g[0]*xy[5]*(g[1]*xy[5] + g[0]*xy[4]) - 
        2*((dc15[2]*pow(g[0],2) + g[1]*((dc14[2] + dc56[2])*g[0] + dc46[2]*g[1]) + ((dc13[2] + dc55[2])*g[0] + (dc36[2] + dc45[2])*g[1])*g[2] + dc35[2]*pow(g[2],2))*xy[5] + 
        (dc56[2]*pow(g[0],2) + g[1]*((dc25[2] + dc46[2])*g[0] + dc24[2]*g[1]) + ((dc36[2] + dc45[2])*g[0] + (dc23[2] + dc44[2])*g[1])*g[2] + dc34[2]*pow(g[2],2))*xy[4])*xy[3] - 
        (dc55[2]*pow(g[0],2) + 2*dc45[2]*g[0]*g[1] + dc44[2]*pow(g[1],2) + 2*dc35[2]*g[0]*g[2] + 2*dc34[2]*g[1]*g[2] + dc33[2]*pow(g[2],2))*pow(xy[3],2))/2;

        rhs[3] /= grd->d3; rhs[4] /= grd->d2; rhs[5] /= grd->d1; // Correct normalization from eno
}

void grid3genani_s2c_rhs(void* par /* grid */, 
		float* xy /* coordinate [6] */,
		float* rhs  /* right-hand side [6] */,
		int* rungecount /* count of runge-kutta solve */,
		float** oldg /* Polarization from storage*/)
/*< right-hand side for the most coupled qS2 ray tracing system >*/
{
    grid3genani grd;
    float x, y, z, n1, n2, n3, g[3];
    float c11,c12,c13,c14,c15,c16,c22,c23,c24,c25,c26,c33,c34,c35,c36,c44,c45,c46,c55,c56,c66,the,phi;
    float dc11[3],dc12[3],dc13[3],dc14[3],dc15[3],dc16[3],dc22[3],dc23[3],dc24[3],dc25[3],dc26[3],dc33[3],dc34[3],dc35[3],dc36[3],dc44[3],dc45[3],dc46[3],dc55[3],dc56[3],dc66[3],dthe[3],dphi[3];
    int i, j, k;
    
    grd = (grid3genani) par;
    
    x = (xy[0]-grd->o1)/grd->d1; i = x; x -= i;
    y = (xy[1]-grd->o2)/grd->d2; j = y; y -= j;
    z = (xy[2]-grd->o3)/grd->d3; k = z; z -= k;

    eno3_apply(grd->p11, i, j, k, x, y, z, &c11, dc11, BOTH);
    eno3_apply(grd->p12, i, j, k, x, y, z, &c12, dc12, BOTH);
    eno3_apply(grd->p13, i, j, k, x, y, z, &c13, dc13, BOTH);
    eno3_apply(grd->p14, i, j, k, x, y, z, &c14, dc14, BOTH);
    eno3_apply(grd->p15, i, j, k, x, y, z, &c15, dc15, BOTH);
    eno3_apply(grd->p16, i, j, k, x, y, z, &c16, dc16, BOTH);
    eno3_apply(grd->p22, i, j, k, x, y, z, &c22, dc22, BOTH);
    eno3_apply(grd->p23, i, j, k, x, y, z, &c23, dc23, BOTH);
    eno3_apply(grd->p24, i, j, k, x, y, z, &c24, dc24, BOTH);
    eno3_apply(grd->p25, i, j, k, x, y, z, &c25, dc25, BOTH);
    eno3_apply(grd->p26, i, j, k, x, y, z, &c26, dc26, BOTH);
    eno3_apply(grd->p33, i, j, k, x, y, z, &c33, dc33, BOTH);
    eno3_apply(grd->p34, i, j, k, x, y, z, &c34, dc34, BOTH);
    eno3_apply(grd->p35, i, j, k, x, y, z, &c35, dc35, BOTH);
    eno3_apply(grd->p36, i, j, k, x, y, z, &c36, dc36, BOTH);
    eno3_apply(grd->p44, i, j, k, x, y, z, &c44, dc44, BOTH);
    eno3_apply(grd->p45, i, j, k, x, y, z, &c45, dc45, BOTH);
    eno3_apply(grd->p46, i, j, k, x, y, z, &c46, dc46, BOTH);
    eno3_apply(grd->p55, i, j, k, x, y, z, &c55, dc55, BOTH);
    eno3_apply(grd->p56, i, j, k, x, y, z, &c56, dc56, BOTH);
    eno3_apply(grd->p66, i, j, k, x, y, z, &c66, dc66, BOTH);
    eno3_apply(grd->pthe, k, j, i, z, y, x, &the, dthe, BOTH);
    eno3_apply(grd->pphi, k, j, i, z, y, x, &phi, dphi, BOTH);
    
    /* LA PACK to solve Christoffel */
    /* definition for LAPACK SVD ROUTINEs */
    char    jobz='V';  // for SVD 
    char    uplo='U';  // for SVD 
    int     M=3;       // for SVD 
    int     LDA=M;     // for SVD 
    int     LWORK=8*M; // for SVD 
    int     INFO;      // for SVD 
    double  Chr[9], ww[3], work[24];  // Lapack SVD array 
    
    float one = sqrt(xy[3]*xy[3] + xy[4]*xy[4] + xy[5]*xy[5]);
    n1 = xy[5]/one; n2 = xy[4]/one; n3 = xy[3]/one;
    
    /*Bond transformation*/
    bond(&phi,&the,&c11,&c12,&c13,&c14,&c15,&c16,&c22,&c23,&c24,&c25,&c26,&c33,&c34,&c35,&c36,&c44,&c45,&c46,&c55,&c56,&c66);
    
    /* Christoffel matrix */
    Chr[0] = c11*n1*n1 + c66*n2*n2 + c55*n3*n3 + 2*c16*n1*n2 + 2*c15*n1*n3 + 2*c56*n2*n3;
    Chr[4] = c66*n1*n1 + c22*n2*n2 + c44*n3*n3 + 2*c26*n1*n2 + 2*c46*n1*n3 + 2*c24*n2*n3;
    Chr[8] = c55*n1*n1 + c44*n2*n2 + c33*n3*n3 + 2*c45*n1*n2 + 2*c35*n1*n3 + 2*c34*n2*n3;
    Chr[1] = Chr[3] = c16*n1*n1 + c26*n2*n2 + c45*n3*n3 + (c12+c66)*n1*n2 + (c14+c56)*n1*n3 + (c25+c46)*n2*n3; 
    Chr[2] = Chr[6] = c15*n1*n1 + c46*n2*n2 + c35*n3*n3 + (c14+c56)*n1*n2 + (c13+c55)*n1*n3 + (c36+c45)*n2*n3; 
    Chr[5] = Chr[7] = c56*n1*n1 + c24*n2*n2 + c34*n3*n3 + (c25+c46)*n1*n2 + (c36+c45)*n1*n3 + (c23+c44)*n2*n3; 

    /* LAPACK's ssyev routine (slow but accurate) */
    dsyev_(&jobz, &uplo, &M, Chr, &LDA, ww, work, &LWORK, &INFO);
    
    int quot = (int) floorf(rungecount[0]/3);
    int rem = (int) rungecount[0]-quot*3;
    
    // qS2
    g[0] = Chr[0]; g[1] = Chr[1]; g[2] = Chr[2];
    xy[5] = n1/sqrt(ww[0]); xy[4] = n2/sqrt(ww[0]); xy[3] = n3/sqrt(ww[0]);
    
    // Switch to qS1 if not coupled
    if (quot > 1 && fabsf(g[0]*oldg[quot-1][0]+g[1]*oldg[quot-1][1]+g[2]*oldg[quot-1][2]) < fabsf(Chr[3]*oldg[quot-1][0]+Chr[4]*oldg[quot-1][1]+Chr[5]*oldg[quot-1][2])) {
        g[0] = Chr[3]; g[1] = Chr[4]; g[2] = Chr[5];
        xy[5] = n1/sqrt(ww[1]); xy[4] = n2/sqrt(ww[1]); xy[3] = n3/sqrt(ww[1]); 
        
    }
    
    if(rem == 0) {
    	oldg[quot][0] = g[0]; oldg[quot][1] = g[1]; oldg[quot][2] = g[2]; // Store for every last solve of RK3 at each timme step
    }
    rungecount[0] += 1;
    

// z
rhs[0] = g[2]*(c13*g[0]*xy[5] + c36*g[1]*xy[5] + c35*g[2]*xy[5] + c36*g[0]*xy[4] + c23*g[1]*xy[4] + c34*g[2]*xy[4] + c35*g[0]*xy[3] + c34*g[1]*xy[3] + c33*g[2]*xy[3]) + 
   g[1]*(c14*g[0]*xy[5] + c46*g[1]*xy[5] + c45*g[2]*xy[5] + c46*g[0]*xy[4] + c24*g[1]*xy[4] + c44*g[2]*xy[4] + c45*g[0]*xy[3] + c44*g[1]*xy[3] + c34*g[2]*xy[3]) + 
   g[0]*(c15*g[0]*xy[5] + c56*g[1]*xy[5] + c55*g[2]*xy[5] + c56*g[0]*xy[4] + c25*g[1]*xy[4] + c45*g[2]*xy[4] + c55*g[0]*xy[3] + c45*g[1]*xy[3] + c35*g[2]*xy[3]);

// y
rhs[1] = g[1]*(c12*g[0]*xy[5] + c26*g[1]*xy[5] + c25*g[2]*xy[5] + c26*g[0]*xy[4] + c22*g[1]*xy[4] + c24*g[2]*xy[4] + c25*g[0]*xy[3] + c24*g[1]*xy[3] + c23*g[2]*xy[3]) + 
   g[2]*(c14*g[0]*xy[5] + c46*g[1]*xy[5] + c45*g[2]*xy[5] + c46*g[0]*xy[4] + c24*g[1]*xy[4] + c44*g[2]*xy[4] + c45*g[0]*xy[3] + c44*g[1]*xy[3] + c34*g[2]*xy[3]) + 
   g[0]*(c16*g[0]*xy[5] + c66*g[1]*xy[5] + c56*g[2]*xy[5] + c66*g[0]*xy[4] + c26*g[1]*xy[4] + c46*g[2]*xy[4] + c56*g[0]*xy[3] + c46*g[1]*xy[3] + c36*g[2]*xy[3]);

// x
rhs[2] = g[0]*(c11*g[0]*xy[5] + c16*g[1]*xy[5] + c15*g[2]*xy[5] + c16*g[0]*xy[4] + c12*g[1]*xy[4] + c14*g[2]*xy[4] + c15*g[0]*xy[3] + c14*g[1]*xy[3] + c13*g[2]*xy[3]) + 
   g[2]*(c15*g[0]*xy[5] + c56*g[1]*xy[5] + c55*g[2]*xy[5] + c56*g[0]*xy[4] + c25*g[1]*xy[4] + c45*g[2]*xy[4] + c55*g[0]*xy[3] + c45*g[1]*xy[3] + c35*g[2]*xy[3]) + 
   g[1]*(c16*g[0]*xy[5] + c66*g[1]*xy[5] + c56*g[2]*xy[5] + c66*g[0]*xy[4] + c26*g[1]*xy[4] + c46*g[2]*xy[4] + c56*g[0]*xy[3] + c46*g[1]*xy[3] + c36*g[2]*xy[3]);

// pz
rhs[3] =  (-(dc11[0]*pow(g[0],2)*pow(xy[5],2)) - dc66[0]*pow(g[1],2)*pow(xy[5],2) - 2*dc15[0]*g[0]*g[2]*pow(xy[5],2) - 2*dc56[0]*g[1]*g[2]*pow(xy[5],2) - dc55[0]*pow(g[2],2)*pow(xy[5],2) - 
        2*dc12[0]*g[0]*g[1]*xy[5]*xy[4] - 2*dc66[0]*g[0]*g[1]*xy[5]*xy[4] - 2*dc26[0]*pow(g[1],2)*xy[5]*xy[4] - 2*dc14[0]*g[0]*g[2]*xy[5]*xy[4] - 2*dc56[0]*g[0]*g[2]*xy[5]*xy[4] - 
        2*dc25[0]*g[1]*g[2]*xy[5]*xy[4] - 2*dc46[0]*g[1]*g[2]*xy[5]*xy[4] - 2*dc45[0]*pow(g[2],2)*xy[5]*xy[4] - dc66[0]*pow(g[0],2)*pow(xy[4],2) - 2*dc26[0]*g[0]*g[1]*pow(xy[4],2) - 
        dc22[0]*pow(g[1],2)*pow(xy[4],2) - 2*dc46[0]*g[0]*g[2]*pow(xy[4],2) - 2*dc24[0]*g[1]*g[2]*pow(xy[4],2) - dc44[0]*pow(g[2],2)*pow(xy[4],2) - 2*dc16[0]*g[0]*xy[5]*(g[1]*xy[5] + g[0]*xy[4]) - 
        2*((dc15[0]*pow(g[0],2) + g[1]*((dc14[0] + dc56[0])*g[0] + dc46[0]*g[1]) + ((dc13[0] + dc55[0])*g[0] + (dc36[0] + dc45[0])*g[1])*g[2] + dc35[0]*pow(g[2],2))*xy[5] + 
        (dc56[0]*pow(g[0],2) + g[1]*((dc25[0] + dc46[0])*g[0] + dc24[0]*g[1]) + ((dc36[0] + dc45[0])*g[0] + (dc23[0] + dc44[0])*g[1])*g[2] + dc34[0]*pow(g[2],2))*xy[4])*xy[3] - 
        (dc55[0]*pow(g[0],2) + 2*dc45[0]*g[0]*g[1] + dc44[0]*pow(g[1],2) + 2*dc35[0]*g[0]*g[2] + 2*dc34[0]*g[1]*g[2] + dc33[0]*pow(g[2],2))*pow(xy[3],2))/2;

// py
rhs[4] =  (-(dc11[1]*pow(g[0],2)*pow(xy[5],2)) - dc66[1]*pow(g[1],2)*pow(xy[5],2) - 2*dc15[1]*g[0]*g[2]*pow(xy[5],2) - 2*dc56[1]*g[1]*g[2]*pow(xy[5],2) - dc55[1]*pow(g[2],2)*pow(xy[5],2) - 
        2*dc12[1]*g[0]*g[1]*xy[5]*xy[4] - 2*dc66[1]*g[0]*g[1]*xy[5]*xy[4] - 2*dc26[1]*pow(g[1],2)*xy[5]*xy[4] - 2*dc14[1]*g[0]*g[2]*xy[5]*xy[4] - 2*dc56[1]*g[0]*g[2]*xy[5]*xy[4] - 
        2*dc25[1]*g[1]*g[2]*xy[5]*xy[4] - 2*dc46[1]*g[1]*g[2]*xy[5]*xy[4] - 2*dc45[1]*pow(g[2],2)*xy[5]*xy[4] - dc66[1]*pow(g[0],2)*pow(xy[4],2) - 2*dc26[1]*g[0]*g[1]*pow(xy[4],2) - 
        dc22[1]*pow(g[1],2)*pow(xy[4],2) - 2*dc46[1]*g[0]*g[2]*pow(xy[4],2) - 2*dc24[1]*g[1]*g[2]*pow(xy[4],2) - dc44[1]*pow(g[2],2)*pow(xy[4],2) - 2*dc16[1]*g[0]*xy[5]*(g[1]*xy[5] + g[0]*xy[4]) - 
        2*((dc15[1]*pow(g[0],2) + g[1]*((dc14[1] + dc56[1])*g[0] + dc46[1]*g[1]) + ((dc13[1] + dc55[1])*g[0] + (dc36[1] + dc45[1])*g[1])*g[2] + dc35[1]*pow(g[2],2))*xy[5] + 
        (dc56[1]*pow(g[0],2) + g[1]*((dc25[1] + dc46[1])*g[0] + dc24[1]*g[1]) + ((dc36[1] + dc45[1])*g[0] + (dc23[1] + dc44[1])*g[1])*g[2] + dc34[1]*pow(g[2],2))*xy[4])*xy[3] - 
        (dc55[1]*pow(g[0],2) + 2*dc45[1]*g[0]*g[1] + dc44[1]*pow(g[1],2) + 2*dc35[1]*g[0]*g[2] + 2*dc34[1]*g[1]*g[2] + dc33[1]*pow(g[2],2))*pow(xy[3],2))/2;

// px
rhs[5] =  (-(dc11[2]*pow(g[0],2)*pow(xy[5],2)) - dc66[2]*pow(g[1],2)*pow(xy[5],2) - 2*dc15[2]*g[0]*g[2]*pow(xy[5],2) - 2*dc56[2]*g[1]*g[2]*pow(xy[5],2) - dc55[2]*pow(g[2],2)*pow(xy[5],2) - 
        2*dc12[2]*g[0]*g[1]*xy[5]*xy[4] - 2*dc66[2]*g[0]*g[1]*xy[5]*xy[4] - 2*dc26[2]*pow(g[1],2)*xy[5]*xy[4] - 2*dc14[2]*g[0]*g[2]*xy[5]*xy[4] - 2*dc56[2]*g[0]*g[2]*xy[5]*xy[4] - 
        2*dc25[2]*g[1]*g[2]*xy[5]*xy[4] - 2*dc46[2]*g[1]*g[2]*xy[5]*xy[4] - 2*dc45[2]*pow(g[2],2)*xy[5]*xy[4] - dc66[2]*pow(g[0],2)*pow(xy[4],2) - 2*dc26[2]*g[0]*g[1]*pow(xy[4],2) - 
        dc22[2]*pow(g[1],2)*pow(xy[4],2) - 2*dc46[2]*g[0]*g[2]*pow(xy[4],2) - 2*dc24[2]*g[1]*g[2]*pow(xy[4],2) - dc44[2]*pow(g[2],2)*pow(xy[4],2) - 2*dc16[2]*g[0]*xy[5]*(g[1]*xy[5] + g[0]*xy[4]) - 
        2*((dc15[2]*pow(g[0],2) + g[1]*((dc14[2] + dc56[2])*g[0] + dc46[2]*g[1]) + ((dc13[2] + dc55[2])*g[0] + (dc36[2] + dc45[2])*g[1])*g[2] + dc35[2]*pow(g[2],2))*xy[5] + 
        (dc56[2]*pow(g[0],2) + g[1]*((dc25[2] + dc46[2])*g[0] + dc24[2]*g[1]) + ((dc36[2] + dc45[2])*g[0] + (dc23[2] + dc44[2])*g[1])*g[2] + dc34[2]*pow(g[2],2))*xy[4])*xy[3] - 
        (dc55[2]*pow(g[0],2) + 2*dc45[2]*g[0]*g[1] + dc44[2]*pow(g[1],2) + 2*dc35[2]*g[0]*g[2] + 2*dc34[2]*g[1]*g[2] + dc33[2]*pow(g[2],2))*pow(xy[3],2))/2;


        rhs[3] /= grd->d3; rhs[4] /= grd->d2; rhs[5] /= grd->d1; // Correct normalization from eno
}


int grid3genani_term (void* par /* grid */, 
		 float* xy /* location [3] */)
/*< Termination criterion. returns 0 if xy (data coordinates)
  are inside the grid >*/
{
    grid3genani grd;
    
    grd = (grid3genani) par;
    return (xy[0] < grd->o1 || xy[0] > grd->o1 + (grd->n1-1)*grd->d1 || 
	    xy[1] < grd->o2 || xy[1] > grd->o2 + (grd->n2-1)*grd->d2 ||
	    xy[2] < grd->o3 || xy[2] > grd->o3 + (grd->n3-1)*grd->d3);
}

void grid3genani_close(grid3genani grd)
/*< Free internal storage >*/
{
    eno3_close (grd->p11);
    eno3_close (grd->p12);
    eno3_close (grd->p13);
    eno3_close (grd->p14);
    eno3_close (grd->p15);
    eno3_close (grd->p16);
    eno3_close (grd->p22);
    eno3_close (grd->p23);
    eno3_close (grd->p24);
    eno3_close (grd->p25);
    eno3_close (grd->p26);
    eno3_close (grd->p33);
    eno3_close (grd->p34);
    eno3_close (grd->p35);
    eno3_close (grd->p36);
    eno3_close (grd->p44);
    eno3_close (grd->p45);
    eno3_close (grd->p46);
    eno3_close (grd->p55);
    eno3_close (grd->p56);
    eno3_close (grd->p66);
    free (grd);
}

/* 	$Id$	 */
