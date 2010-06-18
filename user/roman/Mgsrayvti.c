/* Gauss Seidel iterative solver for phase space escape positions, angle and traveltime */
/*
   Downwind first-order scheme with angular grid ]-180,180]
   Escape angle from vertical for slowness pointing upward

   Equation for escape positions or angle is (with appropriate BC)
   -S.sin(a).dX/dx - S.cos(a).dX/dz - [cos(a).Sx - sin(a).Sz].dX/da = 0

   Equation for escape traveltime T is (with appropriate BC)
   -S.sin(a).dT/dx - S.cos(a).dT/dz - [cos(a).Sx - sin(a).Sz].dT/da = -S*S
*/
/*
  Copyright (C) 2009 University of Texas at Austin
  
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

#include <rsf.h>
#include <assert.h>

#include "gsray.h"

const float epsilon=1e-5f;

struct VTIStruct {
    float 
        c11, c33, c13, c44, c66,
	alfa, beta, eps, gamma, delta,
	c1, c2, c3, c4, c5; 
    float 
    c11x, c11z, c33x, c33z, c13x, c13z, c44x, c44z, c66x, c66z; 
    float 
    c1x, c1z, c2x, c2z, c3x, c3z, c4x, c4z, c5x, c5z; 
};
/* c33=3*c44, c44=gaussian */

void vti_init_greenriver_c(struct VTIStruct *V)
/* float * c11, float *c33, float *c13, float *c44, float *c66, float *c44x, float *c44z) */
{
    const float
	a11 = 15.0638,
	a13 = 1.6381,
	a12 = 6.5616,
	a44 = 3.1258,
	a33 = 3*a44; /*10.8373,*/


    /* (c11-a12)/2 = c66 */
    V->c11 = a11;
    V->c33 = a33;
    V->c13 = a13;
    V->c44 = a44;
    V->c66 = (a11 - a12)/2.f;    
    
    V->c11x = V->c11z =
	V->c44x = V->c44z = 
	V->c33x = V->c33z =
	V->c13x = V->c13z =
	V->c66x = V->c66z = 0.f; 
}
/* rho == 1 */
void vti_thomsen2c(float alfa, float beta, float eps, float gamma, float delta,
		   float * c11, float *c33, float *c13, float *c44, float *c66)
{
    *c33 = alfa*alfa;
    *c44 = beta*beta;
    *c11 = eps * 2.f * (*c33) + (*c33);
    *c66 = gamma * 2.f * (*c44) + (*c44);
    *c13 = 0.5f * (delta * 2.f * (*c33)*(*c33) + ((*c33)-(*c44))*(*c11 + *c33 - 2*(*c44)));
    assert(*c13 > 0);
    *c13 = sqrtf(*c13) -(*c44);
}
void vti_c2thomsen(struct VTIStruct *V)
//float *alfa, float *beta, float *eps, float *gamma, float *delta,
//float c11, float c33, float c13, float c44, float c66)
{
    (V->alfa) = sqrtf(V->c33);
    (V->beta) = sqrtf(V->c44);
    (V->eps) = (V->c11 - V->c33)/(2.f*V->c33);
    (V->gamma) = (V->c66 - V->c44) / (2.f*V->c44);
    (V->delta) = (2*(V->c13+V->c44)*(V->c13+V->c44)-(V->c33-V->c44)*(V->c11+V->c33-2*V->c44))/(2*V->c33*V->c33);
}
void vti_gauss_2_c(struct VTIStruct *V, 
		   float v, float vx, float vz, 
		   float eps, float gamma, float delta)
{        
    const float v2 = v*v;
    float tmp;

    V->c33 = 3*v;
    V->c44 = v;

    V->c11 = (eps * 2.f + 1.f) * V->c33;
    V->c66 = (gamma*2.f + 1.f) * V->c44;

    tmp = delta*9*v2 + v*(V->c11+v); 
    assert(tmp > 0.f); 
    V->c13 = sqrtf(tmp) - V->c44;
    assert(V->c13>0);

    /* V->a12 = V->c11 + 2.f*V->c66; */

    /* grad x, grad z */
    V->c33x=vx*3; 
    V->c33z=vz*3;

    V->c44x=vx;   
    V->c44z=vz;

    V->c11x=(eps*2.f+1.f)*V->c33x;
    V->c11z=(eps*2.f+1.f)*V->c33z;

    V->c66x = (gamma*2.f + 1.f) * V->c44x;
    V->c66z = (gamma*2.f + 1.f) * V->c44z;

    V->c13x = -vx + 0.5f/tmp*(18*delta*v*vx + 2*v*vx + V->c11x*v + V->c11*vx);
    V->c13z = -vz + 0.5f/tmp*(18*delta*v*vz + 2*v*vz + V->c11z*v + V->c11*vz);
}
void vti_c11_2_c1(struct VTIStruct *V)
{
    V->c1 = V->c11 * V->c44;	
    V->c2 = V->c11 * V->c33 + V->c44*V->c44 - (V->c13 + V->c44)*(V->c13+V->c44);
    V->c3 = V->c33 * V->c44;
    V->c4 = -(V->c11 + V->c44);
    V->c5 = -(V->c33 + V->c44);
    /* V->a12 = V->c11 - 2.f*V->c66;*/

    V->c1x = V->c11x * V->c44 + V->c11 * V->c44x;
    V->c1z = V->c11z * V->c44 + V->c11 * V->c44z;

    V->c2x = V->c11x*V->c33 + V->c11*V->c33x + 2*V->c44x*V->c44 - 2*(V->c13 + V->c44)*(V->c13x+V->c44x);
    V->c2z = V->c11z*V->c33 + V->c11*V->c33z + 2*V->c44z*V->c44 - 2*(V->c13 + V->c44)*(V->c13z+V->c44z);

    V->c3x = V->c33x*V->c44+V->c33*V->c44x;
    V->c3z = V->c33z*V->c44+V->c33*V->c44z;

    V->c4x = -(V->c11x + V->c44x);
    V->c4z = -(V->c11z + V->c44z);

    V->c5x = -(V->c33x + V->c44x);
    V->c5z = -(V->c33z + V->c44z);
}
float vti_slowness_c(struct VTIStruct *V, float ca2, float sa2, char is_vti, float * Si) 
{
    double v2, y1, y2;
    const double 
	ca4=ca2*ca2, 
	sa4=sa2*sa2;

    switch (is_vti) {
	case 'p':
	    y1 = V->c4*ca2 + V->c5*sa2;
	    y2 = V->c1*ca4 + V->c2*ca2*sa2 + V->c3*sa4;

	    v2 = 0.5f * (-y1 + sqrt(y1*y1 - 4*y2));
	    break;
	case 'v':
	    y1 = V->c4*ca2 + V->c5*sa2;
	    y2 = V->c1*ca4 + V->c2*ca2*sa2 + V->c3*sa4;

	    v2 = 0.5f * (-y1 - sqrt(y1*y1 - 4*y2));
	    break;
	case 'h':
	    v2 = V->c66*ca2 + V->c44*sa2;
	    break;
	default:
	    assert(0);
    }
    assert(v2 > 0);
    *Si = sqrtf(v2);
    return 1.f/(*Si); /* slowness */
}
void Svti(struct VTIStruct *V, 
	  float S, 
	  float Si,  /* inverse slowness S */
	  char is_vti, 
	  float ca, float sa, 
	  float * Hp, float *Hq, float * dt_dsigma, float * dl_dsigma, float *da_dsigma)
{    
    const double 
	p = -ca*S, 
	p2 = p*p,
	p3 = p2*p,
	q = -sa*S,
	q2 = q*q,
	q3 = q2*q;
    
    const double 
	h = V->c1*p3*p+V->c2*p2*q2+V->c3*q3*q+V->c4*p2+V->c5*q2 + 1.0,
	hh = V->c66*p2 + V->c44*q2 - 1.0;
    float Hx = 0, Hz = 0;
    

    switch (is_vti) {
	case 'p':
	case 'v':
	    *Hp = 4*V->c1*p3 + 2*V->c2*p*q2 + 2*V->c4*p;
	    *Hq = 4*V->c3*q3 + 2*V->c2*q*p2 + 2*V->c5*q;

	    if (da_dsigma) {
		Hx = V->c1x*p3*p+V->c2x*p2*q2+V->c3x*q3*q+V->c4x*p2+V->c5x*q2;
		Hz = V->c1z*p3*p+V->c2z*p2*q2+V->c3z*q3*q+V->c4z*p2+V->c5z*q2;
	    }

	    assert (fabs(h) < epsilon);
	    break;
	case 'h':
	    *Hp = 2*V->c66*p; /* (V->c11 - V->a12)*/
	    *Hq = 2*V->c44*q;
	    /* c44=gaussian anomaly(x,z) and c33=3*c44 */
	    if (da_dsigma) {
		Hx = V->c66x*p2 + V->c44x*q2;
		Hz = V->c66z*p2 + V->c44z*q2;
	    }

	    assert (fabs(hh) < 1e-6f);
	    break;
	default:
	    assert(0);
    }

    if (dt_dsigma)
	*dt_dsigma = fabs (p*(*Hp) + q*(*Hq));

    if (dl_dsigma)
	*dl_dsigma = sqrt( (*Hp)*(*Hp) + (*Hq)*(*Hq) );

    if (da_dsigma)
	*da_dsigma = Si*(Hx*ca - Hz*sa);
}
float vti_slowness_thomsen(struct VTIStruct *V, 
			   float ca, float sa, float ca2, float sa2, char is_vti)
{
    const float 
	a2 = V->alfa* V->alfa,
	b2 =  V->beta* V->beta,
	ba = b2/a2,
	ab = a2/b2;

    float d, v2;

    switch (is_vti) {
	case 'p':
	    d = 4.f* V->delta*sa2*ca2/((1.f-ba)*(1.f-ba));

	    d += 4.f*(1.f-ba+ V->eps)* V->eps*sa2*sa2/((1.f-ba)*(1.f-ba));

	    d *= 0.5f * (1.f - ba);

	    v2 = a2*(1.f +  V->eps*sa2 + d);
	    break;

	case 'v':
	    d = 4.f*V->delta*sa2*ca2/((1.f-ba)*(1.f-ba));

	    d += 4.f*(1.f-ba+ V->eps)* V->eps*sa2*sa2/((1.f-ba)*(1.f-ba));

	    d *= 0.5f * (1.f - ba);


	    v2 = b2*(1.f + ab*(sa2-d));
	    break;

	case 'h':
	    v2 = b2*(1.f + 2.f* V->gamma*sa2);
	    break;
	default:
	    assert(0);
    }

    return (float)(1./sqrt(v2)); /* slowness */
}
/***************************************************************/
int main(int argc, char* argv[])
{
    int nz,nx,na;
    float oz,ox,oa;
    float dz,dx,da;

    int iq;                        /* escape variables switch */  

    int ix;                        /* grid points in x */
    int iz;                        /* grid points in z */
    int ia;                        /* grid points in a */

    int iter, niter;               /* number of iterations */

    int ix0, ix1, ixs;
    int iz0, iz1, izs;
    int ia0, ia1, ias;

    float a;                   /* angle */

    float ***t;                    /* escape variable */
    float **v,**vx,**vz;           /* slowness, gradients */

    float vel, velx, velz;

    float cs,sn;
    float new_val;
    float cvt,tol;

    sf_file in,out,slow,slowz,slowx;
    float sl_vti, sl_vti_inv;
    char  is_P_SH_SV, *whatV;
    float Hp, Hq, p_gradH, dl_dsigma, da_dsigma;
    float 
	vti_eps = 0.20, /*0.2, */
	vti_gamma=0.20, /*0.2, */
	vti_delta=-0.450; /*-0.45;*/

    struct VTIStruct VTI;

    /* bool swap_dir_x=1, swap_dir_z=1; */

    sf_init(argc,argv);

    in = sf_input("in");
    out = sf_output("out");

    if (!sf_getint("iq",&iq)) sf_error("Need iq=");
    /* switch for escape variable 0=x, 1=a, 2=t, 3=z, 4=l */

    whatV = sf_getstring ("vti");
    /* what to compute (p=qP, v=qSV, h=SH) */
    if (NULL == whatV) {
        whatV = "p";
    } else {
        if (whatV[0] != 'p' && whatV[0] != 'v' && whatV[0] != 'h')
            sf_error ("Need vti=p|v|h");
    }
    is_P_SH_SV = whatV[0];


    /* read input file parameters */
    if (SF_FLOAT != sf_gettype(in)) sf_error("Need float");

    if (!sf_histint(in,"n1",&nz)) sf_error("No n1=");
    if (!sf_histfloat(in,"d1",&dz)) sf_error("No d1=");
    if (!sf_histfloat(in,"o1",&oz)) sf_error("No o1=");

    if (!sf_histint(in,"n2",&nx)) sf_error("No n2=");
    if (!sf_histfloat(in,"d2",&dx)) sf_error("No d2=");
    if (!sf_histfloat(in,"o2",&ox)) sf_error("No o2=");

    if (!sf_histint(in,"n3",&na)) sf_error("No n3=");
    if (!sf_histfloat(in,"d3",&da)) sf_error("No d3=");
    if (!sf_histfloat(in,"o3",&oa)) sf_error("No o3=");
    /* angle in degrees */

    if (!sf_getint("niter",&niter)) niter=50;
    /* number of Gauss-Seidel iterations */

    if (!sf_getfloat("tol",&tol)) tol=0.000002*nx*nz;
    /* accuracy tolerance */

    if (!sf_getfloat("vti_eps",&vti_eps)) vti_eps=0.0;
    if (!sf_getfloat("vti_gamma",&vti_gamma)) vti_gamma=0.0;
    if (!sf_getfloat("vti_delta",&vti_delta)) vti_delta=0.0;
    /* VTI constants Thomsen  */


    /* memory allocations */

    /* Solution vector */
    t = sf_floatalloc3(nz,nx,na);

    /* read input escape variable - initial guess */
    sf_floatread(t[0][0],nz*nx*na,in);
    
    /*
    for (int ix=0; ix < nx; ix++)
	    for (int iz=0; iz < nz; iz++)
    		for (int ia=0; ia < na; ia++)
			t[ia][ix][iz] = 1e10f;
    */

    /* read auxiliary slowness file */
    slow = sf_input("vel");
    v = sf_floatalloc2(nz,nx);
    sf_floatread(v[0],nz*nx,slow);

    /* read auxiliary slowness z-gradient file */
    slowz = sf_input("velz");
    vz = sf_floatalloc2(nz,nx);
    sf_floatread(vz[0],nz*nx,slowz);

    /* read auxiliary slowness x-gradient file */
    slowx = sf_input("velx");
    vx = sf_floatalloc2(nz,nx);
    sf_floatread(vx[0],nz*nx,slowx);

    /* convert to radians */
    oa *= SF_PI/180.;
    da *= SF_PI/180.;

    ia0=0; ia1=na; ias=1;

    gsray_init(nz,nx,na,
	       oz,ox,oa,
	       dz,dx,da);

    /* Greenriver shale example 
     alfa=3.29 beta=1.77 eps=0.195 gamma=0.18 delat=-0.45 */
    vti_init_greenriver_c(&VTI);
    vti_c2thomsen(&VTI);


    for (iter=0; iter < niter; iter++) {

	cvt = 0.0; /* X_next - X_prev */
	int num_fails = 0;

        /* Gauss-Seidel iteration on angle */
	for (ia = ia0; ia != ia1; ia += ias) {
			    
	    a = oa + ia*da;

	    cs = cosf(a);
	    sn = sinf(a);
    
	    /* assert (fabs(sn) > 1e-6 && fabs(cs) > 1e-6);
	    Hq = -sn;
	    Hp = -cs; */
	    vti_gauss_2_c(&VTI, v[0][0], vx[0][0], vz[0][0], vti_eps, vti_gamma, vti_delta);

	    vti_c11_2_c1(&VTI);

	    sl_vti = vti_slowness_c(&VTI, cs*cs, sn*sn, is_P_SH_SV, &sl_vti_inv);

	    Svti(&VTI, sl_vti, sl_vti_inv, is_P_SH_SV, cs, sn, &Hp, &Hq, (float*)0, (float*)0, (float*)0);//&da_dsigma);

	    if (1e-6 < Hq) {
		ix0=nx-2; ix1=-1; ixs=-1;
		for (iz = 0; iz < nz; iz++) 
		    boundary_mat(t,iz,nx-1,ia,iq);

	    } else {
		ix0=1; ix1=nx; ixs=1;
		for (iz = 0; iz < nz; iz++) 
		    boundary_mat(t,iz,0,ia,iq);
	    }

	    if (1e-6 < Hp) {
		iz0=nz-2; iz1=-1; izs=-1;
		for (ix = 0; ix < nx; ix++) 
		    boundary_mat(t,nz-1,ix,ia,iq);
	    } else {
		iz0=1; iz1=nz; izs=1;
		for (ix = 0; ix < nx; ix++) 
		    boundary_mat(t,0,ix,ia,iq);
	    }

	    /* loop over grid Z.X */
	    for (ix = ix0; ix != ix1; ix += ixs) {

		for (iz = iz0; iz != iz1; iz += izs) {
		    
		    vel = v[ix][iz]; 
		    velx = vx[ix][iz];
		    velz = vz[ix][iz];

		    vti_gauss_2_c(&VTI, vel, velx, velz, vti_eps, vti_gamma, vti_delta);

		    vti_c11_2_c1(&VTI);

		    sl_vti = vti_slowness_c(&VTI, cs*cs, sn*sn, is_P_SH_SV, &sl_vti_inv);

		    Svti(&VTI, sl_vti, sl_vti_inv, is_P_SH_SV, cs, sn, &Hp, &Hq, &p_gradH, &dl_dsigma, &da_dsigma);

		    vti_c2thomsen(&VTI);
		    assert(fabs(VTI.eps - vti_eps) < 1e-6);
		    assert(fabs(VTI.gamma - vti_gamma) < 1e-6);
		    assert(fabs(VTI.delta - vti_delta) < 1e-6);


		    /* Gauss-Seidel update 
 	            new_val = gs_update(t,iz,ix,ia,ss,ssz,ssx,cs,sn,iq); */

		    new_val = gs_update(t,-Hp, -Hq, -da_dsigma, iz,ix,ia, p_gradH, dl_dsigma, iq);
		    
		    /* new_val = fmin(new_val, t[ia][ix][iz]); */

		    cvt += fabsf(t[ia][ix][iz]-new_val);

		    /* t[ia][ix][iz] = fmin(new_val, t[ia][ix][iz]); */

		    t[ia][ix][iz] = new_val;

		    /* if (t[ia][ix][iz] > 1000) num_fails++; */
	
		} /* ix */
		    
	    } /* iz */

	} /* ia */

	sf_warning("Iter = %d, Norm L1 = %g num fails = %d sweep=%g %g",iter,cvt/(nx*nz*na), num_fails,Hp,Hq);

        /* tol is tolerance for convergence */
	if (cvt < tol) break;

        /* alternate updating direction on angle grid*/
	if (0 == ia0) {
	    ia0 = na-1; ia1 = -1; ias = -1;
	} else {
	    ia0 = 0; ia1 = na; ias = 1;
	}

    } /* end G-S iterations */

    /* output */
    sf_floatwrite(t[0][0],nz*nx*na,out);
    
    exit(0);
}


