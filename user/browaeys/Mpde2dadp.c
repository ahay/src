/* Numerical solution of linear pde 2-d (X-Z-a) for phase space escape positions, angle and traveltime */
/*
  Upwind first order scheme with angular grid ]-180,180]
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


typedef enum {dim_Z, dim_X, dim_A} dim_type;
typedef enum {ESC_UNINIT, ESC_UP, ESC_DOWN, ESC_LEFT, ESC_RIGHT} escape_type;
typedef enum {PADD, PERIODYC} padd_type;

static float update_mat(float ***t, int iz, int ix, int ia, 
			float fzmax, float fzmin, 
			float fxmax, float fxmin, 
			float famax, float famin, int nz, int nx, int na)
{
    float ts;

    ts = 0.0;

    /* z-direction - DOWNWIND SCHEME */
    if (iz == 0) {

	ts +=  /*fzmin*t[-1][ix][ia] + */ fzmax*t[1][ix][ia];
	assert (fzmin < 1e-6);

    } else if (iz == (nz-1)) {

	ts += fzmin*t[nz-2][ix][ia];
	assert (fzmax < 1e-6);

    } else {

	ts += fzmax*t[iz+1][ix][ia] + fzmin*t[iz-1][ix][ia];

    }

    /* x-direction - DOWNWIND SCHEME */
    if (ix == 0) {

	ts += fxmax*t[iz][1][ia];
	assert (fxmin < 1e-6);

    } else if (ix == (nx-1)) {

	ts += fxmin*t[iz][nx-2][ia];
	assert (fxmax < 1e-6);

    } else {

	ts += fxmax*t[iz][ix+1][ia] + fxmin*t[iz][ix-1][ia];

    }

    /* escape variable periodic on a-sides ]-PI,+PI] - DOWNWIND SCHEME */
    if (ia == 0) {

	ts += famin*t[iz][ix][na-1] + famax*t[iz][ix][1];

    } else if (ia == (na-1)) {

	ts += famin*t[iz][ix][na-2] + famax*t[iz][ix][0];

    } else {

	ts += famin*t[iz][ix][ia-1] + famax*t[iz][ix][ia+1];

    }

    return (ts);
}


static void boundary_mat(float ***t, int iz, int ix, int ia, int iq, int *** t_esc, int esc_dir, 
			 float oz, float ox, float oa, float dz, float dx, float da)
{
    /* Outgoing rays i.e. four different boundary conditions */
    /*  each BC consists of two edges: (X=0,Z) and (X,Z=zmax) etc. */
    /* iq defines the function: Tscape, Xescape, Zescape, or Aescape. 
       BC should be set s.t. this type */
    
    switch(iq) {

	case 0:
            /* escape location x - it is just x :) */
	    t[iz][ix][ia] = ox + ix * dx;
	    break;
	case 1:
	    /* escape angle a - it is just a */
	    t[iz][ix][ia] = oa + ia*da;
	    break;
	case 2:
	    /* escape traveltime - it is zero for finite domains
	       , for infinite domains it is well defined */
	    t[iz][ix][ia] = 0.0;
	    t_esc[iz][ix][ia] = esc_dir;
	    break;
	case 3:
	    /* escape depth z - just z */
	    t[iz][ix][ia] = oz + iz*dz;
	    break;
    }

    return;
}
	    
/* z0 = either 0 or NZ-1,  escape either at X=0 or X=NX-1
   therefor t = min (t escape at X=0, t escape at X=NX-1 */
void boundary_init_z0(int ia, int nx, int z0, float dx,
		      float ** sl, float *** t, float ** _sz, float px, float pz, int *** t_esc)
{
    int ix;

    // the escape trajectory is along diagonal tan(pz/px): dx / cos()
    const float dx_diag = dx / fabs(px / sqrt(px*px + pz*pz));
    assert(fabs(px) > 1e-4);

    if (0 > px) { // escape from X=0
    
	t[z0][0][ia] = 0.f;

	for (ix = 1; ix < nx; ix++) {

	    t[z0][ix][ia] = t[z0][ix - 1][ia] + dx_diag * sl[ix - 1][z0];

	    t_esc[z0][ix][ia] = ESC_UP;

	    //if (z0==0) sf_warning("iz=%d (==0) ix = %d (increase) t=%g (nonzero) ",z0, ix, t[z0][ix][ia]);

	    assert(fabs(_sz[ix - 1][z0])<1e-6);
	}
    }
    else { // escape from X=Nx-1
 
	t[z0][nx - 1][ia] = 0.f;

	for (ix = nx-2; ix > 0; ix--) {

	    t[z0][ix][ia] = t[z0][ix + 1][ia] + dx_diag * sl[ix + 1][z0];

	    t_esc[z0][ix][ia] = ESC_DOWN;

	    //sf_warning("iz=%d ix = %d t=%g",z0, ix, t[z0][ix][ia]);
	    assert(fabs(_sz[ix - 1][z0])<1e-6);
	}
    }
}
/* x0 = either 0 or Nx-1,  escape either at X=0 or X=NX-1
   therefor t = min (t escape at X=0, t escape at X=NX-1 */
void boundary_init_x0(int ia, int nz, int ix0, float dz,
		      float ** sl, float *** t, float ** _sx, float px, float pz, 
		      int *** t_esc)
{
    int iz;

    const float cos_a = fabs(pz / sqrt(px*px + pz*pz));

    const float dz_diag = dz / cos_a;

    assert (cos_a > 1e-4);

    if (0 > pz) { // escape through Z = 0 
	assert (0 == ix0 || px > 0);

	t[0][ix0][ia] = 0.f;

	for (iz = 1; iz < nz; iz++) {

	    t[iz][ix0][ia] = t[iz - 1][ix0][ia] + dz_diag * sl[ix0][iz - 1];
      
	    t_esc[iz][ix0][ia] = ESC_LEFT;
    
	    //sf_warning("boundary_init_x0: iz=%d ix0 = %d t=%g diag_dz = %g",iz, ix0, t[iz][ix0][ia],dz_diag);

	    assert(fabs(_sx[ix0][iz-1])<1e-6);
	}
    }
    else { // pz > 0 => escape from Nz-1
	assert(0 == ix0 || px > 0);

	t[nz - 1][ix0][ia] = 0.f;

	for (iz = nz-2; iz > 0; iz--) {
      
	    t[iz][ix0][ia] = t[iz + 1][ix0][ia] + dz_diag * sl[ix0][iz + 1];

	    t_esc[iz][ix0][ia] = ESC_RIGHT;

	    //sf_warning("boundary_init_x0: iz=%d ix = %d t=%g ",iz, ix0, t[iz][ix0][ia]);

	    assert(fabs(_sx[ix0][iz-1])<1e-6);
	}
    }
}

void boundary_sweep_setup(int iq,
			  float px, float pz, int nx, int nz, int ia, 
			  float ** slo, float ** sx, float ** sz,
			  //float dx, float dz,
			  int is_xinf, int is_zinf,
			  int * ix0, int * ix1, int * ixs, 
			  int * iz0, int * iz1, int * izs, 
			  float *** t, int *** t_esc,
			  float oz, float ox, float oa, float dz, float dx, float da)
{
    int iz, ix;

    if (1e-6 < px) { // escape to the right thru X=Xmax

	*ix0=nx-2; *ix1=-1; *ixs=-1;

	if (2 == iq && is_xinf) {
      
	    boundary_init_x0(ia, nz, nx - 1, dz, slo, t, sx, px, pz, t_esc);

	}
	else {
	    for (iz = 0; iz < nz; iz++) 
		boundary_mat(t,iz,nx-1,ia,iq, t_esc, ESC_RIGHT, oz, ox, oa, dz, dx, da);
	}
  
    }
    else {// px < 0 => init left edge
    
	*ix0=1; *ix1=nx; *ixs=1;
  
	if (2 == iq && is_xinf) {
      
	    boundary_init_x0(ia, nz, 0, dz, slo, t, sx, px, pz, t_esc);

	}
	else {
	    for (iz = 0; iz < nz; iz++) 
		boundary_mat(t,iz,0,ia,iq, t_esc, ESC_LEFT, oz, ox, oa, dz, dx, da);
    
	}
    } // end if px><0

    if (1e-6 < pz) {

	*iz0=nz-2; *iz1=-1; *izs=-1;
 
	if (2 == iq && is_zinf) { // infinite in Z=>init at z=Nz
      
	    boundary_init_z0(ia, nx, nz - 1, dx, slo, t, sz, px, pz, t_esc);

	}
	else { // finite in Z=>zero at Z=maxZ -edge

	    for (ix = 0; ix < nx; ix++) 
	
		boundary_mat(t,nz-1,ix,ia,iq, t_esc, ESC_DOWN, oz, ox, oa, dz, dx, da);
	}
    } 
    else {// pz < 0

	*iz0=1; *iz1=nz; *izs=1;    

	if (2 == iq && is_zinf) { // infinite in Z=>init at Z==0

	    boundary_init_z0(ia, nx, 0, dx, slo, t, sz, px, pz, t_esc);

	}
	else { // finite in Z=>zero at Z=0 -edge

	    for (ix = 0; ix < nx; ix++) 
	
		boundary_mat(t,0,ix,ia,iq, t_esc, ESC_UP, oz, ox, oa, dz, dx, da);
	}
    }
}

	
/******************************************************

   Finite volumes: 
   
   S = const => T(z,x, a=a0) - 2-D solution

   example 1: method==0: UPWIND first order 
              1-st order convergence as dx*dz -> 0

   example 2: Courant-Fridrich-Levy necessary stability criteria
              in UPWIND

   example 3: Lax higher order - high resolution 
              (operates with slopes)

   example 4: high-resolution adaptive slope:
              handles singularities and not dessipate as others

   S != const => truly T(z,x,a) - 3-D solution
   
***************************************************/    
float minmod (float a, float b)
{
    if (a*b > 1e-6) {
	if (fabs(a) < fabs(b))
	    return a;
	return b;
    }
    return 0.f; 
}

float maxmod (float a, float b)
{
    if (a*b > 1e-6) {
	if (fabs(a) > fabs(b))
	    return a;
	return b;
    }
    return 0.f; 
}

// float psi(float theta) {  return minmod(1.f, theta); }
/**************
   SMOOTH CORNERS AT X=Z=0
**************/
int ix_smooth (int iz_smooth_edge_pos, int ixsmooth, int izsmooth, int nx)
{
    int	    ix_smooth_edge_pos = 0;

    if (izsmooth >= iz_smooth_edge_pos) {

	// iz = izsmooth (1 - sin (pi/2 * ix / ixsmooth))
	ix_smooth_edge_pos = 
	    2.f * (float)ixsmooth / 3.1415 * asin (1.f - iz_smooth_edge_pos / (float)izsmooth);

	if (ix_smooth_edge_pos >=  nx)
	    ix_smooth_edge_pos = nx-1;

	if (ix_smooth_edge_pos < 0)
	    ix_smooth_edge_pos = 0;
	// circle arc:
	//floor(0.5 + sm_rad - sqrt ((float)(sm_rad*sm_rad - (iz_smooth_edge_pos-sm_rad) * (iz_smooth_edge_pos-sm_rad))));
    }

//if (1 == ia && ix_smooth_edge_pos) 
//sf_warning("ix0 = %d iz_smooth=%d ix_smooth=%d",ix_smooth_edge_pos,iz_smooth_edge_pos,ix_smooth_edge_pos);

    return ix_smooth_edge_pos;
}
/****************
   a(na) := a(0) and a(-1):=a(na)
****************/
int periodic_BC (int in_ia, int na)
{
    int ia = in_ia; 

    if (ia < 0) 
	ia = na + ia;
    else  
	if (ia > na - 1)
	    ia = ia - na;
    assert (0 <= ia && na > ia);
    return ia;
}
/***************************/
int is_beyond_padd(int iz, int nz, int * ret_iz)// = (int*)0)
{    
    if (iz  >= 0 && iz < nz) {
	if (ret_iz)
	    * ret_iz = iz;
	return 0;
    }

    if (ret_iz) {
	if (iz < 0) 
	    *ret_iz = 0;
	if (iz > nz - 1)
	    *ret_iz = nz - 1;
    }

    return 1;
}

/************************
  Init arrays for CTU multi-dim:
   U, V, F, G Nx x 3={ia-1,ia,ia+1}
   fluxes F = G = zero
   coefs Tz + U Tx + V Ta = ...
***********************/
void init_fluctuation_step_dz_dx_da (int nz, int nx, int na, 
				     int iz, int ia, 
				     int izs,
				     float oa, float da, 
				     float ** s, float ** sx, float ** sz, 
				     float ** *U, float ** *V, float ** *F, float ** *G)
{
    int indx, inda;

    *U = sf_floatalloc2(3, nx);
    *V = sf_floatalloc2(3, nx);
    *F = sf_floatalloc2(3, nx);
    *G = sf_floatalloc2(3, nx);

    for (indx=0; indx < nx; indx++) {
	for (inda=0; inda < 3; inda++) {
			  
	    (*F)[indx][inda] = (*G)[indx][inda] = 0.f;
			  
	    const float
		a = oa + (ia + inda - 1) *da,
		cs = cosf(a),
		sn = sinf(a),
		ss = s[indx][iz - izs],
		ssx = sx[indx][iz - izs],
		ssz = sz[indx][iz - izs];
			  
	    // u = fx/fz = sn/cs
	    (*U)[indx][inda] = sn / cs;
	    // fa/fz = (cs * ssx - sn * ssz) / (ss*cs)
	    (*V)[indx][inda] = (cs * ssx - sn * ssz) / (ss*cs);
	}
    }
}
/* CTU employs fluctuation form with A+-, B+- and fluxes F,G 
   Leveque p.442 (19.19) Finite volumes*/
float fluctuation_step_dz_dx_da (int iserles_method,
				 float dz, float dx, float da, 
				 int iz, int ix, int ia,
				 int izs, int ixs, int ias, 
				 int nz, int nx, int na,
				 float ** u, float **v, 
				 float ** F, float ** G,
				 float *** t)
{
    int 
	ia_prev = periodic_BC(ia - ias, na),
	ia_plus = periodic_BC(ia + ias, na),
	ix_plus = ix + ixs;

    if (ix_plus > nx-1) 
	ix_plus = nx - 1;

    assert (ix - ixs >= 0 && iz - izs>=0 && ia >= 0 && ia < na);

    float 
	dQ_i_j_iminus1_j = t[iz - izs][ix][ia] - t[iz - izs][ix-ixs][ia],//Q_i_j - Q_iminus1_j,	 

	dQ_i_j_i_jminus1 = t[iz - izs][ix][ia] - t[iz - izs][ix][ia_prev], //Q_i_j - Q_i_jminus1,

	dQ_i_jplus1_i_j = t[iz - izs][ix][ia_plus] - t[iz-izs][ix][ia];

    float dQ_iplus1_j_i_j;

    if (ix_plus < 0 || ix_plus >= nx) {

	dQ_iplus1_j_i_j = dQ_i_j_iminus1_j;

	if (ix_plus < 0) ix_plus = 0;
	if (ix_plus >= nx) ix_plus = nx-1;
    }
    else

	dQ_iplus1_j_i_j = t[iz - izs][ix_plus][ia] - t[iz-izs][ix][ia];//Q_i+1_j - Q_i_j,
    

    // trick: since dim (u,v,F,G) in direction "a" is 3.
    // here for ia0 = 1, => ia0-1 = 0 and ia0+1=2
    const int ia0 = 1;

    float
	// (19.19 p.442)
	delta_Q_nplus1_i_j =

	-dz / dx * (SF_MAX(0.f, u[ix][ia0]) * dQ_i_j_iminus1_j  + SF_MIN(0.f, u[ix][ia0])* dQ_iplus1_j_i_j)

	-dz / da * (SF_MAX(0.f, v[ix][ia0]) * dQ_i_j_i_jminus1  + SF_MIN(0.f, v[ix][ia0])* dQ_i_jplus1_i_j);


    if (0 == iserles_method) { // UPWIND 

	if (0 == ia && 1 == iz && 1 == ix) sf_warning(" 2-D UPWIND Fluctuation");

	return  delta_Q_nplus1_i_j;
    }
     

    //   Corner Transport Upwind 
    if (0 == ia && 1 == iz && 1 == ix) sf_warning(" 2-D Corner Transport UPWIND Fluctuation"); 
    assert (ix - ixs >= 0 && ix-ixs < nx);

    delta_Q_nplus1_i_j +=  // (19.19 p.442)

	-dz / dx * (F[ix_plus][ia0] - F[ix][ia0])
	 
	-dz / da * (G[ix][ia0 + ias] - G[ix][ia0]);

    const float
	dt_2dx = 0.5f * dz / dx,
	dt_2da = 0.5f * dz / da;

    G[ix - ixs][ ia0 ]       -= dt_2dx * SF_MIN(0.f, v[ix - ixs][ia0])     * SF_MIN(0.f, u[ix][ia0]) * dQ_i_j_iminus1_j;

    G[ix - ixs][ia0 + ias]   -= dt_2dx * SF_MAX(0.f, v[ix - ixs][ia0 + ias]) * SF_MIN(0.f, u[ix][ia0]) * dQ_i_j_iminus1_j;

    G[ix][ia0]               -= dt_2dx * SF_MIN(0.f, v[ix][ia0])           * SF_MAX(0.f, u[ix][ia0]) * dQ_i_j_iminus1_j;

    G[ix][ia0 + ias]         -= dt_2dx * SF_MAX(0.f, v[ix][ia0 + ias])     * SF_MAX(0.f, u[ix][ia0]) * dQ_i_j_iminus1_j;
	   
    F[ix][ia0 - ias]         -= dt_2da * SF_MIN(0.f, u[ix][ia0 - ias])     * SF_MIN(0.f, v[ix][ia0]) * dQ_i_j_i_jminus1;

    F[ix_plus][ia0 - ias]    -= dt_2da * SF_MAX(0.f, u[ix_plus][ia0 - ias]) * SF_MIN(0.f, v[ix][ia0]) * dQ_i_j_i_jminus1;

    F[ix][ia0]               -= dt_2da * SF_MIN(0.f, u[ix][ia0])           * SF_MAX(0.f, v[ix][ia0]) * dQ_i_j_i_jminus1;

    F[ix_plus][ia0]          -= dt_2da * SF_MAX(0.f, u[ix_plus][ia0])     * SF_MAX(0.f, v[ix][ia0]) * dQ_i_j_i_jminus1;


 
    return delta_Q_nplus1_i_j;

}
/* 1-D UPWIND using fluctuation form A+- and F */
/* Fluctuation FORM: Qn+1 =Qn - dt/dx (Aminus dQiplus12 + Aplus dQiminus12) - dt/dx (Fi+12 - Fi-12)
   F_n_iminus12 = F(Q_i,Q_i-1) , F_n_iplus12  =  F(Q_i+1,Q_i) */
float fluctuation_step_dz (int iserles_method, int ia, int iz, int ix, 
			   float u_tag_minus, float u_tag_plus, float u_tag_abs,
			   float Q_n_iminus1, float Q_n_i, float Q_n_iplus1,
			   float dz, float dx, 
			   //float dQ_n_iminus12_minus1, float dQ_n_iminus12, float dQ_n_iplus12,
			   int is_iminus2_padd,
			   float Q_n_iminus2, 
			   int is_iplus1_padd,
			   float u_tag)
	       
{
    const float
	//dQ_n_iminus12_minus1 = Q_n_iminus1 - Q_n_iminus2,
	dQ_n_iminus12 = Q_n_i - Q_n_iminus1,
	dQ_n_iplus12  = Q_n_iplus1 - Q_n_i;

    /* iserles_method==0: UPWIND METHOD FLUX: unstable for tg>>1, dispersing singularity*/			
    float 
	F_n_iminus12 = 0.f, 
	F_n_iplus12  = 0.f; 
  
    if (0 == iserles_method) { /* UPWIND */
	if (0 == ia && 1 == iz && 1 == ix) sf_warning(" UPWIND Fluctuation");
    }
    else
	assert (0);

    float delta_Q_nplus1_i = -dz / dx * (u_tag_minus*dQ_n_iplus12 + u_tag_plus * dQ_n_iminus12)
	- dz/dx * (F_n_iplus12 - F_n_iminus12);
 
    return delta_Q_nplus1_i;
}
/* 1-D step FLUX DIFFERENCING FORM: Qn+1 =Qn = dt/dx (Fi+12 - Fi-12)
   F_n_iminus12 = F(Q_i,Q_i-1) 
   F_n_iplus12  =  F(Q_i+1,Q_i) */
float flux_difference_step_dz (int iserles_method, //int ia, int iz, int ix, 
			       //float u_tag_minus, float u_tag_plus, float _u_tag_abs,
			       float Q_n_iminus1, float Q_n_i, float Q_n_iplus1,
			       //float _dz, float _dx, 
			       int is_iminus2_padd,
			       float Q_n_iminus2, 
			       int is_iplus1_padd,
			       float u_tag, 
			       float dz_divided_dx) //CFL_abs)
	       
{
    const float
	u_tag_plus  = SF_MAX(u_tag,0.f),
	u_tag_minus = SF_MIN(u_tag,0.f),
	u_tag_abs   = fabs(u_tag),
	CFL_abs = fabs (u_tag * dz_divided_dx), // _dz / _dx * _u_tag_abs,
	dQ_n_iminus12_minus1 = Q_n_iminus1 - Q_n_iminus2,
	dQ_n_iminus12 = Q_n_i - Q_n_iminus1,
	dQ_n_iplus12  = Q_n_iplus1 - Q_n_i;

    // iserles_method==0: UPWIND METHOD FLUX: unstable for tg>>1, dispersing singularity
    float 
	F_n_iminus12 = u_tag_minus * Q_n_i      + u_tag_plus * Q_n_iminus1, 
	F_n_iplus12  = u_tag_minus * Q_n_iplus1 + u_tag_plus * Q_n_i; 

    if (0 == iserles_method) { // UPWIND 
	;//if (0 == ia && 1 == iz && 1 == ix) sf_warning(" UPWIND ");
    }
    else { /* ifnot UPWIND */

	if (1 == iserles_method) { // Lax-Wendroff flux: 
	    //if (0 == ia && 1 == iz && 1 == ix) sf_warning(" Lax-Wendroff ");

	    F_n_iminus12 += 0.5f * u_tag_abs*(1.f - CFL_abs )*dQ_n_iminus12;
	
	    F_n_iplus12  += 0.5f * u_tag_abs*(1.f - CFL_abs )*dQ_n_iplus12;
			
	}
	else { // Piecewise Linear with slop-limiter 

	    float sigma_n_iminus1, sigma_n_i;

	    if (2 == iserles_method) { /* minmod */
		//if (0 == ia&& 1 == iz && 1 == ix) sf_warning(" minmod ");
			      			   
		sigma_n_iminus1 = minmod(dQ_n_iminus12_minus1, dQ_n_iminus12);
			      
		// padd with 1str order boundary conditions on the left
		if (is_iminus2_padd) 
		    sigma_n_iminus1 = dQ_n_iminus12;

		sigma_n_i = minmod(dQ_n_iminus12, dQ_n_iplus12);

		// padd with 1str order boundary conditions on the right
		if (is_iplus1_padd) 
		    sigma_n_i = dQ_n_iminus12;
	    }
	    else { /* superbee limiter */
			      
		//if (0 == ia && 1 == iz && 1 == ix) sf_warning(" superbee ");
		assert (3==iserles_method);
			      
		float 
		    sigma1_iminus1 = minmod (dQ_n_iminus12, 2.f*dQ_n_iminus12_minus1),

		    sigma2_iminus1 = minmod (2.f * dQ_n_iminus12, dQ_n_iminus12_minus1),

		    sigma1_i = minmod (dQ_n_iplus12, 2.f*dQ_n_iminus12),

		    sigma2_i = minmod (2.f * dQ_n_iplus12, dQ_n_iminus12);

		sigma_n_iminus1 = maxmod(sigma1_iminus1, sigma2_iminus1);

		// padd with 1str order boundary conditions on the left
		if (is_iminus2_padd) 
		    sigma_n_iminus1 = dQ_n_iminus12;

		sigma_n_i = maxmod(sigma1_i, sigma2_i);

		// padd with 1str order boundary conditions on the right
		if (is_iplus1_padd) 
		    sigma_n_i = dQ_n_iminus12;
	    }

	    const float 
		lamda_n_iminus12 = sigma_n_iminus1,

		lamda_n_iplus12 = sigma_n_i;
			    
	    F_n_iminus12 +=  0.5f * u_tag_abs*(1.f - CFL_abs  ) * lamda_n_iminus12;

	    F_n_iplus12  +=  0.5f * u_tag_abs*(1.f - CFL_abs  ) * lamda_n_iplus12;

	} 
    }

		
    const float 
	//Q_nplus1_i = Q_n_i - dz / dx * (F_n_iplus12 - F_n_iminus12);
	delta_Q_nplus1_i =   - dz_divided_dx /*dz / dx*/  * (F_n_iplus12 - F_n_iminus12);

    return delta_Q_nplus1_i;
}
/* Flux-fifferencing in direction Tz + u Tx = 0 
   float step_dz_dx_flux_difference (int iserles_method, 
   float u_tag, 
   int iz, int ix, int ia,
   int izs, int ixs, int ias, 
   int nz, int nx, int na,
   float *** t)
   {
   const float
   u_tag_plus  = SF_MAX(u_tag,0.f),
   u_tag_minus = SF_MIN(u_tag,0.f),
   u_tag_abs   = fabs(u_tag),

   Q_n_i       = t[iz-izs][ix][ia],
   Q_n_iminus1 = t[iz-izs][ix-ixs][ia];
			
   int is_iplus1_padd = 0;
   float Q_n_iplus1 = Q_n_i;
   if (ix + ixs >= 0 && ix + ixs < nx) 
   Q_n_iplus1  = t[iz-izs][ix+ixs][ia];
   else
   is_iplus1_padd = 1;

   int is_iminus2_padd = 0;
   float Q_n_iminus2 = Q_n_iminus1;
   if (ix - ixs - ixs >= 0 && ix - ixs - ixs < nx)
   Q_n_iminus2 = t[iz-izs][ix - ixs - ixs][ia];
   else
   is_iminus2_padd = 1;

   const float delta_Q_nplus1_i = 
   //fluctuation_step_dz
   flux_difference_step_dz
   (iserles_method,  //ia,  iz,  ix, 
   //u_tag_minus,  u_tag_plus,  u_tag_abs,
   Q_n_iminus1,   Q_n_i, Q_n_iplus1,
   //dz,  dx, 
   //dQ_n_iminus12_minus1,  dQ_n_iminus12,  dQ_n_iplus12,
   is_iminus2_padd,
   Q_n_iminus2, 
   is_iplus1_padd,
   u_tag,
   dz/dx,
   );

   return delta_Q_nplus1_i;
   } */
/* Flux-fifferencing in direction Tz + v Ta = 0 
   float step_dz_da_flux_difference (int iserles_method, 
   float v_tag, 
   int iz, int ix, int ia,
   int izs, int ixs, int ias, 
   int nz, int nx, int na,
   float *** t)
   {
  
   const float 
   v_tag_plus  = SF_MAX(v_tag,0.f),
   v_tag_minus = SF_MIN(v_tag,0.f),
   v_tag_abs   = fabs(v_tag);
			
   const int 
   ia_minus_ias = periodic_BC(ia - ias, na),
   ia_minus2_ias = periodic_BC (ia - ias- ias, na),
   ia_plus_ias = periodic_BC (ia + ias, na),
   is_iplus1_padd = 0,
   is_iminus2_padd = 0,
   izizs = iz - izs; // iz-izs
			
   const float
   Q_n_i       = t[izizs][ix][ia], // new_val for split method 
   Q_n_iminus1 = t[izizs][ix][ia_minus_ias],
   Q_n_iplus1  = t[izizs][ix][ia_plus_ias],
   Q_n_iminus2 = t[izizs][ix][ia_minus2_ias];

   const float delta_G_nplus1_i = 
   //fluctuation_step_dz
   flux_difference_step_dz
   (iserles_method,  //ia,  iz,  ix, 
   //v_tag_minus,  v_tag_plus,  v_tag_abs,
   Q_n_iminus1,  Q_n_i,  Q_n_iplus1,
   //dz,  da, 
   //dQ_n_iminus12_minus1,  dQ_n_iminus12,  dQ_n_iplus12,
   is_iminus2_padd,
   Q_n_iminus2, 
   is_iplus1_padd,
   v_tag,
   dz / da,
   );
   //if (1==iz && 1 == ix) sf_warning("iz = %d ix = %d G_nplus1_i=%g ",iz, ix, delta_G_nplus1_i);

   return delta_G_nplus1_i;
   }*/
//////////////////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////////////////////
//
//		      i_prev = prev_index(iz,ix,ia, is_Z_X_A_fast, nz, nx, na, is_PADD_PERIODYC),
//		      i_curr = next_index(i_prev, izx, ixs, ias, is_Z_X_A_fast, nz, nx, na, is_PADD_PERYODIC),			
//		      i_slow_prev = prev_index(iz,ix,ia, is_Z_X_A_slow, nz, nx, na, is_PADD_PERIODYC);
//		      Q_n_i       = get_t (t, iz, ix, ia, i_slow_prev, is_Z_X_A_slow, i_curr, is_Z_X_A_fast), //t[iz-izs][ix_curr][ia],
//			i_iplus1 = -1,
//			i_iplus1_padd = 
int bc_index(int i, int n, padd_type is_PADD_PERIODYC)
{
  
    if (i < 0) {
	assert (-i < n);

	if (PERIODYC == is_PADD_PERIODYC) 
	    return n + i;      

	if (PADD == is_PADD_PERIODYC) 
	    return 0;
      
	assert (0);
    }

    if (i < n-1) {
	assert (i < n + n - 1);

	if (PERIODYC == is_PADD_PERIODYC) 
	    return i - (n - 1);
      

	if (PADD == is_PADD_PERIODYC) 
	    return n - 1;
      
	assert(0);
    }
//if (i >= 0 && i < n) 
    return i;
}

//////////////////////////////////////////////////////////////////////////////

int is_beyond(int i, int nz, int nx, int na, int * i_ret, dim_type is_Z_X_A_fast, padd_type is_PADD_PERIODYC)
{
    int n = -1, ret_is_beyond = 0;

    * i_ret = i;

    switch (is_Z_X_A_fast) {
  
	case dim_Z: n = nz; assert (is_PADD_PERIODYC == PADD); break;

	case dim_X: n = nx; assert (is_PADD_PERIODYC == PADD); break;

	case dim_A: n = na; assert (is_PADD_PERIODYC == PERIODYC);  break;
    }

    if (i < 0 || i > n - 1) {
    
	ret_is_beyond = 1;

	*i_ret =  bc_index(i, n, is_PADD_PERIODYC); 
    }
 
    return ret_is_beyond;
}

///////////////////////////////////////////////////////////////////////////////

int prev_index(int iz, int ix, int ia, int izs, int ixs, int ias, dim_type is_Z_X_A_fast,int nz, int nx, int na, padd_type is_PADD_PERIODYC)
{
    int ret_i  = -1;

    switch (is_Z_X_A_fast) {

	case dim_Z: ret_i = bc_index(iz - izs, nz, is_PADD_PERIODYC); break;

	case dim_X: ret_i = bc_index(ix - ixs, nx, is_PADD_PERIODYC); break;

	case dim_A: ret_i = bc_index(ia - ias, na, is_PADD_PERIODYC); break;
    }

    return ret_i;
}

/////////////////////////////////////////////////////////////////////////////////

int next_index(int i_curr, int izs, int ixs, int ias, dim_type is_Z_X_A_fast,int nz, int nx, int na, padd_type is_PADD_PERIODYC)
{
    int ret_i  = -1;

    switch (is_Z_X_A_fast) {

	case dim_Z: ret_i = bc_index(i_curr + izs, nz, is_PADD_PERIODYC); break;

	case dim_X: ret_i = bc_index(i_curr + ixs, nx, is_PADD_PERIODYC); break;

	case dim_A: ret_i = bc_index(i_curr + ias, na, is_PADD_PERIODYC); break;
    }

    return ret_i;
}


//////////////////////////////////////////////////////////////////////////////////

float get_t (float *** t, int *** t_esc, int iz, int ix, int ia, int i_slow, dim_type is_Z_X_A_slow, int i_fast, dim_type is_Z_X_A_fast, int * ret_t_esc)//t[iz-izs][ix_curr][ia],
{
    float ret_t = -1.f;
						
    assert ((t && !ret_t_esc) || (t_esc && ret_t_esc));

    switch (is_Z_X_A_slow) {
 
	case dim_Z:
   
	    if (dim_X == is_Z_X_A_fast) {
	
		if (t)
		    ret_t = t[i_slow][i_fast][ia]; 
		else
		    * ret_t_esc = t_esc[i_slow][i_fast][ia]; 
	    }
	    else {

		assert (dim_A == is_Z_X_A_fast);

		if (t)
		    ret_t = t[i_slow][ix][i_fast];
		else 
		    *ret_t_esc = t_esc[i_slow][ix][i_fast];
	    }
  

	    break;

	case dim_X:

	    assert (dim_Z == is_Z_X_A_fast);

	    if (t)
		ret_t = t[i_fast][i_slow][ia]; 
	    else
		* ret_t_esc = t_esc[i_fast][i_slow][ia]; 

	    break; 
      
	case dim_A:

	    assert (dim_Z == is_Z_X_A_fast);

	    if (t)
		ret_t = t[i_fast][ix][i_slow];
	    else
		* ret_t_esc = t_esc[i_fast][ix][i_slow];

	    break; 
    }

    return ret_t;
}

//////////////////////////////////////////////////////////////////////////////////

int get_t_esc (int *** t_esc, int iz, int ix, int ia, int i_slow, dim_type is_Z_X_A_slow, int i_fast, dim_type is_Z_X_A_fast)//t_esc [iz-izs][ix_curr][ia],
{
    int ret_t_esc;

    (void)get_t((float***)0, t_esc, iz, ix, ia, i_slow, is_Z_X_A_slow, i_fast, is_Z_X_A_fast, & ret_t_esc);//t_esc[iz-izs][ix_curr][ia];

    return ret_t_esc;
}

///////////////////////////////////////////////////////////////////////////////////////			

int prev_index_half_lagrangian(int iz, int ix, int ia, int izs, int ixs, int ias, int is_semi_lagrangian, float u_tag, 
			       float dz, float dx, float da, 
			       int nz, int nx, int na, 
			       float coeff, dim_type is_Z_X_A_slow, dim_type is_Z_X_A_fast)
{
    assert (!is_semi_lagrangian);// for CFL>>10 not helping...

    int i_prev = -1;

    //const float coeff = 0.5f;// coeff_i = 1.f / coeff;
    if (!is_semi_lagrangian) {// CFL criteria
	  
	switch (is_Z_X_A_slow) {

	    case dim_Z: 

		if (dim_X == is_Z_X_A_fast) { // Tz + u Tx = 0

		    i_prev = ix - ixs;
		    assert (fabs(u_tag * dz) < coeff * dx);
		    assert (0 <= i_prev && nx > i_prev);
		}
		else {
		    assert (dim_A == is_Z_X_A_fast); // Tz + u Ta = 0

		    i_prev = bc_index (ia - ias, na, PERIODYC);
		    assert (fabs(u_tag * dz) < coeff * da);
		    assert (0 <= i_prev && na > i_prev);
		}
		break;

	    case dim_X:
		assert (dim_Z== is_Z_X_A_fast); // Tx + u Tz = 0

		i_prev = iz - izs;
		assert (fabs(u_tag * dx) < coeff * dz);
		assert (0 <= i_prev && nz > i_prev);
		break;

	    case dim_A: assert(0); break;
	}
    }
    else {
	int delta_is = -999999;
	switch (is_Z_X_A_slow) {

	    case dim_Z: 

		if (dim_X == is_Z_X_A_fast) { // Tz + u Tx = 0
		    // (fabs(u_tag * dz) > 1e-6 + coeff * dx) {
		    delta_is = 1 + floor( fabs(u_tag * dz/dx));

		    //if (delta_ixs > 5) sf_warning("ix_prev = %d ix=%d  u_tag=%g CFL = %g",ix_prev, ix, u_tag, u_tag * dz / dx);
		    assert (fabs(u_tag * dz) < delta_is * dx);
	      
		    i_prev = bc_index(ix - ixs * delta_is, nx, PADD);
		}
		else {
		    assert (dim_A == is_Z_X_A_fast); // Tz + u Ta = 0

		    delta_is = 1 + floor( fabs(u_tag * dz/da));

		    assert (fabs(u_tag * dz) < delta_is * da);

		    i_prev = bc_index(ia - ias * delta_is, na, PERIODYC);
		}
		break;

	    case dim_X:
		if (dim_Z == is_Z_X_A_fast) { // Tx + u Tz = 0

		    // (fabs(u_tag * dx) > 1e-6 + coeff * dz) {

		    delta_is = 1 + floor( fabs(u_tag * dx/dz));

		    //if (delta_ixs > 5) sf_warning("ix_prev = %d ix=%d  u_tag=%g CFL = %g",ix_prev, ix, u_tag, u_tag * dz / dx);
		    assert (fabs(u_tag * dx) < delta_is * dz);
	      		     
		    i_prev = bc_index(iz - izs * delta_is, nz, PADD);

		}
		else { 
		    assert (0);
		}
		break;

	    case dim_A: assert(0); break;

	}
    }
       
    assert (0 <= i_prev && nx > i_prev);

    return i_prev;		      
}
///////////////////////////////////////////////
float val_implicit_rk(int ia, int ix, int iz, 
		      float dai, float dxi, float dzi,
		      float cs, float sn, float ss,
		      float *** t,
		      float ssx, float ssz, int iq, 
		      int nz, int nx, int na)		
{
    float new_val = -1.f;
	
    if (2 == ia && 2 == ix && 2 == iz) sf_warning (" implicit R-K");

    // Downwind scheme factors / coefficient  
    const float 
	fz = -dzi*cs*ss,
	fx = -dxi*sn*ss,
	fa = -dai*(cs*ssx - sn*ssz),

	fzmin = SF_MAX(-fz,0.),
	fzmax = SF_MAX(fz,0.),

	fxmin = SF_MAX(-fx,0.),
	fxmax = SF_MAX(fx,0.),

	famin = SF_MAX(-fa,0.),
	famax = SF_MAX(fa,0.),
				    
	// Diagonal term 
	dd = (fxmax + fxmin + fzmax + fzmin + famax + famin);
	    
    // Sum of non-diagonal terms: all the non-diag terms in the right part of G-S scheme
    float ts = update_mat(t,iz,ix,ia,	fzmax,fzmin,  	fxmax,fxmin,  	famax,famin, nz, nx, na);
	    
    // Right hand side: b = S^2 for travel time,zero for Xescape,Zescape,Aescape 
    if (iq == 2) 
	ts += ss*ss;

    new_val = ts / dd;
	
    return new_val;
}
/////////////////////////////////////////////////////////////

//////////////////////////////////////////////////////////
float propagate_in_z_Tz_u_Tx(float * pQ_n_i, int iserles_method, int is_semi_lagrangian, 
			     int iz, int izs, int nz,
			     int ix, int ixs, int nx, 
			     int ia, int ias, int na,
			     float *** t, 
			     int *** t_esc, 
			     float u_tag, 
			     dim_type is_Z_X_A_slow, dim_type is_Z_X_A_fast, padd_type is_PADD_PERIODYC,
			     int debug_is_Z_slow, 
			     float dz, float dx)
{
    /*const float
      u_tag_plus  = SF_MAX(u_tag,0.f),
      u_tag_minus = SF_MIN(u_tag,0.f);
      u_tag_abs   = fabs(u_tag); 
		    
      if (0) {

      assert (is_Z_X_A_slow != is_Z_X_A_fast);

      const int is      = (is_Z_X_A_fast == dim_X ? ixs : (is_Z_X_A_fast == dim_Z ? izs : ias ));
	
      const int i_curr  = (is_Z_X_A_fast == dim_X ? ix  : (is_Z_X_A_fast == dim_Z ? iz  : ia));
	
      const int 
      i_prev =      prev_index(iz, ix, ia, izs, ixs, ias, is_Z_X_A_fast, nz, nx, na, is_PADD_PERIODYC),
			
      i_slow_prev = prev_index(iz, ix, ia, izs, ixs, ias, is_Z_X_A_slow, nz, nx, na, is_PADD_PERIODYC);

      const float 

      Q_n_i       = get_t (t, (int***)0, iz, ix, ia, i_slow_prev, is_Z_X_A_slow, i_curr, is_Z_X_A_fast, (int*)0), //t[iz-izs][ix_curr][ia],

      Q_n_iminus1 = get_t (t, (int***)0, iz, ix, ia, i_slow_prev, is_Z_X_A_slow, i_prev, is_Z_X_A_fast, (int*)0); //t[iz-izs][ix_prev][ia];


      int 
      i_iplus1 = -1,
      is_iplus1_padd = is_beyond(i_prev + is + is, nz, nx, na, & i_iplus1, is_Z_X_A_fast, is_PADD_PERIODYC);
      float 
      Q_n_iplus1 = Q_n_i;
		    
      if (!is_iplus1_padd)
      Q_n_iplus1  = get_t (t, (int***)0, iz, ix, ia, i_slow_prev, is_Z_X_A_slow, i_iplus1, is_Z_X_A_fast, (int*)0); // t[iz-izs][ix_prev + ixs + ixs][ia];

      int 
      i_iminus2 = -1,
      is_iminus2_padd = is_beyond(i_prev - is, nz, nx, na, & i_iminus2, is_Z_X_A_fast, is_PADD_PERIODYC);
      float 
      Q_n_iminus2 = Q_n_iminus1;
		    
      if (!is_iminus2_padd)
      Q_n_iminus2  =get_t (t, (int***)0, iz, ix, ia, i_slow_prev, is_Z_X_A_slow, i_iminus2, is_Z_X_A_fast, (int*)0); //t[iz-izs][ix_prev - ixs][ia];
			    

      const float ret_delta_Q_nplus1_i = 
      //fluctuation_step_dz
      flux_difference_step_dz
      (iserles_method,  ia,  iz,  ix, 
      u_tag_minus,  u_tag_plus,  u_tag_abs,
      Q_n_iminus1,   Q_n_i, Q_n_iplus1,
      dz,  dx, 
      is_iminus2_padd, Q_n_iminus2, is_iplus1_padd, u_tag
      );


      *pQ_n_i = Q_n_i;		    

      (void)get_t((float***)0, t_esc, iz, ix, ia, i_slow_prev, is_Z_X_A_slow, i_curr, is_Z_X_A_fast, & t_esc[iz][ix][ia]);//t_esc[iz-izs][ix_curr][ia];

      assert(ESC_UNINIT != get_t_esc (t_esc, iz, ix, ia, i_slow_prev, is_Z_X_A_slow, i_prev, is_Z_X_A_fast)); 
      assert(ESC_UNINIT != get_t_esc (t_esc, iz, ix, ia, i_slow_prev, is_Z_X_A_slow, i_curr, is_Z_X_A_fast)); 
      assert(ESC_UNINIT != get_t_esc (t_esc, iz, ix, ia, i_slow_prev, is_Z_X_A_slow, i_iplus1, is_Z_X_A_fast)); 
      assert(ESC_UNINIT != get_t_esc (t_esc, iz, ix, ia, i_slow_prev, is_Z_X_A_slow, i_iminus2, is_Z_X_A_fast)); 

      return ret_delta_Q_nplus1_i;
      }
    */
    /* ------------------------------- */

    if (1) {//ESC_UNINIT == t_esc[iz][ix][ia]) {

	if (debug_is_Z_slow) { 
	    const int 
		ix_prev = ix - ixs, //prev_index(ix, ixs, is_semi_lagrangian, u_tag, dz, dx, da, nz, nx, na, 1.0f),
		ix_curr = ix_prev + ixs;
			    
	    const float 
		Q_n_i       = t[iz-izs][ix_curr][ia],
		Q_n_iminus1 = t[iz-izs][ix_prev][ia];
			    
	    int 
		ix_iplus1 = -1,
		is_iplus1_padd = is_beyond_padd(ix_prev + ixs + ixs, nx,  & ix_iplus1);
	    float 
		Q_n_iplus1 = Q_n_i;

	    if (!is_iplus1_padd)
		Q_n_iplus1  = t[iz-izs][ix_prev + ixs + ixs][ia];
			    
	    int 
		ix_iminus2 = -1,
		is_iminus2_padd = is_beyond_padd(ix_prev - ixs, nx,  & ix_iminus2);
	    float 
		Q_n_iminus2 = Q_n_iminus1;
			    
	    if (!is_iminus2_padd)
		Q_n_iminus2 = t[iz-izs][ix_prev - ixs][ia];
			    
	    const float delta_Q_nplus1_i = 
		//fluctuation_step_dz
		flux_difference_step_dz
		(iserles_method,  //ia,  iz,  ix, 
		 //u_tag_minus,  u_tag_plus,  u_tag_abs,
		 Q_n_iminus1,   Q_n_i, Q_n_iplus1,
		 //dz,  dx, 
		 is_iminus2_padd, Q_n_iminus2, is_iplus1_padd, 
		 u_tag,
		 dz / dx
		    );

	    *pQ_n_i = Q_n_i;

	    t_esc[iz][ix][ia] = t_esc[iz-izs][ix_curr][ia];

	    //assert(ESC_UNINIT != t_esc[iz-izs][ix_prev][ia]);
	    //assert(ESC_UNINIT != t_esc[iz-izs][ix_curr][ia]);
	    //assert(ESC_UNINIT != t_esc[iz-izs][ix_iplus1][ia]);
			    
	    return delta_Q_nplus1_i;
	} // end if Tz + u Tx = 0
	else {
	    const int
		//iz_prev = prev_index_half_lagrangian(iz, izs, is_semi_lagrangian, u_tag, dx, dz, nz, 1.0f),
		iz_prev = iz - izs, //prev_index(iz, ix, ia, izs, ixs, ias, dim_Z , nz,  nx,  na, PADD),
		iz_curr = iz_prev + izs,
		ix_ixs = ix - ixs;

	    const float 
		Q_n_i       = t[iz_curr][ix_ixs][ia],
		Q_n_iminus1 = t[iz_prev][ix_ixs][ia];

	    int 
		iz_iplus1 = -1,
		is_iplus1_padd = is_beyond_padd(iz_prev + izs + izs, nz,  & iz_iplus1);
	    float 
		Q_n_iplus1 = Q_n_i;

	    if (!is_iplus1_padd)
		Q_n_iplus1  = t[iz_prev + izs + izs][ix_ixs][ia];
			    
	    int 
		iz_iminus2 = -1,
		is_iminus2_padd = is_beyond_padd(iz_prev - izs, nz,  & iz_iminus2);
	    float 
		Q_n_iminus2 = Q_n_iminus1;
 
	    if (!is_iminus2_padd)
		Q_n_iminus2 = t[iz_prev - izs][ix_ixs][ia];

	    const float delta_Q_nplus1_i = 
		//fluctuation_step_dz
		flux_difference_step_dz
		(iserles_method,  //ia,  iz,  ix, 
		 //u_tag_minus,  u_tag_plus,  u_tag_abs,
		 Q_n_iminus1,   Q_n_i, Q_n_iplus1,
		 //dz,  dx, 
		 is_iminus2_padd, Q_n_iminus2, is_iplus1_padd, 
		 u_tag,
		 dx / dz
		    );

	    *pQ_n_i = Q_n_i;
	    t_esc[iz][ix][ia] = t_esc[iz_curr][ix-ixs][ia];

	    //assert(ESC_UNINIT != t_esc[iz_prev] [ix - ixs][ia]);
	    //assert(ESC_UNINIT != t_esc[iz_curr] [ix - ixs][ia]);
	    //assert(ESC_UNINIT != t_esc[iz_iplus1][ix - ixs][ia]);			

	    //assert(fabs(ret_delta_Q_nplus1_i - delta_Q_nplus1_i) < 1e-20);
	    return delta_Q_nplus1_i;
	}
    }
}


//////////////////////////////////////////////////////////

const int is_semi_lagrangian = 0;
int main(int argc, char* argv[])
{
    int iz, ix, ia;
    int icv, iter;

    /* Phase space grid */
    static int nz,nx,na;
    static float oz,ox,oa;
    static float dz,dx,da;

    int iq;                        /* escape variables switch */  
    int  iserles_method_2d = 0, iserles_method = 0;
    
    int niter;               /* number of iterations */
    int  ncv;                  /* convergence loop  */

    //int ix0, ix1, ixs;    
    //int iz0, iz1, izs;
    //int ia0, ia1, ias;
    //int i;

    float dxi;                     /* inverse dx */
    float dzi;                     /* inverse dz */
    float dai;                   /* inverse da, angle */

    float ***t, ***t_a0;                    /* escape variable */
    float ***tsol=NULL;                 /* escape variable (converged) */ 
    int *** t_esc; /* flags: ESC_UP or DOWN or LEFT or RIGHT */

    float **s,**sx,**sz;           /* slowness, gradients */
    float *cvg=NULL;               /* convergence values */

    float cvt=0,tol,cvs=0;   



    float xsmooth, zsmooth;
    int   ixsmooth, izsmooth, is_xinf, is_zinf;

    char *cvgcefile = NULL;

    sf_file in,out,slow,slowz,slowx,cvgce=NULL;

    sf_init(argc,argv);

    in = sf_input("in");
    out = sf_output("out");

    cvgcefile = sf_getstring ("cvgce");
    /* output file for convergence */

    ncv = 1;
    if (cvgcefile) {
	ncv = 2;
        cvgce = sf_output (cvgcefile);
    }


    if (!sf_getint("iq",&iq)) sf_error("Need iq=");
    /* switch for escape variable 0=x, 1=a, 2=t, 3=z */

    if (!sf_getint("method",&iserles_method)) sf_error("Need method=");
    assert (iserles_method >=-1 && iserles_method < 5);

    if (!sf_getint("method_2d",&iserles_method_2d)) sf_error("Need method_2d=0 for split, 1 for full 2-D");
    assert (iserles_method_2d >=0 && iserles_method_2d <= 2);

    /* read input file parameters */
    if (SF_FLOAT != sf_gettype(in)) sf_error("Need float");

    /* angle in degrees */
    if (!sf_histint(in,"n1",&na)) sf_error("No n1=");
    if (!sf_histfloat(in,"d1",&da)) sf_error("No d1=");
    if (!sf_histfloat(in,"o1",&oa)) sf_error("No o1=");

    if (!sf_histint(in,"n2",&nx)) sf_error("No n2=");
    if (!sf_histfloat(in,"d2",&dx)) sf_error("No d2=");
    if (!sf_histfloat(in,"o2",&ox)) sf_error("No o2=");

    if (!sf_histint(in,"n3",&nz)) sf_error("No n3=");
    if (!sf_histfloat(in,"d3",&dz)) sf_error("No d3=");
    if (!sf_histfloat(in,"o3",&oz)) sf_error("No o3=");

    if (!sf_getint("niter",&niter)) niter=100;
    /* number of Gauss-Seidel iterations */

    //if (!sf_getint("ixsmooth",&ixsmooth)) sf_error("Need radius (int) ixsmooth=");
    //if (!sf_getint("izsmooth",&izsmooth)) sf_error("Need radius (int) izsmooth=");
    if (!sf_getfloat("xsmooth",&xsmooth)) sf_error("Need radius (float) xsmooth=");
    if (!sf_getfloat("zsmooth",&zsmooth)) sf_error("Need radius (float) zsmooth=");
    ixsmooth = floor (0.5 + xsmooth / dx);
    izsmooth = floor (0.5 + zsmooth / dz);

    if (!sf_getint("is_xinf",&is_xinf)) sf_error("Need is_xinf (int) is_xinf=");
    if (!sf_getint("is_zinf",&is_zinf)) sf_error("Need is_zinf (int) is_zinf=");
    
    sf_warning("ixsmooth = %d izsmooth=%d xsmooth=%g zsmooth=%g",ixsmooth,izsmooth,xsmooth,zsmooth);

    sf_warning("nz=%d nx=%d na=%d dz=%g dx=%g da=%g",nz,nx,na,dz,dx,da);

    //sf_warning("gx = %g gz = %g sc = %g ",gx,gz,sc);


    assert(!is_xinf || !is_zinf); /* infinite either in X or in Z */

    if (cvgcefile) {

        sf_putint (cvgce, "n1", niter);
        sf_putint (cvgce, "n2", 1);
        sf_putint (cvgce, "n3", 1);
        sf_putfloat (cvgce, "d1", 1.0);
        sf_putfloat (cvgce, "d2", 1.0);
        sf_putfloat (cvgce, "d3", 1.0);
        sf_putfloat (cvgce, "o1", 0.0);
        sf_putfloat (cvgce, "o2", 0.0);
        sf_putfloat (cvgce, "o3", 0.0);

        cvg = sf_floatalloc (niter);
	tsol = sf_floatalloc3(na,nx,nz);

	for (iz = 0; iz < nz; iz ++) 
	    for (ix = 0; ix < nx; ix ++) 
		for (ia = 0; ia < na; ia++) 
		    tsol[iz][ix][ia] = 0.0f;

    }

    if (!sf_getfloat("tol",&tol)) tol=0.000002*nx*nz;
    /* accuracy tolerance */

    /* memory allocations */

    /* T - THE SOLUTION VECTOR for U(t (z,x,a) ) = S(x,z)^2 */
    t = sf_floatalloc3(na,nx,nz);
    t_esc = sf_intalloc3(na, nx, nz);
    t_a0 = sf_floatalloc3(na,nx,nz);

    for (iz = 0; iz < nz; iz ++) {
	for (ix = 0; ix < nx; ix ++) {
	    for (ia = 0; ia < na; ia++) {
		t_esc[iz][ix][ia] = ESC_UNINIT;
		t    [iz][ix][ia] = 0.; //1e10f;
		t_a0 [iz][ix][ia] = 0.; //1e10f;
	    }
	}		
    }

    /* read input escape variable - it is the inistability tial guess */
    sf_floatread(t[0][0],na*nx*nz,in);

    /* read auxiliary slowness file */
    slow = sf_input("slow");
    s = sf_floatalloc2(nz,nx);
    sf_floatread(s[0],nz*nx,slow);

    /* read auxiliary slowness z-gradient file */
    slowz = sf_input("slowz");
    sz = sf_floatalloc2(nz,nx);
    sf_floatread(sz[0],nz*nx,slowz);

    /* read auxiliary slowness x-gradient file */
    slowx = sf_input("slowx");
    sx = sf_floatalloc2(nz,nx);
    sf_floatread(sx[0],nz*nx,slowx);

    /* convert to radians */
    oa *= SF_PI/180.;
    da *= SF_PI/180.;

    dzi = 1.0/dz;
    dxi = 1.0/dx;
    dai = 1.0/da;

    for(icv=0; icv < ncv; icv++) {
	/* first - find the solution
	   after that - compute the intricateerrors |X_iter - X_solution|
	*/

        int ia0=0, ia1=na, ias=1;

	for (iter=0; iter < niter; iter++) {

	    for (ia = ia0; ia != ia1; ia += ias) {

		int ix0, ix1, ixs, iz0, iz1, izs;

		const float 
		    a = oa + ia*da,
		    cs = cosf(a),
		    sn = sinf(a),
		    px = -sn,
		    pz = -cs;

		assert (fabs(px) > 1e-6 && fabs(pz) > 1e-6);
		
		boundary_sweep_setup( iq,
				      px,  pz,  nx,  nz,  ia,  
				      s, sx, sz, 
				      is_xinf,  is_zinf,
				      &ix0, &ix1, &ixs, 
				      &iz0, &iz1, &izs, 
				      t_a0, 
				      t_esc,
				      oz,  ox,  oa,  dz,  dx,  da);

		t_esc[iz0][ix0][ia] = ESC_UP; // TBD
		
		cvt = 0.0; /* X_next - X_prev */
		cvs = 0.0; /* X_exact - X_iter */

		if (fabs(sn/cs * dz)  >= dx) { 
		    // propagate in X(slow) Z (fast) T_x + 1/u_tag T_z = S^2 / u_tag

		    int ix_smooth_edge_pos = 0,
			iz_smooth_edge_pos = 0;

		    /* the iteration: loop over the grid Z x X x A */
		    for (ix = ix0; ix != ix1; ix += ixs) {	   
 
			iz_smooth_edge_pos = ix_smooth (ix_smooth_edge_pos, izsmooth, ixsmooth, nz);

			ix_smooth_edge_pos++;		    
			      
			for (iz = iz0  + izs * iz_smooth_edge_pos; iz != iz1; iz += izs) {

			    float new_val = 0.f;

			    const float 		
				ss = s[ix - ixs][iz],
				ssx = sx[ix - ixs][iz],
				ssz = sz[ix - ixs][iz],
				u_tag = cs/sn * (float)(ixs * izs);

			    if (-1 == iserles_method) { // implicit R-K 1st order
				  
				new_val = val_implicit_rk(ia, ix, iz, dai, dxi, dzi, cs, sn, ss, t_a0, ssx, ssz, iq,
							  nz,  nx,  na);
				  
			    }
			    else { // explicit
						
				assert (0 == iserles_method_2d);
  			  
				assert (u_tag > 1e-6 && fabs(u_tag * dx / dz) < 1.f + 1e-6);
			  
				if (u_tag > 1e-6) { //Tx + 1/u_tag Tz = s^2 / u_tag

				    float Q_n_i = -1.f;

				    const float delta_Q_nplus1_i = 

					propagate_in_z_Tz_u_Tx (& Q_n_i, iserles_method, is_semi_lagrangian, 
								iz, izs, nz,
								ix, ixs, nx, 
								ia, ias, na,
								t_a0, 
								t_esc, 
								u_tag, 
								dim_X, dim_Z, PADD,
								0,
								dz, dx
					    ); // is_Z_slow = 0 => X_slow
				    //is_Z_X_A_slow, is_Z_X_A_fast, is_PADD_PERIODYC);

				    new_val = Q_n_i + delta_Q_nplus1_i;
				}
				if (iq == 2) {
				    new_val += ss / sn * dx * (float)ixs;
				}						
			    } // end explicit
	
			      //cvt += fabsf(t[iz][ix][ia]-new_val);

			    t_a0[iz][ix][ia] = new_val;
			      
			    //if (icv == 1) 
			    //cvs += fabsf(tsol[iz][ix][ia]-new_val);

			} // end iz
		    } // end ix
		} // end T_x + ui T_z = s^2/ ui
	
		else { // propagate Z (slow), X (fast)

		    int ix_smooth_edge_pos,
			iz_smooth_edge_pos = 0;

		    /* the iteration: loop over the grid Z x X x A */
		    for (iz = iz0; iz != iz1; iz += izs) {	   
 
			ix_smooth_edge_pos = ix_smooth (iz_smooth_edge_pos, ixsmooth, izsmooth, nx);

			iz_smooth_edge_pos++;		    
		      
			for (ix = ix0  + ixs * ix_smooth_edge_pos; ix != ix1; ix += ixs) {

			    float new_val = 0.f;

			    const float 		
				ss = s[ix][iz - izs],
				ssx = sx[ix][iz - izs],
				ssz = sz[ix][iz - izs],
				u_tag = sn/cs * (float)(ixs * izs);

			    if (-1 == iserles_method)// || fabs(u_tag * dz / dx) > 1.0f || fabs(v_tag*dz/da)>1.0f ) 
			    { // implicit R-K 1st order
			  
				new_val = val_implicit_rk(ia, ix, iz, dai, dxi, dzi, cs, sn, ss, t_a0, ssx, ssz, iq,
							  nz,  nx,  na);

			    }
			    else { // explicit
						
				assert (0 == iserles_method_2d);
			
				assert (u_tag > 1e-6 && u_tag * dz/dx  < 1e-6 + 1.0f);
			  
				if (u_tag > 1e-6) { //Tz + u_tag Tx = s^2

				    float Q_n_i = -1.f;
				    const float delta_Q_nplus1_i = 
					propagate_in_z_Tz_u_Tx (& Q_n_i, iserles_method, is_semi_lagrangian, 
								iz, izs, nz,
								ix, ixs, nx, 
								ia, ias, na,
								t_a0, 
								t_esc, u_tag, 
								dim_Z, dim_X, PADD,
								1, 
								dz, dx);

//if (5 == iz) sf_warning(" ix = %d,  plus =%d, curr =%d, prev=%d, minus2=%d", ix, ix_iplus1, ix_curr, ix_prev, ix_iminus2);
//if (5 == iz) sf_warning(" t_ix = %d,  t_plus =%d, t_curr =%d, t_prev=%d, t_minus2=%d", t_esc[iz-izs][ix][ia], t_esc[iz-izs][ix_iplus1][ia], t_esc[iz-izs][ix_curr][ia], t_esc[iz-izs][ix_prev][ia], t_esc[iz-izs][ix_iminus2][ia]);

				    new_val = Q_n_i + delta_Q_nplus1_i;

				    if (iq == 2) {
					new_val += ss / cs * dz * (float)izs;
				    }

				} // end if u != 0
			  
			    } // end explicit
			    /////////////////////////////////////////////
			    /* Right hand side: b = S^2 for travel time,zero for Xescape,Zescape,Aescape */
			
			    //cvt += fabsf(t[iz][ix][ia]-new_val);

			    t_a0[iz][ix][ia] = new_val;
		      
			    //if (icv == 1) 
			    //	  cvs += fabsf(tsol[iz][ix][ia]-new_val);
			
			} /* ix */
		    
		    } /* iz */

		} // end of propagate Xslow-Zfast vs. Zslow-Xfast
	    } /* ia */

	    /* the iteration: loop over the grid Z x X x A */
	    const int ixs = 1;
	    for (ix = 0; ix < nx; ix++) {	   
		const int izs = 1;
		for (iz = 0; iz < nz; iz++) {

		    const int ias = 1;
		    for (ia = 0; ia < na; ia++) {
			    
			const float 
			    a = oa + ia*da,
			    cs = cosf(a),
			    sn = sinf(a),
			    ss = s[ix][iz],
			    ssx = sx[ix][iz],
			    ssz = sz[ix][iz],
			    fa = (cs*ssx - sn*ssz),
			    v_tag_z = fa / (ss * cs) * (float)(ias * izs),
			    v_tag_x = fa / (ss * sn) * (float)(ias * ixs);

			//assert (fabs(v_tag_z * dz / da) < 1.f || fabs(v_tag_x *dx/da) < 1.f);
			if (fabs(v_tag_z * dz / da) > 1.f)
			    if (fabs(v_tag_x *dx/da) > 1.f)
				;//sf_warning(" cfl1 = %g,  cfl2 = %g", fabs(v_tag_z * dz / da) , fabs(v_tag_x *dx/da));

			if (-1 == iserles_method) {

			    cvt += fabsf(t[iz][ix][ia]-t_a0[iz][ix][ia]);

			    t[iz][ix][ia] = t_a0[iz][ix][ia];

			}
			else { // explicit , NOT implicit R-K 1st order
  			  							
			    if (fabs(fa) < 1e-6) {

				cvt += fabsf(t[iz][ix][ia]-t_a0[iz][ix][ia]);

				t[iz][ix][ia] = t_a0[iz][ix][ia];
			    }
			    else { // fa != 0
				// Now Tz dz + fa/fz Ta da = 0 
				float 
				    dzx_da = dz / da, 
				    v_tag = v_tag_z;

				if (fabs(v_tag_z * dz / da) >= 1.f) {

				    //assert(fabs(v_tag_x * dx / da) <=1.f);

				    // Now Tx dx + fa/fx Ta da = 0 
				    dzx_da = dx / da;

				    v_tag = v_tag_x;
				}

				{
				    int ia_prev = bc_index(ia - ias, na, PERIODYC); 
	    
				    const int 
					ia_minus2_ias = periodic_BC (ia_prev - ias, na),
					ia_plus_ias = periodic_BC (ia_prev + ias + ias, na),
					is_iplus1_padd = 0,
					is_iminus2_padd = 0;

				    const float
					Q_n_i       = t_a0[iz][ix][ia], 
					Q_n_iminus1 = t_a0[iz][ix][ia_prev],
					Q_n_iplus1  = t_a0[iz][ix][ia_plus_ias],
					Q_n_iminus2 = t_a0[iz][ix][ia_minus2_ias];
				    
				    const float delta_G_nplus1_i = 
					//fluctuation_step_dz
					flux_difference_step_dz
					(iserles_method,
					 Q_n_iminus1,  Q_n_i,  Q_n_iplus1,
					 is_iminus2_padd,
					 Q_n_iminus2, 
					 is_iplus1_padd,
					 v_tag,
					 dzx_da
					    );

				    const float new_val = Q_n_i + delta_G_nplus1_i;
					
				    //if (1e-6 < fabs(delta_G_nplus1_i) ) //if (1==iz && 1 == ix)  sf_warning("iz = %d ix = %d G_nplus1_i=%g ",iz, ix, delta_G_nplus1_i);
					
				    cvt += fabsf(t[iz][ix][ia]-new_val);
				
				    t[iz][ix][ia] = new_val;
					
				} // end CFS < 1 in Tz + v_tag Ta = 0	
				/*else { 
				// TBD: bad results for  Ta da + fz/fa Tz dz = 0 
				const float v = 1.f / v_tag;
					
				assert (fabs(v * da / dz) < 1.f);
			    
				const int 
				iz_prev = bc_index(iz - izs, nz, PADD),
				iaias = ia; // bc_index(ia - ias, na, PERIODYC); 
			    
				const float
				Q_n_i       = t_a0[iz][ix][iaias],
				Q_n_iminus1 = t_a0[iz_prev][ix][iaias];

				int 
				iz_iplus1 = -1,
				is_iplus1_padd = is_beyond_padd(iz_prev + izs + izs, nz,  & iz_iplus1);
				float 
				Q_n_iplus1 = Q_n_i;

				if (!is_iplus1_padd)
				Q_n_iplus1  = t_a0[iz_prev + izs + izs][ix][iaias];
			    
				int 
				iz_iminus2 = -1,
				is_iminus2_padd = is_beyond_padd(iz_prev - izs, nz,  & iz_iminus2);
				float 
				Q_n_iminus2 = Q_n_iminus1;
 
				if (!is_iminus2_padd)
				Q_n_iminus2 = t_a0[iz_prev - izs][ix][iaias];
				    
				const float delta_G_nplus1_i = 
				flux_difference_step_dz(iserles_method,  
				Q_n_iminus1,  Q_n_i,  Q_n_iplus1,
				is_iminus2_padd,
				Q_n_iminus2, 
				is_iplus1_padd,
				v,
				da / dz);			

				const float new_val = Q_n_i + delta_G_nplus1_i;

				cvt += fabsf(t[iz][ix][ia] - new_val);

				t[iz][ix][ia] = new_val;

				} // end CFL > 1 
				*/
			    } // end if fa != 0 

			} // end if explicit	
			    
			if (icv == 1) 
			    cvs += fabsf(tsol[iz][ix][ia] - t[iz][ix][ia]);
			    
		    } // end for a

		} // end for z

	    } // ix

	    sf_warning("Iter = %d, Norm L1 (tprev - t) = %g (%g normlzed)",iter,cvt, cvt / (float)(nz*nx*na));

            /* tol - tolerance for convergence */
	    //if (cvt < tol) 
	    //break;  

	    if (0 == ia0) { 
		ia0 = na-1; ia1 = -1; ias = -1; 
	    }
	    else {
		ia0 = 0; ia1 = na; ias = 1;
	    }


	    if (1 == icv)
		cvg[iter] = cvs;	

	} // end for iter
	
	if (cvgcefile && icv == 0) {
	    memcpy (tsol, t, nz*nx*na*sizeof (float) );
	    /* read input escape variable - it is the initial guess */
	    sf_floatread(t[0][0],na*nx*nz,in);
	}
	
    } // end for icv

    /* output */
    sf_floatwrite(t[0][0],na*nx*nz,out);

    /* output convergence */
    if (cvgcefile) 
	sf_floatwrite (cvg,niter,cvgce);
    
    exit(0);
}


