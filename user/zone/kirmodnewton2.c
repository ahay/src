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

#include <stdlib.h>
#include <float.h>
#include <rsf.h>

#include "kirmod2.h"
#include "kirmodnewton.h"
/*^*/

#ifndef _kirmodnewton2_h


typedef struct Velocity2 {
	float *v, *gx, *gz, *xref, *zref, *thick, *sumthick, **aniso;
} velocity2;
/*^*/

#endif

struct Surface {
    int is, ih;
    float x;
    ktable **ta;
};


static int ny, nc, nx, **map, vstatus;
static float x0, dx, r0, dr;
static velocity2 v;

/* Reflector function----------------------------------------------------------------------------------*/

static sf_eno *eno, *deno; /* Interpolation structure */	

static float z(int k,float x) 
/*Function*/
{
	int i;
	float f, f1;
	
	x = (x-r0)/dr; 
	i = floorf(x);
	
	sf_eno_apply (eno[k],i,x-i,&f,&f1,FUNC);
	return f;
}

static float zder(int k,float x) 
/*First derivative*/
{
	int i;
	float f, f1;
	
	x = (x-r0)/dr; 
	i = floorf(x);
	
	sf_eno_apply (deno[k],i,x-i,&f,&f1,FUNC);
	return f;
}

static float zder2(int k,float x) 
/*Second derivative*/
{
	int i;
	float f, f1;
	
	x = (x-r0)/dr; 
	i = floorf(x);
	
	sf_eno_apply (deno[k],i,x-i,&f,&f1,DER);
	return f1/dr;	
}

void kirmodnewton_init(float **temp_rr /* Reflectors data of dimension N2xN1 */,
					   float **temp_rd /* Dip data*/,
					   int *updown /* Direction of the ray */,
					   float r01 /* Origin of reflectors */,
					   float dr1 /* Increment between elements in reflectors */,
					   int N1 /* Number of elements in each reflector */,
					   int n /* Number of reflections */,
					   int order /* Order of eno interpolation */,
					   int N2 /* Number of reflectors including surface*/,
					   int vstatus1 /* Type of model (vconstant(0) or vgradient(1) or VTI(2))*/,
					   float *xref_inp /* x-coordinate reference points */,
					   float *zref_inp /* z-coorditnate reference points */,
					   float *v_inp /* Velocities at reference points */,
					   float *gx_inp /* x-gradient at the reference points */,
					   float *gz_inp /* z-gradeint at the reference points */,
					   float **aniso_input /* anisotropy parameters*/)
/*<Initialize reflectors for kirmodnewton>*/

{
	float **rr, **rd;
	int ir2;
	int d1,d2,d2n,d3,d4,d5,d6,p3; /*counter*/
	/*int p1,p2; Temp value*/
	float p4=0; /*Temp value*/
	int ithick;
	
	rr = sf_floatalloc2(N1,n+2); /* Reflector values according to updown*/
	rd = sf_floatalloc2(N1,n+2); /* Slope values according to updown*/
	
	/*At each point*/
	v.v = sf_floatalloc(n+1);  /* Velocity array used in calculation generated according to where the ray travels*/
	v.gx = sf_floatalloc(n+1); /* Velocity gradient in x-direction used in calculation generated according to where the ray travels*/
	v.gz = sf_floatalloc(n+1); /* Velocity gradient in z-direction used in calculation generated according to where the ray travels*/
	v.xref = sf_floatalloc(n+1); /* Reference point x-coordinate used in calculation generated according to where the ray travels*/
	v.zref = sf_floatalloc(n+1); /* Reference point z-coordinate used in calculation generated according to where the ray travels*/
	v.thick = sf_floatalloc(n+1); /*Avg thickness of each layer for xintial*/
	v.sumthick = sf_floatalloc(n+1); /*Avg thickness of each layer for xintial*/
	v.aniso = sf_floatalloc2(4,n+1); /* Anisotropy parameters*/

	r0 = r01;
	dr = dr1;
	vstatus = vstatus1;

	/* Check the array, consecutive two inputs must not differ by more than 1-----------------------------*/
	
	/* Allow skipping of layer if the thickness is zero*/
	/*for (p1=0; p1<n-1; p1++) {
		p2 = updown[p1]-updown[p1+1];
		if (p2>1) {
			sf_warning("The layer number array indicates skipping of layers. Please reenter the array.\n");
			exit(0);
		}
	}*/
	if (vstatus == 0){
			int index;
		for(index=0;index<N2-1;index++) {
			gx_inp[index] = 0.0;
			gz_inp[index] = 0.0;
		}
	}
	/* Check whether the gradient and vstatus match-------------------------------------------------------*/
	if (vstatus != 2) {
		for (p3=0; p3<N2-1; p3++) {
			p4 = p4 + fabs(gx_inp[p3])+fabs(gz_inp[p3]);
		
			if (p3==N2-2 && p4/(2*N2-2)!=0 && vstatus==0) {
				sf_warning("The gradients are not zero. Please reenter nonzero vstatus.\n");
				exit(0);
			}
			if (p3==N2-2 && p4/(2*N2-2)==0 && vstatus!=0) {
				sf_warning("The gradients are zero. Please enter vstatus=0 for constant velocity model.\n");
				exit(0);
			}
		}
	}
	/* Generate layer input according to the reflection sequence-----------------------------------------------*/
	for (d2=0; d2<n+1; d2++) {
		
		/* Set velocity, gradient, and reference points arrays-------------------------------------------*/
		if (d2<1) { 
			v.v[d2] = v_inp[0];
			v.gx[d2] = gx_inp[0];
			v.gz[d2] = gz_inp[0];
			v.xref[d2] = xref_inp[0];
			v.zref[d2] = zref_inp[0];
			for(d6=0; d6<4; d6++){
				v.aniso[d2][d6] = aniso_input[0][d6];
			}
		}
		else {
			d3 = updown[d2-1]; /* Need d3, d4, and d5 because array argument needs to be an interger*/
			d4 = updown[d2];
			
			if (d4-d3>0) {
				v.v[d2] = v_inp[d3];
				v.gx[d2] = gx_inp[d3];
				v.gz[d2] = gz_inp[d3];
				v.xref[d2] = xref_inp[d3];
				v.zref[d2] = zref_inp[d3];
				for(d6=0; d6<4; d6++){
					v.aniso[d2][d6] = aniso_input[d3][d6];
				}
			}	
			if(d4-d3<0){
				v.v[d2] = v_inp[d4];
				v.gx[d2] = gx_inp[d4];
				v.gz[d2] = gz_inp[d4];
				v.xref[d2] = xref_inp[d4];
				v.zref[d2] = zref_inp[d4];
				for(d6=0; d6<4; d6++){
					v.aniso[d2][d6] = aniso_input[d4][d6];
				}
			}
		}
	}
	/* Generate point input according to the reflection sequence-----------------------------------------------*/
	for (d2n=0; d2n<n+2; d2n++) {	
		for (d1=0; d1<N1; d1++) { /* Set layers according to updown*/
			
			if (d2n == 0) {
				rr[d2n][d1] = temp_rr[0][d1];
				rd[d2n][d1] = temp_rd[0][d1];
			}
			else {
			d5 = updown[d2n-1];
			rr[d2n][d1] = temp_rr[d5][d1];
			rd[d2n][d1] = temp_rd[d5][d1];
			}	
		}
		
	}
	for(ithick = 0; ithick < n+1; ithick++){ /*To calculate the average thickness of each layer measured from both ends for xinitial*/
		v.thick[ithick] = ((rr[ithick+1][0] - rr[ithick][0]) + (rr[ithick+1][N1-1] - rr[ithick][N1-1]))/2;
		if (ithick==0){
			v.sumthick[ithick] = v.thick[ithick];
		}
		else {
			v.sumthick[ithick] = v.sumthick[ithick-1] + v.thick[ithick] ;
		}
	}
	
	/* Initialize interpolation------------------------------------------------------------------------*/
	
	eno = (sf_eno*) sf_alloc(n+2,sizeof(*eno)); /* Allocation for eno array*/
	deno = (sf_eno*) sf_alloc(n+2,sizeof(*deno));
	
	/* Compute reflector slope---------------------------------------------------------------------*/
	
	for (ir2=0; ir2 < n+2; ir2++) { /* Loop through eno*/
		eno[ir2]  = sf_eno_init(order,N1); /* Function values*/
		deno[ir2]  = sf_eno_init(order,N1); /* Dip values*/
		sf_eno_set (eno[ir2],rr[ir2]); 
		sf_eno_set (deno[ir2],rd[ir2]);
		
	}
	
}

static int surface_comp(const void *a, const void *b)
/* compare by the surface coordinate same as in kirmod2.c*/
{
    float aa, bb;
	
    aa = ((surface) a)->x;
    bb = ((surface) b)->x;
	
    if (aa <  bb) return (-1);
    if (aa == bb) return 0;
    return 1;
}

surface kirmodnewton2_init(int ns,  float s0,  float ds  /* source/midpoint axis */,
					 int nh,  float h0,  float dh  /* offset axis */,
					 int nx1, float x01, float dx1 /* reflector axis */,
					 int nc1                       /* number of reflectors */,
                     bool cmp                      /* if CMP instead of shot gather */,
                     bool absoff                   /* use absolute offset */)
/*< Initialize surface locations same as in kirmod2 >*/ 
{
    int is, ih, iy;
    float s, h;
    surface yi, y;
	
    nx = nx1;
    dx = dx1;
    x0 = x01;
    nc = nc1;
	
    if (cmp) {
		ny = 2*ns*nh; 
		map = sf_intalloc2(2*nh,ns);
    } else {
		ny = ns*(nh+1); 
		map = sf_intalloc2(nh+1,ns);
    }	
	
    y = (surface) sf_alloc(ny,sizeof(*y));
	
    yi = y;
	
    for (is=0; is < ns; is++) {
		s = s0 + is*ds;
		if (cmp) {
			for (ih=0; ih < nh; ih++, yi++) {
				h = 0.5*(h0 + ih*dh);
				yi->x = s - h;
				yi->is = is;
				yi->ih = 2*ih;
				yi++;
				yi->x = s + h;
				yi->is = is;
				yi->ih = 2*ih+1;
			}
		} else {
			for (ih=0; ih < nh; ih++, yi++) {
				yi->x = absoff ? h0 + ih*dh : s + h0 + ih*dh;
				yi->is = is;
				yi->ih = ih;
			}
			yi->x = s;
			yi->is = is;
			yi->ih = nh;
			yi++;
		}
    }
	
    qsort(y,ny,sizeof(*y),surface_comp);
	
    for (iy=0; iy < ny; iy++) {
		yi = y+iy;
		map[yi->is][yi->ih] = iy;
    }
	
    return y;
}

void kirmodnewton2_table(surface y /* Surface structure*/,
				   bool debug /* Debug Newton */,
				   bool fwdxini /* Use the result of the previous iteration for the next one*/,
				   int niter /* Number of iteration for Newton */,
				   double tolerance /* Tolerance level for Newton */)
/*<Compute traveltime map>*/
{
	
	int ix, iy, ic, iv, num;
	float x2, x1, xp=0.;
	float *xinitial, **oldans;
	bool skip=false;
	ktable **ta=NULL;
	
	if (fwdxini && nc==1) {
		sf_warning("Cannot fwdxini in the case of 1 reflector (nc=1). Please enter fwdxini=n");
		exit(0);
	}
	
	if (fwdxini) {
	    oldans = sf_floatalloc2(nc-1,nc-1); /* To store old ans for the case of fwdxini*/
	} else {
	    oldans = NULL;
	}
	
    for (iy=0; iy < ny; iy++) {	/* source/midpoint and offset axes */
		x1 = y[iy].x; /* x1 is on the surface */
		if (0==iy || x1 != xp) { /* new point */
			ta = (ktable**) sf_alloc(nx,sizeof(ktable*));
			
			for (ix=0; ix < nx; ix++) { /* reflector axis */
				ta[ix] = (ktable*) sf_alloc(nc,sizeof(ktable));
				
				x2 = x0 + ix*dx; /* x2 is on the reflector */
				
				for (ic=0; ic < nc; ic++) { 
					ta[ix][ic] = (ktable) sf_alloc(1,sizeof(ta[ix][ic][0]));
					
					/* Additional loop from kirmod2.c for generating extra parameters*/
					num = ic;
					if (num!=0) {
						xinitial = sf_floatalloc(num);
						
						for (iv=0; iv<num; iv++) { /* How many reflection*/
							if (fwdxini) {
									
								if (ix==0 || skip) { /* Initialize this old result array for the very first calculation or the previous ray can't be traced*/
									/*oldans[iv][num-1] = x1+(x2-x1)*(iv+1)/(num+1);*/
									oldans[iv][num-1] = x1+(x2-x1)*v.sumthick[iv]/(v.sumthick[num]);
								}
								
								xinitial[iv] = oldans[iv][num-1]; /* 2nd axis is initial points axis*/
							
							}
							else {
								/*xinitial[iv] = x1+(x2-x1)*(iv+1)/(num+1);*/
								xinitial[iv] = x1+(x2-x1)*v.sumthick[iv]/(v.sumthick[num]);
							}
						}
					}
					else {
						xinitial = sf_floatalloc(1);
						xinitial[0] = 0; /* The case where there is NO reflection (to avoid warning)*/
					}

					
					kirmodnewton_table(vstatus, debug, x1, x2, x1, x2, niter, tolerance, num, xinitial, v.xref, v.zref,v.v, v.gx, v.gz,v.aniso,z, zder, zder2, oldans, skip, ta[ix][ic]);
				
				} 
			} 
		} 
		y[iy].ta = ta;
		xp = x1;
    }

}

ktable kirmodnewton2_map(surface y, int is, int ih, int ix, int ic) 
/*< Extract from traveltime/amplitude map >*/
{
    return y[map[is][ih]].ta[ix][ic];
}

