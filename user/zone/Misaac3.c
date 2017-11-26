/* 3D Bending ray tracing in Multi-layered media*/
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
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "vectorops.h"
#include "general_traveltime_3D.h"
#include "ml_traveltime_vconstant_3D.h"
#include "ml_traveltime_vgradient_3D.h"
#include "ml_traveltime_vti_3D.h"
#include "setvelocity_3D.h"
#include "matrixops_2D.h"


/* Reflector function--------------------------------------------------------------------------------------------------*/

static sf_eno2 *eno, *d1eno, *d2eno; /* Interpolation structure */	
static float r0, dr,r1,dr1,r2,dr2;

static float z(int k,float x, float y) 
/* Function */
{
	int i,j;
	float f, f1[2];
	
	x = (x-r0)/dr; 
	i = floorf(x);
	
	y = (y-r1)/dr1; 
	j = floorf(y);
	
	sf_eno2_apply (eno[k],i,j,x-i,y-j,&f,f1,FUNC);
	return f;
}

static float zder(int k,float x, float y, int m /* m=0 for x-direction and m=1 for y-direction*/) 
/* First derivative with respect to x_k or y_k*/
{
	int i,j;
	float f,f1[2];
	
	x = (x-r0)/dr; 
	i = floorf(x);
	
	y = (y-r1)/dr1; 
	j = floorf(y);
	
	if (m==1) {
		dr=dr1; /*Change dr to dr1 if we want the derivative with respect to y_k*/
	}
	
	sf_eno2_apply (eno[k],i,j,x-i,y-j,&f,f1,DER);
	return f1[m]/dr;
}

static float zder2_1(int k,float x, float y,int m /* m=0 for x-direction and m=1 for y-direction*/) 
/* Second derivative with respect to x_k*/
{
	int i,j;
	float f, f1[2];
	
	x = (x-r0)/dr; 
	i = floorf(x);
	
	y = (y-r1)/dr1; 
	j = floorf(y);
	
	if (m==1) {
		dr=dr1; /* Change dr to dr1 if we want the derivative with respect to y_k*/
	}
	
	sf_eno2_apply (d1eno[k],i,j,x-i,y-j,&f,f1,DER);
	return f1[m]/dr;	
}

static float zder2_2(int k,float x, float y, int m /* m=0 for x-direction and m=1 for y-direction*/) 
/* Second derivative with respect to y_k*/
{
	int i,j;
	float f, f1[2];
	
	x = (x-r0)/dr; 
	i = floorf(x);
	
	y = (y-r1)/dr1; 
	j = floorf(y);
	
	if (m==1) {
		dr=dr1; /* Change dr to dr1 if we want the derivative with respect to y_k*/
	}
	
	sf_eno2_apply (d2eno[k],i,j,x-i,y-j,&f,f1,DER);
	return f1[m]/dr;	
}


/* Main program--------------------------------------------------------------------------------------------------------*/
int main(int argc, char* argv[])
{
	int nr1, nr2, nr3, N, ir1, ir2, ir3, nt, nt2,nt3, order, niter, vstatus, count, q=0/* Counter for the bog loop*/;
	float x, y, dt, t0, bmin, bmax, xbmin, xbmax, ybmin, ybmax, tt;
	float ***rr, ***temp_rr, ***rd1, ***rd2, **ans, **xx, **xxnew, *xinitial, *yinitial, *updown, *v_inp, *gx_inp, *gy_inp, *gz_inp, *xref_inp, *yref_inp, *zref_inp ,*v, *gx, *gy, *gz, *xref, *yref, *zref, **F, **dk, **xxtem, **zk,***ck_inv, **aniso, **aniso_inp; 
	float **t1_temp,**t2_temp,**t3_temp,**t4_temp,*t5_temp,*t6_temp,**t7_temp,*t8_temp,*t9_temp;
	double tol;
	bool debug;
	sf_file refl, xrefl, vti;
	
	sf_init(argc,argv); /* Initialize - always call first */
	
	/* Set input------------------------------------------------------------------------------*/
	refl = sf_input("in"); /* reflector */
	if (!sf_histint(refl,"n1",&nr1)) sf_error("No n1= in input");/* Number of points in x-direction*/
	if (!sf_histint(refl,"n2",&nr2)) sf_error("No n2= in input");/* Number of points in y-direction*/
	if (!sf_histint(refl,"n3",&N)) sf_error("No n3= in input"); /* Number of reflector including the surface at z=0*/
	
	if (!sf_histfloat(refl,"o1",&r0)) r0=0.;
	if (!sf_histfloat(refl,"o2",&r1)) r1=0.;
	if (!sf_histfloat(refl,"o3",&r2)) r2=0.;
	
	if (!sf_histfloat(refl,"d1",&dr)) dr=1.;
	if (!sf_histfloat(refl,"d2",&dr1)) dr1=1.;
	if (!sf_histfloat(refl,"d3",&dr2)) dr2=1.;
	
	
	if (!sf_getint("number",&nr3)) sf_error("Please enter the number of reflections [nr3]");
	/* Number of reflectors*/
	
	/* Allocate space--------------------------------------------------------------------------*/
	rr = sf_floatalloc3(nr1,nr2,nr3+2); /* Reflector according to the updown*/
	temp_rr = sf_floatalloc3(nr1,nr2,N); /* Reflector input*/
	rd1 = sf_floatalloc3(nr1,nr2,nr3+2); /* Reflector slope for x-direction*/
	rd2 = sf_floatalloc3(nr1,nr2,nr3+2); /* Reflector slope for y-direction*/
	ans = sf_floatalloc2(3,nr3+2); /* Answer reflection point for 2 axes (x,y,z)*/
	aniso_inp = sf_floatalloc2(4,N-1); /* VTI parameters of the model*/
	aniso = sf_floatalloc2(4,nr3+2); /* VTI parameters of the model*/
	
	xx = sf_floatalloc2(2,nr3+2); /* Reflection position at all interfaces*/
	xxnew = sf_floatalloc2(2,nr3+2);
	
	if (nr3!=0) {
		xinitial = sf_floatalloc(nr3); /* Initial guess of the relfection position*/
		yinitial = sf_floatalloc(nr3);
	}
	updown = sf_floatalloc(nr3+1);
	
	v_inp = sf_floatalloc(N-1); /* Input velocity array*/
	gx_inp = sf_floatalloc(N-1); /* Input velocity gradient in x-direction*/
	gy_inp = sf_floatalloc(N-1); /* Input velocity gradient in x-direction*/
	gz_inp = sf_floatalloc(N-1); /* Input velocity gradient in z-direction*/
	xref_inp = sf_floatalloc(N-1); /* Input reference point x-coordinate*/
	yref_inp = sf_floatalloc(N-1); /* Input reference point x-coordinate*/
	zref_inp = sf_floatalloc(N-1); /* Input reference point z-coordinate*/
	
	v = sf_floatalloc(nr3+2);  /* Velocity array used in calculation generated according to where the ray travels*/
	gx = sf_floatalloc(nr3+2); /* Velocity gradient in x-direction used in calculation generated according to where the ray travels*/
	gy = sf_floatalloc(nr3+2); /* Velocity gradient in y-direction used in calculation generated according to where the ray travels*/
	gz = sf_floatalloc(nr3+2); /* Velocity gradient in z-direction used in calculation generated according to where the ray travels*/
	xref = sf_floatalloc(nr3+2); /* Reference point x-coordinate used in calculation generated according to where the ray travels*/
	yref = sf_floatalloc(nr3+2); /* Reference point y-coordinate used in calculation generated according to where the ray travels*/
	zref = sf_floatalloc(nr3+2); /* Reference point z-coordinate used in calculation generated according to where the ray travels*/

	F = sf_floatalloc2(2,nr3+2);
	dk = sf_floatalloc2(2,nr3+2);
	xxtem = sf_floatalloc2(2,nr3+2);
	ck_inv = sf_floatalloc3(2,2,nr3+2);
	zk = sf_floatalloc2(2,nr3+2);
	t1_temp = sf_floatalloc2(2,2); /*Temporary value*/
	t2_temp = sf_floatalloc2(2,2);
	t3_temp = sf_floatalloc2(2,2);
	t4_temp = sf_floatalloc2(2,2);
	t5_temp = sf_floatalloc(2);
	t6_temp = sf_floatalloc(2);
	t7_temp = sf_floatalloc2(2,2);
	t8_temp = sf_floatalloc(2);
	t9_temp = sf_floatalloc(2);
	
	/*Set input------------------------------------------------------------------------------*/
	
	if (!sf_getfloat("xs",&xx[0][0])) sf_error("Please enter the source x-coordinate position");
	/* x-Source*/
	
	if (!sf_getfloat("ys",&xx[0][1])) sf_error("Please enter the source y-coordinate position");
	/* y-Source*/
	
	if (!sf_getfloat("xr",&xx[nr3+1][0])) sf_error("Please enter the receiver x-coordinate position");
	/* x-Receiver*/
	
	if (!sf_getfloat("yr",&xx[nr3+1][1])) sf_error("Please enter the receiver y-coordinate position");
	/* y-Receiver*/
	
	if (!sf_getfloats("layer",updown,nr3+1)) sf_error("Please enter the layer number array [nr3+1]");
	/* Layer sequence*/
	
	if (!sf_getint("vstatus",&vstatus)) sf_error("Please enter the status of velocity (0 for constant v,1 for gradient v, and 2 for VTI)");
	/* Velocity status (0 for constant v, 1 for gradient v, and 2 for VTI)*/
	
	if (vstatus != 2){ /*Don't need all of these if  consider vti case*/
		if (!sf_getfloats("velocity",v_inp,N-1)) sf_error("Please enter the velocity array [N-1]");
		/* Assign velocity km/s*/
		
		if (!sf_getfloats("xgradient",gx_inp,N-1)) sf_error("Please enter the x-gradient array [N-1]");
		/* Assign x-gradient*/
		
		if (!sf_getfloats("ygradient",gy_inp,N-1)) sf_error("Please enter the y-gradient array [N-1]");
		/* Assign y-gradient*/
		
		if (!sf_getfloats("zgradient",gz_inp,N-1)) sf_error("Please enter the z-gradient array [N-1]");
		/* Assign z-gradient */
		
		
		if (!sf_getfloats("xref",xref_inp,N-1)) sf_error("Please enter the x-reference points array [N-1]");
		/* Assign x-reference point*/
		
		if (!sf_getfloats("yref",yref_inp,N-1)) sf_error("Please enter the y-reference points array [N-1]");
		/* Assign y-reference point*/
		
		if (!sf_getfloats("zref",zref_inp,N-1)) sf_error("Please enter the z-reference points array [N-1]");
		/* Assign z-reference point*/
	}
	
	if (!sf_getint("order",&order)) order=3;
	/* Interpolation order*/
	
	if (!sf_getfloat("xmin",&xbmin)) {
		xbmin= (xx[0][0]<xx[nr3+1][0])? xx[0][0]:xx[nr3+1][0];
	}
	/* The x-minimum boundary if not entered, set to min(xs,xr)*/
	
	if (!sf_getfloat("ymin",&ybmin)) {
		ybmin=(xx[0][1]<xx[nr3+1][1])? xx[0][1]:xx[nr3+1][1];
	}
	/* The y-minimum boundary if not entered, set to min(ys,yr)*/
	
	if (!sf_getfloat("xmax",&xbmax)) {
		xbmax=(xx[0][0]>xx[nr3+1][0])? xx[0][0]:xx[nr3+1][0];
	}
	/* The x-maximum boundary if not entered, set to max(xr,xr)*/
	
	if (!sf_getfloat("ymax",&ybmax)) {
		ybmax= (xx[0][1]>xx[nr3+1][1])? xx[0][1]:xx[nr3+1][1];
	}
	/* The y-maximum boundary if not entered, set to max(ys,yr)*/
	
	if (!sf_getint("niter",&niter)) sf_error("Please enter the number of iterations");
	/* The number of iterations*/
	
	if (!sf_getbool("debug",&debug)) debug=false ;
	/* Debug flag*/
	
	if (!sf_getdouble("tol",&tol)) tol=0.000001/v_inp[0];
	/* Assign a default value for tolerance*/
	
	/* Set output 2D array reflection point--------------------------------------------------*/
	xrefl = sf_output("out"); /* Output reflection points*/
	if (!sf_getint("ns",&nt)) nt=3 ;
	/* Dimension of output reflection points (x,y,z) */
	
	if (!sf_getint("ns2",&nt2)) nt2=nr3+2; 
	/* Dimension of output reflection points (the number of points)*/
	
	if (!sf_getfloat("ds",&dt)) dt=1; 
	/* Step increment*/
	
	if (!sf_getfloat("s0",&t0)) t0=0;
	/* Staring position*/
	
	nt3=1; /*Force the 3rd dim to be 1*/
	
	sf_putint(xrefl,"n1",nt);
	sf_putint(xrefl,"n2",nt2);
	sf_putint(xrefl,"n3",nt3);
	
	sf_putfloat(xrefl,"d1",dt);
	sf_putfloat(xrefl,"d2",dt);
	sf_putfloat(xrefl,"o1",t0);
		
	/* Read input------------------------------------------------------------------------------*/
	sf_floatread(temp_rr[0][0],nr1*nr2*N,refl);
	
	if (vstatus == 2){
		vti = sf_input("aniso"); /* anisotropy*/
		sf_floatread(aniso_inp[0],4*(N-1),vti);
	}
	
	/* Check the array, consecutive two inputs must not differ by more than 1------------------*/
	int d1,d2,d3,d4,d5,d6,p1,p3; /* Counter*/
	int p2; /* Temp value*/
	float p4=0; /* Temp value*/
	
	for (p1=0; p1<nr3-1; p1++) {
		p2 = updown[p1]-updown[p1+1];
		if (p2>1) {
			sf_warning("The layer number array indicates skipping of layers. Please reenter the array.\n");
			exit(0);
		}
	}
	
	/* Check whether the gradient and vstatus match--------------------------------------------*/
	if (vstatus ==0 || vstatus ==1){ /*Skip this loop for vti*/
		for (p3=0; p3<N-1; p3++) {
			p4 = p4 + gx_inp[p3]+gy_inp[p3]+gz_inp[p3];
			
			if (p3==N-2 && p4/(3*N-2)!=0 && vstatus==0) {
				sf_warning("The gradients are not zero. Please reenter nonzero vstatus.\n");
				exit(0);
			}
			if (p3==N-2 && p4/(3*N-2)==0 && vstatus!=0) {
				sf_warning("The gradients are zero. Please enter vstatus=0 for constant velocity model.\n");
				exit(0);
			}
		}
	}
	/* Generate input according to the reflection sequence-------------------------------------*/
	for (d3=0; d3<nr3+2; d3++) {
		/* Set velocity, gradient, and reference points arrays---------------------------------*/
		if (d3<1) {
			v[d3] = v_inp[0];
			gx[d3] = gx_inp[0];
			gy[d3] = gy_inp[0];
			gz[d3] = gz_inp[0];
			xref[d3] = xref_inp[0];
			yref[d3] = yref_inp[0];
			zref[d3] = zref_inp[0];
			for(d6=0; d6<4; d6++){
				aniso[d3][d6] = aniso_inp[0][d6];
			}
		}
		else {
			d4 = updown[d3-1]; /* Need d3, d4, and d5 because array argument needs to be an interger*/
			d5 = updown[d3];
			
			if (d5-d4>0) {
				v[d3] = v_inp[d4];
				gx[d3] = gx_inp[d4];
				gy[d3] = gy_inp[d4];
				gz[d3] = gz_inp[d4];
				xref[d3] = xref_inp[d4];
				yref[d3] = yref_inp[d4];
				zref[d3] = zref_inp[d4];
				for(d6=0; d6<4; d6++){
					aniso[d3][d6] = aniso_inp[d4][d6];
				}
			}	
			
			if(d5-d4<0){
				v[d3] = v_inp[d5];
				gx[d3] = gx_inp[d5];
				gy[d3] = gy_inp[d5];
				gz[d3] = gz_inp[d5];
				xref[d3] = xref_inp[d5];
				yref[d3] = yref_inp[d5];
				zref[d3] = zref_inp[d5];
				for(d6=0; d6<4; d6++){
					aniso[d3][d6] = aniso_inp[d5][d6];
				}
			}
		}
		
		for (d2=0; d2<nr2; d2++) {
			
			for (d1=0; d1<nr1; d1++) { /* Set layers according to updown*/
				
				if (d3 == 0) {
					rr[d3][d2][d1] = temp_rr[0][d2][d1];
				}
				else {
					d5 = updown[d3-1];
					
					rr[d3][d2][d1] = temp_rr[d5][d2][d1];
				}

			}
		}
	}
	
	/* Initialize interpolation------------------------------------------------------------------*/
	
	if (!sf_getint("order",&order)) order=3;/* Interpolation order*/
	eno = (sf_eno2*) sf_alloc(nr3+2,sizeof(*eno)); /* Allocation for eno array*/
	d1eno = (sf_eno2*) sf_alloc(nr3+2,sizeof(*d1eno));
	d2eno = (sf_eno2*) sf_alloc(nr3+2,sizeof(*d2eno));
	
	/* Compute reflector slope-------------------------------------------------------------------*/
	
	for (ir3=0; ir3 < nr3+2; ir3++) { /* Loop through eno*/
		
		eno[ir3]  = sf_eno2_init(order,nr1,nr2); /* Function values*/
		sf_eno2_set (eno[ir3],rr[ir3]); 
		
		for (ir2=0; ir2<nr2; ir2++) {
			for (ir1=0; ir1 < nr1; ir1++) {
				
				x = r0+ir1*dr; /* Distance in x-direction */
				y = r1+ir2*dr1; /* Distance in y-direction */
				rd1[ir3][ir2][ir1] = zder(ir3,x,y,0);
				rd2[ir3][ir2][ir1] = zder(ir3,x,y,1);
			}
		}
		d1eno[ir3] = sf_eno2_init(order,nr1,nr2);	/* Derivatives*/
		sf_eno2_set (d1eno[ir3],rd1[ir3]);
		d2eno[ir3] = sf_eno2_init(order,nr1,nr2);	/* Derivatives*/
		sf_eno2_set (d2eno[ir3],rd2[ir3]);
	}
	
	/* Set vconstant or vgradient----------------------------------------------------------------*/
	
	func3 f;
	
	f.T_k = 0; /* Initialize structure f to prevent warning*/
	f.T_k_k_1 = 0;
	f.T_k_k_2 = 0;
	f.T_k_k1_1 = 0;
	f.T_k_k1_2 = 0;
	f.T_k_k_k_1 = 0;
	f.T_k_k_k_2 = 0;
	f.T_k_k_k_12 = 0;
	f.T_k_k1_k1_1 = 0;
	f.T_k_k1_k1_2 = 0;
	f.T_k_k1_k1_12 = 0;
	f.T_k_k_k1_1 = 0;
	f.T_k_k_k1_2 = 0;
	f.T_k_k_k1_12 = 0;
	f.T_k_k_k1_21 = 0;
	f.T_k_zk = 0;
	f.T_k_zk1 = 0;
	f.T_k_zk_zk = 0;
	f.T_k_zk1_zk1 = 0;
	f.T_k_zk_zk1 = 0;
	f.T_k_k_zk_1 = 0;
	f.T_k_k_zk_2 = 0;
	f.T_k_k1_zk1_1 = 0;
	f.T_k_k1_zk1_2 = 0;
	f.T_k_k_zk1_1 = 0;
	f.T_k_k_zk1_2 = 0;
	f.T_k_k1_zk_1 = 0;
	f.T_k_k1_zk_2 = 0;
	
	setfunc(vstatus,&f); /* Set value of structure f*/
	
	
	/* To avoid the case where there is NO reflection*/
	if (nr3==0) {
		goto mark; /* If there is no reflection*/
	}
	
	/* If no initial points specified*/
	/*if (!sf_getfloats("xinitial",xinitial,nr3)) {
		for (count=0; count<nr3; count++) {
			xinitial[count] = xx[0][0]+(count+1)*(xx[nr3+1][0]-xx[0][0])/(nr3+1);*/
			/* Divide the distance from s to r equally and set the initial points accordingly*/
	/*	}	
	}*/
	/* x-initial position*/
	
	/*if (!sf_getfloats("yinitial",yinitial,nr3)) {
		for (count=0; count<nr3; count++) {
			yinitial[count] = xx[0][1]+(count+1)*(xx[nr3+1][1]-xx[0][1])/(nr3+1);*/
			/* Divide the distance from s to r equally and set the initial points accordingly*/
	/*	}	
	}*/
	/* y-initial position*/
	
	int ithick;
	float *thick, *sumthick; 
	
	if (!sf_getfloats("xinitial",xinitial,nr3)) {
		thick = sf_floatalloc(nr3+1); /*Avg thickness of each layer for xintial*/
		sumthick = sf_floatalloc(nr3+1); /*Avg thickness of each layer for xintial*/
		for(ithick = 0; ithick < nr3+1; ithick++){ /*To calculate the average thickness of each layer measured from both ends for xinitial*/
			thick[ithick] = ((rr[ithick+1][0][0] - rr[ithick][0][0]) + (rr[ithick+1][nr2-1][nr1-1] - rr[ithick][nr2-1][nr1-1]))/2;
			if (ithick==0){
				sumthick[ithick] = fabsf(thick[ithick]);
			}
			else {
				sumthick[ithick] = sumthick[ithick-1] + fabs(thick[ithick]) ;
			}
		}
		for (count=0; count<nr3; count++) {
			xinitial[count] = xx[0][0]+(xx[nr3+1][0]-xx[0][0])*sumthick[count]/(sumthick[nr3]);
			yinitial[count] = xx[0][1]+(xx[nr3+1][1]-xx[0][1])*sumthick[count]/(sumthick[nr3]);
		}
	}
	
	
	
	
	
/* Step 1: Calculate F(y) to see if it is sufficiently close to zero-------------------------------*/
	
	int i,j,i1,j1,i2,j2,i3; /* Counter*/
	float Ftem=0;
	for (i=0; i<nr3; i++) {
		for (j=0; j<2; j++) {
			if (j==0) {
				xx[i+1][j] = xinitial[i]; /* Create an array of points of intersection from soure to receiver*/
			}
			if (j==1) {
				xx[i+1][j] = yinitial[i]; /* Create an array of points of intersection from soure to receiver*/
			}
		}
	}
	
	for (i1=0; i1<nr3; i1++) {
		initialize(i1+1,nr3,xx,v,xref,yref,zref,gx,gy,gz,aniso,z,zder,zder2_1,zder2_2); /*Initialize y_k and y_k1*/
		for (j1=0; j1<2; j1++) {
			if (j1==0) {
				F[i1+1][j1] = T_hat_1k_k_1(f.T_k_k1_1,f.T_k_zk1) + T_hat_k_k_1(f.T_k_k_1,f.T_k_zk);
			}
			if (j1==1) {
				F[i1+1][j1] = T_hat_1k_k_2(f.T_k_k1_2,f.T_k_zk1) + T_hat_k_k_2(f.T_k_k_2,f.T_k_zk);
			}
		}
	}	
	
	for (i2=0; i2<nr3; i2++) { /* Check the tolerance*/
		for (j2=0; j2<2; j2++) {
			
			Ftem = Ftem+fabsf(F[i2+1][j2]);
			
			if (Ftem<2*nr3*tol && i2 == nr3-1) {
				for (i3=0; i3<nr3; i3++) {
					sf_warning("F(%d) is sufficiently close to zero. x[%d] = %g and y[%d] = %g \n",i3+1,i3+1,xx[i3+1][0],i3+1,xx[i3+1][1]);
				}				
				goto mark; /* Exit the loop to the part for writing the result*/
			}	
		}
	}

/* MAIN LOOP through the output for repeating yk=yk-dk-----------------------------------------------*/
	int i4,j4,i5,j5,i6; /* Counter*/
	int w = niter; /* Number of loops for yk=yk-dk*/
	
	for (q=0; q<w; q++) {
		if (q!=0){
			Ftem=0; /* Reset Ftem to zero*/
			for (i4=0; i4<nr3; i4++) { /* Recalculate F for new y*/
				initialize(i4+1,nr3,xx,v,xref,yref,zref,gx,gy,gz,aniso,z,zder,zder2_1,zder2_2); /* Initialize y_k and y_k1*/
				for (j4=0; j4<2; j4++) {
					if (j4==0) {
						F[i4+1][j4] = T_hat_1k_k_1(f.T_k_k1_1,f.T_k_zk1) + T_hat_k_k_1(f.T_k_k_1,f.T_k_zk);
					}
					if (j4==1) {
						F[i4+1][j4] = T_hat_1k_k_2(f.T_k_k1_2,f.T_k_zk1) + T_hat_k_k_2(f.T_k_k_2,f.T_k_zk);
					}
				}
			}
			
			for (i5=0; i5<nr3; i5++) { /* Check the tolerance*/
				for (j5=0; j5<2; j5++) {
					Ftem = Ftem+fabsf(F[i5+1][j5]);
					if (Ftem<2*nr3*tol && i5 == nr3-1) {
						for (i6=0; i6<nr3; i6++) {
							sf_warning("F(%d) is sufficeintly close to zero. x[%d] = %g and y[%d] = %g \n",i6+1,i6+1,xx[i6+1][0],i6+1,xx[i6+1][1]);
						}				
						goto mark; /* Exit the loop to the part for writing the result*/
					}	
				}
			}
		}
/* Step 2: Forward recursion-------------------------------------------------------------------------*/

		int l,l1,l2; /* Counter*/
		for (l=0; l<nr3; l++) {
			initialize(l+1,nr3,xx,v,xref,yref,zref,gx,gy,gz,aniso,z,zder,zder2_1,zder2_2); /* Initialize y_k and y_k1*/
			for (l1=0; l1<2; l1++) {
				if (l==0) {
					if (l1==0) { /* First column of the 1st element T0,11 + T1,11*/
						ck_inv[1][l1][0]= T_hat_1k_k_k_1(f.T_k_k1_k1_1,f.T_k_k1_zk1_1,f.T_k_zk1,f.T_k_zk1_zk1) + T_hat_k_k_k_1(f.T_k_k_k_1,f.T_k_k_zk_1,f.T_k_zk,f.T_k_zk_zk);
						ck_inv[1][l1][1]= T_hat_1k_k_k_12(f.T_k_k1_k1_12,f.T_k_k1_zk1_1,f.T_k_k1_zk1_2,f.T_k_zk1,f.T_k_zk1_zk1) + T_hat_k_k_k_12(f.T_k_k_k_12,f.T_k_k_zk_1,f.T_k_k_zk_2,f.T_k_zk,f.T_k_zk_zk);
						zk[1][l1] = T_hat_1k_k_1(f.T_k_k1_1,f.T_k_zk1) +T_hat_k_k_1(f.T_k_k_1,f.T_k_zk);
					}
					if (l1==1) { /* First row of the 1st element T0,11 + T1,11*/
						ck_inv[1][l1][0]= T_hat_1k_k_k_12(f.T_k_k1_k1_12,f.T_k_k1_zk1_1,f.T_k_k1_zk1_2,f.T_k_zk1,f.T_k_zk1_zk1) + T_hat_k_k_k_12(f.T_k_k_k_12,f.T_k_k_zk_1,f.T_k_k_zk_2,f.T_k_zk,f.T_k_zk_zk);
						ck_inv[1][l1][1]= T_hat_1k_k_k_2(f.T_k_k1_k1_2,f.T_k_k1_zk1_2,f.T_k_zk1,f.T_k_zk1_zk1) + T_hat_k_k_k_2(f.T_k_k_k_2,f.T_k_k_zk_2,f.T_k_zk,f.T_k_zk_zk);
						zk[1][l1] = T_hat_1k_k_2(f.T_k_k1_2,f.T_k_zk1) +T_hat_k_k_2(f.T_k_k_2,f.T_k_zk);
						mat_inverse(ck_inv[1]);
					}
				}
				
				if(l!=0){
					for (l2=0; l2<2; l2++) { /* ck_inv*/
						if (l2==0) {
							t1_temp[l2][0] = T_hat_1k_1k_k_1(f.T_k_k_k1_1,f.T_k_k1_zk_1,f.T_k_k_zk1_1, f.T_k_zk_zk1);
							t1_temp[l2][1] = T_hat_1k_1k_k_21(f.T_k_k_k1_21,f.T_k_k1_zk_1,f.T_k_k_zk1_2, f.T_k_zk_zk1);
						}
						if (l2==1) {
							t1_temp[l2][0] = T_hat_1k_1k_k_12(f.T_k_k_k1_12,f.T_k_k1_zk_2,f.T_k_k_zk1_1, f.T_k_zk_zk1);
							t1_temp[l2][1] = T_hat_1k_1k_k_2(f.T_k_k_k1_2,f.T_k_k1_zk_2,f.T_k_k_zk1_2, f.T_k_zk_zk1);
						}
					}
					mat_transp(t1_temp,t2_temp);
					mat_mul(ck_inv[l],t1_temp,t3_temp);
					mat_mul(t2_temp,t3_temp,t4_temp);
					matv_mul(ck_inv[l],zk[l],t5_temp);
					matv_mul(t2_temp,t5_temp,t6_temp);
					
					if (l1==0) { /*First column*/
						ck_inv[l+1][l1][0]= T_hat_1k_k_k_1(f.T_k_k1_k1_1,f.T_k_k1_zk1_1,f.T_k_zk1,f.T_k_zk1_zk1) + T_hat_k_k_k_1(f.T_k_k_k_1,f.T_k_k_zk_1,f.T_k_zk,f.T_k_zk_zk) - t4_temp[l1][0];
						ck_inv[l+1][l1][1]= T_hat_1k_k_k_12(f.T_k_k1_k1_12,f.T_k_k1_zk1_1,f.T_k_k1_zk1_2,f.T_k_zk1,f.T_k_zk1_zk1) + T_hat_k_k_k_12(f.T_k_k_k_12,f.T_k_k_zk_1,f.T_k_k_zk_2,f.T_k_zk,f.T_k_zk_zk) - t4_temp[l1][1];
						zk[l+1][l1] = T_hat_1k_k_1(f.T_k_k1_1,f.T_k_zk1) +T_hat_k_k_1(f.T_k_k_1,f.T_k_zk) - t6_temp[l1];
					}
					if (l1==1) { /*First row*/
						ck_inv[l+1][l1][0]= T_hat_1k_k_k_12(f.T_k_k1_k1_12,f.T_k_k1_zk1_1,f.T_k_k1_zk1_2,f.T_k_zk1,f.T_k_zk1_zk1) + T_hat_k_k_k_12(f.T_k_k_k_12,f.T_k_k_zk_1,f.T_k_k_zk_2,f.T_k_zk,f.T_k_zk_zk) - t4_temp[l1][0];
						ck_inv[l+1][l1][1]= T_hat_1k_k_k_2(f.T_k_k1_k1_2,f.T_k_k1_zk1_2,f.T_k_zk1,f.T_k_zk1_zk1) + T_hat_k_k_k_2(f.T_k_k_k_2,f.T_k_k_zk_2,f.T_k_zk,f.T_k_zk_zk)- t4_temp[l1][1];
						zk[l+1][l1] = T_hat_1k_k_2(f.T_k_k1_2,f.T_k_zk1) +T_hat_k_k_2(f.T_k_k_2,f.T_k_zk) - t6_temp[l1];
						mat_inverse(ck_inv[l+1]);
					}
				}
				
				if (isnan(ck_inv[l+1][l1][0]) != 0 || isinf(ck_inv[l+1][l1][0]) != 0) {
					sf_warning("ck_inv doesn't exist. The solutions do not converge.\n");
					exit(0);
				}
			}
		}	

/* Step 3: Backward recursion-----------------------------------------------------------------------*/

		int m,m1; /* Counter*/
		for (m=nr3-1; m>=0; m--) { 
			initialize(m+1,nr3,xx,v,xref,yref,zref,gx,gy,gz,aniso,z,zder,zder2_1,zder2_2); /* Initialize y_k and y_k1*/
			if (m==nr3-1) {
				matv_mul(ck_inv[m+1],zk[m+1],dk[m+1]);
			}
			if (m!=nr3-1) {
				for (m1=0; m1<2; m1++) { /*ck_inv*/
					if (m1==0) {
						t7_temp[m1][0] = T_hat_k_k_k1_1(f.T_k_k_k1_1,f.T_k_k1_zk_1,f.T_k_k_zk1_1, f.T_k_zk_zk1);
						t7_temp[m1][1] = T_hat_k_k_k1_21(f.T_k_k_k1_21,f.T_k_k1_zk_1,f.T_k_k_zk1_2, f.T_k_zk_zk1);
					}
					if (m1==1) {
						t7_temp[m1][0] = T_hat_k_k_k1_12(f.T_k_k_k1_12,f.T_k_k1_zk_2,f.T_k_k_zk1_1, f.T_k_zk_zk1);
						t7_temp[m1][1] = T_hat_k_k_k1_2(f.T_k_k_k1_2,f.T_k_k1_zk_2,f.T_k_k_zk1_2, f.T_k_zk_zk1);
					}
				}
				matv_mul(t7_temp,dk[m+2],t8_temp);
				vector_sub3d_v(2,zk[m+1],2,t8_temp,t9_temp,0);
				matv_mul(ck_inv[m+1],t9_temp,dk[m+1]);
			}
		}
		
		/*Apply boundary---------------------------------------------------------------------------*/
		
		int a,a1,a2,a3,b3,b4; /* Counter*/
		for (a=0; a<nr3; a++) {
			a2=0;
			a3=0;
			for (a1=0; a1<2; a1++) {
				xxtem[a+1][a1] = xx[a+1][a1]-dk[a+1][a1];
				
				if (a1==0) { /* Switch bmin and bmax for x- and y- coordinates*/
					bmin = xbmin;
					bmax = xbmax;
				}
				else if (a1==1){
					bmin = ybmin;
					bmax = ybmax;
				}
				
				while (xxtem[a+1][a1]<bmin && a2<20) {/* Maximum times to multiply is 20*/
					
					dk[a+1][a1]=0.5*dk[a+1][a1];
					if (debug) {
						if (a1==0) {
							sf_warning("The new x[%d] value exceeds the minimum boundary. dk[%d] is reduced to %g\n",a+1,a+1,dk[a+1][a1]);
						}
						if (a1==1) {
							sf_warning("The new y[%d] value exceeds the minimum boundary. dk[%d] is reduced to %g\n",a+1,a+1,dk[a+1][a1]);
						}
					}
					xxtem[a+1][a1] = xx[a+1][a1]-dk[a+1][a1];
					a2++;
				}
				while(xxtem[a+1][a1]>bmax && a3<20) {/* Maximum times to multiply is 20*/
					
					dk[a+1][a1]=0.5*dk[a+1][a1];
					if (debug) {
						if (a1==0) {
							sf_warning("The new x[%d] value exceeds the maximum boundary. dk[%d] is reduced to %g\n",a+1,a+1,dk[a+1][a1]);
						}
						if (a1==1) {
							sf_warning("The new y[%d] value exceeds the maximum boundary. dk[%d] is reduced to %g\n",a+1,a+1,dk[a+1][a1]);
						}
					}
					xxtem[a+1][a1] = xx[a+1][a1]-dk[a+1][a1];
					a3++;
				}
				
				if (a2>=20) {
					if (a1==0) {
						sf_warning("The position x[%d] still exceed the minimum boundary after being halved for 20 times. Please reenter a more appropriate set of initial\n", a+1);
					}
					if (a1==1) {
						sf_warning("The position y[%d] still exceed the minimum boundary after being halved for 20 times. Please reenter a more appropriate set of initial\n", a+1);
					}
					
					exit(0);
				}
				if (a3>=20) {
					if (a1==0) {
						sf_warning("The position x[%d] still exceed the maximum boundary after being halved for 20 times. Please reenter a more appropriate set of initial\n", a+1);
					}
					if (a1==1) {
						sf_warning("The position y[%d] still exceed the maximum boundary after being halved for 20 times. Please reenter a more appropriate set of initial\n", a+1);
					}
					
					exit(0);
				}
			}
		}
		
		xxtem[0][0] = xx[0][0]; /*Fixed source*/
		xxtem[0][1] = xx[0][1]; /*Fixed source*/
		xxtem[nr3+1][0] = xx[nr3+1][0]; /*Fixed receiver*/
		xxtem[nr3+1][1] = xx[nr3+1][1]; /*Fixed receiver*/
		float *dktemp;
		dktemp = sf_floatalloc(2);
		/*Zero thickness---we modify the dk from newton to converge to the same point*/
		for(b3=0; b3<nr3+1; b3++) {
			for (b4=0; b4<2; b4++) {
				if (fabsf(z(b3,xxtem[b3][0],xxtem[b3][1])-z(b3+1,xxtem[b3+1][0],xxtem[b3+1][1]))<0.00001){
					if(fabsf(dk[b3][b4])>fabsf(dk[b3+1][b4])) dktemp[b4] = dk[b3+1][b4];
					if(fabsf(dk[b3][b4])<fabsf(dk[b3+1][b4])) dktemp[b4] = dk[b3][b4];
					
					if (xx[b3][b4]<xx[b3+1][b4]){
						if(b3==0){ /*First layer = zero thickness*/
							dk[b3+1][b4] = fabsf(xx[b3][b4]-xx[b3+1][b4])/2 +dktemp[b4];
							dk[b3][b4] = dktemp[b4];
						}
						else if(b3==nr3) { /*Last layer = zero thickness*/
							dk[b3+1][b4] = dktemp[b4];
							dk[b3][b4] = (-1)*fabsf(xx[b3][b4]-xx[b3+1][b4])/2 +dktemp[b4];
						}
						else { /*Any other layer*/
							dk[b3+1][b4] = fabsf(xx[b3][b4]-xx[b3+1][b4])/2 +dktemp[b4];
							dk[b3][b4] = (-1)*fabsf(xx[b3][b4]-xx[b3+1][b4])/2 +dktemp[b4];
						}
					}
					else {
						if(b3==0){ /*First layer = zero thickness*/
							dk[b3+1][b4] = (-1)*fabsf(xx[b3][b4]-xx[b3+1][b4])/2 +dktemp[b4];
							dk[b3][b4] = dktemp[b4];
						}
						else if(b3==nr3) { /*Last layer = zero thickness*/
							dk[b3+1][b4] = dktemp[b4];
							dk[b3][b4] = fabsf(xx[b3][b4]-xx[b3+1][b4])/2 +dktemp[b4];
						}
						else { /*Any other layer*/
							dk[b3][b4] = fabsf(xx[b3][b4]-xx[b3+1][b4])/2 +dktemp[b4];
							dk[b3+1][b4] = (-1)*fabsf(xx[b3][b4]-xx[b3+1][b4])/2 +dktemp[b4];
						}
						dk[b3][b4] = fabsf(xx[b3][b4]-xx[b3+1][b4])/2 +dktemp[b4];
						dk[b3+1][b4] = (-1)*fabsf(xx[b3][b4]-xx[b3+1][b4])/2 +dktemp[b4];
					}
				}
			}
		}
		
		
/* Step 4: Update xx------------------------------------------------------------------------------*/
		
		vector_sub3d(nr3+2,xx,nr3+2,dk,xxnew,0);
		int t,t1;/*counter*/
		for (t=0; t<nr3+2;t++) {
			for (t1=0; t1<2; t1++) {
				if (debug) {
					if (t1==0) {
						sf_warning("The original value of x[%d] is %g, d[%d] is %g and the new value of x[%d] is %g\n",t,xx[t][t1],t,dk[t][t1],t,xxnew[t][t1]);
					}
					if (t1==1) {
						sf_warning("The original value of y[%d] is %g, d[%d] is %g and the new value of y[%d] is %g\n",t,xx[t][t1],t,dk[t][t1],t,xxnew[t][t1]);
					}
					if (t==nr3+1 && t1==1) {
						sf_warning("Iteration:%d\n\n",q+1);
					}					
				}
				
				xx[t][t1] = xxnew[t][t1]; /* Update xx values*/
			}
		}
	}
	
/* END OF MAIN LOOP--------------------------------------------------------------------------------*/
	
	/* Write result in 2D & Compute traveltime-----------------------------------------------------*/
	int c,c1,c2,c3; /* Counter*/
	
mark: /* Mark point for goto*/
	
	for (c2=0; c2<nr3+2; c2++) {
		for (c1=0; c1<3; c1++) {
			if (c1==2) {
				ans[c2][2] = z(c2,xx[c2][0],xx[c2][1]);
			}
			else {
				ans[c2][0] = xx[c2][0]; /*x*/
				ans[c2][1] = xx[c2][1]; /*y*/
			}
		}
	}
	
	if (q == niter){
		for (c3=0; c3<nr3; c3++) {
			sf_warning("F(%d) is sufficeintly close to zero. x[%d] = %g and y[%d] = %g \n",c3+1,c3+1,xx[c3+1][0],c3+1,xx[c3+1][1]);
		}
	}
	
	tt=0; /* Initialize traveltime tt*/
	
	for (c=0; c<nr3+1; c++) {
		half_initialize(c,nr3,xx,v,xref,yref,zref,gx,gy,gz,aniso,z,zder,zder2_1,zder2_2);
		tt = tt + T_hat_k(f.T_k);
		if (c==nr3) {
			sf_warning("Traveltime is %g and the total number of iterations is %d",tt,q);
		}
	}
	
	/* Write output */
	sf_floatwrite(ans[0],3*(nr3+2),xrefl);
	
	exit(0);
}
