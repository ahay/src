/* 2D Bending ray tracing in multi-layered media
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
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "vectorops.h"
#include "general_traveltime.h"
#include "ml_traveltime_vconstant.h"
#include "ml_traveltime_vgradient.h"
#include "ml_traveltime_vti.h"
#include "setvelocity.h"

/* Reflector function----------------------------------------------------------------------------------*/

static sf_eno *eno, *deno; /* Interpolation structure */	
static float r0, dr,r1,dr1;

static float z(int k,float x) 
/* Function */
{
	int i;
	float f, f1;
	
	x = (x-r0)/dr; 
	i = floorf(x);
	
	sf_eno_apply (eno[k],i,x-i,&f,&f1,FUNC);
	return f;
}

static float zder(int k,float x) 
/* First derivative */
{
	int i;
	float f, f1;
	
	x = (x-r0)/dr; 
	i = floorf(x);
	
	sf_eno_apply (eno[k],i,x-i,&f,&f1,DER);
	return f1/dr;
}

static float zder2(int k,float x) 
/* Second derivative */
{
	int i;
	float f, f1;
	
	x = (x-r0)/dr; 
	i = floorf(x);
	
	sf_eno_apply (deno[k],i,x-i,&f,&f1,DER);
	return f1/dr;	
}


/* Main program------------------------------------------------------------------------------------------*/
int main(int argc, char* argv[])
{
	int nr1,nr2, N, ir1,ir2, nt, nt2, order, niter, vstatus, count, q=0/* Counter for the bog loop*/;
	float x, dt, t0, bmin, bmax, tt;
	float **rr, **temp_rr, **rd, **ans, *xx, *xxnew, *xinitial, *updown, *v_inp, *gx_inp, *gz_inp, *xref_inp,*zref_inp ,*v, *gx, *gz, *xref, *zref, *F, *dk, *xxtem, *zk,*ck_inv, **aniso, **aniso_inp; 
	double tol;
	bool  debug;
	sf_file refl, xrefl, vti;
	
	sf_init(argc,argv); /* initialize - always call first */
	
	/* Set input-----------------------------------------------------------------------------------------*/
	refl = sf_input("in"); /* reflector */
	if (!sf_histint(refl,"n1",&nr1)) sf_error("No n1= in input");
	if (!sf_histint(refl,"n2",&N)) sf_error("No n2= in input");
	
	if (!sf_histfloat(refl,"o1",&r0)) r0=0.;
	if (!sf_histfloat(refl,"o2",&r1)) r1=0.;
	
	if (!sf_histfloat(refl,"d1",&dr)) dr=1.;
	if (!sf_histfloat(refl,"d2",&dr1)) dr1=1.;
	
	
	if (!sf_getint("number",&nr2)) sf_error("Please enter the number of reflections [nr2]");
	/* Number of intersecting points [nr2]*/
	
	/* Allocate space-------------------------------------------------------------------------------------*/
	temp_rr = sf_floatalloc2(nr1,N); /* Input reflector values*/
	rr = sf_floatalloc2(nr1,nr2+2); /* Reflector values according to updown*/
	rd = sf_floatalloc2(nr1,nr2+2); /* Slope values according to updown*/
	ans = sf_floatalloc2(2,nr2+2); /* Final answer x- and z- coordinates*/
	aniso_inp = sf_floatalloc2(4,N-1); /* VTI parameters of the model*/
	aniso = sf_floatalloc2(4,nr2+2); /* VTI parameters of the model*/
	
	xx = sf_floatalloc(nr2+2); /* Positions of intersection points to be iteratively perturbed*/
	xxnew = sf_floatalloc(nr2+2); /* Temporary array for xx*/
	if (nr2!=0) {
		xinitial = sf_floatalloc(nr2); /* Initial guess of intersection points*/
	}
	updown = sf_floatalloc(nr2+1); /* Array indicating direction of the ray*/
	
	v_inp = sf_floatalloc(N-1); /* Input velocity array*/
	gx_inp = sf_floatalloc(N-1); /* Input velocity gradient in x-direction*/
	gz_inp = sf_floatalloc(N-1); /* Input velocity gradient in z-direction*/
	xref_inp = sf_floatalloc(N-1); /* Input reference point x-coordinate*/
	zref_inp = sf_floatalloc(N-1); /* Input reference point z-coordinate*/
	
	v = sf_floatalloc(nr2+2);  /* Velocity array used in calculation generated according to where the ray travels*/
	gx = sf_floatalloc(nr2+2); /* Velocity gradient in x-direction used in calculation generated according to where the ray travels*/
	gz = sf_floatalloc(nr2+2); /* Velocity gradient in z-direction used in calculation generated according to where the ray travels*/
	xref = sf_floatalloc(nr2+2); /* Reference point x-coordinate used in calculation generated according to where the ray travels*/
	zref = sf_floatalloc(nr2+2); /* Reference point z-coordinate used in calculation generated according to where the ray travels*/
	
	F = sf_floatalloc(nr2+2); /* Array of first derivative of traveltime (Snell's law at each interface)*/
	dk = sf_floatalloc(nr2+2); /* Changes to xx*/
	xxtem = sf_floatalloc(nr2+2); /* Pre-calculated xx to check the boundary condition*/
	ck_inv = sf_floatalloc(nr2+2);
	zk = sf_floatalloc(nr2+2);
	
	
	
	/* Set input------------------------------------------------------------------------------------------*/
	
	if (!sf_getfloat("xs",&xx[0])) sf_error("Please enter the source position");
	/* Source*/
	
	if (!sf_getfloat("xr",&xx[nr2+1])) sf_error("Please enter the receiver position");
	/* Receiver*/
	
	if (!sf_getfloats("layer",updown,nr2+1)) sf_error("Please enter the layer number array [nr2+1]");
	/* Layer sequence*/
	
	if (!sf_getint("vstatus",&vstatus)) sf_error("Please enter the status of velocity (0 for constant v,1 for gradient v, and 2 for VTI)");
	/* Velocity status (0 for constant v, 1 for gradient v, and 2 for VTI)*/
	
	if (vstatus == 1){ /*Don't need all of these if  consider vti case*/
		if (!sf_getfloats("velocity",v_inp,N-1)) sf_error("Please enter the velocity array [N-1]");
		/* Assign velocity km/s*/
		
		if (!sf_getfloats("xgradient",gx_inp,N-1)) sf_error("Please enter the x-gradient array [N-1]");
		/* Assign x-gradient*/
		
		if (!sf_getfloats("zgradient",gz_inp,N-1)) sf_error("Please enter the z-gradient array [N-1]");
		/* Assign z-gradient */
		
		if (!sf_getfloats("xref",xref_inp,N-1)) sf_error("Please enter the x-reference points array [N-1]");
		/* Assign x-reference point*/
		
		if (!sf_getfloats("zref",zref_inp,N-1)) sf_error("Please enter the z-reference points array [N-1]");
		/* Assign z-reference point*/
	}
	else {
		if (!sf_getfloats("velocity",v_inp,N-1)) sf_error("Please enter the velocity array [N-1]");
		/* Assign velocity km/s*/
		int index;
		for(index=0;index<N-1;index++) {
			gx_inp[index] = 0.0;
			gz_inp[index] = 0.0;
		}
	}
	
	if (!sf_getfloat("min",&bmin)) bmin=xx[0];
	/* The minimum boundary if not entered, set to xs*/
	
	
	if (!sf_getfloat("max",&bmax)) bmax=xx[nr2+1];
	/* The maximum boundary if not entered, set to xr*/
	
	if (!sf_getint("niter",&niter)) niter=100;
	/* The number of iterations*/
	
	if (!sf_getbool("debug",&debug)) debug=false;
	/* Debug flag*/
	
	if (!sf_getdouble("tol",&tol)) tol=0.000001/v_inp[0];
	/* Assign a default value for tolerance*/
	
	
	/* Set output 2D array reflection point----------------------------------------------------------------*/
	xrefl = sf_output("out"); /* Output reflection points*/
	
	if (!sf_getint("ns",&nt)) nt=2;
	/* Dimension of output reflection points (x,z)*/
	
	if (!sf_getint("ns2",&nt2)) nt2=nr2+2; 
	/* Dimension of output reflection points (the number of points)*/
	
	
	if (!sf_getfloat("ds",&dt)) dt=1; 
	/* Step increment*/
	
	if (!sf_getfloat("s0",&t0)) t0=0;
	/* Staring position*/
	
	sf_putint(xrefl,"n1",nt);
	sf_putint(xrefl,"n2",nt2);
	
	sf_putfloat(xrefl,"d1",dt);
	sf_putfloat(xrefl,"o1",t0);
	
	
	/* Read input-----------------------------------------------------------------------------------------*/	
	sf_floatread(temp_rr[0],nr1*N,refl);
	
	if (vstatus == 2){
		vti = sf_input("aniso"); /* anisotropy*/
		sf_floatread(aniso_inp[0],4*(N-1),vti);
	}
	/* Check the array, consecutive two inputs must not differ by more than 1-----------------------------*/
	
	int d1,d2,d3,d4,d5,d6,p1,p3; /*counter*/
	int p2; /*Temp value*/
	float p4=0; /*Temp value*/
	
	for (p1=0; p1<nr2-1; p1++) {
		p2 = updown[p1]-updown[p1+1];
		if (p2>1) {
			sf_warning("The layer number array indicates skipping of layers. Please reenter the array.\n");
			exit(0);
		}
	}
	
	/* Check whether the gradient and vstatus match-------------------------------------------------------*/
	if (vstatus ==0 || vstatus ==1){ /*Skip this loop for vti*/
		for (p3=0; p3<N-1; p3++) {
			p4 = p4 + gx_inp[p3]+gz_inp[p3];
			
			if (p3==N-2 && p4/(2*N-2)!=0 && vstatus==0) {
				sf_warning("The gradients are not zero. Please reenter nonzero vstatus.\n");
				exit(0);
			}
			if (p3==N-2 && p4/(2*N-2)==0 && vstatus!=0) {
				sf_warning("The gradients are zero. Please enter vstatus=0 for constant velocity model.\n");
				exit(0);
			}
		}
	}
	
	/* Generate input according to the reflection sequence-----------------------------------------------*/
	for (d2=0; d2<nr2+2; d2++) {
		
		/* Set velocity, gradient, and reference points arrays-------------------------------------------*/
			if (d2<1) {
				v[d2] = v_inp[0];
				gx[d2] = gx_inp[0];
				gz[d2] = gz_inp[0];
				xref[d2] = xref_inp[0];
				zref[d2] = zref_inp[0];
				for(d6=0; d6<4; d6++){
					aniso[d2][d6] = aniso_inp[0][d6];
				}
			}
			else {
				d3 = updown[d2-1]; /* Need d3, d4, and d5 because array argument needs to be an interger*/
				d4 = updown[d2];
				
				if (d4-d3>0) {
					v[d2] = v_inp[d3];
					gx[d2] = gx_inp[d3];
					gz[d2] = gz_inp[d3];
					xref[d2] = xref_inp[d3];
					zref[d2] = zref_inp[d3];
					for(d6=0; d6<4; d6++){
						aniso[d2][d6] = aniso_inp[d3][d6];
					}
				}	
				if(d4-d3<0){
					v[d2] = v_inp[d4];
					gx[d2] = gx_inp[d4];
					gz[d2] = gz_inp[d4];
					xref[d2] = xref_inp[d4];
					zref[d2] = zref_inp[d4];
					for(d6=0; d6<4; d6++){
						aniso[d2][d6] = aniso_inp[d4][d6];
					}
				}
			}
		
		for (d1=0; d1<nr1; d1++) { /* Set layers according to updown*/
			
			if (d2 == 0) {
				rr[d2][d1] = temp_rr[0][d1];
			}
			else {
			d5 = updown[d2-1];
			rr[d2][d1] = temp_rr[d5][d1];
			}	
		}
	}
	
	/* Initialize interpolation------------------------------------------------------------------------*/
	
	if (!sf_getint("order",&order)) order=3;/* Interpolation order*/
	eno = (sf_eno*) sf_alloc(nr2+2,sizeof(*eno)); /* Allocation for eno array*/
	deno = (sf_eno*) sf_alloc(nr2+2,sizeof(*deno));

		/* Compute reflector slope---------------------------------------------------------------------*/
	
		for (ir2=0; ir2 < nr2+2; ir2++) { /* Loop through eno*/
			eno[ir2]  = sf_eno_init(order,nr1); /* Function values*/
			sf_eno_set (eno[ir2],rr[ir2]); 
			
			for (ir1=0; ir1 < nr1; ir1++) {
				x = r0+ir1*dr; /* Distance */
				rd[ir2][ir1] = zder(ir2,x);
			}
			deno[ir2] = sf_eno_init(order,nr1);	/* Derivatives*/
			sf_eno_set (deno[ir2],rd[ir2]);
		}
	
	/* Set vconstant or vgradient----------------------------------------------------------------------*/
	
	func3 f;
	
	f.T_k = 0; /* Initialize structure f*/
	f.T_k_k = 0;
	f.T_k_k1 = 0;
	f.T_k_k_k = 0;
	f.T_k_k1_k1 = 0;
	f.T_k_k_k1 = 0;
	f.T_k_zk = 0;
	f.T_k_zk1 = 0;
	f.T_k_zk_zk = 0;
	f.T_k_zk1_zk1 = 0;
	f.T_k_zk_zk1 = 0;
	f.T_k_k_zk = 0;
	f.T_k_k1_zk1 = 0;
	f.T_k_k_zk1 = 0;
	f.T_k_k1_zk = 0;
	
	setfunc(vstatus,&f); /* Set value of structure f*/
	
	/* To avoid the case where there is NO reflection*/
	if (nr2==0) {
		goto mark; /* If there is no reflection*/
	}
	
	int ithick;
	float *thick, *sumthick; 
	
	if (!sf_getfloats("xinitial",xinitial,nr2)) {
		thick = sf_floatalloc(nr2+1); /*Avg thickness of each layer for xintial*/
		sumthick = sf_floatalloc(nr2+1); /*Avg thickness of each layer for xintial*/
		for(ithick = 0; ithick < nr2+1; ithick++){ /*To calculate the average thickness of each layer measured from both ends for xinitial*/
			thick[ithick] = ((rr[ithick+1][0] - rr[ithick][0]) + (rr[ithick+1][nr1-1] - rr[ithick][nr1-1]))/2;
			if (ithick==0){
				sumthick[ithick] = fabsf(thick[ithick]);
			}
			else {
				sumthick[ithick] = sumthick[ithick-1] + fabs(thick[ithick]) ;
			}
		}
		for (count=0; count<nr2; count++) {
			xinitial[count] = xx[0]+(xx[nr2+1]-xx[0])*sumthick[count]/(sumthick[nr2]);			
		}
		/*for (count=0; count<nr2; count++) {
			xinitial[count] = xx[0]+(count+1)*(xx[nr2+1]-xx[0])/(nr2+1); Divide the distance from s to r equally and set the initial points accordingly
		} */	
	}
	/* Initial position*/
	
	
/* Step 1: Calculate F(y) to see if it is sufficiently close to zero----------------------------------*/
	
	int i,j1,i3,i4; /* Counter*/
	float Ftem=0;
	
	for (i3=0; i3<nr2; i3++) {
		xx[i3+1] = xinitial[i3]; /* Create an array of points of intersection from source to receiver*/
	}
	
	for (i=0; i<nr2; i++) {
		initialize(i+1,nr2,xx,v,xref,zref,gx,gz,aniso,z,zder,zder2); /*Initialize y_1k, y_k and y_k1*/
		F[i+1] = T_hat_1k_k(f.T_k_k1,f.T_k_zk1) + T_hat_k_k(f.T_k_k,f.T_k_zk);
	}	
	
	for (i4=0; i4<nr2; i4++) { /* Check the tolerance*/
		Ftem = Ftem+fabsf(F[i4+1]);
		if (Ftem<nr2*tol && i4 == nr2-1) {
			for (j1=0; j1<nr2; j1++) {
				sf_warning("F(%d) is sufficiently close to zero. y[%d] = %g \n",j1+1,j1+1,xx[j1+1]);
			}
			goto mark; /* Exit the loop to the part for writing the result*/
		}	
	}
/* MAIN LOOP through the output for repeating yk=yk-dk-----------------------------------------------*/
	
	int i2,j2,i5; /* Counter*/
	int w = niter; /* Number of loops for yk=yk-dk*/
	
	for (q=0; q<w; q++) {
		Ftem=0; /* Reset Ftem to zero*/
		
		for (i2=0; i2<nr2; i2++) { /* Recalculate F for new y (Repeat Step 1)*/
			initialize(i2+1,nr2,xx,v,xref,zref,gx,gz,aniso,z,zder,zder2);
			F[i2+1] = T_hat_1k_k(f.T_k_k1,f.T_k_zk1) + T_hat_k_k(f.T_k_k,f.T_k_zk);
			
			/*Zero thickness layer case*/
			if (fabsf(z(i2+1,xx[i2+1])-z(i2+2,xx[i2+2]))<0.00001 || fabsf(z(i2,xx[i2])-z(i2+1,xx[i2+1]))<0.00001){ 
			/*Cheap trick to get it to exit. Not respect Snell's law. Traveltime derivative is 0/0*/
				F[i2+1] = 0.0;
			}
		}
		
		for (i5=0; i5<nr2; i5++) { /* Check the tolerance*/
			Ftem = Ftem+fabsf(F[i5+1]);
			if (Ftem<nr2*tol && i5 == nr2-1) {
				for (j2=0; j2<nr2; j2++) {
					sf_warning("F(%d) is sufficiently close to zero. y[%d] = %g \n",j2+1,j2+1,xx[j2+1]);
				}
				goto mark; /* Exit the loop to the part for writing the result*/
			}
		}
		
/* Step 2: Forward recursion-------------------------------------------------------------------------*/
		
		int l; /* Counter*/
		for (l=0; l<nr2; l++) {
			initialize(l+1,nr2,xx,v,xref,zref,gx,gz,aniso,z,zder,zder2);
			if (l==0) {
				ck_inv[1]= 1/(T_hat_1k_k_k(f.T_k_k1_k1,f.T_k_k1_zk1,f.T_k_zk1,f.T_k_zk1_zk1) + T_hat_k_k_k(f.T_k_k_k,f.T_k_k_zk,f.T_k_zk,f.T_k_zk_zk));
				zk[1] = T_hat_1k_k(f.T_k_k1,f.T_k_zk1) +T_hat_k_k(f.T_k_k,f.T_k_zk);
			}
			else {
				ck_inv[l+1]= 1/(T_hat_1k_k_k(f.T_k_k1_k1,f.T_k_k1_zk1,f.T_k_zk1,f.T_k_zk1_zk1) + T_hat_k_k_k(f.T_k_k_k,f.T_k_k_zk,f.T_k_zk,f.T_k_zk_zk) - T_hat_1k_1k_k(f.T_k_k_k1,f.T_k_k1_zk,f.T_k_k_zk1,f.T_k_zk_zk1)*ck_inv[l]*T_hat_1k_1k_k(f.T_k_k_k1,f.T_k_k1_zk,f.T_k_k_zk1,f.T_k_zk_zk1));
				zk[l+1] = T_hat_1k_k(f.T_k_k1,f.T_k_zk1) + T_hat_k_k(f.T_k_k,f.T_k_zk) - T_hat_1k_1k_k(f.T_k_k_k1,f.T_k_k1_zk,f.T_k_k_zk1,f.T_k_zk_zk1)*ck_inv[l]*zk[l];
			}

			if (isnan(1/ck_inv[l+1]) != 0 || isinf(1/ck_inv[l+1]) != 0) {
				sf_warning("ck_inv doesn't exist. The solutions do not converge.\n");
				exit(0);
			}
		}	
		
/* Step 3: Backward recursion-----------------------------------------------------------------------*/
		
		int m; /* Counter*/
		for (m=nr2-1; m>=0; m--) { 
			initialize(m+1,nr2,xx,v,xref,zref,gx,gz,aniso,z,zder,zder2);
			if (m==nr2-1) {
				dk[m+1] = ck_inv[m+1]*zk[m+1];
			}
			else {
				dk[m+1] = ck_inv[m+1]*(zk[m+1]-T_hat_k_k_k1(f.T_k_k_k1,f.T_k_k1_zk,f.T_k_k_zk1,f.T_k_zk_zk1)*dk[m+2]);
			}
		}
	
		/*Apply boundary & check zero thickness---------------------------------------------------------------------------*/
		
		dk[0] = 0;
		dk[nr2+1] = 0;
		xxtem[0] = xx[0]; /*Fixed source*/
		xxtem[nr2+1] = xx[nr2+1]; /*Fixed receiver*/
		int t,a,b1,b2,b3; /* Counter*/
		float dktemp;
		for (a=0; a<nr2; a++) {
			b1=0;
			b2=0;
			xxtem[a+1] = xx[a+1]-dk[a+1];
			while (xxtem[a+1]<bmin && b1<20) {/* Maximum times to multiply is 20*/
				
				dk[a+1]=0.5*dk[a+1]; /* Decrease the change by half*/
				if (debug) {
				sf_warning("The new x value exceeds the minimum boundary. dk[%d] is reduced to %g\n",a+1,dk[a+1]);
				}
				xxtem[a+1] = xx[a+1]-dk[a+1]; /* Recompute xxtem to see if it is still exceed the boundary*/
				b1++;
			}
			while(xxtem[a+1]>bmax && b2<20) {/* Maximum times to multiply is 20*/
				
				dk[a+1]=0.5*dk[a+1];
				if (debug) {
				sf_warning("The new x value exceeds the maximum boundary. dk[%d] is reduced to %g\n",a+1,dk[a+1]);
				}
				xxtem[a+1] = xx[a+1]-dk[a+1];
				b2++;
			}
			if (b1>=20) {
				sf_warning("The position x[%d] still exceed the minimum boundary after being halved for 20 times. Please reenter a more appropriate set of xinitial\n", a+1);
				exit(0);
			}
			if (b2>=20) {
				sf_warning("The position x[%d] still exceed the maximum boundary after being halved for 20 times. Please reenter a more appropriate set of xinitial\n", a+1);
				exit(0);
			}
		}
		
		/*Zero thickness---we modify the dk from newton to converge to the same point*/
		for(b3=0; b3<nr2+1; b3++) {
			if (fabsf(z(b3,xxtem[b3])-z(b3+1,xxtem[b3+1]))<0.00001){
				if(fabsf(dk[b3])>fabsf(dk[b3+1])) dktemp = dk[b3+1];
				if(fabsf(dk[b3])<fabsf(dk[b3+1])) dktemp = dk[b3];
				
				if (xx[b3]<xx[b3+1]){
					if(b3==0){ /*First layer = zero thickness*/
						dk[b3+1] = fabsf(xx[b3]-xx[b3+1])/2 +dktemp;
						dk[b3] = dktemp;
					}
					else if(b3==nr2) { /*Last layer = zero thickness*/
						dk[b3+1] = dktemp;
						dk[b3] = (-1)*fabsf(xx[b3]-xx[b3+1])/2 +dktemp;
					}
					else { /*Any other layer*/
						dk[b3+1] = fabsf(xx[b3]-xx[b3+1])/2 +dktemp;
						dk[b3] = (-1)*fabsf(xx[b3]-xx[b3+1])/2 +dktemp;
					}
				}
				else {
					if(b3==0){ /*First layer = zero thickness*/
						dk[b3+1] = (-1)*fabsf(xx[b3]-xx[b3+1])/2 +dktemp;
						dk[b3] = dktemp;
					}
					else if(b3==nr2) { /*Last layer = zero thickness*/
						dk[b3+1] = dktemp;
						dk[b3] = fabsf(xx[b3]-xx[b3+1])/2 +dktemp;
					}
					else { /*Any other layer*/
						dk[b3] = fabsf(xx[b3]-xx[b3+1])/2 +dktemp;
						dk[b3+1] = (-1)*fabsf(xx[b3]-xx[b3+1])/2 +dktemp;
					}
					dk[b3] = fabsf(xx[b3]-xx[b3+1])/2 +dktemp;
					dk[b3+1] = (-1)*fabsf(xx[b3]-xx[b3+1])/2 +dktemp;
				}
			}
		}
		
/* Step 4: Update xx------------------------------------------------------------------------------*/
		
		vector_sub(nr2+2,xx,nr2+2,dk,xxnew,0); /* Update xx (Newton)*/
		for (t=0; t<nr2+2;t++) {
			if (debug) {
				sf_warning("The original value of y[%d] is %g, d[%d] is %g and the new value of y[%d] is %g\n",t,*(xx+t),t,*(dk+t),t,*(xxnew+t));
				if (t==nr2+1) {
					sf_warning("Iteration:%d\n",q+1);
				}					
			}
			*(xx+t) = *(xxnew+t); /* Update xx values*/
		}
	}
	
/* END OF MAIN LOOP--------------------------------------------------------------------------------*/	
	
	/* Write result in 2D & Compute traveltime-----------------------------------------------------*/
	int c1,c2,c3;
	int c; /* Counter*/
	
mark: /* Mark point for goto*/
	
	for (c2=0; c2<nr2+2; c2++) {
		for (c1=0; c1<2; c1++) {
			if (c1==0) {
				ans[c2][c1] = xx[c2];
			}
			else {
				ans[c2][c1] = z(c2,xx[c2]);
			}
		}
	}
	if (q == niter){
		for (c3=0; c3<nr2; c3++) {
			sf_warning("F(%d) is sufficiently close to zero. y[%d] = %g \n",c3+1,c3+1,xx[c3+1]);
		}
	}
	tt=0; /* Initialize traveltime tt*/
	
	for (c=0; c<nr2+1; c++) {
		half_initialize(c,nr2,xx,v,xref,zref,gx,gz,aniso,z,zder,zder2);
		tt = tt + T_hat_k(f.T_k);
		if (c==nr2) {
			sf_warning("Traveltime is %g and the total number of iterations is %d\n",tt,q);
		}
	}
	
	/* Write output*/
	sf_floatwrite(ans[0],2*(nr2+2),xrefl);
	
	exit(0);
}
