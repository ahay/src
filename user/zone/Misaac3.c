/*
 *  Misaac3.c
 *  
 *
 *  Created by Yanadet Sripanich on 11/11/12.
 *  Copyright 2012 __MyCompanyName__. All rights reserved.
 *
 */


/* Reflection traveltime for Multi-layered at any specified source and receiver location */
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
#include "general_traveltime.h"
#include "vectorsub.h"
#include "ml_traveltime_vconstant.h"

/*Reflector function--------------------------------------------------------------------------------*/

static sf_eno *eno, *deno; /* interpolation structure */	
static float r0, dr,r1,dr1;

static float z(int k,float x) 
/* function */
{
	int i;
	float f, f1;
	
	x = (x-r0)/dr; 
	i = floorf(x);
	
	sf_eno_apply (eno[k],i,x-i,&f,&f1,FUNC);
	return f;
}

static float zder(int k,float x) 
/* first derivative */
{
	int i;
	float f, f1;
	
	x = (x-r0)/dr; 
	i = floorf(x);
	
	sf_eno_apply (eno[k],i,x-i,&f,&f1,DER);
	return f1/dr;
}

static float zder2(int k,float x) 
/* second derivative */
{
	int i;
	float f, f1;
	
	x = (x-r0)/dr; 
	i = floorf(x);
	
	sf_eno_apply (deno[k],i,x-i,&f,&f1,DER);
	return f1/dr;	
}


/*Main program-------------------------------------------------------------------------------------------*/
int main(int argc, char* argv[])
{
	int nr1,nr2, N, ir1,ir2, nt, nt2, order, niter;
	float x, dt, t0, bmin, bmax, tt;
	float **rr, **temp_rr, **rd, **ans, *xx, *xxnew, *xinitial, *updown, *v_inp, *gx_inp, *gz_inp, *xref_inp,*zref_inp ,*v, *gx, *gz, *xref, *zref, *F, *dk, *dk_old, *xxtem, *zk,*ck_inv; 
	double tol;
	sf_file refl, xrefl;
	
	sf_init(argc,argv); /* initialize - always call first */
	
	
	
	/*Set input------------------------------------------------------------------------------*/
	refl = sf_input("in"); /* reflector */
	if (!sf_histint(refl,"n1",&nr1)) sf_error("No n1= in input");
	if (!sf_histint(refl,"n2",&N)) sf_error("No n2= in input");
	
	if (!sf_histfloat(refl,"o1",&r0)) r0=0.;
	if (!sf_histfloat(refl,"o2",&r1)) r1=0.;
	
	if (!sf_histfloat(refl,"d1",&dr)) dr=1.;
	if (!sf_histfloat(refl,"d2",&dr1)) dr1=1.;
	
	
	if (!sf_getint("number",&nr2)) sf_error("Please enter the number of reflections [nr2]");
	/* Number of reflectors*/
	
	/*Allocate space--------------------------------------------------------------------------*/
	rr = sf_floatalloc2(nr1,nr2+2);
	temp_rr = sf_floatalloc2(nr1,N);
	rd = sf_floatalloc2(nr1,nr2+2);
	ans = sf_floatalloc2(nr2+2,2);
	
	
	
	xx = sf_floatalloc(nr2+2);
	xxnew = sf_floatalloc(nr2+2);
	xinitial = sf_floatalloc(nr2);
	updown = sf_floatalloc(nr2+1);
	
	v_inp = sf_floatalloc(N-1); /* Input velocity array*/
	gx_inp = sf_floatalloc(N-1); /* Input velocity gradient in x-direction*/
	gz_inp = sf_floatalloc(N-1); /* Input velocity gradient in z-direction*/
	xref_inp = sf_floatalloc(N-1); /* Input reference point x-coordinate*/
	zref_inp = sf_floatalloc(N-1); /* Input reference point z-coordinate*/
	
	v = sf_floatalloc(nr2+1);  /* Velocity array used in calculation generated according to where the ray travels*/
	gx = sf_floatalloc(nr2+1); /* Velocity gradient in x-direction used in calculation generated according to where the ray travels*/
	gz = sf_floatalloc(nr2+1); /* Velocity gradient in z-direction used in calculation generated according to where the ray travels*/
	xref = sf_floatalloc(nr2+1); /* Reference point x-coordinate used in calculation generated according to where the ray travels*/
	zref = sf_floatalloc(nr2+1); /* Reference point z-coordinate used in calculation generated according to where the ray travels*/
	
	F = sf_floatalloc(nr2+2);
	dk = sf_floatalloc(nr2+2);
	dk_old = sf_floatalloc(nr2+2);
	xxtem = sf_floatalloc(nr2+2);
	ck_inv = sf_floatalloc(nr2+2);
	zk = sf_floatalloc(nr2+2);

	
	
	/*Set input------------------------------------------------------------------------------*/
	
	
	if (!sf_getfloat("xs",&xx[0])) sf_error("Please enter the source position");
	/* Source*/
	
	if (!sf_getfloat("xr",&xx[nr2+1])) sf_error("Please enter the receiver position");
	/* Receiver*/
	
	if (!sf_getfloats("layer",updown,nr2+1)) sf_error("Please enter the layer number array [nr2+1]");
	/* layer sequence*/
	
	if (!sf_getfloats("xinitial",xinitial,nr2)) sf_error("Please enter the initial position array [nr2]");
	/* initial position*/
	
	if (!sf_getfloats("velocity",v_inp,N-1)) sf_error("Please enter the velocity array [N-1]");
	/* assign velocity km/s*/
	
	if (!sf_getfloats("xgradient",gx_inp,N-1)) sf_error("Please enter the x-gradient array [N-1]");
	/* assign x-gradient*/
	
	if (!sf_getfloats("zgradient",gz_inp,N-1)) sf_error("Please enter the z-gradient array [N-1]");
	/* assign z-gradient */
	
	if (!sf_getfloats("xref",xref_inp,N-1)) sf_error("Please enter the x-reference points array [N-1]");
	/* assign x-reference point*/
	
	if (!sf_getfloats("zref",zref_inp,N-1)) sf_error("Please enter the z-reference points array [N-1]");
	/* assign z-reference point*/
	
	if (!sf_getfloat("min",&bmin)) sf_error("Please enter the minimum boundary");
	/* The minimum boundary*/
	
	if (!sf_getfloat("max",&bmax)) sf_error("Please enter the maximum boundary");
	/* The maximum boundary*/
	
	if (!sf_getint("niter",&niter)) sf_error("Please enter the number of iterations");
	/* The number of iterations*/
	
	if (!sf_getdouble("tol",&tol)) tol=0.000001/v_inp[0];
	/* Assign a default value for tolerance*/
	
	
	/*Set output 2D array reflection point---------------------------------------*/
	xrefl = sf_output("out"); /* Output reflection points*/
	
	if (!sf_getint("ns",&nt)) nt=nr2+2;
	/* Dimension of output reflection points (the number of points)*/
	
	if (!sf_getint("ns2",&nt2)) nt2=2; 
	/* Dimension of output reflection points (2 dim)*/
	
	if (!sf_getfloat("ds",&dt)) dt=1; 
	/* Step increment*/
	
	if (!sf_getfloat("s0",&t0)) t0=0;
	/* Staring position*/
	
	sf_putint(xrefl,"n1",nt);
	sf_putint(xrefl,"n2",nt2);
	
	sf_putfloat(xrefl,"d1",dt);
	sf_putfloat(xrefl,"o1",t0);
	
	
	/*read input------------------------------------------------------------------------------*/
	int d1,d2,d3,d4,d5,p1; /*counter*/
	int p2; /*Temp value*/
	
	sf_floatread(temp_rr[0],nr1*N,refl);
	
	/*Check the array, consecutive two inputs must not differ by more than 1*/
	
	for (p1=0; p1<nr2-1; p1++) {
		p2 = updown[p1]-updown[p1+1];
		if (p2>1) {
			sf_warning("The layer number array indicates skipping of layers. Please reenter the array.\n");
			exit(0);
		}
	}
	
	/*generate input according to the reflection sequence-----------------------------*/
	for (d2=0; d2<nr2+2; d2++) {
		
				
		if (d2<nr2+1) { /*Set velocity, gradient, and reference points arrays*/
			
			if (d2 == 0) {
				v[d2] = v_inp[0];
				gx[d2] = gx_inp[0];
				gz[d2] = gz_inp[0];
				xref[d2] = xref_inp[0];
				zref[d2] = zref_inp[0];
 			}
			d3 = updown[d2-1]; /*need d3, d4, and d5 because array argument needs to be an interger*/
			d4 = updown[d2];
			
			if (d4-d3>0) {
				
				v[d2] = v_inp[d3];
				gx[d2] = gx_inp[d3];
				gz[d2] = gz_inp[d3];
				xref[d2] = xref_inp[d3];
				zref[d2] = zref_inp[d3];
			}	
				
			if(d4-d3<0){
			
				v[d2] = v_inp[d4];
				gx[d2] = gx_inp[d4];
				gz[d2] = gz_inp[d4];
				xref[d2] = xref_inp[d4];
				zref[d2] = zref_inp[d4];
				
			}
		}
		
		for (d1=0; d1<nr1; d1++) {
			
			if (d2 == 0) {
				rr[d2][d1] = temp_rr[0][d1];
			}
			d5 = updown[d2-1];
			
			rr[d2][d1] = temp_rr[d5][d1];
		}
		
	}
	
	/*Initialize interpolation-------------------------------------------------------------------*/
	if (!sf_getint("order",&order)) order=3;/*interpolation order*/
	
	
	/*Compute reflector slope--------------------------------------------------------------------*/
	
	eno = (sf_eno*) sf_alloc(nr2+2,sizeof(*eno)); /*allocation for eno array*/
	deno = (sf_eno*) sf_alloc(nr2+2,sizeof(*deno));
	
	for (ir2=0; ir2 < nr2+2; ir2++) {
		
		
		/*Loop through eno*/
		
		eno[ir2]  = sf_eno_init(order,nr1); /*function values*/
		sf_eno_set (eno[ir2],rr[ir2]); 
		
		for (ir1=0; ir1 < nr1; ir1++) {
			
			x = r0+ir1*dr; /* distance */
			rd[ir2][ir1] = zder(ir2,x);
		}
		
		deno[ir2] = sf_eno_init(order,nr1);	/*derivatives*/	
		sf_eno_set (deno[ir2],rr[ir2]);
	}
	
	
	/*Step 1: Calculate F(y) to see if it is sufficiently close to zero-------------------*/
	
	int i,j1,i3,i4; /*counter*/
	float Ftem=0;
	
	for (i3=0; i3<nr2; i3++) {
		xx[i3+1] = xinitial[i3];
	}
	
	
	/*Set initial structure variables' values---------------------------------------------------*/
	
	for (i=0; i<nr2; i++) {
		
		initialize(i+1,xx,v,xref,zref,gx,gz,z,zder,zder2); /*Initialize y_k and y_k1*/
		
		F[i+1] = T_hat_1k_k(T_k_k1,T_k_zk1) + T_hat_k_k(T_k_k,T_k_zk);
		
		/*F[i+1] = traveltime_1k_k(i+1,v[i],xx[i],xx[i+1],xx[i+2],z,zder,zder2) + traveltime_k_k(i+1,v[i+1],xx[i],xx[i+1],xx[i+2],z,zder,zder2);*/
		
		
		for (i4=0; i4<nr2; i4++) { /*check the tolerance*/
			
			Ftem = Ftem+fabsf(F[i4+1]);
			
			if (Ftem<nr2*tol) {
				for (j1=0; j1<nr2; j1++) {
					sf_warning("F(%d) is sufficeintly close to zero. y[%d] = %g \n",j1+1,j1+1,xx[j1+1]);
				}
				goto mark; /*Exit the loop to the writing the result part*/
			}	
			
			if (i4 == nr2-1) {
				Ftem=0;
			}
			
			
		}
		
	}	
	
	/*Loop through the output for repeating yk=yk-dk----------------------------------------------------------------------*/
	
	int q; /*counter for big loop*/
	int i2,j2,i5; /*counter*/
	int w = niter; /*number of loops for yk=yk-dk*/
	
	Ftem=0; /*reset Ftem to zero*/
	
	for (q=0; q<w; q++) {
		
		for (i2=0; i2<nr2; i2++) { /*Recalculate F for new y*/
			
			initialize(i2+1,xx,v,xref,zref,gx,gz,z,zder,zder2);
			
			F[i2+1] = T_hat_1k_k(T_k_k1,T_k_zk1) + T_hat_k_k(T_k_k,T_k_zk);
			
			/*F[i2+1] = traveltime_1k_k(i2+1,v[i2],xx[i2],xx[i2+1],xx[i2+2],z,zder,zder2) + traveltime_k_k(i2+1,v[i2+1],xx[i2],xx[i2+1],xx[i2+2],z,zder,zder2);*/
			
			
			for (i5=0; i5<nr2; i5++) { /*check the tolerance*/
				
				Ftem = Ftem+fabsf(F[i5+1]);
				
				if (Ftem<nr2*tol) {
					for (j2=0; j2<nr2; j2++) {
						sf_warning("F(%d) is sufficeintly close to zero. y[%d] = %g \n",j2+1,j2+1,xx[j2+1]);
					}
					goto mark; /*Exit the loop to the writing the result part*/
				}
				
				if (i5 == nr2-1) {
					Ftem=0;
				}
				
				
			}
			
		}
		
		/*Step 2: Forward recursion-----------------------------------------------------------*/
		
		int l; /*counter*/
		for (l=0; l<nr2; l++) {
			
		initialize(l+1,xx,v,xref,zref,gx,gz,z,zder,zder2);
			
			if (l==0) {
				
				ck_inv[1]= 1/(T_hat_1k_k_k(T_k_k1_k1,T_k_k1_zk1,T_k_zk1,T_k_zk1_zk1) + T_hat_k_k_k(T_k_k_k,T_k_k_zk,T_k_zk,T_k_zk_zk));
				zk[1] = T_hat_1k_k(T_k_k1,T_k_zk1) +T_hat_k_k(T_k_k,T_k_zk);
				
				/*ck_inv[1] = 1/(traveltime_1k_k_k(l+1,v[l],xx[0],xx[1],xx[2],z,zder,zder2) + traveltime_k_k_k(l+1,v[l+1],xx[0],xx[1],xx[2],z,zder,zder2)); */
				/*zk[1] = traveltime_1k_k(l+1,v[l],xx[0],xx[1],xx[2],z,zder,zder2) + traveltime_k_k(l+1,v[l+1],xx[0],xx[1],xx[2],z,zder,zder2); */
			}
			else {
				
				ck_inv[l+1]= 1/(T_hat_1k_k_k(T_k_k1_k1,T_k_k1_zk1,T_k_zk1,T_k_zk1_zk1) + T_hat_k_k_k(T_k_k_k,T_k_k_zk,T_k_zk,T_k_zk_zk) - T_hat_1k_1k_k(T_k_k_k1,T_k_k1_zk,T_k_k_zk1,T_k_zk_zk1)*ck_inv[l]*T_hat_1k_1k_k(T_k_k_k1,T_k_k1_zk,T_k_k_zk1,T_k_zk_zk1));
				zk[l+1] = T_hat_1k_k(T_k_k1,T_k_zk1) + T_hat_k_k(T_k_k,T_k_zk) - T_hat_1k_1k_k(T_k_k_k1,T_k_k1_zk,T_k_k_zk1,T_k_zk_zk1)*ck_inv[l]*zk[l];
				
				/*ck_inv[l+1] = 1/(traveltime_1k_k_k(l+1,v[l],xx[l],xx[l+1],xx[l+2],z,zder,zder2) + traveltime_k_k_k(l+1,v[l+1],xx[l],xx[l+1],xx[l+2],z,zder,zder2)-traveltime_1k_1k_k(l+1,v[l],xx[l],xx[l+1],xx[l+2],z,zder,zder2)*ck_inv[l]*traveltime_1k_1k_k(l+1,v[l],xx[l],xx[l+1],xx[l+2],z,zder,zder2));*/
				/*zk[l+1] = traveltime_1k_k(l+1,v[l],xx[l],xx[l+1],xx[l+2],z,zder,zder2) + traveltime_k_k(l+1,v[l+1],xx[l],xx[l+1],xx[l+2],z,zder,zder2)-traveltime_1k_1k_k(l+1,v[l],xx[l],xx[l+1],xx[l+2],z,zder,zder2)*ck_inv[l]*zk[l];*/	
			}
			
			if (isnan(1/ck_inv[l+1]) != 0 || isinf(1/ck_inv[l+1]) != 0) {
				sf_warning("ck_inv doesn't exist. The solutions do not converge.\n");
				exit(0);
			}
			
		}	
		
		/*Step 3: Backward recursion----------------------------------------------------------*/
		int m,u,w; /*counter*/
		for (m=nr2-1; m>=0; m--) { 
			
			initialize(m+1,xx,v,xref,zref,gx,gz,z,zder,zder2);
			
			if (m==nr2-1) {
				dk[m+1] = ck_inv[m+1]*zk[m+1];
			}
			else {
				dk[m+1] = ck_inv[m+1]*(zk[m+1]-T_hat_k_k_k1(T_k_k_k1,T_k_k1_zk,T_k_k_zk1,T_k_zk_zk1)*dk[m+2]);
				/*dk[m+1] = ck_inv[m+1]*(zk[m+1]-traveltime_k_k_k1(m+1,v[m+1],xx[m],xx[m+1],xx[m+2],z,zder,zder2)*dk[m+2]);*/
			}
			
		}
		
		/*Check whether the increment dk is diverging or not and multiply by 0.5 if it is*/
		
		if (q==0) {
			for (w=0; w<nr2+2; w++) {
				*(dk_old+w) = *(dk+w); /*The first loop >> set dk_old = dk*/
			}
		}
		
		for (u=0; u<nr2; u++) {
			while (fabsf(dk[u+1])>fabsf(dk_old[u+1])) {
				dk[u+1] = 0.5*dk[u+1];
				sf_warning("The increment dk[%d] diverges and is halved to be %g\n",u+1,dk[u+1]);
				if (fabsf(dk[u+1])<tol) {
					break;
				}
			}
		} 
		
		
		/*Step 4: Update xx--------------------------------------------------------------------*/
		int t,a,b1,b2; /*counter*/
		
		/*Apply boundary------------------------------------------*/
		
		for (a=0; a<nr2; a++) {
			b1=0;
			b2=0;
			xxtem[a+1] = xx[a+1]-dk[a+1];
			
			if (xxtem[a+1]<bmin) {
				switch (b1) {
					case 20: /*maximum times to multiply is 20*/
						break;
					default:
						dk[a+1]=0.5*dk[a+1];
						sf_warning("The new y value exceeds the minimum boundary. dk[%d] is reduced to %g\n",a+1,dk[a+1]);
						xxtem[a+1] = xx[a+1]-dk[a+1];
						b1++;
						break;
				}
			}
			else if (xxtem[a+1]>bmax) {
				switch (b2) {
					case 20: /*maximum times to multiply is 20*/
						break;
					default:
						dk[a+1]=0.5*dk[a+1];
						sf_warning("The new y value exceeds the minimum boundary. dk[%d] is reduced to %g\n",a+1,dk[a+1]);
						xxtem[a+1] = xx[a+1]-dk[a+1];
						b2++;
						break;
				}
			}
			
		}
		
		/*--------------------------------------------------------*/
		
		vector_sub(nr2+2,xx,nr2+2,dk,xxnew,0);
		
		for (t=0; t<nr2+2;t++) {
			sf_warning("The original value of y[%d] is %g, d[%d] is %g and the new value of y[%d] is %g\n",t,*(xx+t),t,*(dk+t),t,*(xxnew+t));
			if (t==nr2+1) {
				printf("\n\n");
			}	
			*(xx+t) = *(xxnew+t); /*update xx values*/
			*(dk_old+t) = *(dk+t); /*store the working dk*/
			
			
		}
		
		
	}
	
	
	/*Write result in 2D & Compute traveltime-----------------------------------------------------------------------*/
	int c1,c2;
	int c; /*counter*/
	
mark: /*mark point for goto*/
	
	for (c2=0; c2<nr2+2; c2++) {
		for (c1=0; c1<2; c1++) {
			
			if (c1==0) {
				ans[c1][c2] = xx[c2];
			}
			else {
				ans[c1][c2] = z(c2,xx[c2]);
			}
			
		}
	}
	
	tt=0; /*initialize traveltime tt*/
	
	for (c=0; c<nr2+1; c++) {
		
		half_initialize(c,xx,v,xref,zref,gx,gz,z,zder,zder2);
		
		tt = tt + T_hat_k(T_k);
		/*tt = tt + traveltime_k(c, v[c],xx[c],xx[c+1],z,zder,zder2);*/
		
		if (c==nr2) {
			sf_warning("Traveltime is %g",tt);
		}
	}
	
	/* write output */
	sf_floatwrite(ans[0],2*(nr2+2),xrefl);
	
	exit(0);
}
