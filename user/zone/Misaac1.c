/* Reflection traveltime for Pre-stack */
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
#include <math.h>
#include "ps_traveltime.h"
#include "newton.h"

#define pi 3.14159265358979323846

static sf_eno eno, deno; /* interpolation structure */	
static float r0, dr,*rr;
static int nr;

static float z(float x) 
/* function */
{
	int i;
	float f, f1;
	
	x = (x-r0)/dr; 
	i = floorf(x);
	
	sf_eno_apply (eno,i,x-i,&f,&f1,FUNC);
	return f;
}

static float zder(float x) 
/* first derivative */
{
	int i;
	float f, f1;
	
	x = (x-r0)/dr; 
	i = floorf(x);
	
	sf_eno_apply (eno,i,x-i,&f,&f1,DER);
	return f1/dr;
}

static float zder2(float x) 
/* second derivative */
{
	int i;
	float f, f1;
	
	x = (x-r0)/dr; 
	i = floorf(x);
	
	sf_eno_apply (deno,i,x-i,&f,&f1,DER);
	return f1/dr;	
}

static float zder22(float x) 
/* second derivative from finite difference */
{
	int i;
	float f;
	
	x = (x-r0)/dr; 
	i = floorf(x);

	if (i<=1) {
		
		f = (2*rr[i]-5*rr[i+1]+4*rr[i+2]-rr[i+3])/(dr*dr);/*2nd order Right-sided finite difference for left end*/
		
	} else if (i>=nr-2){
		
		f = (2*rr[i]-5*rr[i-1]+4*rr[i-2]-rr[i-3])/(dr*dr);/*2nd order Left-sided finite difference for right end*/
		
	} else {
		f = ((-1)*rr[i+2]+16*rr[i+1]-30*rr[i]+16*rr[i-1]-rr[i-2])/(12*dr*dr); /*4th order Central finite difference for middle points*/
	}

	return f;
}

static float d1(float x,float h,float alpha)
/*Second derivative for hyperbolic reflector*/
{
	float f;
	
	f = x*pow(tan(alpha),2)/hypot(h,x*tan(alpha));
	
	return f;
}

static float d2(float x,float h,float alpha)
/*Second derivative for hyperbolic reflector*/
{
	float f;
	
	f = pow(h,2)*pow(tan(alpha),2)/pow(hypot(h,x*tan(alpha)),3);
	
	return f;
}

int main(int argc, char* argv[])
{
	int nt, it, ir, order, niter, nt2, it2;
	float x, t0, dt, velocity, tol, xinitial,dt2,t02; /*,max_extent, scale,max_x,xmax,hmax;*/
	float zd1,zd1_real, zd2, zd2_real; /*temp second derivative value*/
	
	float *rd, **tt, **xx; 
	double xs,xr;
	sf_file refl, ttime, xrefl;
	
	sf_init(argc,argv); /* initialize - always call first */
	
	/* Set input */
	refl = sf_input("in"); /* reflector */
	if (!sf_histint(refl,"n1",&nr)) sf_error("No n1= in input");
	if (!sf_histfloat(refl,"o1",&r0)) r0=0.;
	if (!sf_histfloat(refl,"d1",&dr)) dr=1.;
	
	/* Set output 2D traveltime and 1D reflection point*/
	ttime = sf_output("out"); /* Output traveltime */
	
	if (!sf_getint("ns",&nt)) nt=nr; /* number of sources for midpoint*/
	if (!sf_getint("ns2",&nt2)) nt2=nr; /* number of sources for offset*/

	if (!sf_getfloat("ds",&dt)) dt=dr; 	/* source sampling for midpoint*/
	if (!sf_getfloat("ds2",&dt2)) dt2=dr; 	/* source sampling for offset*/
	
	if (!sf_getfloat("s0",&t0)) t0=r0;/* origin for midpoint*/
	if (!sf_getfloat("s02",&t02)) t02=r0;/* origin for offset*/
	
	
	sf_putint(ttime,"n1",nt2); /*S&R increment axis (Faster axis)*/
	sf_putfloat(ttime,"d1",dt2);
	sf_putfloat(ttime,"o1",t02);
	
	sf_putint(ttime,"n2",nt); /* Midpoint axis*/
	sf_putfloat(ttime,"d2",dt);
	sf_putfloat(ttime,"o2",t0);
	
	xrefl = sf_output("xrefl"); /* Output reflection point */
	sf_putint(xrefl,"n1",nt);
	sf_putfloat(xrefl,"d1",dt);
	sf_putfloat(xrefl,"o1",t0);	
	
	/* Allocate space */
	rr = sf_floatalloc(nr);
	rd = sf_floatalloc(nr);
	
	tt = sf_floatalloc2(nt2,nt);/*Row is always for midpoint and column is always for offset*/
	xx = sf_floatalloc2(nt2,nt);
	
	/* read input */
	sf_floatread(rr,nr,refl);
	
	/* Initialize interpolation */
	if (!sf_getint("order",&order)) order=3;/*interpolation order*/
	
	
	if (!sf_getfloat("velocity",&velocity)) velocity=2.0;/*assign velocity km/s*/
	
	
	if (!sf_getfloat("tol",&tol)) tol=1/(1000000*velocity);/* assign a default value for tolerance*/
	
	
	/*if (!sf_getfloat("max_extent",&max_extent)) hmax=500;
	assign max extent from the center of the xs and xr in meter*/
	
	/*if (!sf_getfloat("max_x",&max_x)) xmax=500;
	assign how much x0 go in meter*/
	
	/*if (!sf_getfloat("scale",&scale)) scale=5.0;
	assign the scale of h and x0 for each increment in meter*/
	
	eno  = sf_eno_init(order,nr);
	sf_eno_set (eno,rr);
	
	/* compute reflector slope */
	for (ir=0; ir < nr; ir++) {
		x = r0+ir*dr; /* distance */
		rd[ir] = zder(x);
	}
	deno = sf_eno_init(order,nr);		
	sf_eno_set (deno,rd);
	
	
	niter = 60000; /* number of iterations for Newton's method*/
	
	
	/* Loop through the output */
	for (it=0; it < nt; it++){ /*How many times we move the location (Midpoint) =>Each Row*/
		
		
		if (it >= 200) {
			sf_warning("Line = %d\n",it);
		}
		
		for (it2=0; it2 < nt2; it2++) { /*How many times we compute xs&xr at each location =>Each Column*/
			xs = t0+it*dt-it2*dt2; /* source location */
			xr = t0+it*dt+it2*dt2; /* receiver location */
			
			traveltime_init(z,zder,zder22,xs,xr,velocity);
			
			if (it2==0) {
				xinitial = t0+it*dt; /*The first calculation starting from the origin*/
				x = newton(dtdx,d2tdx2,xinitial,niter,tol);
				xinitial = x; /*Use the answer of the previous iteration to be the starting point of the next iteration*/
				
				zd1 = zder(xinitial);
				zd1_real = d1(xinitial,1, 45*pi/180.0);
				zd2 = zder22(xinitial);
				zd2_real = d2(xinitial,1, 45*pi/180.0);
				
				sf_warning("x=%g,zder1=%g,zd1_real=%g \n\tzder2=%g,zd2_real=%g",xinitial,zd1,zd1_real,zd2,zd2_real);
			}
			
			tt[it][it2] = traveltime(x);
			xx[it][it2] = x;
		}
		
	}	
	
	
	/* write output */
	sf_floatwrite(*tt,nt*nt2,ttime); 
	sf_floatwrite(*xx,nt*nt2,xrefl);
	
	exit(0);
}
