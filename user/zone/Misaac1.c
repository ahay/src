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
#include "ps_traveltime.h"
#include "newton.h"

static sf_eno eno, deno; /* interpolation structure */	
static float r0, dr;

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

int main(int argc, char* argv[])
{
	int nr, nt, it, ir, order, niter, nt2, it2;
	float x, t0, dt, velocity, tol, xinitial,dt2,t02; /*,max_extent, scale,max_x,xmax,hmax;*/
	float *rr, *rd, **tt, **xx; 
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
	
	sf_putint(ttime,"n1",nt); /* Midpoint axis*/
	sf_putfloat(ttime,"d1",dt);
	sf_putfloat(ttime,"o1",t0);
	
	sf_putint(ttime,"n2",nt2); /*S&R increment axis*/
	sf_putfloat(ttime,"d2",dt2);
	sf_putfloat(ttime,"o2",t02);
	
	xrefl = sf_output("xrefl"); /* Output reflection point */
	sf_putint(xrefl,"n1",nt);
	sf_putfloat(xrefl,"d1",dt);
	sf_putfloat(xrefl,"o1",t0);	
	
	/* Allocate space */
	rr = sf_floatalloc(nr);
	rd = sf_floatalloc(nr);
	
	tt = sf_floatalloc2(nt,nt2);
	xx = sf_floatalloc2(nt,nt2);
	
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
	for (it2=0; it2<nt2; it2++){ /*How many times we move the location => column*/
		
		
		for (it=0; it < nt; it++) { /*How many times we compute xs&xr at each location => row*/
			xs = t0+it2*dt-it*dt2; /* source location */
			xr = t0+it2*dt+it*dt2; /* receiver location */
			
			traveltime_init(z,zder,zder2,xs,xr,velocity);
			
			if (it==0) {
				xinitial = t0+it2*dt; /*The first calculation starting from the origin*/
				x = newton(dtdx,d2tdx2,xinitial,niter,tol);
				xinitial = x; /*Use the answer of the previous iteration to be the starting point of the next iteration*/
			}
			
			tt[it][it2] = traveltime(x);
			xx[it][it2] = x;
		}
		
	}	
	
	/* write output */
	sf_floatwrite(tt[0],nt*nt2,ttime); 
	sf_floatwrite(xx[0],nt*nt2,xrefl);
	
	exit(0);
}
