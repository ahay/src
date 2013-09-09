/* Zero-offset bending ray tracing in one-layered media
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
#include "traveltime.h"
#include "newton.h"

static sf_eno eno, deno; /* interpolation structure */	
static float r0, dr;

static float z(float x) 
/* first derivative */
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
	int nr, nt, it, ir, order, niter;
	float x, t0, dt, velocity, tol, xinitial;
	float *rr, *rd, *tt, *xx; 
	sf_file refl, ttime, xrefl;
	
	sf_init(argc,argv); /* initialize - always call first */
	
	/* Set input */
	refl = sf_input("in"); /* reflector */
	if (!sf_histint(refl,"n1",&nr)) sf_error("No n1= in input");
	if (!sf_histfloat(refl,"o1",&r0)) r0=0.;
	if (!sf_histfloat(refl,"d1",&dr)) dr=1.;
	

	if (!sf_getint("ns",&nt)) nt=nr;/* Number of sources */
	
	
	if (!sf_getfloat("ds",&dt)) dt=dr;/* source sampling */
	

	if (!sf_getfloat("s0",&t0)) t0=r0;/* source origin */
	
	/* Set output */
	ttime = sf_output("out"); /* Output traveltime */	
	sf_putint(ttime,"n1",nt);
	sf_putfloat(ttime,"d1",dt);
	sf_putfloat(ttime,"o1",t0);
	
	xrefl = sf_output("xrefl"); /* Output reflection point */
	sf_putint(xrefl,"n1",nt);
	sf_putfloat(xrefl,"d1",dt);
	sf_putfloat(xrefl,"o1",t0);	
	
	/* Allocate space */
	rr = sf_floatalloc(nr);
	rd = sf_floatalloc(nr);
	
	tt = sf_floatalloc(nt);
	xx = sf_floatalloc(nt);
	
	/* read input */
	sf_floatread(rr,nr,refl);
	
	/* Initialize interpolation */
	if (!sf_getint("order",&order)) order=3;/* interpolation order */

	if (!sf_getfloat("velocity",&velocity)) velocity=2.0;/* assign velocity km/s*/
	
	if (!sf_getfloat("tol",&tol)) tol=0.0001/velocity;/* assign a default value for tolerance*/
	
	eno  = sf_eno_init(order,nr);
	sf_eno_set (eno,rr);

	/* Compute reflector slope */
	for (ir=0; ir < nr; ir++) {
		x = r0+ir*dr; /* distance */
		rd[ir] = zder(x);
	}
	deno = sf_eno_init(order,nr);		
	sf_eno_set (deno,rd);
	
	
	niter = 6000; /* number of iterations for Newton's method*/
	
	
	/* Loop through the output */
	for (it=0; it < nt; it++) {
		x = t0+it*dt; /* source  and receiver location */
		
		traveltime_init(z,zder,zder2, x,velocity);
		
		if (it==0) {
			xinitial = t0; /*The first calculation starting from the origin*/
		}
		x = newton(dtdx,d2tdx2,xinitial,niter,tol);
		xinitial = x; /*Use the answer of the previous iteration to be the starting point of the next iteration*/
		tt[it] = traveltime(x);
		xx[it] = x;
	}
	
	
	/* write output */
	sf_floatwrite(tt,nt,ttime);
	sf_floatwrite(xx,nt,xrefl);
	
	exit(0);
}
