/* Pre-stack bending ray tracing in one-layered media
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

static float zder2_eno(float x) 
/* second derivative */
{
	int i;
	float f, f1;
	
	x = (x-r0)/dr; 
	i = floorf(x);
	
	sf_eno_apply (deno,i,x-i,&f,&f1,DER);
	return f1/dr;	
}

static float zder2_cen(float x) 
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

int main(int argc, char* argv[])
{
	int nt, it, ir, order, niter, nt2, it2, type, print;
	float x, t0, dt, velocity, tol, xinitial=0, xmid,dt2,t02; /*,max_extent, scale,max_x,xmax,hmax;*/
	float zd1, zd2_eno,zd2_cen /*zd2_real,rel_cen,rel_eno,zd1_real*/; /*temp second derivative value*/
	bool stop;
	func1 zder2=0;
	
	float *rd, **tt, **xx; 
	double xs,xr;
	sf_file refl, ttime, xrefl;
	
	sf_init(argc,argv); /* initialize - always call first */
	
	/* Set input */
	refl = sf_input("in"); /* reflector */
	if (!sf_histint(refl,"n1",&nr)) sf_error("No n1= in input");
	if (!sf_histfloat(refl,"o1",&r0)) r0=0.;
	if (!sf_histfloat(refl,"d1",&dr)) dr=1.;
	
	/* Set output 2D traveltime and 2D reflection point*/
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
	
	xrefl = sf_output("xrefl"); /* Output reflection point*/
	
	sf_putint(xrefl,"n1",nt2); /*S&R increment axis (Faster axis)*/
	sf_putfloat(xrefl,"d1",dt2);
	sf_putfloat(xrefl,"o1",t02);	
	
	sf_putint(xrefl,"n2",nt); /* Midpoint axis*/
	sf_putfloat(xrefl,"d2",dt);
	sf_putfloat(xrefl,"o2",t0);
	
	/* Allocate space */
	rr = sf_floatalloc(nr);
	rd = sf_floatalloc(nr);
	
	tt = sf_floatalloc2(nt2,nt);/* Row is always for midpoint and column is always for offset*/
	xx = sf_floatalloc2(nt2,nt);
	
	/* read input */
	sf_floatread(rr,nr,refl);
	
	/* Initialize interpolation */
	if (!sf_getint("print",&type)) print=0;/* Print the actual calculated values by sf_eno and central finite difference 0=No and 1=Yes*/
	
	if (!sf_getint("type",&type)) type=1;/* Interpolation type 0=sf_eno and 1=central finite difference*/
	
	if (!sf_getint("order",&order)) order=4;/* Interpolation order if choose to use sf_eno*/
	
	if (!sf_getfloat("velocity",&velocity)) velocity=2.0;/* Assign velocity km/s*/
	
	if (!sf_getfloat("tol",&tol)) tol=1/(1000000*velocity);/* Assign a default value for tolerance*/
	
	if (!sf_getbool("break",&stop)) stop=false;/* Go beyond zero or not*/
	
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
	
	
	if (type==0) { /*set the type of interpolation to use*/
		zder2 = zder2_eno;
	} else if (type==1) {
		zder2 = zder2_cen;
	}

	
	
	/* Loop through the output */
	for (it=0; it < nt; it++){ /*How many times we move the location (Midpoint) =>Each Row*/
		
		xmid = t0+it*dt;
		
		for (it2=0; it2 < nt2; it2++) { /*How many times we compute xs&xr at each location =>Each Column*/
			xs = xmid-it2*dt2; /* source location */
			xr = xmid+it2*dt2; /* receiver location */
			
			if ((xs<0 || xr<0) && stop) { /* To stop thr program from going beyond zero*/
				break;
			}
			
			traveltime_init(z,zder,zder2,xs,xr,velocity);
			
				if (it==0) xinitial = xmid;  /*The first calculation starting from the origin*/
				
				x = newton(dtdx,d2tdx2,xinitial,niter,tol);
				xinitial = x; /*Use the answer of the previous iteration to be the starting point of the next iteration*/
				
				if (print==1) {
					zd1 = zder(xinitial);
					zd2_eno = zder2_eno(xinitial);
					zd2_cen = zder2_cen(xinitial);
				
					sf_warning("x=%g, zder1=%g, zder2_eno=%g, zder2_cen=%g",xinitial,zd1,zd2_eno,zd2_cen);
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
