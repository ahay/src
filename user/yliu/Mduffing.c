/* Duffing differential equation solved by 4th order Runge-Kutta method. 
Duffing equation: x'' + d x' - a x + b x^3 = gamma H(2 PI freq t) + alpha input(t)
*/
/*
  Copyright (C) 2013 Jilin University
  
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
#include <math.h>

static float stochasforce(int i, int n, int m, float deltat, float omega, bool ricker)
/* stochastic force function */
{
    float h;
    int j;
    if(ricker) {
	h = 0.;
	for(j=0; j < n; j++) {
	    h += (1.-0.5*omega*omega*(i-0.5*m-m*j)*(i-0.5*m-m*j)*deltat*deltat)
		*expf(-0.25*omega*omega*(i-0.5*m-m*j)*(i-0.5*m-m*j)
		      *deltat*deltat);
	}
    } else {
	h = cos(omega*i*deltat);
    }
    return h;
}

int main (int argc, char* argv[]) 
{
    int nt, n2, i2, i1, p1, p2, m, n;
    float dt, ot, x0, y0, d, a, b, gamma, alpha, freq, omega;
    float x1, y1, f1, x2, y2, f2, x3, y3, f3, x4, y4, f4;
    bool ricker;
    float *x, *y, *input;
    sf_complex *result;

    sf_file in, out;
    
    sf_init (argc, argv); 
    in = sf_input("in");
    out = sf_output("out");

    sf_setformat(out,"native_complex");
    
    if (!sf_histint(in,"n1",&nt)) sf_error("No n1= in input");
    if (!sf_histfloat(in,"d1",&dt)) sf_error("No d1= in input");
    if (!sf_histfloat(in,"o1",&ot)) sf_error("No o1= in input");
    n2 = sf_leftsize(in,1);

    if (!sf_getfloat("x0",&x0)) x0=0.;
    /* intial x coordinate */
    if (!sf_getfloat("y0",&y0)) y0=0.;
    /* intial y coordinate */

    if (!sf_getfloat("d",&d)) d=1.;
    /* Duffing equation coefficients d */
    if (!sf_getfloat("a",&a)) a=1.;
    /* Duffing equation coefficients a */
    if (!sf_getfloat("b",&b)) b=1.;
    /* Duffing equation coefficients b */
    if (!sf_getfloat("gamma",&gamma)) gamma=1.;
    /* Duffing equation coefficients gamma */
    if (!sf_getfloat("alpha",&alpha)) alpha=0.;
    /* Duffing equation coefficients alpha */

    if (!sf_getbool("ricker",&ricker)) ricker=false;
    /* if y, stochastic force function is Ricker, otherwise, cosine */

    if(ricker) {
	if (!sf_getint("n",&n)) n=40;
	if (!sf_getint("m",&m)) m=80;
	/* Parameters for Ricer wavelet force function */
    } else {
	n = 1;
	m = 1;
    }

    if (!sf_getfloat("freq",&freq)) freq=1.;
    /* stochastic force coefficients freq */
    omega = 2.*SF_PI*freq;

    if (!sf_getint("p1",&p1)) p1=1.;
    /* coefficients p1 */
    if (!sf_getint("p2",&p2)) p2=3.;
    /* coefficients p2 */

    x = sf_floatalloc(nt);
    y = sf_floatalloc(nt);
    input = sf_floatalloc(nt);

    result = sf_complexalloc(nt);

    x[0] = x0;
    y[0] = y0;

    result[0] = sf_cmplx(x[0],y[0]);

    for(i2=0; i2 < n2; i2++) {
	sf_floatread(input,nt,in);

	for(i1=0; i1 < nt-1; i1++) {
	    x1 = x[i1];
	    y1 = y[i1];
	    f1 = -1.*d*y1+a*powf(x1,p1)-b*powf(x1,p2)
		+gamma*stochasforce(i1,n,m,dt,omega,ricker)
		+alpha*input[i1];
	    x2 = x[i1]+y1*dt/2.;
	    y2 = y[i1]+f1*dt/2.;
	    f2 = -1.*d*y2+a*powf(x2,p1)-b*powf(x2,p2)
		+gamma*stochasforce(i1,n,m,dt,omega,ricker)
		+alpha*input[i1];
	    x3 = x[i1]+y2*dt/2.;
	    y3 = y[i1]+f2*dt/2.;
	    f3 = -1.*d*y3+a*powf(x3,p1)-b*powf(x3,p2)
		+gamma*stochasforce(i1+1,n,m,dt,omega,ricker)
		+alpha*input[i1+1];
	    x4 = x[i1]+y3*dt;
	    y4 = y[i1]+f3*dt;
	    f4 = -1.*d*y4+a*powf(x4,p1)-b*powf(x4,p2)
		+gamma*stochasforce(i1+1,n,m,dt,omega,ricker)
		+alpha*input[i1+1];
	    x[i1+1] = x[i1]+dt/6.*(y1+2*y2+2*y3+y4);
	    y[i1+1] = y[i1]+dt/6.*(f1+2*f2+2*f3+f4);
	    result[i1+1] = sf_cmplx(x[i1+1],y[i1+1]);
	}
	sf_complexwrite(result,nt,out);

    }

    exit (0);
}

/* 	$Id$	 */
