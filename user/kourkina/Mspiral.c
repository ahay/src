/* Spiral function */
/*
  Copyright (C) 2008 New York University
  
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
#include <math.h>

#include <rsf.h>

static float v( float x, float y ); /* slowness function */
static float solve_newton(float (*myf)(float,float),float (*dmyf)(float),
			   float alpha,float x0,float meps);
static float newton_th0(float rad,float theta,float x0,float meps); 
static float fun(float,float);
static float dfun(float);
static float rho(float);

static float xc, yc, eps, v0, v1, b, r0, r1, fac, sp_r, sp_t;

static float v( float x, float y ) {
    float cx ,cy, ang, tau, r, d, th0;
    cx=x-xc;
    cy=y-yc;
    ang=atan2f(cx,cy);
    d=hypotf(x-xc,y-yc);
    th0=newton_th0(d,ang,ang,eps);
    tau=solve_newton(fun,dfun,th0,th0,eps);
    r=rho(tau);
    return v0+v1*expf(-b*d*d/(r*r));
}

/*--------------------------------------------*/
float rho(float tau) {
    return r0-r1*cosf(3.0*tau);
}

/*--------------------------------------------*/
float fun(float alpha,float tau) {
    return alpha-tau+fac*cosf(3.0*tau);
}

/*--------------------------------------------*/
float dfun(float tau) {
  return -1.0-3.0*fac*sin(3.0*tau);
}

/*--------------------------------------------*/
static float solve_newton(float (*myf)(float,float),float (*dmyf)(float),
			  float alpha,float x0,float meps){
    float delta=1.0e+4, f0,fp,x;
    while( delta>meps ) {
	f0=fun(alpha,x0);
	fp=dfun(x0);
	x=x0-f0/fp;
	delta=fabs(x-x0);
	x0=x;
    }
    return x0;
}

/*---------------------------------------------*/
static float newton_th0(float rad,float theta,float x0,float meps)
{
    float delta=1.0e+4, f0, fp, x, rad0, tau, aux;
    
    x=x0;
    while( delta>meps ) {
	tau=solve_newton(fun,dfun,x0,x0,meps);
	rad0=rho(tau);
	f0=rad-rad0-sp_r*(theta-x0)/sp_t;
	aux=3.0*sin(3.0*tau);
	fp=r1*aux/(1.0+fac*aux)+sp_r/sp_t;
	x=x0-f0/fp;
	delta=fabs(x-x0);
	x0=x;
    }
    return x;
}

/************* F U N C T I O N   M A I N **************/
int main(int argc, char* argv[]) 
{
    int nx, ny, ix, iy;
    float dx, dy, xmin, ymin, xmax, ymax, *vmesh;
    sf_file spiral;

    sf_init(argc,argv);
    spiral = sf_output("out");

    sf_setformat(spiral,"native_float");
    if (!sf_getint("nx",&nx)) nx=500;
    if (!sf_getint("ny",&ny)) ny=200;

    if (!sf_getfloat("xmax",&xmax)) xmax=20.;
    if (!sf_getfloat("ymax",&ymax)) ymax=3.;

    if (!sf_getfloat("xmin",&xmin)) xmin=0.;
    if (!sf_getfloat("ymin",&ymin)) ymin=0.;

    if (!sf_getfloat("dx",&dx)) dx=(xmax-xmin)/(nx-1);
    if (!sf_getfloat("dy",&dy)) dy=(ymax-ymin)/(ny-1);

    if (!sf_getfloat("xc",&xc)) xc=10.;
    if (!sf_getfloat("yc",&yc)) yc=5.;

    if (!sf_getfloat("eps",&eps)) eps=1.0e-6;
    if (!sf_getfloat("v0",&v0)) v0=2.0;
    if (!sf_getfloat("v1",&v1)) v1=2.0;
    if (!sf_getfloat("r0",&r0)) r0=1.0;
    if (!sf_getfloat("r1",&r1)) r1=0.4;
    /* paramters of original shape */
    if (!sf_getfloat("b",&b)) b=0.1;
    /* exponential decay factor */
    if (!sf_getfloat("fac",&fac)) fac=0.2;
    if (!sf_getfloat("sp_r",&sp_r)) sp_r=1.;
    /* speed in radius */
    if (!sf_getfloat("sp_t",&sp_t)) sp_t=0.05;
    /* speed in angle */

    sf_putint(spiral,"n1",nx);
    sf_putint(spiral,"n2",ny);

    sf_putfloat(spiral,"d1",dx);
    sf_putfloat(spiral,"d2",dy);

    sf_putfloat(spiral,"o1",xmin);
    sf_putfloat(spiral,"o2",ymin);

    vmesh = sf_floatalloc(nx);

    for( iy=0; iy<ny; iy++ ) {
	for( ix=0; ix<nx; ix++ ) {
	    vmesh[ix] = v(xmin+ix*dx,ymin+iy*dy);
	}
	sf_floatwrite(vmesh,nx,spiral);
    }
    
    exit(0);
}




