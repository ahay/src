/* Convert interval velocity to Dix velocity */
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

struct darr {
  double x;
  double y;
  double a;
  double t;
};

static void shooting(void);
static struct darr rkm( float x0, float y0, float a0, float t0, float hs ); /* 4th order Runge-Kutta method */
static void rhs( float x, float y, float a, float *kx, float *ky, float *ka );
static void s(float x, float y);

static float xmax, ymax, ht, hx, hy, v, gradv[2];
static int nt, nx, ny;
static float **sc;
static float *w,*xr,*yr;
static sf_eno2 pnt;

static void shooting(void) 
{
    int i, j, ind;
    float xs, ys, as, ts;
    struct darr arr;

    for( i=0; i<nx; i++ ) {
	*(w+i)=sc[0][i];
	*(xr+i)=i*hx;
	*(yr+i)=0.0;
	for( j=1; j<nt; j++ ) {
	    ind=i+nx*j;
	    *(w+ind)=0.0;
	    *(xr+ind)=0.0;
	    *(yr+ind)=0.0;
	}
    }

    for( i=0; i<nx; i++ ) {
	xs=i*hx;
	ys=0.0;
	as=0.0;
	ts=0.0;

	for( j=1; j<nt; j++) {
	    arr=rkm(xs,ys,as,ts,ht);
	    xs=arr.x;
	    ys=arr.y;
	    as=arr.a;
	    ts=arr.t;

	    ind=i+nx*j; 

	    if( xs>=0.0 && xs<=xmax && 
		ys>=0.0 && ys<=ymax ) {
		xr[ind]=xs;
		yr[ind]=ys;
		s(xs,ys);
		w[ind]=v;
	    } else {
		for(; j<nt; j++) {
		    ind=i+nx*j;
		    xr[ind]=xs;
		    yr[ind]=ys;
		    w[ind]=v;
		}
	    }
	}
    }
}

static struct darr rkm( float x0, float y0, float a0, 
			float t0, float hs ) 
/* Runge-Kutta Method */
{
    float k1x,k2x,k3x,k4x,k1y,k2y,k3y,k4y,k1a,k2a,k3a,k4a;
    struct darr tmp;

    rhs(x0,y0,a0,&k1x,&k1y,&k1a);
    rhs(x0+0.5*hs*k1x,
	y0+0.5*hs*k1y,
	a0+0.5*hs*k1a,&k2x,&k2y,&k2a);
    rhs(x0+0.5*hs*k2x,y0+0.5*hs*k2y,a0+0.5*hs*k2a,&k3x,&k3y,&k3a);
    rhs(x0+hs*k3x,y0+hs*k3y,a0+hs*k3a,&k4x,&k4y,&k4a);

    tmp.x=x0+hs*(k1x+2.0*k2x+2.0*k3x+k4x)/6.;
    tmp.y=y0+hs*(k1y+2.0*k2y+2.0*k3y+k4y)/6.;
    tmp.a=a0+hs*(k1a+2.0*k2a+2.0*k3a+k4a)/6.;
    tmp.t=t0+hs;

    return tmp;
}

/************* R I G H T   H A N D   S I D E **************/

static void rhs( float x, float y, float a, float *kx, float *ky, float *ka) 
/* R I G H T   H A N D   S I D E */
{
    s(x,y);
    *kx = sinf(a)*v;
    *ky = cosf(a)*v;
    *ka = gradv[1]*sin(a)/hy-gradv[0]*cos(a)/hx;
}

/*******************************************************/
 
static void s( float xs, float ys ) 
/* compute velocity with gradient */
{
    int ix, iy;

    xs /= hx; ix = floorf(xs); xs -= ix;
    ys /= hy; iy = floorf(ys); ys -= iy;

    sf_eno2_apply (pnt, ix, iy, xs, ys, &v, gradv, BOTH);
}

/******************************************************/

int main(int argc, char* argv[]) 
{
    sf_file fid=NULL, fx=NULL, fy=NULL, fsc=NULL;
    int i, j, ind, order;
    int nx1, ny1, ntx;
    float q;

    sf_init(argc,argv);
    fsc = sf_input("in");
    fid = sf_output("out");
    fx = sf_output("x");
    fy = sf_output("z");

    if (!sf_histint(fsc,"n1",&nx)) sf_error("No n1= in input");
    if (!sf_histint(fsc,"n2",&ny)) sf_error("No n2= in input");
    if (!sf_histfloat(fsc,"d1",&hx)) sf_error("No d1= in input");
    if (!sf_histfloat(fsc,"d2",&hy)) sf_error("No d2= in input");

    nx1=nx-1;
    ny1=ny-1;
    xmax = nx1*hx;
    ymax = ny1*hy;

    sc = sf_floatalloc2(nx,ny);
    sf_floatread(sc[0],nx*ny,fsc);

    if (!sf_getint("nt",&nt)) sf_error("Need nt=");
    if (!sf_getfloat("dt",&ht)) sf_error("Need dt=");

    if (!sf_getint("order",&order)) order=4;
    /* interpolation order */

    pnt = sf_eno2_init (order,nx,ny);
    sf_eno2_set(pnt,sc);

    sf_putint(fid,"n2",nt);
    sf_putfloat(fid,"d2",ht);

    sf_putint(fx,"n2",nt);
    sf_putfloat(fx,"d2",ht);

    sf_putint(fy,"n2",nt);
    sf_putfloat(fy,"d2",ht);

    ht /= 2;

    ntx = nt*nx;
    w=sf_floatalloc(ntx);
    xr=sf_floatalloc(ntx);
    yr=sf_floatalloc(ntx);

    /* do the work */
    shooting();

    nx1=nx-1;
    for( j=0; j<nt; j++ ) {
	for( i=0; i<nx; i++ ) {
	    ind=i+nx*j;
	    if( i==0 ) {
		q=hypotf(xr[ind+1]-xr[ind],yr[ind+1]-yr[ind])/hx;
	    }
	    else if( i==nx1 ) {
		q=hypotf(xr[ind-1]-xr[ind],yr[ind-1]-yr[ind])/hx;
	    } else {
		q=hypotf(xr[ind+1]-xr[ind-1],yr[ind+1]-yr[ind-1])*0.5/hx;
	    }

	    if (w[ind] != 0.0) w[ind] /= q;
	}
    }

    sf_floatwrite(w,ntx,fid);
    sf_floatwrite(xr,ntx,fx);
    sf_floatwrite(yr,ntx,fy);


    exit(0);
}

