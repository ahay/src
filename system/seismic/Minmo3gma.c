/* 3-D Inverse generalized normal moveout.
   Velocity file contains slowness squared with n2=16 (wx,wxy,wy,A1,A2,A3,A4,A5,B1,B2,B3,C1,C2,C3,C4,C5)
   following Sripanich and Fomel (2015).
*/
/*
  Copyright (C) 2004 University of Texas at Austin

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

#include <math.h>
#include <rsf.h>

int main (int argc, char* argv[])
{
    sf_map4 nmo; /* using cubic spline interpolation */
    bool half;
    int it,ix,iy, nt, nx, ny, nw, n;
    float dt, t0, x, y, x0, y0, f, dx, dy, eps, A, B, C;
    float *trace=NULL, *wx=NULL, *wy=NULL, *wxy=NULL, *str=NULL, *out=NULL;
    float *a1,*a2,*a3,*a4,*a5,*b1,*b2,*b3,*c1,*c2,*c3,*c4,*c5;

    sf_file cmp=NULL, nmod=NULL, vel=NULL;

    sf_init (argc,argv);
    cmp = sf_input("in");
    nmod = sf_output("out");

    if (SF_FLOAT != sf_gettype(cmp)) sf_error("Need float input");
    if (!sf_histint(cmp,"n1",&nt)) sf_error("No n1= in input");
    if (!sf_histfloat(cmp,"d1",&dt)) sf_error("No d1= in input");
    if (!sf_histfloat(cmp,"o1",&t0)) sf_error("No o1= in input");

    if (!sf_histint(cmp,"n2",&nx)) sf_error("No n2= in input");
    if (!sf_histint(cmp,"n3",&ny)) sf_error("No n3= in input");

    vel = sf_input("velocity");
    if (SF_FLOAT != sf_gettype(vel)) sf_error("Need float velocity");
    if (!sf_histint(vel,"n1",&n) || nt != n) sf_error("Need n1=%d in velocity",nt);
    if (!sf_histint(vel,"n2",&n) || 16 != n) sf_error("Need n2=16 in velocity (for 16 coefficients)",nt);

    if (!sf_histfloat(cmp,"d2",&dx)) sf_error("No d2= in input");
    if (!sf_histfloat(cmp,"o2",&x0)) sf_error("No o2= in input");
    if (!sf_histfloat(cmp,"d3",&dy)) sf_error("No d3= in input");
    if (!sf_histfloat(cmp,"o3",&y0)) sf_error("No o3= in input");

    if (!sf_getbool("half",&half)) half=false;
    /* if x,y , the second and third axes are half-offset instead of full offset */

    if (half) { /* get full offset - multiply by 2 */
	dx *= 2.;
	dy *= 2.;
	x0 *= 2.;
	y0 *= 2.;
	sf_warning("Since half=y, offsets were doubled.");
    }else{
	sf_warning("Since half=n, offsets not doubled.");
    }

    if (!sf_getfloat("eps",&eps)) eps=0.01;
    /* stretch regularization */

    trace = sf_floatalloc(nt);
    wx = sf_floatalloc(nt);
    wxy = sf_floatalloc(nt);
	wy = sf_floatalloc(nt);
    a1 = sf_floatalloc(nt);
    a2 = sf_floatalloc(nt);
    a2 = sf_floatalloc(nt);
    a3 = sf_floatalloc(nt);
    a4 = sf_floatalloc(nt);
    a5 = sf_floatalloc(nt);
    b1 = sf_floatalloc(nt);
    b2 = sf_floatalloc(nt);
    b3 = sf_floatalloc(nt);
    c1 = sf_floatalloc(nt);
    c2 = sf_floatalloc(nt);
    c3 = sf_floatalloc(nt);
    c4 = sf_floatalloc(nt);
    c5 = sf_floatalloc(nt);

    str = sf_floatalloc(nt);
    out = sf_floatalloc(nt);

    if (!sf_getint("extend",&nw)) nw=8;
    /* trace extension */

    nmo = sf_stretch4_init (nt, t0, dt, nt, eps);

    sf_floatread (wx,nt,vel);
    sf_floatread (wxy,nt,vel);
    sf_floatread (wy,nt,vel);
    sf_floatread (a1,nt,vel);
    sf_floatread (a2,nt,vel);
    sf_floatread (a3,nt,vel);
    sf_floatread (a4,nt,vel);
    sf_floatread (a5,nt,vel);
    sf_floatread (b1,nt,vel);
    sf_floatread (b2,nt,vel);
    sf_floatread (b3,nt,vel);
    sf_floatread (c1,nt,vel);
    sf_floatread (c2,nt,vel);
    sf_floatread (c3,nt,vel);
    sf_floatread (c4,nt,vel);
    sf_floatread (c5,nt,vel);

    for (iy = 0; iy < ny; iy++) {
	y = y0+iy*dy;

	for (ix = 0; ix < nx; ix++) {
	    x = x0 + ix*dx;

	    sf_floatread (trace,nt,cmp);

	    for (it=0; it < nt; it++) {
		f = t0 + it*dt;
		A = a1[it]*pow(x,4) + a2[it]*pow(x,3)*y +a3[it]*x*x*y*y +a4[it]*pow(y,3)*x+a5[it]*pow(y,4);
		B = b1[it]*x*x + b2[it]*x*y + b3[it]*y*y;
		C = c1[it]*pow(x,4) + c2[it]*pow(x,3)*y +c3[it]*x*x*y*y +c4[it]*pow(y,3)*x+c5[it]*pow(y,4);
		f = f*f + x*x*wx[it]+y*y*wy[it]+ x*y*wxy[it] + A/(f*f + B + sqrt(pow(f,4) + 2*f*f*B + C) );

		if (f < 0.) {
		    str[it]=t0-10.*dt;
		} else {
		    str[it] = sqrtf(f);
		}
	    }

	    sf_stretch4_define (nmo,str,false);
	    sf_stretch4_apply (false,nmo,trace,out);

	    sf_floatwrite (out,nt,nmod);
	} /* ix */
    } /* iy */

    exit (0);
}
