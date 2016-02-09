/* 3-D Inverse normal moveout.

velocity file contains slowness squared with n2=3 (wx,wy,wxy)
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
    int it,ix,iy, nt, nx, ny, nw, i4, n4, n;
    float dt, t0, x, y, x0, y0, f, dx, dy, eps;
    float *trace=NULL, *vx=NULL, *vy=NULL, *vxy=NULL, *str=NULL, *out=NULL;
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
    n4 = sf_leftsize(cmp,3);

    vel = sf_input("velocity");
    if (SF_FLOAT != sf_gettype(vel)) sf_error("Need float velocity");
    if (!sf_histint(vel,"n1",&n) || nt != n) sf_error("Need n1=%d in velocity",nt);
    if (!sf_histint(vel,"n2",&n) || 3 != n) sf_error("Need n2=3 in velocity",nt);
    if (n4 != sf_leftsize(vel,2)) sf_error("Wrong dimensions in velocity");

    if (!sf_histfloat(cmp,"d2",&dx)) sf_error("No d2= in input");
    if (!sf_histfloat(cmp,"o2",&x0)) sf_error("No o2= in input");
    if (!sf_histfloat(cmp,"d3",&dy)) sf_error("No d3= in input");
    if (!sf_histfloat(cmp,"o3",&y0)) sf_error("No o3= in input");

    if (!sf_getbool("half",&half)) half=true;
    /* if y, the second and third axes are half-offset instead of full offset */

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
    vx = sf_floatalloc(nt);
    vy = sf_floatalloc(nt);
    vxy = sf_floatalloc(nt);

    str = sf_floatalloc(nt);
    out = sf_floatalloc(nt);

    if (!sf_getint("extend",&nw)) nw=8;
    /* trace extension */

    nmo = sf_stretch4_init (nt, t0, dt, nt, eps);

    for (i4=0; i4 < n4; i4++) { /* loop over cmps */
      sf_floatread (vx,nt,vel);
      sf_floatread (vy,nt,vel);
      sf_floatread (vxy,nt,vel);

      for (iy = 0; iy < ny; iy++) {
	y = y0+iy*dy;

	for (ix = 0; ix < nx; ix++) {
	  x = x0 + ix*dx;

	  sf_floatread (trace,nt,cmp);

	  for (it=0; it < nt; it++) {
	    f = t0 + it*dt;
	    f = f*f + x*x*vx[it]+y*y*vy[it]+2*x*y*vxy[it];

	    if (f < 0.) {
	      str[it]=t0-10.*dt;
	    } else {
	      str[it] = sqrtf(f);
	    }
	  }

	  sf_stretch4_define (nmo,str);
	  sf_stretch4_apply (false,nmo,trace,out);

	  sf_floatwrite (out,nt,nmod);
	} /* ix */
      } /* iy */
    } /* i4 */

    exit (0);
}
