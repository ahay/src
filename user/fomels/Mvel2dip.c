/* Velocity-dependent dip */
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

int main(int argc, char* argv[])
{
    int it, nt, ix, nx;
    float ot, dt, ox, dx, t0, ti, x, *t, *p, v, *vt;
    sf_file vel, tp;

    sf_init(argc,argv);
    vel = sf_input("in");
    tp = sf_output("out");

    if (!sf_histint(vel,"n1",&nt)) sf_error("Need n1= in input");
    if (!sf_histfloat(vel,"d1",&dt)) sf_error("Need d1= in input");
    if (!sf_histfloat(vel,"o1",&ot)) sf_error("Need o1= in input");

    if (!sf_getint("nx",&nx)) nx=1;
    if (!sf_getfloat("dx",&dx)) dx=1.0f;
    if (!sf_getfloat("x0",&ox)) ox=0.0f;
    /* offset sampling */

    sf_putint(tp,"n1",nt+1);
    sf_putint(tp,"n2",2);
    sf_putint(tp,"n3",nx);
    sf_putfloat(tp,"d3",dx);
    sf_putfloat(tp,"o3",ox);
    
    vt = sf_floatalloc(nt);
    sf_floatread(vt,nt,vel);

    t = sf_floatalloc(nt+1);
    p = sf_floatalloc(nt+1);

    for (ix=0; ix < nx; ix++) {
	x = ox+ix*dx;
	for (it=0; it < nt; it++) {
	    t0 = ot+it*dt;
	    v = vt[it];
	    
	    ti = sqrtf(t0*t0+x*x/(v*v));
	    t[it] = ti;
	    if (ti < dt) ti=dt; 
	    p[it] = x/(ti*v*v);
	}
	t[nt] = ot-dt;
	p[nt] = 0.0f;
	sf_floatwrite(t,nt+1,tp);
	sf_floatwrite(p,nt+1,tp);
    }

    exit(0);
}
