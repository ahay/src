/* Finite-difference 3-D heat-flow equation using helix */
/*
  Copyright (C) 2006 University of Texas at Austin
   
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

#include "helify.h"

#define NXY 10000

int main(int argc, char* argv[])
{
    const int nx = 100, ny = 100, a = 2, nf = 17, n1 = 8;
    const float gamma = 0.666666;
    float q[NXY], d[NXY];
    sf_filter aa;
    float alpha, scale, middle;
    int it, ix, iy, nt, nh;
    sf_file out;

    sf_init (argc,argv);
    out = sf_output("out");

    if (!sf_getint("n3",&nt)) nt=10;
    if (!sf_getint("nh",&nh)) nh=5;

    sf_setformat(out,"native_float");
    sf_putint(out,"n1",nx);
    sf_putint(out,"n2",ny);
    sf_putint(out,"n3",nt);

    alpha = 0.5*a; 
    middle = 4.*(gamma + (1.-gamma)*0.5);

    /* Initial temperature step */
    for (ix=0; ix < NXY; ix++) {
	q[ix] = 0.;
    }
    for (ix = ny/2-nh-1; ix < ny/2+nh; ix++) {
	for (iy = nx/2-nh-1+nx*ix; iy < nx/2+nh + nx*ix; iy++) {
	    q[iy] = 1.;
	}
    }

    aa = sf_allocatehelix(nf);
    scale = helify(1.,alpha,n1,nf,aa->flt); 

    for (it=0; it < n1; it++) {
	aa->lag[it] = it+1;
	sf_warning("aa[%d]=%g",aa->lag[it],aa->flt[it]);
    }
    for (it=n1; it < nf; it++) {
	aa->lag[it] = nx + 2 - nf + it;
	sf_warning("aa[%d]=%g",aa->lag[it],aa->flt[it]);
    }
    sf_polydiv_init (NXY, aa);

    for (it=0; it < nt; it++) { 
	sf_floatwrite (q,NXY,out);

	for (ix=0; ix < NXY; ix++) {
	    d[ix] = 0.;
	}
	for (ix = nx+1; ix < NXY-nx-1; ix++) {
	    d[ix] = (1. - alpha*middle) * q[ix] + 
		alpha*gamma * (q[ix-1] + q[ix+1] + q[ix-nx] + q[ix+nx]) + 
		0.5*alpha*(1. - gamma) * 
		(q[ix-nx-1] + q[ix-nx+1] + q[ix+nx-1] + q[ix+nx+1]); 
	}
	sf_polydiv_lop (true, false,  NXY, NXY, q, d);
	sf_polydiv_lop (false, false, NXY, NXY, q, d);
	for (ix=0; ix < NXY; ix++) {   
	    q[ix] = d[ix]*scale;
	}
    }

    exit(0);
}
