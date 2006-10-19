/* Finite-difference solution of 2-D heat-flow equation */
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

int main(int argc, char* argv[])
{
    const int nx = 100, nt = 10;
    const float beta = 0.122996;
    float q[nx], d[nx], alpha, b1, b2;
    int it, ix;
    bool impl;
    sf_tris slv;
    sf_file out;

    sf_init(argc,argv);
    out = sf_output("out");
    sf_setformat(out,"native_float");
    sf_putint(out,"n1",nx);
    sf_putint(out,"n2",nt);

    if (!sf_getbool("impl",&impl)) impl=false;
    /* if y, use implicit scheme */
    if (!sf_getfloat("alpha",&alpha)) alpha=1.;
    alpha *= 0.5;

    /* initial temperature */
    for (ix=0; ix < nx/2-1; ix++) {
	q[ix]=0.;
    }
    q[nx/2-1] = 0.5;
    for (ix=nx/2; ix < nx; ix++) {
	q[ix] = 1.;
    }

    if (impl) {
	b1 = alpha + 0.5*beta;
	b2 = alpha - 0.5*beta;
	alpha = b1;

	slv = sf_tridiagonal_init (nx);
	sf_tridiagonal_const_define (slv,1.+2.*b2,-b2,true);
    } else {
	slv = NULL;
    }

    for (it=0; it < nt; it++) { 
	sf_floatwrite(q,nx,out);
   
	d[0] = q[0] + alpha*(q[1]-q[0]);
	for (ix=1; ix < nx-1; ix++) {
	    d[ix] = q[ix] + alpha*(q[ix-1]-2.*q[ix]+q[ix+1]);
	}
	d[nx-1] = q[nx-1] + alpha*(q[nx-2]-q[nx-1]);
	for (ix=0; ix < nx; ix++) {
	    q[ix] = d[ix];
	} 

	if (impl) sf_tridiagonal_solve (slv,q);
    }
    
    exit(0);
}
