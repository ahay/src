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

#define NX 100
#define NT 10

int main(int argc, char* argv[])
{
    const float beta = 0.122996;
    float q[NX], d[NX], alpha, b1, b2;
    int it, ix;
    bool impl;
    sf_tris slv;
    sf_file out;

    sf_init(argc,argv);
    out = sf_output("out");
    sf_setformat(out,"native_float");
    sf_putint(out,"n1",NX);
    sf_putint(out,"n2",NT);

    if (!sf_getbool("impl",&impl)) impl=false;
    /* if y, use implicit scheme */
    if (!sf_getfloat("alpha",&alpha)) alpha=1.;
    alpha *= 0.5;

    /* initial temperature */
    for (ix=0; ix < NX/2-1; ix++) {
	q[ix]=0.;
    }
    q[NX/2-1] = 0.5;
    for (ix=NX/2; ix < NX; ix++) {
	q[ix] = 1.;
    }

    if (impl) {
	b1 = alpha + 0.5*beta;
	b2 = alpha - 0.5*beta;
	alpha = b1;

	slv = sf_tridiagonal_init (NX);
	sf_tridiagonal_const_define (slv,1.+2.*b2,-b2,true);
    } else {
	slv = NULL;
    }

    for (it=0; it < NT; it++) { 
	sf_floatwrite(q,NX,out);
   
	d[0] = q[0] + alpha*(q[1]-q[0]);
	for (ix=1; ix < NX-1; ix++) {
	    d[ix] = q[ix] + alpha*(q[ix-1]-2.*q[ix]+q[ix+1]);
	}
	d[NX-1] = q[NX-1] + alpha*(q[NX-2]-q[NX-1]);
	for (ix=0; ix < NX; ix++) {
	    q[ix] = d[ix];
	} 

	if (impl) sf_tridiagonal_solve (slv,q);
    }
    
    exit(0);
}
