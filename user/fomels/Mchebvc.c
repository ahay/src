/* Post-stack 2-D velocity continuation by Chebyshev-tau method. */
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

#include "tricheb.h"

int main (int argc, char* argv[])
{
    int n, n1, n2, i1, i2, iv, nv;
    float dv, vel, v2, k, dk, dt;
    float *t, *d;
    sf_file in, out;

    sf_init(argc,argv);
    in = sf_input("in");
    out = sf_output("out");

    if (SF_FLOAT != sf_gettype(in)) sf_error("Need float input");
    if (!sf_histint(in,"n1",&n1)) sf_error("No n1= in input");
    if (!sf_histint(in,"n2",&n2)) sf_error("No n2= in input");
    n = n1-1;

    if (!sf_getint("nv",&nv)) sf_error("Need nv=");
    if (!sf_getfloat("vel",&vel)) sf_error("Need vel=");

    v2 = 0.5*vel*vel*SF_SIG(vel);
    dv = v2/(nv-1);

    if (!sf_histfloat(in,"dt",&dt)) sf_error("No dt= in input");
    if (!sf_histfloat(in,"d2",&dk)) sf_error("No dk= in input");
    dk *= 2*SF_PI;

    t = sf_floatalloc(n1);
    d = sf_floatalloc(n);
    tricheb_init(n);

    for (i2=0; i2 < n2; i2++) {
	k = i2*dk;
	k *= k*dt*dv*0.25;

	sf_floatread(t,n1,in);
	tricheb_step(-dv,k);

	for (iv=0; iv < nv; iv++) {
	    tricheb_apply (d,t);
	    tricheb_solve (d,t+1);
	    t[0] = 0.;
	    for (i1=1; i1 < n; i1 +=2) {
		t[0] += t[i1]*SF_SIG(dv) - t[i1+1];
	    }
	}

	sf_floatwrite(t,n1,out);
    }

    exit(0);
}
