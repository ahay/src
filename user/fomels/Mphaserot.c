/* Non-stationary phase rotation. */
/*
  Copyright (C) 2010 University of Texas at Austin
  
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
    int i1, n1, i2, n2, ia, na, n;
    float *trace, *hilbt, *phase, **rotat, a, a0, da, c, deg2rad;
    sf_file inp, pha, out;

    sf_init(argc,argv);
    inp = sf_input("in");
    pha = sf_input("phase");
    out = sf_output("out");

    if (!sf_histint(inp,"n1",&n1)) sf_error("No n1= in input");
    n2 = sf_leftsize(inp,1);

    if (!sf_getint("na",&na)) na=721; /* number of angles */
    if (!sf_getfloat("da",&da)) da=1.0; /* angle increment */
    if (!sf_getfloat("a0",&a0)) a0=-360.; /* first angle */

    /* convert to radians */
    deg2rad = SF_PI/180.;

    trace = sf_floatalloc(n1);
    hilbt = sf_floatalloc(n1);
    phase = sf_floatalloc(n1);
    rotat = sf_floatalloc2(n1,na);

    if (!sf_getint("order",&n)) n=100;
    /* Hilbert transformer order */
    if (!sf_getfloat("ref",&c)) c=1.;
    /* Hilbert transformer reference (0.5 < ref <= 1) */

    sf_hilbert_init(n1, n, c);

    for (i2=0; i2 < n2; i2++) {
	sf_floatread(trace,n1,inp);

	sf_hilbert(trace,hilbt);

	/* loop over angles */
	for (ia=0; ia < na; ia++) {
	    a = deg2rad*(a0 + ia*da);

	    for (i1=0; i1 < n1; i1++) {
		rotat[ia][i1] = trace[i1]*cosf(a) + hilbt[i1]*sinf(a);
	    }
	}

	sf_floatread(phase,n1,pha);

	/* linear interpolation */
	for (i1=0; i1 < n1; i1++) {
	    a = (phase[i1]-a0)/da;
	    ia = floorf(a); a -= ia;
	    if (ia >=0 && ia < na-1) {
		trace[i1] = a * rotat[ia+1][i1] + (1.-a) * rotat[ia][i1];
	    } else {
		trace[i1] = 0.;
	    }
	}

	sf_floatwrite(trace,n1,out);
    }

    exit(0);
}
