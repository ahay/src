/* Rotate phase and compute similarity with enevelope. */
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
    bool inv, verb;
    int i1, n1, ia, na, i2, n2, rect, niter, n;
    float *trace, *sim1, *sim2, *envel, *hilbt, *rotat, a0, da, a, c, norm;
    sf_file inp, sim;

    sf_init(argc,argv);
    inp = sf_input("in");
    sim = sf_output("out");

    if (!sf_histint(inp,"n1",&n1)) sf_error("No n1= in input");
    n2 = sf_leftsize(inp,1);
    
    if (!sf_getint("na",&na)) na=360; /* number of angles */
    if (!sf_getfloat("da",&da)) da=1.0; /* angle increment */
    if (!sf_getfloat("a0",&a0)) a0=-180.; /* first angle */

    sf_shiftdim(inp, sim, 2);
    
    sf_putint(sim,"n2",na);
    sf_putfloat(sim,"d2",da);
    sf_putfloat(sim,"o2",a0);
    sf_putstring(sim,"label2","Angle");
    sf_putstring(sim,"unit2","\\^o\\_");

    /* convert to radians */
    da *= SF_PI/180.;
    a0 *= SF_PI/180.;

    trace = sf_floatalloc(n1);
    sim1 = sf_floatalloc(n1);
    sim2 = sf_floatalloc(n1);
    hilbt = sf_floatalloc(n1);
    envel = sf_floatalloc(n1);
    rotat = sf_floatalloc(n1);

    if (!sf_getint("order",&n)) n=100;
    /* Hilbert transformer order */
    if (!sf_getfloat("ref",&c)) c=1.;
    /* Hilbert transformer reference (0.5 < ref <= 1) */

    if (!sf_getbool("verb",&verb)) verb=true;
    /* verbosity */

    if (!sf_getbool("inv",&inv)) inv=true;
    /* inverse similarity */

    if (!sf_getint("niter",&niter)) niter=20;
    /* maximum number of iterations */
    
    if (!sf_getint("rect",&rect)) rect=3;
    /* smoothing radius */

    sf_hilbert_init(n1, n, c);
    sf_divn_init(1, n1, &n1, &rect, niter, verb);

    for (i2=0; i2 < n2; i2++) {
	if (!verb) sf_warning("trace %d of %d;",i2+1,n2);

	sf_floatread(trace,n1,inp);

	sf_hilbert(trace,hilbt);

	/* compute envelope and normalize */

	norm = 0.;
	for (i1=0; i1 < n1; i1++) {
	    envel[i1] = hypotf(trace[i1],hilbt[i1]);
	    norm += envel[i1]*envel[i1];
	}
	norm = sqrtf(n1/norm);	
	for (i1=0; i1 < n1; i1++) {
	    envel[i1] *= norm;
	}

	/* loop over angles */
	for (ia=0; ia < na; ia++) {
	    a = a0 + ia*da;

	    /* rotate phase */
	    norm = 0.;
	    for (i1=0; i1 < n1; i1++) {
		rotat[i1] = trace[i1]*cosf(a) + hilbt[i1]*sinf(a);
		norm += rotat[i1]*rotat[i1];
	    }
	    norm = sqrtf(n1/norm);	
	    for (i1=0; i1 < n1; i1++) {
		rotat[i1] *= norm;
	    }

	    /* compute similarity */
	    sf_divn(rotat,envel,sim1);
	    sf_divn(envel,rotat,sim2);

	    sf_divn_combine(sim1,sim2,sim1);

	    if (inv) {
		for (i1=0; i1 < n1; i1++) {
		    sim1[i1] = 1.0/sim1[i1];
		}
	    }
	    
	    sf_floatwrite(sim1,n1,sim);
	}
    }
    if (!verb) sf_warning(".");

    exit(0);
}
