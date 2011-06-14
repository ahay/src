/* Zero crossings with sub-pixel resolution. */
/*
  Copyright (C) 2011 University of Texas at Austin
  
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

static sf_eno eno;
static int it;

static float func_eno(float t);

int main (int argc, char* argv[])
{
    bool down;
    int nt, i2, n2, *n0, nz, nw;
    float *trace, dt, t0, a, b, t;
    sf_file inp, nzero, out;

    sf_init (argc,argv);
    inp = sf_input("in");
    out = sf_output("out");

    nzero = sf_output("nzero");

    if (!sf_histint(inp,"n1",&nt)) sf_error("No n1= in input");
    if (!sf_histfloat(inp,"d1",&dt)) dt=1.;
    if (!sf_histfloat(inp,"o1",&t0)) t0=0.;

    n2 = sf_leftsize(inp,1);

    trace = sf_floatalloc(nt);

    n0 = sf_intalloc(n2);
    sf_putint(nzero,"n1",1);
    sf_settype(nzero,SF_INT);

    if (!sf_getint ("nw",&nw)) nw=4;
    /* Interpolation accuracy */

    if (!sf_getbool("down",&down)) down=false;
    /* only zeros on the way down */

    eno = sf_eno_init (nw, nt);

    for (i2=0; i2 < n2; i2++) {
	sf_floatread(trace,nt,inp);

	nz = 0;
	for (it = 0; it < nt-1; it++) {
	    a = trace[it];
	    b = trace[it+1];
	
	    if ((a <= 0. && b > 0. && !down) ||
		(a >= 0. && b < 0.)) {	  
		t = sf_zero(func_eno,0.,1.,a,b,1.e-3,false);
	  
		trace[nz] = t0+(it+t)*dt;
		nz++;
	    }
	}

	sf_floatwrite(trace,nt,out);

	n0[i2] = nz;
    }

    sf_intwrite(n0,n2,nzero);

    exit(0);
}

static float func_eno(float t)
/* interpolation function */
{
    float f, g;
    sf_eno_apply (eno,it,t,&f,&g,FUNC);
    return f;
}
