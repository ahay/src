/* Picking local maxima on the first axis */
/*
  Copyright (C) 2007 University of Texas at Austin
  
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
    int i1, n1, i2, n2, ip, np, *npick;
    float o1, d1, t0, t1, t2, t, a, *trace, *pick, *ampl;
    sf_file in, picks, npicks, ampls;

    sf_init(argc, argv);
    in = sf_input("in");
    picks = sf_output("out");
    ampls = sf_output("ampl");
    npicks = sf_output("picks");
    
    /* number of picks at each trace */

    if (!sf_histint(in,"n1",&n1)) sf_error("No n1= in input");
    n2 = sf_leftsize(in,1);

    if (!sf_histfloat(in,"d1",&d1)) d1=1.;
    if (!sf_histfloat(in,"o1",&o1)) o1=0.;

    if (!sf_getint("np",&np)) np=n1;
    /* maximum number of picks */

    sf_putint(picks,"n1",np);
    sf_putint(ampls,"n1",np);

    sf_putint(npicks,"n1",1);
    sf_settype(npicks,SF_INT);

    trace = sf_floatalloc(n1);
    pick = sf_floatalloc(np);
    ampl = sf_floatalloc(np);
    npick = sf_intalloc(n2);

    for (i2=0; i2 < n2; i2++) {
	sf_floatread(trace,n1,in);

	t0 = trace[0];
	t1 = trace[1];
	ip = 0;
	for (i1=2; i1 < n1; i1++) {
	    t2 = trace[i1];

	    if (ip < np && t1 > t0 && t1 > t2) {
		/* parabolic approximation */
		t = 0.5*(t2-t0)/(2*t1-t0-t2);
		a = t1+0.25*(t2-t0)*t;

		if (t < -1.) {
		    t=-1;
		    a=t0;
		} else if (t > 1.) {
		    t=1.;
		    a=t2;
		} 

		pick[ip] = o1+(i1-1+t)*d1;
		ampl[ip] = a;
		ip++;
	    }

	    t0 = t1;
	    t1 = t2;
	}
	npick[i2] = ip;

	if (0==ip) {
	    pick[0] = o1-d1;
	    ip++;
	}
	
	for (i1=ip; i1 < np; i1++) {
	    pick[i1] = pick[ip-1];
	    ampl[i1] = 0.;
	}

	sf_floatwrite(pick,np,picks);
	sf_floatwrite(ampl,np,ampls);
    }
    sf_intwrite(npick,n2,npicks);

    exit(0);
}
