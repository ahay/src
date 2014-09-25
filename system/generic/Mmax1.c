/* Picking local maxima on the first axis. 

Outputs complex numbers (time,amplitude) sorted by amplitude.

September 2014 program of the month:
http://ahay.org/rsflog/index.php?/archives/403-Program-of-the-month-sfmax1.html
*/
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

static int pick_compare (const void *p1, const void *p2)
{
    float f1 = cimagf(* (sf_complex*) p1);
    float f2 = cimagf(* (sf_complex*) p2);
    return (f1 < f2)? 1: (f1 > f2)? -1: 0;
}

int main(int argc, char* argv[])
{
    int i1, n1, i2, n2, ip, np;
    float o1, d1, t0, t1, t2, t, a, *trace=NULL;
    float min, max, x;
    sf_complex *pick=NULL;
    sf_file in=NULL, out=NULL;

    sf_init(argc, argv);
    in = sf_input("in");
    out = sf_output("out");

    if (!sf_histint(in,"n1",&n1)) sf_error("No n1= in input");
    n2 = sf_leftsize(in,1);

    if (!sf_histfloat(in,"d1",&d1)) d1=1.;
    if (!sf_histfloat(in,"o1",&o1)) o1=0.;

    if (!sf_getfloat("min",&min)) min=o1;
    /* minimum value of time */

    if (!sf_getfloat("max",&max)) max=o1+(n1-1)*d1; 
    /* maximum value of time */ 

    if (!sf_getint("np",&np)) np=n1;
    /* maximum number of picks */

    sf_putint(out,"n1",np);
    sf_settype(out,SF_COMPLEX);

    trace = sf_floatalloc(n1);
    pick = sf_complexalloc(np);

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

		x = o1+(i1-1+t)*d1;

		if (x >= min && x <= max) {
			pick[ip] = sf_cmplx(x,a);	
			ip++;
		}
	    }

	    t0 = t1;
	    t1 = t2;
	}

	if (0==ip) {
	    pick[0] = sf_cmplx(o1-d1,0.);
	    ip++;
	}

	qsort(pick,ip,sizeof(sf_complex),pick_compare);
	
	for (i1=ip; i1 < np; i1++) {
	    pick[i1] = sf_cmplx(crealf(pick[ip-1]),0.);
	}

	sf_complexwrite(pick,np,out);
    }

    exit(0);
}
