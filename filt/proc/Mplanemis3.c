/* Missing data interpolation in 3-D using plane-wave destruction. */
/*
  Copyright (C) 2004 University of Texas at Austin
  
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
#include "allp3.h"

int main(int argc, char* argv[])
{
    int i, niter, nw, n1, n2, n3, n123;
    float *mm, *dd, ***pp, ***qq;
    bool *known, verb;
    sf_file in, out, dip, mask;

    sf_init (argc,argv);
    in = sf_input("in");
    out = sf_output("out");
    dip = sf_input("dip");

    if (!sf_histint(in,"n1",&n1)) sf_error("No n1= in input");
    if (!sf_histint(in,"n2",&n2)) sf_error("No n2= in input");
    if (!sf_histint(in,"n3",&n3)) sf_error("No n3= in input");
    n123 = n1*n2*n3;

    if (!sf_getint("niter",&niter)) niter=100;
    /* number of iterations */

    if (!sf_getint("order",&nw)) nw=1;
    /* [1,2,3] accuracy order */
    if (nw < 1 || nw > 3) 
	sf_error ("Unsupported nw=%d, choose between 1 and 3",nw);

    if (!sf_getbool("verb",&verb)) verb = false;
    /* verbosity flag */

    pp = sf_floatalloc3(n1,n2,n3);
    qq = sf_floatalloc3(n1,n2,n3);

    sf_floatread(pp[0][0],n123,dip);
    sf_floatread(qq[0][0],n123,dip);

    mm = sf_floatalloc(n123);
    dd = sf_floatalloc(2*n123);
    known = sf_boolalloc(n123);

    sf_floatread(mm,n123,in);
    
    if (NULL != sf_getstring ("mask")) {
	mask = sf_input("mask");
	sf_floatread(dd,n123,mask);

	for (i=0; i < n123; i++) {
	    known[i] = (dd[i] != 0.);
	}
    } else {
	for (i=0; i < n123; i++) {
	    known[i] = (mm[i] != 0.);
	}
    }

    for (i=0; i < 2*n123; i++) {
	dd[i] = 0.;
    }

    allpass3_init(allpass_init(nw, 1, n1,n2,n3, pp),
		  allpass_init(nw, 1, n1,n2,n3, qq));
    sf_solver(allpass3_lop, sf_cgstep, n123, 2*n123, mm, dd, niter,
	      "known", known, "x0", mm, "verb", verb, "end");

    sf_floatwrite (mm,n123,out);

    exit(0);
}
