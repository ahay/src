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
#include <time.h>

#include <rsf.h>
#include "allp3.h"

int main(int argc, char* argv[])
{
    int i, niter, nw, n1, n2, n3, n123, nj1, nj2, seed, i4, n4;
    float *mm, *dd, *pp, *qq, a, var;
    bool *known, verb, drift;
    sf_file in, out, dip, mask;

    sf_init (argc,argv);
    in = sf_input("in");
    out = sf_output("out");
    dip = sf_input("dip");

    if (!sf_histint(in,"n1",&n1)) sf_error("No n1= in input");
    if (!sf_histint(in,"n2",&n2)) sf_error("No n2= in input");
    if (!sf_histint(in,"n3",&n3)) sf_error("No n3= in input");
    n123 = n1*n2*n3;
    n4 = sf_leftsize(in,3);

    if (!sf_getint("niter",&niter)) niter=100;
    /* number of iterations */

    if (!sf_getint("order",&nw)) nw=1;
    /* [1,2,3] accuracy order */
    if (nw < 1 || nw > 3) 
	sf_error ("Unsupported nw=%d, choose between 1 and 3",nw);

    if (!sf_getint("nj1",&nj1)) nj1=1;
    if (!sf_getint("nj2",&nj2)) nj2=1;
    /* antialiasing */

    if (!sf_getbool("drift",&drift)) drift=false;
    /* if shift filter */

    if (!sf_getbool("verb",&verb)) verb = false;
    /* verbosity flag */
    
    if (!sf_getint("seed",&seed)) seed = time(NULL);
    /* random seed */
    init_genrand((unsigned long) seed);

    if (!sf_getfloat("var",&var)) var=0.;
    /* noise variance */
    a = sqrtf(var);

    pp = sf_floatalloc(n123);
    qq = sf_floatalloc(n123);

    mm = sf_floatalloc(n123);
    dd = sf_floatalloc(2*n123);
    known = sf_boolalloc(n123);

    if (NULL != sf_getstring ("mask")) {
	mask = sf_input("mask");
    } else {
	mask = NULL;
    }

    allpass3_init(allpass_init(nw, nj1, n1,n2,n3, drift, pp),
		  allpass_init(nw, nj2, n1,n2,n3, drift, qq));

    for (i4=0; i4 < n4; i4++) {
	sf_warning("slice %d of %d",i4+1,n4);
	sf_floatread(pp,n123,dip);
	sf_floatread(qq,n123,dip);
	
	sf_floatread(mm,n123,in);
	
	if (NULL != mask) {
	    sf_floatread(dd,n123,mask);	    
	    for (i=0; i < n123; i++) {
		known[i] = (bool) (dd[i] != 0.);
	    }
	} else {
	    for (i=0; i < n123; i++) {
		known[i] = (bool) (mm[i] != 0.);
	    }
	}
	
	for (i=0; i < 2*n123; i++) {
	    dd[i] = a*sf_randn_one_bm();
	}
	
	sf_solver(allpass3_lop, sf_cgstep, n123, 2*n123, mm, dd, niter,
		  "known", known, "x0", mm, "verb", verb, "end");
	sf_cgstep_close();

	sf_floatwrite (mm,n123,out);
    }

    exit(0);
}
