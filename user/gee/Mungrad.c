/* Phase unwrapping by least squares. */
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

#include  "unwrap.h"

int main(int argc, char* argv[])
{
    int n1, n2, i1, i2, niter;
    float **b, *h, ***rt, d, maxd;
    sf_complex **data;
    sf_file inp, out, bad;

    sf_init (argc,argv);
    inp = sf_input("in");
    out = sf_output("out");
    
    if (SF_COMPLEX != sf_gettype(inp)) sf_error("Need complex input");
    if (!sf_histint(inp,"n1",&n1)) sf_error("No n1= in input");
    if (!sf_histint(inp,"n2",&n2)) sf_error("No n2= in input");

    sf_settype(out,SF_FLOAT);

    if (!sf_getint("niter",&niter)) niter=0;
    /* number of iterations */

    if (NULL != sf_getstring("badness")) {
	bad = sf_output("badness");
	/* (optional) badness attribute file */
	sf_settype(bad,SF_FLOAT);
    } else {
	bad = NULL;
    }

    data = sf_complexalloc2(n1,n2);
    sf_complexread(data[0],n1*n2,inp);

    /* normalize */
    maxd = 0;
    for (i2=0; i2 < n2; i2++) {
	for (i1=0; i1 < n1; i1++) {
	    d = cabsf(data[i2][i1]);
	    if (maxd < d) maxd=d;
	}
    }
    for (i2=0; i2 < n2; i2++) {
	for (i1=0; i1 < n1; i1++) {
#ifdef SF_HAS_COMPLEX_H
	    data[i2][i1] = data[i2][i1]/maxd;
#else
	    data[i2][i1] = sf_crmul(data[i2][i1],1.0f/maxd);
#endif
	}
    }

    if (NULL != bad) {
	b = sf_floatalloc2(n1,n2);
	rt = sf_floatalloc3(n1,n2,2);

	grad2init (n1, n2, data, rt);

	for (i2=0; i2 < n2; i2++) {
	    b[i2][n1-1] = 0.;
	    b[i2][n1-2] = 0.;
	}

	for (i1=0; i1 < n1; i1++) {
	    b[n2-1][i1] = 0.;
	    b[n2-2][i1] = 0.;
	}
	
	for (i2=0; i2 < n2-2; i2++) {
	    for (i1=0; i1 < n1-2; i1++) {
		b[i2][i1] = 
		    (rt[0][i2+1][i1] - rt[0][i2][i1]) - 
		    (rt[1][i2][i1+1] - rt[1][i2][i1]);
	    }
	}

	sf_floatwrite(b[0],n1*n2,bad);
    }

    if (niter > 0) {
	h = sf_floatalloc(n1*n2);
	unwraper (n1,n2, data, h, niter);
	sf_floatwrite(h,n1*n2,out);
    }

    exit(0);
}




