/* Estimate complex PEF on the first axis. */
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

#include "ctcaf1.h"

int main(int argc, char* argv[])
{
    bool *mask;
    int i1, n1, i2, n2, nf, nr, niter;
    sf_complex *trace, *a, *resid;
    sf_file in, out;

    sf_init(argc,argv);
    in = sf_input("in");
    out = sf_output("out");

    if (SF_COMPLEX != sf_gettype(in)) sf_error("Need complex input");
    if (!sf_histint(in,"n1",&n1)) sf_error("No n1= in input");
    n2 = sf_leftsize(in,1);
    
    if (!sf_getint("nf",&nf)) sf_error("Need nf=");
    /* filter length */

    if (!sf_getint("niter",&niter)) niter=2*nf;
    /* number of iterations */

    sf_putint(out,"n1",nf);
    nr = n1+nf-1;

    a = sf_complexalloc(nf);
    mask = sf_boolalloc(nf);

    a[0] = sf_cmplx(1.0f,0.0f);
    mask[0] = true;

    for (i1=1; i1 < nf; i1++) {
	a[i1] = sf_cmplx(0.0f,0.0f);
	mask[i1] = false;
    }

    trace = sf_complexalloc(n1);
    resid = sf_complexalloc(nr);

    for (i1=0; i1 < nr; i1++) {
	resid[i1] = sf_cmplx(0.0f,0.0f);
    }

    ctcaf1_init(n1,trace);

    for (i2=0; i2 < n2; i2++) {
	sf_complexread(trace,n1,in);

	sf_csolver(ctcaf1_lop, sf_ccgstep, nf, nr, a, 
		   resid, niter, "x0", a, "known", mask, "end");
	sf_ccgstep_close();

	sf_complexwrite(a,nf,out);
    }

    exit(0);
}
