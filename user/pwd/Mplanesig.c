/* Signal separation using plane-wave destruction filters. */
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
    
#include "allp2.h"
#include "planesig.h"

int main(int argc, char* argv[])
{
    int n1, n2, n12, nj, niter, nw, ic, nc, m1, m2, ns;
    float *eps, *d, *s, *dd, ***ss, norm;
    bool verb, drift;
    sf_file in, out, dips;

    sf_init (argc,argv);
    in = sf_input("in");
    dips = sf_input("dips");
    out = sf_output("out");

    if (!sf_histint(in,"n1",&n1)) sf_error("No n1= in input");
    if (!sf_histint(in,"n2",&n2)) sf_error("No n2= in input");
    n12 = n1*n2;

    if (!sf_histint(dips,"n1",&m1)) sf_error("No n1= in dips");
    if (!sf_histint(dips,"n2",&m2)) sf_error("No n2= in dips");
    if (!sf_histint(dips,"n3",&nc)) nc=1;
    if (m1 != n1 || m2 != n2) sf_error("Wrong dimensions in dips");
    
    sf_putint (out,"n3",nc);

    if (!sf_getint ("niter",&niter)) niter=50;
    /* maximum number of iterations */

    eps = sf_floatalloc(nc);
    if (!sf_getfloats ("eps",eps,nc)) 
	/* regularization parameter */
    {
	for (ic=0; ic < nc; ic++) {
	    eps[ic] = 1.0f;
	}
    } else {
	norm = 1.0f;
	for (ic=0; ic < nc; ic++) {
	    norm *= eps[ic];
	}
	norm = powf(norm,-1.0f/nc);
	for (ic=0; ic < nc; ic++) {
	    eps[ic] *= norm;
	}
    }
		
    if (!sf_getint("order",&nw)) nw=1;
    if (!sf_getint("nj",&nj)) nj=1;
    /* antialiasing */

    if (!sf_getbool("drift",&drift)) drift=false;
    /* if shift filter */

    if (!sf_getbool("verb",&verb)) verb = false;
    /* verbosity flag */

    ns = nc*n12;
    dd = sf_floatalloc(2*ns);
    s = sf_floatalloc(ns);
    d = sf_floatalloc(n12);
    ss = sf_floatalloc3(n1,n2,nc);

    planesig_init (nc, nw, nj, n1,n2, drift, ss, eps);

    sf_floatread (d,n12,in);
    sf_floatread (ss[0][0],ns,dips);

    planesig_rhs(d,dd);

    sf_solver (planesig_lop, sf_cgstep, ns, 2*ns, s, dd, niter, 
	       "verb", verb, "end");
    sf_cgstep_close();

    sf_floatwrite(s,ns,out);

    exit(0);
}
