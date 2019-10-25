/* Signal and noise separation using plane-wave destruction filters.  

If n3=1 in the output, outputs both signal and noise. Otherwise, only signal.
*/
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
#include "planesignoi.h"

    int main(int argc, char* argv[])
{
    int i, n1, n2, n12, nj1, nj2, niter, nw, n3, i3;
    float eps, *d, *s, *nd, **nn, **ss;
    bool verb, drift;
    sf_file in, out, ndip, sdip;

    sf_init (argc,argv);
    in = sf_input("in");
    ndip = sf_input("ndip");
    sdip = sf_input("sdip");
    out = sf_output("out");

    if (!sf_histint(in,"n1",&n1)) sf_error("No n1= in input");
    if (!sf_histint(in,"n2",&n2)) sf_error("No n2= in input");
    n12 = n1*n2;

    n3 = sf_leftsize(in,2);
    if (1==n3) sf_putint (out,"n3",2);

    if (!sf_getint ("niter",&niter)) niter=50;
    /* maximum number of iterations */

    if (!sf_getfloat ("eps",&eps)) eps=1.;
    /* regularization parameter */

    if (!sf_getint("order",&nw)) nw=1;
    /* [1,2,3] accuracy order */
    if (nw < 1 || nw > 3) 
	sf_error ("Unsupported nw=%d, choose between 1 and 3",nw);
    if (!sf_getint("nj1",&nj1)) nj1=1;
    /* antialiasing for noise dip */
    if (!sf_getint("nj2",&nj2)) nj2=1;
    /* antialiasing for signal dip */

    if (!sf_getbool("drift",&drift)) drift=false;
    /* if shift filter */

    if (!sf_getbool("verb",&verb)) verb = false;
    /* verbosity flag */

    nd = sf_floatalloc(2*n12);
    s = sf_floatalloc(n12);
    d = sf_floatalloc(n12);
    nn = sf_floatalloc2(n1,n2);
    ss = sf_floatalloc2(n1,n2);

    for (i=0; i < n12; i++) {
	nd[n12+i] = 0.;
    }

    planesignoi_init (nw, nj1,nj2, n1,n2, drift, nn, ss, eps);

    for (i3=0; i3 < n3; i3++) {
	sf_floatread (d,n12,in);
	sf_floatread (nn[0],n12,ndip);
	sf_floatread (ss[0],n12,sdip);

	allpass22_init(allpass2_init(nw,nj1,n1,n2,drift,nn));
	allpass21_lop (false,false,n12,n12,d,nd);
	

	sf_solver (planesignoi_lop, sf_cgstep, n12, n12*2, s, nd, niter, 
		   "verb", verb, "end");

	sf_cgstep_close();

	sf_floatwrite(s,n12,out);

	if (1==n3) {
	    for (i=0; i < n12; i++) {
		d[i] -= s[i];
	    }
	    sf_floatwrite(d,n12,out);
	}
    }

    exit(0);
}
