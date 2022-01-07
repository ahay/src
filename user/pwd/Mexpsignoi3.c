/* Signal and noise separation using frequency components. 
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

#include "expont4.h"
#include "expsignoi3.h"

int main(int argc, char* argv[])
{
    int i, n1, n2,n3, n123, niter;
    float eps, *d, *s, *nd, **nnss;
    bool verb;
    sf_file in, out, freq;

    sf_init (argc,argv);
    in = sf_input("in");
    freq = sf_input("freq");
    out = sf_output("out");

    if (!sf_histint(in,"n1",&n1)) sf_error("No n1= in input");
    if (!sf_histint(in,"n2",&n2)) sf_error("No n2= in input");
    if (!sf_histint(in,"n3",&n3)) sf_error("No n3= in input");
    n123 = n1*n2*n3;

    if (!sf_getint ("niter",&niter)) niter=50;
    /* maximum number of iterations */

    if (!sf_getfloat ("eps",&eps)) eps=1.;
    /* regularization parameter */

    sf_putint (out,"n4",2);

    if (!sf_getbool("verb",&verb)) verb = false;
    /* verbosity flag */

    nd = sf_floatalloc(4*n123);
    s = sf_floatalloc(n123);
    d = sf_floatalloc(n123);
    nnss = sf_floatalloc2(n123,8);

    sf_floatread (d,n123,in);
    sf_floatread (nnss[0],n123*8,freq);

    expont4_init(n1,n2,n3,nnss[0],nnss[1],nnss[4],nnss[5]);
    expont4_lop (false,false,n123,2*n123,d,nd);
    for (i=n123; i < 3*n123; i++) {
	nd[n123+i] = 0.;
    }

    expsignoi3_init (n1,n2,n3, eps, nnss);
    sf_solver (expsignoi3_lop, sf_cgstep, n123, n123*4, s, nd, niter, 
	       "verb", verb, "end");

    sf_floatwrite(s,n123,out);
    for (i=0; i < n123; i++) {
	d[i] -= s[i];
    }
    sf_floatwrite(d,n123,out);

    exit(0);
}
