/* Smoothing in 1-D by simple regularization.
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

int main(int argc, char* argv[])
{ 
    int i1, n1, i2, n2, ir, nr;
    float *trace, eps,*diag, *offd;
    sf_tris slv;
    sf_file in, smooth;

    sf_init (argc, argv);
    in = sf_input("in");
    smooth = sf_output("out");

    if(!sf_histint(in,"n1",&n1)) sf_error ("No n1= in input");

    if (!sf_getfloat("eps",&eps)) eps=1.;
    /* smoothness parameter */
    eps = (eps*eps-1.)/12.;

    if (!sf_getint("repeat",&nr)) nr=1;
    /* repeat smoothing */

    n2 = sf_leftsize(in,1);

    trace = sf_floatalloc (n1);

    diag = sf_floatalloc (n1);
    offd = sf_floatalloc (n1);

    for (i1=0; i1 < n1; i1++) {
	diag[i1] = 1.+ 2.*eps;
	offd[i1] = -eps;
    }
    diag[0] = diag[n1-1] = 1.+eps;

    slv = sf_tridiagonal_init(n1);
    sf_tridiagonal_define (slv,diag,offd);

    for (i2=0; i2 < n2; i2++) {
	sf_floatread(trace,n1,in);

	for (ir=0; ir < nr; ir++){
	    sf_tridiagonal_solve(slv,trace);
	}

	sf_floatwrite(trace,n1,smooth);
    }

    exit (0);
}

/* 	$Id$	 */
