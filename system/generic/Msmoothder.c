/* Smooth first derivative on the first axis.

Applies D/(I + eps*D'D)
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
    int i1, n1, i2, n2;
    float *trace=NULL, *dtrace=NULL, d1, eps, *diag=NULL, **offd=NULL;
    sf_bands slv;
    sf_file in, der;

    sf_init (argc, argv);
    in = sf_input("in");
    der = sf_output("out");

    if(!sf_histint(in,"n1",&n1)) sf_error ("No n1= in input");
    if(!sf_histfloat(in,"d1",&d1)) d1=1.;

    if (!sf_getfloat("eps",&eps)) eps=0.2;
    /* smoothness parameter */

    n2 = sf_leftsize(in,1);

    trace = sf_floatalloc (n1);
    dtrace = sf_floatalloc (n1);
    diag = sf_floatalloc (n1);
    offd = sf_floatalloc2 (n1,2);

    for (i1=0; i1 < n1; i1++) {
	diag[i1] = 1.+ 6.*eps;
	offd[0][i1] = -4.*eps;
	offd[1][i1] = eps;
    }
    diag[0] = diag[n1-1] = 1.+eps;
    diag[1] = diag[n1-2] = 1.+5.*eps;
    offd[0][0] = offd[0][n1-2] = -2.*eps;

    slv = sf_banded_init(n1,2);
    sf_banded_define (slv,diag,offd);

    for (i2=0; i2 < n2; i2++) {
	sf_floatread(trace,n1,in);

	/* smooth */
	sf_banded_solve(slv,trace);

	/* differentiate */
	for (i1=0; i1 < n1-1; i1++) {
	    dtrace[i1] = trace[i1+1]-trace[i1];
	    dtrace[i1] /= d1;
	}
	dtrace[n1-1] = dtrace[n1-2];

	sf_floatwrite(dtrace,n1,der);
    }

    exit (0);
}
