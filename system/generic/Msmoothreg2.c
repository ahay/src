/* Smoothing in 2-D by simple regularization.
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

#include <math.h>
#include <rsf.h>

int main(int argc, char* argv[])
{
    int i1, n1, n2, n12, ir, nr;
    float *trace=NULL, eps, *out=NULL;
    sf_file in=NULL, smooth=NULL;

    sf_init (argc, argv);
    in = sf_input("in");
    smooth = sf_output("out");

    if(!sf_histint(in,"n1",&n1)) sf_error ("No n1= in input");

    if (!sf_getfloat("eps",&eps)) eps=1.;
    /* smoothness parameter */
    eps = sqrtf((eps*eps-1.)/12.);

    if (!sf_getint("repeat",&nr)) nr=1;
    /* repeat smoothing */

    n2 = sf_leftsize(in,1);
    n12 = n1*n2;

    trace = sf_floatalloc (n12);
    out = sf_floatalloc (n12);

    sf_floatread(trace,n12,in);

    sf_igrad2_init(n1,n2);

    for (ir=0; ir < nr; ir++){
	sf_solver_reg (sf_copy_lop, sf_cgstep, sf_igrad2_lop, 
		       2*n12, n12, n12, out, trace, 
		       2*n12, eps, "end");
	sf_cgstep_close();

	for (i1=0; i1 < n12; i1++) {
	    trace[i1] = out[i1];
	}
    }

    sf_floatwrite(trace,n12,smooth);

    exit (0);
}
