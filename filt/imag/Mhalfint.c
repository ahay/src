/* Half-order integration or differentiation.
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

#include "halfint.h"

int main(int argc, char* argv[])
{
    bool adj, inv;
    int n1,n2,i2;
    float rho, *pp;
    sf_file in, out;

    sf_init(argc,argv);
    in = sf_input("in");
    out = sf_output("out");

    if (!sf_histint(in,"n1",&n1)) sf_error("No n1= in input");
    n2 = sf_leftsize(in,1);

    if (!sf_getbool("adj",&adj)) adj=false;
    /* If y, apply adjoint */
    if (!sf_getbool("inv",&inv)) inv=false;
    /* If y, do differentiation instead of integration */
    if (!sf_getfloat("rho",&rho)) rho = 1.-1./n1;
    /* Leaky integration constant */

    if (adj) {
	sf_warning("%s half-order differentiation",inv? "anticausal":"causal");
    } else {
	sf_warning("%s half-order integration",inv? "anticausal":"causal");
    }

    pp = sf_floatalloc(n1);

    halfint_init (adj, inv, n1, rho);

    for (i2=0; i2 < n2; i2++) {
	sf_floatread (pp,n1,in);
	halfint (pp);
	sf_floatwrite (pp,n1,out);
    }

    exit(0);
}

/* 	$Id: Mhalfint.c,v 1.6 2004/07/02 11:54:20 fomels Exp $	 */
