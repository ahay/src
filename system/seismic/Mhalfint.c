/* Half-order integration or differentiation. 

December 2012 program of the month:
http://ahay.org/rsflog/index.php?/archives/317-Program-of-the-month-sfhalfint.html
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
    bool adj, inv;
    int nn,n1,n2,i1,i2;
    float rho, *pp=NULL;
    sf_file in=NULL, out=NULL;

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

    if (inv) {
	sf_warning("%s half-order differentiation",adj? "anticausal":"causal");
    } else {
	sf_warning("%s half-order integration",adj? "anticausal":"causal");
    }

    nn = 2*kiss_fft_next_fast_size((n1+1)/2);
    pp = sf_floatalloc(nn);

    sf_halfint_init (inv, nn, rho);

    for (i2=0; i2 < n2; i2++) {
	sf_floatread (pp,n1,in);
	for (i1=n1; i1 < nn; i1++) {
	    pp[i1]=0.;
	}
	sf_halfint (adj, pp);
	sf_floatwrite (pp,n1,out);
    }

    exit(0);
}
