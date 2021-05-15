/* Gaussian beam forming for complex data. */
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

#include "cbeamform1.h"

int main(int argc, char* argv[])
{
    bool adj;
    int n1, nc, nd, n3, i3, nb, n1c, n1d, nrep;
    float d2;
    sf_complex *dense, *coarse; 
    sf_file in, out;

    sf_init(argc,argv);
    in = sf_input("in");
    out = sf_output("out");

    if (SF_COMPLEX != sf_gettype(in)) sf_error("Need complex input");
    if (!sf_histint(in,"n1",&n1)) sf_error("No n1= in input");
    if (!sf_histfloat(in,"d2",&d2)) d2=1.;

    if (!sf_getint("rect",&nb)) nb=3;
    /* smoothing radius */

    if (!sf_getbool("adj",&adj)) adj=true;
    /* adjoint flag */
    
    if (!sf_getint("repeat",&nrep)) nrep=1;
    /* repeat triangle smoothing */

    if (adj) {
	if (!sf_histint(in,"n2",&nd)) sf_error("No n2= in input");
	nc = (nd-1)/nb+1;
	sf_putint(out,"n2",nc);
	sf_putfloat(out,"d2",d2*nb);
    } else {
	if (!sf_histint(in,"n2",&nc)) sf_error("No n2= in input");
	nd = (nc-1)*nb+1;
	sf_putint(out,"n2",nd);
	sf_putfloat(out,"d2",d2/nb);
    } 
    n3 = sf_leftsize(in,2);

    n1d = n1*nd;
    n1c = n1*nc;

    cbeamform1_init(nrep, n1, nd, nb);

    dense = sf_complexalloc(n1d);
    coarse = sf_complexalloc(n1c);
 
    for (i3=0; i3 < n3; i3++) {
	if (adj) {
	    sf_complexread(dense,n1d,in);
	} else {
	    sf_complexread(coarse,n1c,in);
	}

	cbeamform1_lop(adj,false,n1c,n1d,coarse,dense);

	if (adj) {
	    sf_complexwrite(coarse,n1c,out);
	} else {
	    sf_complexwrite(dense,n1d,out);
	}
    }

    exit(0);
}

