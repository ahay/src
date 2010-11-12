/* 2-D beam forming. */
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

#include "pbeamform.h"

int main(int argc, char* argv[])
{
    bool adj, gauss;
    int n1, nc, nd, n3, i3, nb, n1c, n1d, order;
    float *dense, *coarse, **slope, d2;
    sf_file in, out, dip;

    sf_init(argc,argv);
    in = sf_input("in");
    dip = sf_input("dip");
    out = sf_output("out");

    if (!sf_histint(in,"n1",&n1)) sf_error("No n1= in input");
    if (!sf_histfloat(in,"d2",&d2)) d2=1.;

    if (!sf_getint("rect",&nb)) nb=3;
    /* smoothing radius */

    if (!sf_getbool("adj",&adj)) adj=true;
    /* adjoint flag */
    
    if (!sf_getbool("gauss",&gauss)) gauss=false;
    /* use pseudo-gaussian */

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

    if (!sf_getint("order",&order)) order=1;
    /* PWD accuracy order */

    pbeamform_init(gauss, n1, nd, order, nb);

    dense = sf_floatalloc(n1d);
    coarse = sf_floatalloc(n1c);
    slope = sf_floatalloc2(n1,nd);
 
    for (i3=0; i3 < n3; i3++) {
	sf_floatread(slope[0],n1d,dip);
	pbeamform_set(slope);

	if (adj) {
	    sf_floatread(dense,n1d,in);
	} else {
	    sf_floatread(coarse,n1c,in);
	}

	pbeamform_lop(adj,false,n1c,n1d,coarse,dense);

	if (adj) {
	    sf_floatwrite(coarse,n1c,out);
	} else {
	    sf_floatwrite(dense,n1d,out);
	}
    }

    exit(0);
}

