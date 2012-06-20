/* 1-D beam forming. */
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

#include "pbeamform1.h"

int main(int argc, char* argv[])
{
    bool adj, gauss;
    int nc, nd, n3, i3, nb;
    float *dense, *coarse, d1;
    sf_file in, out;

    sf_init(argc,argv);
    in = sf_input("in");
    out = sf_output("out");

    if (!sf_histfloat(in,"d1",&d1)) d1=1.;

    if (!sf_getint("rect",&nb)) nb=3;
    /* smoothing radius */

    if (!sf_getbool("adj",&adj)) adj=true;
    /* adjoint flag */
    
    if (!sf_getbool("gauss",&gauss)) gauss=false;
    /* use pseudo-gaussian */

    if (adj) {
	if (!sf_histint(in,"n1",&nd)) sf_error("No n1= in input");
	nc = (nd-1)/nb+1;
	sf_putint(out,"n1",nc);
	sf_putfloat(out,"d1",d1*nb);
    } else {
	if (!sf_histint(in,"n1",&nc)) sf_error("No n1= in input");
	nd = (nc-1)*nb+1;
	sf_putint(out,"n1",nd);
	sf_putfloat(out,"d1",d1/nb);
    } 
    n3 = sf_leftsize(in,1);

    pbeamform1_init(gauss, nd, nb);

    dense = sf_floatalloc(nd);
    coarse = sf_floatalloc(nc);
 
    for (i3=0; i3 < n3; i3++) {
	if (adj) {
	    sf_floatread(dense,nd,in);
	} else {
	    sf_floatread(coarse,nc,in);
	}

	pbeamform1_lop(adj,false,nc,nd,coarse,dense);

	if (adj) {
	    sf_floatwrite(coarse,nc,out);
	} else {
	    sf_floatwrite(dense,nd,out);
	}
    }

    exit(0);
}

