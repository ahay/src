/* Resampling with triangle weights. */
/*
  Copyright (C) 2007 University of Texas at Austin
  
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

#include "tristack.h"

int main(int argc, char* argv[])
{
    bool adj, gauss;
    int nc, nd, nb, i2, n2;
    float *c, *d, d1;
    sf_file in, out;

    sf_init(argc,argv);
    in = sf_input("in");
    out = sf_output("out");

    if (!sf_getbool("adj",&adj)) adj=true;
    /* adjoint flag */

    if (!sf_getint("rect",&nb)) nb=1;
    /* smoothing radius */

    if (!sf_getbool("gauss",&gauss)) gauss=false;
    /* use pseudo-gaussian */

    if (!sf_histfloat(in,"d1",&d1)) d1=1.;

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
    n2 = sf_leftsize(in,1);

    c = sf_floatalloc(nc);
    d = sf_floatalloc(nd);

    tristack_init(gauss, nd, nb);

    for (i2=0; i2 < n2; i2++) {
	if (adj) {
	    sf_floatread(d,nd,in);
	} else {
	    sf_floatread(c,nc,in);
	} 

	tristack(adj,false,nc,nd,c,d);

	if (adj) {
	    sf_floatwrite(c,nc,out);
	} else {
	    sf_floatwrite(d,nd,out);
	} 
    }

    exit(0);
}
