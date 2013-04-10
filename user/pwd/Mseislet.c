/* Seislet transform */
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

#include "seislet.h"

int main(int argc, char *argv[])
{
    int n1, n2, i3, n3, n12, order;
    bool inv, adj, unit, verb;
    char *type;
    float *pp, *qq, **dd, eps;
    sf_file in, out, dip;

    sf_init(argc,argv);

    in = sf_input("in");
    out = sf_output("out");
    dip = sf_input("dip");

    if (!sf_histint(in,"n1",&n1)) sf_error("No n1= in input");
    if (!sf_histint(in,"n2",&n2)) sf_error("No n2= in input");
    n3 = sf_leftsize(in,2);
    n12 = n1*n2;

    pp = sf_floatalloc(n12);
    qq = sf_floatalloc(n12);
    dd = sf_floatalloc2(n1,n2);

    if (!sf_getbool("inv",&inv)) inv=false;
    /* if y, do inverse transform */

    if (!sf_getbool("adj",&adj)) adj=false;
    /* if y, do adjoint transform */

    if (!sf_getfloat("eps",&eps)) eps=0.01;
    /* regularization */

    if (!sf_getbool("unit",&unit)) unit=false;
    /* if y, use unitary scaling */

    if (!sf_getint("order",&order)) order=1;
    /* accuracy order */

    if (NULL == (type=sf_getstring("type"))) type="linear";
    /* [haar,linear,biorthogonal] wavelet type, the default is linear  */

    if (!sf_getbool("verb",&verb)) verb = false;
    /* verbosity flag */

    seislet_init(n1,n2,inv,unit,eps,order,type[0]);
    seislet_set(dd);

    for (i3=0; i3 < n3; i3++) {
	if (verb) sf_warning("slice %d of %d;",i3+1,n3);

	sf_floatread(pp,n12,in);
	sf_floatread(dd[0],n12,dip);

	if (adj) {
	    seislet_lop(adj,false,n12,n12,qq,pp);
	} else {
	    seislet_lop(adj,false,n12,n12,pp,qq);
	}
	sf_floatwrite(qq,n12,out);
    }
    sf_warning(".");
    exit(0);
}
