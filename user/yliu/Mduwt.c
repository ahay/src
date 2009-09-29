/* 1-D digital undecimated (stationary) wavelet transform by lifting scheme */
/*
  Copyright (C) 2008 University of Texas at Austin
   
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

#include "duwt.h"

int main(int argc, char *argv[])
{
    int n1, i2, n2, i1, max, scale;
    bool inv, adj, unit;
    char *type;
    float *pp, *qq;
    sf_file in, out;

    sf_init(argc,argv);

    in = sf_input("in");
    out = sf_output("out");

    if (!sf_histint(in,"n1",&n1)) sf_error("No n1= in input");
    n2 = sf_leftsize(in,1);


    if (!sf_getbool("inv",&inv)) inv=false;
    /* if y, do inverse transform */

    if (!sf_getbool("adj",&adj)) adj=false;
    /* if y, do adjoint transform */

    if (!sf_getbool("unit",&unit)) unit=false;
    /* if y, use unitary scaling */

    if (NULL == (type=sf_getstring("type"))) type="haar";
    /* [haar,linear,biorthogonal] wavelet type, the default is linear  */

    max=0;
    for (i1=1; i1 < n1; i1 *= 2) max++;

    if (!sf_getint("scale",&scale)) scale=max;
    /* decomposition scale */
    if (scale > max) sf_error("scale is over max");
    max = scale + 1;


    if (adj) { 
	n2 = sf_leftsize(in,2);
	sf_unshiftdim(in, out, 2);
	sf_putint(out,"n3",1);
    } else {
	n2 = sf_leftsize(in,1);
	sf_putint(out,"n2",max);
	sf_putint(out,"d2",1);
	(void) sf_shiftdim(in, out, 2);
    }

    pp = sf_floatalloc(n1);
    qq = sf_floatalloc(max*n1);

    wavelet_init(n1,inv,unit,type[0],max);

    for (i2=0; i2 < n2; i2++) {
	if (adj) {
	    sf_floatread(qq,n1*max,in);
	    wavelet_lop(adj,false,n1,n1*max,pp,qq);
	    sf_floatwrite(pp,n1,out);
	} else {
	    sf_floatread(pp,n1,in);
	    wavelet_lop(adj,false,n1,n1*max,pp,qq);
	    sf_floatwrite(qq,n1*max,out);
	}
    }

    exit(0);
}

/* 	$Id$	 */
