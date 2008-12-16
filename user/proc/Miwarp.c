/* Inverse 1-D warping */
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

#include "stretch4.h"

int main(int argc, char* argv[])
{
    map4 mo;
    bool inv;
    int nt, n1, i2, n2;
    float o1, d1, eps;
    float *trace, *str, *trace2;
    sf_file in, out, warp;

    sf_init(argc,argv);
    in = sf_input("in");
    out = sf_output("out");
    warp = sf_input("warp");

    if (!sf_histint(in,"n1",&nt)) sf_error("No n1= in input");
    n2 = sf_leftsize(in,1);
    
    if (!sf_getint("n1",&n1)) n1=nt;
    if (!sf_getfloat("d1",&d1) && !sf_histfloat(in,"d1",&d1)) d1=1.;
    if (!sf_getfloat("o1",&o1) && !sf_histfloat(in,"o1",&o1)) o1=0.;

    if (!sf_getbool("inv",&inv)) inv=true;
    /* inversion flag */

    sf_putint(out,"n1",n1);
    sf_putfloat(out,"d1",d1);
    sf_putfloat(out,"o1",o1);
    
    if (!sf_getfloat("eps",&eps)) eps=0.01;
    /* stretch regularization */

    trace = sf_floatalloc(nt);
    str = sf_floatalloc(nt);
    trace2 = sf_floatalloc(n1);

    mo = stretch4_init (n1, o1, d1, nt, eps);

    for (i2=0; i2 < n2; i2++) {
	sf_floatread(trace,nt,in);
	sf_floatread(str,nt,warp);

	stretch4_define (mo,str);
	if (inv) {
	    stretch4_apply (mo,trace,trace2);
	} else {
	    stretch4_invert (mo,trace2,trace);
	}
	sf_floatwrite (trace2,n1,out);
    }

    exit(0);
}
