/* Inverse 1-D warping. with shifted-linear interpolation

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
#include "stretchsh.h"
#include "shprefilter.h"

int main(int argc, char* argv[])
{
    map4 mo;
    bool inv, each=true;
    int nt, n1, i2, n2, nw;
    float o1, d1, t0, dt, eps;
    float *trace, *str, *trace2;
    sf_file in, out, warp;

    sf_init(argc,argv);
    in = sf_input("in");
    out = sf_output("out");
    warp = sf_input("warp");

    if (!sf_getbool("inv",&inv)) inv=true;
    /* inversion flag */

    if (inv) {
	if (!sf_histint(in,"n1",&nt)) sf_error("No n1= in input");
	
	if (!sf_getint("n1",&n1)) n1=nt; /* output samples - for inv=y */
	if (!sf_getfloat("d1",&d1) && !sf_histfloat(in,"d1",&d1)) d1=1.;
	/*( d1=1 output sampling - for inv=y )*/
	if (!sf_getfloat("o1",&o1) && !sf_histfloat(in,"o1",&o1)) o1=0.;
	/*( o1=0 output origin - for inv=y )*/ 

	sf_putint(out,"n1",n1);
	sf_putfloat(out,"d1",d1);
	sf_putfloat(out,"o1",o1);
    } else {
	if (!sf_histint(in,"n1",&n1)) sf_error("No n1= in input");
	if (!sf_histfloat(in,"d1",&d1)) d1=1.;
	if (!sf_histfloat(in,"o1",&o1)) o1=0.;

	if (!sf_histint(warp,"n1",&nt)) sf_error("No n1= in warp");
	if (!sf_histfloat(warp,"d1",&dt)) dt=d1;
	if (!sf_histfloat(warp,"o1",&t0)) t0=o1;
	
	sf_putint(out,"n1",nt);
	sf_putfloat(out,"d1",dt);
	sf_putfloat(out,"o1",t0);
    }

    n2 = sf_leftsize(in,1);
    nw = sf_leftsize(warp,1);
    if (1 == nw) {
	each = false;
    } else if (n2 != nw) {
	sf_error("Need %d traces in warp, got %d",n2,nw);
    } 
    
    if (!sf_getfloat("eps",&eps)) eps=0.01;
    /* stretch regularization */

    trace = sf_floatalloc(nt);
    str = sf_floatalloc(nt);
    trace2 = sf_floatalloc(n1);

    mo = stretchsh_init (n1, o1, d1, nt, eps);

    for (i2=0; i2 < n2; i2++) {
	if (each || 0==i2) {
	    sf_floatread(str,nt,warp);
	    stretchsh_define (mo,str);
	}

	if (inv) {
	    sf_floatread(trace,nt,in);
	    stretchsh_apply (mo,trace,trace2);
	    sf_floatwrite (trace2,n1,out);
	} else {
	    sf_floatread(trace2,n1,in);
	    sf_int1sh_init (str, o1, d1, n1, sf_lin_int, 2, nt);
	    shprefilter(n1,trace2);
	    sf_int1_lop (false,false,n1,nt,trace2,trace);
	    sf_floatwrite (trace,nt,out);
	}
    }

    exit(0);
}
