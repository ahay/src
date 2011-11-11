/* Calculating frequency attenuation gradient. */
/*
  Copyright (C) 2011 Jilin University
  
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

#include <math.h>

int main(int argc, char* argv[])
{
    bool sign;
    int iw, nw, i1, n1, i2, n2, n1w;
    float d1, dw, w0, etotal, f1=0., f2=0., e1=0., e2=0., ecum, *fdg, *e;
    sf_complex *inp=NULL;
    sf_file in, out;
   
    sf_init(argc,argv);
    in = sf_input("in");
    out = sf_output("out");

    if (!sf_histint(in,"n1",&n1)) sf_error("No n1= in input");
    if (!sf_histfloat(in,"d1",&d1)) d1=1.;
    
    n2 = sf_leftsize(in,2);
    if (!sf_histint(in,"n2",&nw)) sf_error("No n2= in input");
    if (!sf_histfloat(in,"d2",&dw)) sf_error("No d2= in input");
    if (!sf_histfloat(in,"o2",&w0)) sf_error("No o2= in input");
    sf_unshiftdim(in, out, 2);
    n1w = n1*nw;
    
    if (SF_COMPLEX == sf_gettype(in)) {
	inp = sf_complexalloc(n1w);
	sf_settype(out,SF_FLOAT);
    }

    fdg = sf_floatalloc(n1);
    e   = sf_floatalloc(n1w);
 
    for (i2=0; i2 < n2; i2++) {
	sf_warning("slice %d of %d;",i2+1,n2);

	for (i1=0; i1 < n1; i1++) {
	    fdg[i1] = 0.;
	}
	if (SF_COMPLEX != sf_gettype(in)) {
	    sf_floatread(e,n1w,in);
	} else {
	    sf_complexread(inp,n1w,in);
	    for (i1=0; i1 < n1; i1++) {
		for (iw=0; iw < nw; iw++) {
		    e[iw*n1+i1] = sqrtf(crealf(inp[iw*n1+i1])*
					crealf(inp[iw*n1+i1])+
					cimagf(inp[iw*n1+i1])*
					cimagf(inp[iw*n1+i1]));
		}
	    }
	}
	for (i1=0; i1 < n1; i1++) {
	    etotal = 0.;
	    for (iw=0; iw < nw; iw++) {
		etotal += e[iw*n1+i1];
	    }
	    ecum = 0.;
	    sign = false;
	    for (iw=0; iw < nw; iw++) {
		ecum += e[iw*n1+i1];
		if (!sign && ecum >= 0.65*etotal) {
		    f1 = w0+iw*dw;
		    e1 = e[iw*n1+i1];
		    sign = true;
		}
		if (ecum >= 0.85*etotal) {
		    f2 = w0+iw*dw;
		    e2 = e[iw*n1+i1];
		    break;
		}
	    }
	    if (f1 == f2) {
		fdg[i1] = 0.;
	    } else {
		fdg[i1] = (e2-e1)/(f2-f1);
	    }
	}
	sf_floatwrite(fdg,n1,out);
    }
    sf_warning(".");
    
    exit(0);
}

/* 	$Id$	 */
