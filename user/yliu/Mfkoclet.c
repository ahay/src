/* 1-D seislet transform using omega-wavenumber offset continuation 
Forward transform (adj=n inv=y/n) m=T[d]
Inverse transform (adj=y inv=y)   d=T^(-1)[d]
Adjoint transform (adj=y inv=n)   d=T'[d]
*/
/*
  Copyright (C) 2009 University of Texas at Austin
   
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

#include "fkoclet.h"

int main(int argc, char *argv[])
{
    int nk, nh, iw, nw, i4, n4, ik;
    float k0, dk, h0, dh, w0, dw, w, k, eps; 
    bool inv, verb, adj, dwt;

    char *type;
    sf_complex *pp, *qq;
    sf_file in, out;

    sf_init(argc,argv);

    in = sf_input("in");
    out = sf_output("out");

    if (SF_COMPLEX != sf_gettype(in)) sf_error("Need complex input");

    if (!sf_histint(in,"n1",&nh)) sf_error("No n1= in input");
    if (!sf_histint(in,"n2",&nk)) sf_error("No n2= in input");
    if (!sf_histint(in,"n3",&nw)) sf_error("No n3= in input");

    if (!sf_histfloat(in,"d1",&dh)) sf_error("No d1= in input");
    if (!sf_histfloat(in,"d2",&dk)) sf_error("No d2= in input");
    if (!sf_histfloat(in,"d3",&dw)) sf_error("No d3= in input");

    if (!sf_histfloat(in,"o1",&h0)) sf_error("No o1= in input");
    if (!sf_histfloat(in,"o2",&k0)) sf_error("No o2= in input");
    if (!sf_histfloat(in,"o3",&w0)) sf_error("No o3= in input");

    n4 = sf_leftsize(in,3);

    if (!sf_getbool("inv",&inv)) inv=false;
    /* if y, do inverse transform */

    if (!sf_getbool("adj",&adj)) adj=false;
    /* if y, do adjoint transform */

    if (!sf_getbool("dwt",&dwt)) dwt=false;
    /* if y, do wavelet transform */

    if (!sf_getbool("verb",&verb)) verb=true;
    /* verbosity flag */

    if (!sf_getfloat("eps",&eps)) eps=0.01;
    /* regularization */

    pp = sf_complexalloc(nh);   /* data space */
    qq = sf_complexalloc(nh);   /* model space */

    if (NULL == (type=sf_getstring("type"))) type="biorthogonal";
    /* [haar,linear,biorthogonal] wavelet type, the default is biorthogonal  */

    fkoclet_init(nh,nk,dh,dk,dw,h0,k0,inv,false,dwt,eps*eps,type[0]);

    /* loop over n4 */
    for (i4=0; i4 < n4; i4++) {
	for (iw=0; iw < nw; iw++) { /* loop over frequency */
	    if (verb) sf_warning("frequency %d of %d;",iw+1,nw);
	    w = w0 + iw*dw;
	    for (ik=0; ik < nk; ik++) { /* loop over wavenumber */
		k = k0 + ik*dk;
		if (adj) {
		    sf_complexread(qq,nh,in);
		} else {
		    sf_complexread(pp,nh,in);
		} 

		if (adj) {
		    fkoclet_lop(false,false,nh,nh,qq,pp,w,k);
		    sf_complexwrite(pp,nh,out);
		} else {
		    fkoclet_lop(true,false,nh,nh,qq,pp,w,k);
		    sf_complexwrite(qq,nh,out);
		} 
	    }
	}
    }
    sf_warning(".");
    exit(0);
}
/* 	$Id: Mfkoclet.c 9673 2012-12-13 04:51:41Z yang_liu $	 */
