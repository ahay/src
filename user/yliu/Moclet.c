/* Seislet transform in log-stretched frequency-offset-midpoint domain 
Forward transform (adj=y inv=y)   m=T[d]
Adjoint transform (adj=y inv=n)   m=T^(-1)'[d]
Inverse transform (adj=n inv=y/n) d=T^(-1)[m]
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

#include "oclet.h"

int main(int argc, char *argv[])
{
    int nx, nh, iw, nw, nxh, i4, n4;
    float x0, dx, h0, dh, w0, dw, w; 
    bool inv, adj, unit, verb;
    char *type;
    sf_complex *pp, *qq;
    sf_file in, out;

    sf_init(argc,argv);

    in = sf_input("in");
    out = sf_output("out");

    if (SF_COMPLEX != sf_gettype(in)) sf_error("Need complex input");

    if (!sf_histint(in,"n1",&nx)) sf_error("No n1= in input");
    if (!sf_histint(in,"n2",&nh)) sf_error("No n2= in input");
    if (!sf_histint(in,"n3",&nw)) sf_error("No n3= in input");

    if (!sf_histfloat(in,"d1",&dx)) sf_error("No d1= in input");
    if (!sf_histfloat(in,"d2",&dh)) sf_error("No d2= in input");
    if (!sf_histfloat(in,"d3",&dw)) sf_error("No d3= in input");

    if (!sf_histfloat(in,"o1",&x0)) sf_error("No o1= in input");
    if (!sf_histfloat(in,"o2",&h0)) sf_error("No o2= in input");
    if (!sf_histfloat(in,"o3",&w0)) sf_error("No o3= in input");

    n4 = sf_leftsize(in,3);
    nxh = nx*nh;
    dh /= dx;
    h0 /= dx;
    w0 *= 2.*SF_PI;
    dw *= 2.*SF_PI;

    pp = sf_complexalloc(nxh);
    qq = sf_complexalloc(nxh);

    if (!sf_getbool("inv",&inv)) inv=true;
    /* if y, do inverse transform */

    if (!sf_getbool("adj",&adj)) adj=true;
    /* if y, do adjoint transform */

    if (!sf_getbool("unit",&unit)) unit=false;
    /* if y, use unitary scaling */

    if (!sf_getbool("verb",&verb)) verb = false;
    /* verbosity flag */

    if (NULL == (type=sf_getstring("type"))) type="linear";
    /* [haar,linear,biorthogonal] wavelet type, the default is linear  */

    oclet_init(nx,nh,dh,dw,h0,inv,unit,type[0]);

    for (i4=0; i4 < n4; i4++) {
	for (iw=0; iw < nw; iw++) {
	    if (verb) sf_warning("frequency %d of %d;",iw+1,nw);
	    sf_complexread(pp,nxh,in);
	    w = w0+iw*dw;
	    if (adj) {
		oclet_lop(adj,false,nxh,nxh,qq,pp,w);
	    } else {
		oclet_lop(adj,false,nxh,nxh,pp,qq,w);
	    }
	    sf_complexwrite(qq,nxh,out);
	}
    }

    sf_warning(".");
    exit(0);
}
/* 	$Id$	 */
