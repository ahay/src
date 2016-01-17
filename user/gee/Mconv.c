/* 1-D convolution. */
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

#include "tcai1.h"
#include "icai1.h"

int main(int argc, char* argv[])
{
    int nx, nf, ny, i2, n2, lag;
    bool each, trans, adj;
    float *xx, *yy, *ff;
    sf_file in, out, filt;

    sf_init(argc,argv);
    in = sf_input("in");
    out = sf_output("out");
    filt = sf_input("filt");

    if (SF_FLOAT != sf_gettype(in) ||
	SF_FLOAT != sf_gettype(filt)) sf_error("Need float input");
    if (!sf_histint(in,"n1",&nx)) sf_error("No n1= in input");
    n2 = sf_leftsize(in,1);

    if (!sf_histint(filt,"n1",&nf)) sf_error("No n1= in filtin");

    xx = sf_floatalloc(nx);
    ff = sf_floatalloc(nf);

    if (!sf_getbool("trans",&trans)) trans=false;
    /* if y, transient convolution; if n, internal */
    if (!sf_getbool("each",&each)) each=false;
    /* if y, new filter for each trace */
    if (!sf_getbool("adj",&adj)) adj=false;
    /* adjoint flag */
    if (!sf_getint("lag",&lag)) lag=1;
    /* lag for internal convolution */

    if (trans) {
	if (adj) {
	    ny = nx-nf+1;
	} else {
	    ny = nx+nf-1;
	}
    } else {
	ny = nx;
    }
    yy = sf_floatalloc(ny);
    if (trans) sf_putint(out,"n1",ny);

    if (!each) sf_floatread (ff,nf,filt);

    for (i2=0; i2 < n2; i2++) {
        if (each) sf_floatread (ff,nf,filt);
        sf_floatread (xx,nx,in);
        if (trans) {
	    tcai1_init(nf,ff);
	    if (adj) {
		tcai1_lop (true, false,ny,nx,yy,xx);
	    } else {
		tcai1_lop (false, false,nx,ny,xx,yy);
	    } 
	} else {
	    icai1_init(nf,ff,lag);
	    if (adj) {
		icai1_lop (true, false,ny,nx,yy,xx);
	    } else {
		icai1_lop (false, false,nx,ny,xx,yy);
	    }
	}
	sf_floatwrite (yy,ny,out);
    }

    exit(0);
}


