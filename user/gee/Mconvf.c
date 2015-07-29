/* 1-D convolution, adjoint is the filter. */
/*
  Copyright (C) 2013 University of Texas at Austin
  
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

#include "icaf1.h"

int main(int argc, char* argv[])
{
    int n1, nf, i2, n2, lag;
    bool adj;
    float *xx, *yy, *ff;
    sf_file inp, out, oth;

    sf_init(argc,argv);

    inp = sf_input("in");
    out = sf_output("out");
    oth = sf_input("other");

    if (!sf_getbool("adj",&adj)) adj=false;
    /* adjoint flag */

    if (SF_FLOAT != sf_gettype(inp) ||
	SF_FLOAT != sf_gettype(oth)) sf_error("Need float input");

    if (adj) {
	/* input data, output filter */
	if (!sf_histint(inp,"n1",&n1)) sf_error("No n1=");
	if (!sf_getint("nf",&nf)) sf_error("Need nf=");
	/* filter size */

	sf_putint(out,"n1",nf);
    } else {
	/* input filter, output data */
	if (!sf_histint(inp,"n1",&nf)) sf_error("No n1=");
	if (!sf_histint(oth,"n1",&n1)) sf_error("No n1=");

	sf_fileflush(out,oth);  /* copy data dimensions */
    }
    n2 = sf_leftsize(inp,1);

    xx = sf_floatalloc(n1);
    yy = sf_floatalloc(n1);
    ff = sf_floatalloc(nf);

    if (!sf_getint("lag",&lag)) lag=1;
    /* lag for internal convolution */

    icaf1_init(n1,xx,lag);

    for (i2=0; i2 < n2; i2++) {
        sf_floatread (xx,n1,oth);

	if (adj) {
	    sf_floatread (yy,n1,inp);
	    icaf1_lop (true, false,nf,n1,ff,yy);
	    sf_floatwrite (ff,nf,out);
	} else {
	    sf_floatread (ff,nf,inp);
	    icaf1_lop (false, false,nf,n1,ff,yy);
	    sf_floatwrite (yy,n1,out);
	}
    }

    exit(0);
}

