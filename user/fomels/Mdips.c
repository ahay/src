/* Estimate a number of constant dips using plane-wave destruction. */
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
#include "dips.h"
#include "matmult.h"

int main(int argc, char* argv[])
{
    int n1, n2, n12, n3, i3, nd, id, niter, iter, nw, nj;
    float **dat, *dip, *dip0, *inc, **aa, *b;
    bool verb;
    sf_file in, out;

    sf_init(argc,argv);
    in = sf_input("in");
    out = sf_output("out");

    if (!sf_histint(in,"n1",&n1)) sf_error("No n1= in input");
    if (!sf_histint(in,"n2",&n2)) sf_error("No n2= in input");
    n12 = n1*n2;
    n3 = sf_leftsize(in,2);

    if (!sf_getint("nd",&nd)) sf_error("Need nd=");
    /* number of dips */
    sf_putint(out,"n1",nd);
    sf_putint(out,"n2",1);

    dat = sf_floatalloc2(n1,n2);
    dip = sf_floatalloc(nd);
    dip0 = sf_floatalloc(nd);
    inc = sf_floatalloc(nd);
    aa = sf_floatalloc2(nd,n12);
    b = sf_floatalloc(n12);

    if (!sf_getint("order",&nw)) nw=1;
    /* [1,2,3] accuracy order */
    if (nw < 1 || nw > 3) 
	sf_error ("Unsupported nw=%d, choose between 1 and 3",nw);
    if (!sf_getint("nj",&nj)) nj=1;
    /* antialiasing */

    if (!sf_getbool("verb",&verb)) verb=false;
    /* verbosity flag */
    
    sf_floatread(dat[0],n12,in);
    dips_init(nd, nw, nj, n1, n2, dat);
    matmult_init (aa);

    if (!sf_getint("niter",&niter)) niter=10;
    /* number of iterations */
    
    if (!sf_getfloats("dips",dip0,nd)) { 
	/* initial dips */
	for (id=0; id < nd; id++) {
	    dip0[id] = (nd > 1)? -1.+2.*id/(nd-1): 0.;
	}
    }

    for (i3=0; i3 < n3; i3++) {
	/* initialize dips */
	for (id=0; id < nd; id++) {
	    dip[id] = dip0[id];
	}
	
	for (iter=0; iter < niter; iter++) {
	    if (verb) sf_warning("------------\niteration %d",iter+1);
	    dips(dip, b, aa);
	    sf_solver(matmult_lop,sf_cgstep,nd,n12,inc,b,nd+1,"end");
	    sf_cgstep_close();
	    for (id=0; id < nd; id++) {
		if (verb) sf_warning("dip%d=%g",id+1,dip[id]);
		dip[id] += inc[id];
	    }
	}
	sf_floatwrite(dip,nd,out);
    }

    exit(0);
}
