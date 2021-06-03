/* Linear regression */
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

#include "l1.h"
#include "l1step.h"

static float **cc;

static void fit(bool adj, bool add, int nm, int nd, float *m, float *d)
{
    int im, id;

    sf_adjnull(adj, add, nm, nd, m, d);

    for (im=0; im < nm; im++) {
	for (id=0; id < nd; id++) {
	    if (adj) {
		m[im] += cc[im][id]*d[id];
	    } else {
		d[id] += cc[im][id]*m[im];
	    }
	}
    }
}

int main(int argc, char* argv[])
{
    bool verb;
    char ni[5];
    int i, nc, nd, n1, niter, n1iter, method, dim;
    float *c, *d, perc, fact;
    char *type;
    sf_file inp, reg, out;

    sf_init(argc,argv);
    inp = sf_input("in");
    reg = sf_input("reg");
    out = sf_output("out");

    if (!sf_getint("dim",&dim)) dim=1;
    /* number of dimensions */
    nc = sf_leftsize(reg,dim);
    nd = sf_filesize(inp);
    n1 = sf_filesize(reg)/nc;

    if (n1 != nd) sf_error("Dimension mismatch");
    sf_putint(out,"n1",nc);
    for (i=2; i <= dim; i++) {
	sprintf(ni,"n%d",i);
	sf_putint(out,ni,1);
    }
    
    d = sf_floatalloc(nd);
    c = sf_floatalloc(nc);
    cc = sf_floatalloc2(nd,nc);

    sf_floatread(d,nd,inp);
    sf_floatread(cc[0],nd*nc,reg);
    sf_fileclose(reg);

    if (!sf_getint("niter",&niter)) niter=10;
    /* number of CG iterations */

    if (!sf_getint("method",&method)) method=2;
    /* method (L1-like or L2) */

    if (!sf_getbool("verb",&verb)) verb=false;
    /* verbosity flag */

    switch (method) {
	case 1:
	    if (!sf_getint("n1iter",&n1iter)) n1iter=10;
	    /* number of POCS iterations */
	    
	    if (!sf_getfloat("perc",&perc)) perc=90.0;
	    /* percentage for sharpening */

	    if (!sf_getfloat("fact",&fact)) fact=1.5;
	    /* percentage for sharpening */

	    if (NULL == (type = sf_getstring("type"))) type="threshold";
    /* thresholding type */

	    l1_init(nd,n1iter,perc,fact,type,verb);
	    
	    sf_solver(fit,l1step,nc,nd,c,d,niter,"verb",verb,"end");
	    break;
	case 2:
	    sf_solver(fit,sf_cgstep,nc,nd,c,d,niter,"verb",verb,"end");
	    break;
	default:
	    sf_error("Unknown method %d",method);
    }
    
    sf_floatwrite(c,nc,out);

    exit(0);
}
