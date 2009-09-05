/* L1 regression */
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
#include <assert.h>

#include <rsf.h>

static int nc;
static float **cc;

static void fit(bool adj, bool add, int nm, int nd, float *m, float *d)
{
    int ic, id;

    sf_adjnull(adj, add, nm, nd, m, d);
    sf_copy_lop(adj, true, nd, nd, m, d);

    for (ic=0; ic < nc; ic++) {
	for (id=0; id < nd; id++) {
	    if (adj) {
		m[nd+ic] += cc[ic][id]*d[id];
	    } else {
		d[id] += cc[ic][id]*m[nd+ic];
	    }
	}
    }
}

int main(int argc, char* argv[])
{
    int nd, n1, niter, liter, iter, i1, nm;
    float *m, *d, perc;
    sf_file inp, reg, out;

    sf_init(argc,argv);
    inp = sf_input("in");
    reg = sf_input("reg");
    out = sf_output("out");

    if (!sf_histint(inp,"n1",&nd)) sf_error("No n1= in input");
    if (!sf_histint(reg,"n1",&n1) || n1 != nd)
	sf_error("Need n1=%d in reg",nd);
    if (!sf_histint(reg,"n2",&nc)) nc=1;
    nm = nd+nc;

    sf_putint(out,"n1",nc);
    
    d = sf_floatalloc(nd);
    m = sf_floatalloc(nm);
    cc = sf_floatalloc2(nd,nc);

    sf_floatread(d,nd,inp);
    sf_floatread(cc[0],nd*nc,reg);
    sf_fileclose(reg);

    if (!sf_getint("niter",&niter)) niter=100;
    /* number of POCS iterations */
    if (!sf_getint("liter",&liter)) liter=10;
    /* number of CG iterations */

    if (!sf_getfloat("perc",&perc)) perc=90.0;
    /* percentage for sharpening */

    sf_sharpen_init(nd,perc);

    /* initialize with zero */
    for (i1=0; i1 < nd+nc; i1++) {
	m[i1]=0.0;
    }

    for (iter=0; iter < niter; iter++) {
	sf_solver(fit,sf_cgstep,nm,nd,m,d,liter,"x0",m,"end");
	sf_sharpen(m);
	sf_weight_apply(nd,m);
    }
    
    sf_floatwrite(m+nd,nc,out);

    sf_fileclose(inp);
    exit(0);
}
