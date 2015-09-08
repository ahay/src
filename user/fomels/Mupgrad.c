/* Causal gradient operator */
/*
  Copyright (C) 2015 University of Texas at Austin
  
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

int main(int argc, char* argv[])
{
    bool adj, grad;
    char key[4];
    int i, dim, n[SF_MAX_DIM], nt, i3, n3;
    float *x, **g, d[SF_MAX_DIM];
    sf_upgrad upg;
    sf_file inp, out, ref;

    sf_init(argc,argv);
    inp = sf_input("in");
    ref = sf_input("ref"); /* reference time */
    out = sf_output("out");

    if (!sf_getbool("adj",&adj)) adj=false;
    /* adjoint flag */

    if (!sf_getbool("grad",&grad)) grad=true;
    /* if y, gradient; if n, Laplacian */

    dim = sf_filedims(inp,n);
    if (grad && adj) dim--;

    sf_getint("dim",&dim);
    
    if (grad) {
	if (adj) {
	    sprintf(key,"n%d",dim+1);
	    if (n[dim] != dim) sf_error("Need %s=%d in input",key,dim);
	    sf_putint(out,key,1);
	    n3 = sf_unshiftdim(inp,out, dim+1);
	} else {
	    sprintf(key,"n%d",dim+1);
	    sf_putint(out,key,dim);
	    n3 = sf_shiftdim(inp, out, dim+1);
	}
    } else {
	n3 = sf_leftsize(inp,dim);
    }

    nt = 1;
    for (i=0; i < dim; i++) {
	sprintf(key,"d%d",i+1);
	if (!sf_histfloat(inp,key,d+i)) d[i]=1.0;
	nt *= n[i];
    }

    upg = sf_upgrad_init(dim,n,d);

    x = sf_floatalloc(nt);
    g = sf_floatalloc2(nt,dim);

    for (i3=0; i3 < n3; i3++) {
	sf_floatread(x,nt,ref);
	sf_upgrad_set(upg,x);
	
	if (grad) {
	    if (adj) {
		sf_floatread(g[0],nt*dim,inp);
		sf_upgrad_grad_adj(upg,x,g);
		sf_floatwrite(x,nt,out);
	    } else {
		sf_floatread(x,nt,inp);
		sf_upgrad_grad(upg,x,g);
		sf_floatwrite(g[0],nt*dim,out);
	    }
	} else { /* Laplacian */
	    sf_floatread(x,nt,inp);
	    sf_upgrad_grad(upg,x,g);
	    sf_upgrad_grad_adj(upg,x,g);
	    sf_floatwrite(x,nt,out);
	}
    }

    exit(0);
}
