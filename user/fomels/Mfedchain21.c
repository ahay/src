/* Fast explicit diffusion using chains - linear operator (2-D) */
/*
  Copyright (C) 2022 University of Texas at Austin
  
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

#include "fedchain2.h"

int main(int argc, char* argv[])
{
    bool adj;
    int n, n1, n2, nc, n3;
    float *xn, *x1, *dx, *r;
    sf_file inp, out, wht, sig;

    sf_init(argc,argv);
    inp = sf_input("in");
    wht = sf_input("weight");
    sig = sf_input("sig");
    out = sf_output("out");

    if (SF_FLOAT != sf_gettype(inp)) sf_error("Need float input");
    if (!sf_histint(inp,"n1",&n1)) sf_error("No n1=");
    if (!sf_histint(inp,"n2",&n2)) sf_error("No n2=");
    if (!sf_histint(inp,"n3",&nc)) sf_error("No n3=");

    n = n1*n2;
    n3 = n*nc;

    x1 = sf_floatalloc(n);
    sf_floatread(x1,n,sig);

    xn = sf_floatalloc(n3);
    sf_floatread(xn,n3,wht);

    fedchain2_init(n1,n2,nc,x1,xn,xn+n);
    
    dx = sf_floatalloc(n3);
    r =  sf_floatalloc(n3);

    if (!sf_getbool("adj",&adj)) adj=false;
    /* adjoint flag */

    if (adj) {
	sf_floatread(r,n3,inp);
    } else {
	sf_floatread(dx,n3,inp);
    } 

    fedchain2_lop(adj,false,n3,n3,dx,r);

    if (adj) {
	sf_floatwrite(dx,n3,out);
    } else {
	sf_floatwrite(r,n3,out);
    } 

    exit(0);
}
