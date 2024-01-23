/* Fast explicit diffusion using chains - linear operator */
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

#include "fedchain.h"

int main(int argc, char* argv[])
{
    bool adj;
    int n, n1, n2;
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

    n = n1*n2;

    x1 = sf_floatalloc(n1);
    sf_floatread(x1,n1,sig);

    xn = sf_floatalloc(n);
    sf_floatread(xn,n,wht);

    fedchain_init(n1,n2,x1,xn,xn+n1);
    
    dx = sf_floatalloc(n);
    r =  sf_floatalloc(n);

    if (!sf_getbool("adj",&adj)) adj=false;
    /* adjoint flag */

    if (adj) {
	sf_floatread(r,n,inp);
    } else {
	sf_floatread(dx,n,inp);
    } 

    fedchain_lop(adj,false,n,n,dx,r);

    if (adj) {
	sf_floatwrite(dx,n,out);
    } else {
	sf_floatwrite(r,n,out);
    } 

    exit(0);
}
