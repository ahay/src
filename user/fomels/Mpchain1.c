/* Nonstationary Prony by chain of PEFs - linear operator */
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

#include "pchain.h"
#include "csmooth1.h"

int main(int argc, char* argv[])
{
    bool adj;
    int n, nc, n2, n1;
    sf_complex *xn, *x1, *dx, *r;
    sf_file inp, out, pef, sig;

    sf_init(argc,argv);
    inp = sf_input("in");
    pef = sf_input("pef");
    sig = sf_input("sig");
    out = sf_output("out");

    if (SF_COMPLEX != sf_gettype(inp)) sf_error("Need complex input");
    if (!sf_histint(pef,"n1",&n))  sf_error("No n1= in pef");
    if (!sf_histint(pef,"n2",&n2)) sf_error("No n2= in pef");

    nc = (n2+1)/2;
    n2 *= n;
    n1 = (n2+n)/2;

    sf_warning("nc=%d n=%d",nc,n);

    x1 = sf_complexalloc(n);
    sf_complexread(x1,n,sig);

    xn = sf_complexalloc(n2);
    sf_complexread(xn,n2,pef);

    pchain_init(n,nc,x1,xn,xn+n*nc);
    
    dx = sf_complexalloc(n2);
    r =  sf_complexalloc(n1);

    if (!sf_getbool("adj",&adj)) adj=false;
    /* adjoint flag */

    if (adj) {
	sf_complexread(r,n1,inp);
    } else {
	sf_complexread(dx,n2,inp);
    } 

    pchain_lop(adj,false,n2,n1,dx,r);

    if (adj) {
	sf_complexwrite(dx,n2,out);
    } else {
	sf_complexwrite(r,n1,out);
    } 

    exit(0);
}
