/* Chain of PWDs */
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

#include "pwdchain.h"
#include "smooth1.h"

int main(int argc, char* argv[])
{
    int m1, m2, i, n, nc, n2, n1;
    float *xn, *x1, *x2, *r;
    sf_file inp, out, sig;

    sf_init(argc,argv);
    inp = sf_input("in");
    sig = sf_input("sig");
    out = sf_output("out");

    if (SF_FLOAT != sf_gettype(inp)) sf_error("Need float input");
    if (!sf_histint(inp,"n1",&m1)) sf_error("No n1= in dip");
    if (!sf_histint(inp,"n2",&m2)) sf_error("No n2= in dip");
    if (!sf_histint(inp,"n3",&n2)) sf_error("No n3= in dip");

    n = m1*m2;

    nc = (n2+1)/2;

    sf_putint(out,"n3",nc);
 
    n2 *= n; 
    n1 = (n2+n)/2;

    sf_warning("nc=%d n=%d",nc,n);
    
    x1 = sf_floatalloc(n);
    sf_floatread(x1,n,sig);

    x2 = sf_floatalloc(n);
    for (i=0; i < n; i++) {
	x2[i] = 0.0f;
    }
	
    xn = sf_floatalloc(n2);
    sf_floatread(xn,n2,inp);

    pwdchain_init(m1,m2,1,nc,x1,xn,xn+n*nc);
 
    r =  sf_floatalloc(n1);
 
    pwdchain_apply(x2,r);

    sf_floatwrite(r,n1,out);
 
    exit(0);
}
