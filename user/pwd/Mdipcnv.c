/* Convert dip to convection */
/*
  Copyright (C) 2025 University of Texas at Austin
  
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

#include "apfilt.h"

int main (int argc, char* argv[])
{
    int order, n1, n2, i1, i2, i, nw;
    float *p, *a, **b;
    sf_file cnv, dip;

    sf_init(argc, argv);
    dip = sf_input("in");
    cnv = sf_output("out");

    if (SF_FLOAT != sf_gettype(dip)) sf_error("Need float input");
    if (!sf_histint(dip,"n1",&n1)) sf_error("Need n1= in input");
    if (!sf_histint(dip,"n2",&n2)) sf_error("Need n2= in input");

    if (!sf_getint("order",&order)) order=1;
    nw = 2*order;

    sf_putint(cnv,"n2",nw);
    sf_shiftdim(dip, cnv, 2);

    p = sf_floatalloc(n1);
    a = sf_floatalloc(nw+1);
    b = sf_floatalloc2(n1,nw);

    apfilt_init(order);
    
    for (i2=0; i2 < n2; i2++) {
      sf_floatread(p,n1,dip);
      for (i1=0; i1 < n1; i1++) {
	passfilter(p[i1],a);
	for (i=1; i <= order; i++) {
	  b[2*i-2][i1] = a[order+i];
	  b[2*i-1][i1] = a[order-i];
	}
      }
      sf_floatwrite(b[0],n1*nw,cnv);
    }

    exit(0);
}
