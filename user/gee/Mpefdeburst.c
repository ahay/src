/* Burst noise removal using PEF. */
/*
  Copyright (C) 2006 University of Texas at Austin
  
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

#include "pefest.h"
#include "fixbad.h"
#include "bound.h"
#include "printfilter.h"

int main(int argc, char* argv[])
{
    int n1, i1, n2, i2, na, tempna, ia, niter, center=0;
    float *data;
    sf_filter aa;
    sf_file in, out;

    sf_init(argc,argv);
    in = sf_input("in");
    out = sf_output("out");

    if (!sf_histint(in,"n1",&n1)) sf_error("No n1= in input");
    n2 = sf_leftsize(in,1);
    
    if (!sf_getint("na",&na)) na=3; 
    /* PEF length */
    na--;

    if (!sf_getint("niter",&niter)) niter=10;
    /* number of iterations */

    data = sf_floatalloc(n1);
    aa = sf_allocatehelix(na);
    tempna = na;

    for (i2=0; i2 < n2; i2++) {
	sf_floatread(data,n1,in);
	
	aa->flt[0]=-2;
	aa->flt[1]=1;
	
	for (ia=0; ia < na; ia++) {
	    aa->lag[ia] = ia+1;
	}
	
	bound (1, &n1, &n1, &na, aa);
	pefest (na * 2, n1, data, aa);
	
	for (i1=0; i1 < n1; i1++) {
	    aa->mis[i1] = false;
	}
	
	na++;
	print (1, &n1, &center, &na, aa);
	
	fixbad (niter, aa, n1, data);
	
	sf_floatwrite(data,n1,out);
	na = tempna;
    }

    exit(0);
}




