/* 2-D finite-difference gradient */
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

#include "_bool.h"
/*^*/

#include "igrad2.h" 

#include "error.h"
#include "adjnull.h"

static int n1, n2, n12; 

void sf_igrad2_init (int n1_in, int n2_in)
/*< initialize with data dimensions >*/
{
    n1 = n1_in; 
    n2 = n2_in;
    n12 = n1*n2;
}

void sf_igrad2_lop (bool adj, bool add, int np, int nr, float* p, float* r)
/*< linear operator, r[n1*n2*2] is the gradient of p[n1*n2] >*/
{
    int i1,i2,i;

    if (np != n12) sf_error("%s: %d != %d",__FILE__,np,n12);
    if (nr != 2*n12) sf_error("%s: %d != 2*%d",__FILE__,nr,n12);

    sf_adjnull (adj,add,np,nr,p,r);

    for (i2=0; i2 < n2-1; i2++) {  
	for (i1=0; i1 < n1-1; i1++) {
	    i = i1+i2*n1;
	    if (adj) {
		p[i+1]  += r[i]; 
		p[i+n1] += r[i+n12];
		p[i]    -= (r[i] + r[i+n12]);
	    } else {
		r[i]     += (p[i+1]  - p[i]); 
		r[i+n12] += (p[i+n1] - p[i]);
	    }
	}
    }
}
