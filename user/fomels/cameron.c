/* Weighted gradient for time to depth conversion */
/*
  Copyright (C) 2010 University of Texas at Austin
  
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

static int n1, n2, n12;
static float d1, d2, *sqv;

void cameron_init(int m1, int m2,     /* dimensions */ 
		  float c1, float c2, /* steps */
		  float* v            /* Dix velocity [n1*n2] */)
/*< initialize >*/
{
    int i;

    n1=m1;
    n2=m2;
    n12=m1*m2;
    d1=1.0/c1;
    d2=0.5/c2;
    
    sqv = sf_floatalloc(n12);

    for (i=0; i < n12; i++) {
	sqv[i] = sqrtf(v[i]);
    }
}

void cameron_close(void)
/*< free allocated storage >*/
{
    free(sqv);
}


void cameron_lop (bool adj, bool add, int np, int nr, float* p, float* r)
/*< linear operator, r[n1*n2*2] is the weighted gradient of p[n1*n2] >*/
{
    int i1,i2,i;

    if (np != n12) sf_error("%s: %d != %d",__FILE__,np,n12);
    if (nr != 2*n12) sf_error("%s: %d != 2*%d",__FILE__,nr,n12);

    sf_adjnull (adj,add,np,nr,p,r);

    for (i2=1; i2 < n2-1; i2++) {  
	for (i1=1; i1 < n1; i1++) {
	    i = i1+i2*n1;
	    if (adj == true) {
		p[i]    += d1*r[i]/sqv[i-1];
		p[i-1]  -= d1*r[i]/sqv[i-1];
		p[i+n1] += d2*sqv[i]*r[i+n12];
		p[i-n1] -= d2*sqv[i]*r[i+n12];
	    } else {
		r[i]     += d1*(p[i]  - p[i-1])/sqv[i-1]; 
		r[i+n12] += d2*sqv[i]*(p[i+n1] - p[i-n1]);
	    }
	}
    }
}
