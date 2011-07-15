/* Phase unwrapping. */
/*
  Copyright (C) 2011 University of Texas at Austin
  
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

void grad2init(int n1, int n2, 
	       sf_complex **z, float ***rt) 
/*< initialize phase gradient >*/
{
    int  i1, i2;
    sf_complex a, a0, a1;

    for (i2=0; i2 < n2; i2++) {
	rt[0][i2][n1-1] = 0.;
	rt[1][i2][n1-1] = 0.;
    }

    for (i1=0; i1 < n1; i1++) {
	rt[0][n2-1][i1] = 0.;
	rt[1][n2-1][i1] = 0.;
    }

    for (i2=0; i2 < n2-1; i2++) {
	for (i1=0; i1 < n1-1; i1++) {
	    a = z[i2][i1];
	    a0 = z[i2][i1+1];
	    a1 = z[i2+1][i1];
#ifdef SF_HAS_COMPLEX_H	    
	    rt[0][i2][i1] = cimagf(clogf(conjf(a0)*a));
	    rt[1][i2][i1] = cimagf(clogf(conjf(a1)*a));
#else
	    rt[0][i2][i1] = cimagf(clogf(sf_cmul(conjf(a0),a)));
	    rt[1][i2][i1] = cimagf(clogf(sf_cmul(conjf(a1),a)));
#endif
        }
    }
}

void unwraper(int n1, int n2, /* dimensions */
	      sf_complex **zz /* input data */, 
	      float *hh       /* phase */, 
	      int niter       /* number of iterations */) 
/*< Phase unwraper.   Starting from phase hh, improve it. >*/
{
    float ***rt;

    rt = sf_floatalloc3(n1,n2,2);
    grad2init(n1,n2, zz,rt);

    sf_igrad2_init(n1,n2);
    sf_tinysolver(sf_igrad2_lop,sf_cgstep,
		  n1*n2,n1*n2*2,hh,NULL,rt[0][0],niter);
    sf_cgstep_close();

    free(**rt);
    free(*rt);
    free(rt);
}

