/* Operations with complex sinusoids */
/*
  Copyright (C) 2008 University of Texas at Austin
  
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

#include "sin.h"

static sf_complex z0;
static int n2, k2;

void sin_init(sf_complex z1) 
/*< initialize >*/
{
    z0 = conjf(z1);
    
}

void sinpred_init(sf_complex z1, 
		  int n /* data size */,
		  int k /* radius */)
/*< initialize >*/
{
    z0 = conjf(z1);
    n2 = n;
    k2 = k;
    if (k2 > n2-1) sf_error("%s: k2=%d > n2-1=%d",__FILE__,k2,n2-1);
}

void sin_destruct(bool adj, bool add, int nx, int ny, 
		  sf_complex *xx, sf_complex *yy)
/*< destruction operator >*/
{
    int i;

    if (nx != ny) sf_error("%s: wrong dimensions %d != %d",__FILE__,nx,ny);

    sf_ccopy_lop(adj,add,nx,nx,xx,yy);

    for (i = 1; i < nx; i++) {
	if(adj) {
#ifdef SF_HAS_COMPLEX_H	 
	    xx[i-1] -= yy[i] * conjf(z0);
#else
	    xx[i-1] = sf_cadd(xx[i-1],sf_cmul(yy[i],sf_cneg(sf_conjf(z0))));
#endif
	} else {
#ifdef SF_HAS_COMPLEX_H	 
	    yy[i] -= xx[i-1] * z0;
#else
	    yy[i] = sf_cadd(yy[i],sf_cmul(xx[i-1],sf_cneg(z0)));
#endif
	}
    }
}

void sin_construct(bool adj, bool add, int nx, int ny, 
		   sf_complex *xx, sf_complex *yy)
/*< construction operator >*/
{
    int i;
    sf_complex t;

    if (nx != ny) sf_error("%s: wrong dimensions %d != %d",__FILE__,nx,ny);

    sf_cadjnull(adj, add, nx, nx, xx, yy);

    t = sf_cmplx(0.0,0.0);
    if(adj) {	
	for (i = nx-1; i >= 0; i--) {
#ifdef SF_HAS_COMPLEX_H
	    t = yy[i] + t * conjf(z0);
	    xx[i] += t;  
#else
	    t = sf_cadd(yy[i],sf_cmul(t,sf_conjf(z0)));
	    xx[i] = sf_cadd(xx[i],t);
#endif
	}
    } else {
	for (i = 0; i < nx; i++) {
#ifdef SF_HAS_COMPLEX_H	    
	    t = xx[i] + t * z0;
	    yy[i] += t;
#else
	    t = sf_cadd(xx[i],sf_cmul(t,z0));
	    yy[i] = sf_cadd(yy[i],t);
#endif
	}
    }
}

void sinpredicter_lop(bool adj, bool add, int nx, int ny, 
		      sf_complex *xx, sf_complex *yy)
/*< linear operator >*/
{
    int i2;
    sf_complex t;

    if (nx != ny || nx != n2+2*k2) 
	sf_error("%s: Wrong dimensions",__FILE__);

    sf_cadjnull(adj,add,nx,ny,xx,yy);

    t = sf_cmplx(0.0,0.0);
    if (adj) {
	for (i2=n2+2*k2-1; i2 >= 0; i2--) {
#ifdef SF_HAS_COMPLEX_H
	    t = yy[i2] + t * conjf(z0);
	    xx[i2] += t;  
#else
	    t = sf_cadd(yy[i2],sf_cmul(t,sf_conjf(z0)));
	    xx[i2] = sf_cadd(xx[i2],t);
#endif
	}
    } else {
	for (i2=0; i2 < n2+2*k2; i2++) {
#ifdef SF_HAS_COMPLEX_H	    
	    t = xx[i2] + t * z0;
	    yy[i2] += t;
#else
	    t = sf_cadd(xx[i2],sf_cmul(t,z0));
	    yy[i2] = sf_cadd(yy[i2],t);
#endif
	}
    }
}

void sinsubtracter_lop(bool adj, bool add, int nx, int ny, 
		       sf_complex *xx, sf_complex *yy)
/*< linear operator >*/
{
    int i2, j2, m2;
    sf_complex c;

    if (nx != ny || nx != n2+2*k2) 
	sf_error("%s: Wrong dimensions",__FILE__);

    sf_cadjnull(adj,add,nx,ny,xx,yy);

    if (adj) {
	for (j2=0; j2 < n2+2*k2; j2++) {
	    i2=j2+k2;
	    if (i2 < n2+2*k2) {
		c = sf_cmplx(1.,0.);
		for (m2=i2-1; m2 >= j2; m2--) {
#ifdef SF_HAS_COMPLEX_H	 
		    c *= conjf(z0);
#else
		    c = sf_cmul(c,sf_conjf(z0));
#endif
		}
#ifdef SF_HAS_COMPLEX_H	 
		xx[j2] += yy[j2] - yy[i2] * c;
#else 
		xx[j2] = sf_cadd(xx[j2],
				 sf_cadd(yy[j2],
					 sf_cneg(sf_cmul(yy[i2],c))));
#endif
	    } else {
#ifdef SF_HAS_COMPLEX_H
		xx[j2] += yy[j2];
#else	
		xx[j2] = sf_cadd(xx[j2],yy[j2]);
#endif
	    }
	}
    } else {
	for (i2=0; i2 < n2+2*k2; i2++) { 
	    j2=i2-k2;
	    if (j2 >=0) {
		c = sf_cmplx(1.,0.);
		for (m2=j2; m2 < i2; m2++) {
#ifdef SF_HAS_COMPLEX_H	 
		    c *= z0;
#else
		    c = sf_cmul(c,z0);
#endif
		}
#ifdef SF_HAS_COMPLEX_H	 
		yy[i2] += xx[i2] - xx[j2] * c;
#else 
		yy[i2] = sf_cadd(yy[i2],
				 sf_cadd(xx[i2],
					 sf_cneg(sf_cmul(xx[j2],c))));
#endif
	    } else {
#ifdef SF_HAS_COMPLEX_H	 
		yy[i2] += xx[i2];
#else 
		yy[i2] = sf_cadd(yy[i2],xx[i2]);
#endif
	    }
	}
    }
}
