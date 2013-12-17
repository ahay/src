/* Nonstationary convolution operator */
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
/*^*/
#include "nfconv.h"

static sf_complex* filter;
static int n1, na2;

void nfconv_init (sf_complex* bb, int nn1, int naa2) 
/*< initialize with a pointer to a matrix >*/
{

    filter = bb;
    n1     = nn1;
    na2    = naa2;
}

void nfconva_lop (bool adj, bool add, 
		  int nx, int ny, sf_complex *xx, sf_complex *yy) 
/*< linear operator nonstationary convolution >*/
{    
    int itau, it,lag;
   
    lag=(na2+1)/2;
    sf_cadjnull (adj, add, nx, ny, xx, yy);

    /* nonstationary convolution*/
	for (itau=0; itau < na2; itau++) {
	    
	    for (it = na2-lag; it <= n1-lag; it++) {
		    if (adj) {
#ifdef SF_HAS_COMPLEX_H
			    xx[it-itau+lag-1] += yy[it]*conjf(filter[it+n1*itau]);
#else
		        xx[it-itau+lag-1] = sfcadd(xx[it-itau+lag-1],sf_cmul(yy[it],conjf(filter[it+n1*itau])));
#endif
		    } else {
#ifdef SF_HAS_COMPLEX_H
		        yy[it] += xx[it-itau+lag-1]*filter[it+n1*itau];
#else  
		        yy[it] = sfcadd(yy[it],sf_cmul(xx[it-itau+lag-1],filter[it+n1*itau]));
#endif
		    } 
	    }
	}
}

void nfconvf_lop (bool adj, bool add, 
		  int nx, int ny, sf_complex *xx, sf_complex *yy) 
/*< linear operator nonstationary convolution >*/
{    
    int itau, it,lag;
    
    lag=na2;
    sf_cadjnull (adj, add, nx, ny, xx, yy);

    /* nonstationary convolution*/
	for (itau=0; itau < na2; itau++) {
	    
	    for (it = na2; it < n1; it++) {
		    if (adj) {
#ifdef SF_HAS_COMPLEX_H
			    xx[it-na2+itau] += yy[it]*conjf(filter[it+n1*itau]);
#else
		        xx[it-na2+itau] = sfcadd(xx[it-na2+itau],sf_cmul(yy[it],conjf(filter[it+n1*itau])));
#endif
		    } else {
#ifdef SF_HAS_COMPLEX_H
		        yy[it] += xx[it-na2+itau]*filter[it+n1*itau];
#else  
		        yy[it] = sfcadd(yy[it],sf_cmul(xx[it-na2+itau],filter[it+n1*itau]));
#endif
		    } 
	    }
	}
}

void nfconvb_lop (bool adj, bool add, 
		  int nx, int ny, sf_complex *xx, sf_complex *yy) 
/*< linear operator nonstationary convolution >*/
{    
    int itau, it,lag;
   
    lag=na2;
    sf_cadjnull (adj, add, nx, ny, xx, yy);

    /* nonstationary convolution*/
	for (itau=0; itau < na2; itau++) {
	    
	    for (it = 0; it < n1-itau; it++) {
		    if (adj) {
#ifdef SF_HAS_COMPLEX_H
			    xx[it+itau] += yy[it]*conjf(filter[it+n1*itau]);
#else
		        xx[it+itau] = sfcadd(xx[it+itau],sf_cmul(yy[it],conjf(filter[it+n1*itau])));
#endif
		    } else {
#ifdef SF_HAS_COMPLEX_H
		        yy[it] += xx[it+itau]*filter[it+n1*itau];
#else  
		        yy[it] = sfcadd(yy[it],sf_cmul(xx[it+itau],filter[it+n1*itau]));
#endif
		    } 
	    }
	}
}

void nfconv_close () 
/*< free filter memory >*/
{
    free (filter);
}

/* 	$Id: mmmult.c 7107 2011-04-10 02:04:14Z ivlad $	 */
