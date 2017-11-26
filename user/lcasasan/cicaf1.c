/* 1-D internal convolution complex convolution, adjoint is filter */
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
/*^*/

#include "cicaf1.h"


static int nx, lag;
static sf_complex* xx;

void cicaf1_init(int ny           /* data size */, 
		 sf_complex* xx1 /* data [ny] */,
		 int lag1  /* filter lag (lag=1 is causal) */) 
/*< initialize >*/
{

    nx  = ny;
    lag = lag1;
    xx = xx1;
}

void cicaf1_lop(bool adj, bool add, int nb, int ny, sf_complex* bb, sf_complex* yy)
/*< linear operator >*/
{
    		

    int x, b, y;
    if(ny != nx) sf_error("%s: size problem: %d != %d",__FILE__,ny,nx);

    sf_cadjnull(adj, add, nb, ny, bb, yy);
    /* sf_warning("############cicaf############"); */
    for (b=0; b < nb; b++) {
	for (y=nb; y <= ny; y++) { 

	    x = y - b - 1;	
#ifdef SF_HAS_COMPLEX_H
	    if( adj) {		/* CONJUGATE TRANSPOSE OPERATOR*/		
				
				
		bb[b] = bb[b]+ (yy[y-lag] * conjf(xx[x]) );

		/* sf_warning("bb[b]=y[%d]-x[%d]=",y-lag,x); cprint(bb[b]); */
	    }
	    else     { 				

		yy[y-lag]  =yy[y-lag]+( bb[b] * xx[x]);
		/* sf_warning("y[%d]-x[%d]=",y-lag,x); cprint(yy[y]); cprint(xx[x]); cprint(bb[b]); */
	    }
#else
	    if( adj) {		/* CONJUGATE TRANSPOSE OPERATOR*/
		bb[b] = sf_cadd(bb[b], sf_cmul(yy[y-lag], sf_conjf(xx[x])) ); 
		/* sf_warning("bb[b]=y[%d]-x[%d]=",y-lag,x); cprint(bb[b]); */
	    }
	    else     { 				

		yy[y-lag] = sf_cadd(yy[y-lag], sf_cmul(bb[b], xx[x]) ); 
		/* sf_warning("y[%d]-x[%d]=",y-lag,x); cprint(yy[y]); */
	    }
#endif
	}
    }
}

void cicaf1_close(void)
/*<close>*/
{
}
