/* 1-D internal convolution complex convolution, adjoint is input */
/*
  Copyright (C) 2010 Politecnico di Milano
  
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

#include "cicai1.h"


static int nb, lag;
static sf_complex* bb;

void cicai1_init(int nb1           /* filter size */, 
		 sf_complex* bb1  /* filter taps [nb] */,
		 int lag1  /* filter lag (lag=1 is causal) */) 
/*< initialize >*/
{

    nb  = nb1;
    lag = lag1;
    bb = bb1;
}

void cicai1_lop(bool adj, bool add, int nx, int ny, sf_complex* xx, sf_complex* yy)
/*< linear operator >*/
{
    		
    int x, b, y;
    if(ny != nx) sf_error("%s: size problem: %d != %d",__FILE__,ny,nx);

    sf_cadjnull(adj, add, nx, ny, xx, yy);
    /* sf_warning("############cicaf############"); */
    for (b=0; b < nb; b++) {
	for (y=nb; y <= ny; y++) { 

	    x = y - b - 1;	
#ifdef SF_HAS_COMPLEX_H
	    if( adj) {		/* CONJUGATE TRANSPOSE OPERATOR*/		
				
		xx[x] = xx[x] + (yy[y-lag] * conjf(bb[b]) );

		/* sf_warning("bb[b]=y[%d]-x[%d]=",y-lag,x); cprint(bb[b]); */
	    }
	    else     { 				

		yy[y-lag]  =yy[y-lag]+( bb[b] * xx[x]);
		/* sf_warning("y[%d]-x[%d]=",y-lag,x); cprint(yy[y]); cprint(xx[x]); cprint(bb[b]); */
	    }
#else
	    if( adj) {		/* CONJUGATE TRANSPOSE OPERATOR*/
		xx[x] = sf_cadd(xx[x], sf_cmul(yy[y-lag], sf_conjf(bb[b])) ); 
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

void cicai1_close(void)
/*<close>*/
{
}
