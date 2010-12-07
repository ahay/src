/* Weighting with multiple components and complex numbers */
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

#include "cweight2.h"

static int nw;
static sf_complex **w;

void cweight2_init(int nw1        /* number of components */, 
		   int n          /* model size */, 
		   sf_complex *ww /* weight [nw*n] */)
/*< initialize >*/
{
    int iw;

    nw = nw1;
    w = (sf_complex**) sf_alloc(nw,sizeof(sf_complex*));

    for (iw=0; iw < nw; iw++) {
	w[iw] = ww+iw*n;
    }
}

void cweight2_close(void)
/*< free allocated storage >*/
{
    free(w);
}

void cweight2_lop (bool adj, bool add, int nx, int ny, sf_complex* xx, sf_complex* yy)
/*< linear operator >*/
{
    int i, iw;

    if (nw*ny != nx) sf_error("%s: size mismatch: %d*%d != %d",
			      __FILE__,nw,ny,nx);

    sf_cadjnull (adj, add, nx, ny, xx, yy);
  
    for (iw=0; iw < nw; iw++) {
	for (i=0; i < ny; i++) {
	    if (adj) {
#ifdef SF_HAS_COMPLEX_H
		xx[i+iw*ny] += yy[i] * conjf(w[iw][i]);
#else
		xx[i+iw*ny] = sf_cadd(xx[i+iw*ny],
				      sf_cmul(yy[i],conjf(w[iw][i])));
#endif
	    } else {
#ifdef SF_HAS_COMPLEX_H
		yy[i] += xx[i+iw*ny] * w[iw][i];
#else
		yy[i] = sf_cadd(yy[i],sf_cmul(xx[i+iw*ny],w[iw][i]));
#endif
	    }
	}
    }
}

