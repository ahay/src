/* Digital freqlet frame */
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

#include "freqlets.h"
#include "freqlet.h"

static int n, nw;
static float d, *w;
static sf_complex *z;

 
void freqlets_init(int n1 /* data size */, 
		   float d1 /* sampling */,
		   bool inv, bool unit, char type,
		   int nw1 /* number of frequencies */, 
		   float *w1      /* [nw] frequencies */,
		   sf_complex *z1 /* [nw] Z factors */)
/*< allocate space >*/
{
    n = n1;
    d = d1 * 2 * SF_PI;
    nw = nw1;
    w = w1;
    z = z1;
    freqlet_init(n,inv,unit,type);
}

void freqlets_close(void)
/*< deallocate space >*/
{
    freqlet_close();
}

void freqlets_set(float *w1, sf_complex *z1)
/*< set frequency >*/
{
    w = w1;
    z = z1;
}

void freqlets_lop(bool adj, bool add, int nx, int ny, 
		  sf_complex *x, sf_complex *y)
/*< linear operator >*/
{
    int iw;

    if (n*nw != nx && ny != n) sf_error ("%s: wrong size",__FILE__);

    sf_cadjnull(adj,add,nx,ny,x,y);

    /* loop over frequencies */
    for (iw=0; iw < nw; iw++) {
	if (NULL != w) {
	    freqlet_set(w[iw]*d);
	} else {
	    freqlet_setz(z[iw]);
	}
	freqlet_lop(adj,true,n,n,x+iw*n,y);
    }
}
