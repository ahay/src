/* Freqlet frame plus interpolation */
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

#include "freqint.h"
#include "freqlets.h"
 
static int nt;
static sf_complex *t;

void freqint_init (int nd /* irregular data size */,
		   float *coord /* [nd] irregular coordinates */,
		   int n /* regular data size */, 
		   float d /* sampling */,
		   float o /* origin */,
		   sf_interpolator interp     /* interpolation function */, 
		   int nf                  /* interpolator length */,
		   bool inv, bool unit, char type,
		   int nw /* number of frequencies */, 
		   float *w /* [nw] frequencies */,
		   sf_complex *z /* [nw] Z factors */)
/*< allocate space >*/
{
    nt = n;
    t = sf_complexalloc(nt);

    sf_int1_init (coord,o,d,n,interp,nf,nd, 0.0);
    freqlets_init(n,d,inv,unit,type,nw,w,z);
}

void freqint_close(void)
/*< deallocate space >*/
{
    free(t);
    freqlets_close();
}

void freqint_lop(bool adj, bool add, int nx, int ny, 
		 sf_complex *x, sf_complex *y)
/*< linear operator >*/
{
    sf_cadjnull(adj,add,nx,ny,x,y);
    sf_cchain(sf_cint1_lop,freqlets_lop,adj,true,nx,ny,nt,x,y,t);
}
