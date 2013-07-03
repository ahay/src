/* Transient convolution, adjoint is the input */
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

#include "tcaih.h"

static int nb, nt, nm, nh;
static const float* bb;

void tcaih_init (int na          /* filter size */, 
		 const float* aa /* filter [na] */,
                 int n1, int n2, int n3 /* data dimension */)
/*< initialize >*/
{
    nb = na;
    bb = aa;
    nt = n1;
    nm = n2;
    nh = n3;
}

void tcaih_lop (bool adj, bool add, int nx, int ny, float* xx, float* yy) 
/*< linear operator >*/
{
    int b, x, y;
    int i1,i2;

    if(ny < nx+nb-1) sf_error("%s: size problem: %d < %d+%d-1",
			      __FILE__,ny,nx,nb);
    sf_adjnull (adj, add, nx, ny, xx, yy);
    
    for (i2=0; i2<nm; i2++) {
    for (i1=0; i1<nt; i1++) {
    for( b=0; b < nb; b++) {
	for( x=0; x < nh; x++) { y = x + b;
	    if( adj) xx[x*nm*nt+i2*nt+i1] += yy[y*nm*nt+i2*nt+i1] * bb[b];
	    else     yy[y*nm*nt+i2*nt+i1] += xx[x*nm*nt+i2*nt+i1] * bb[b];
	}
    }
    }
    }
}

/* 	$Id: tcai1.c 2139 2006-08-09 21:24:23Z sfomel $	 */
