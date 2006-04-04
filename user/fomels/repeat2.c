/* A clone of repeat */
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

#include "repeat2.h"

static int n1, n2;
static sf_operator oper;

void repeat2_init(int m1            /* trace length */, 
		  int m2            /* number of traces */, 
		  sf_operator oper1 /* operator */)
/*< initialize >*/
{
    n1 = m1;
    n2 = m2;
    oper = oper1;
}

void repeat2_lop (bool adj, bool add, int nx, int ny, float *xx, float *yy)
/*< combined linear operator >*/
{
    int i2;       
    
    if (nx != ny || nx != n1*n2) sf_error("%s: Wrong size (nx=%d ny=%d n1=%d n2=%d)",__FILE__,nx,ny,n1,n2);

    sf_adjnull (adj, add, nx, ny, xx, yy);

    for (i2=0; i2 < n2; i2++) {
	oper(adj,true,n1,n1,xx+i2*n1,yy+i2*n1);
    }
}

/* 	$Id: repeat.c 1571 2005-11-21 05:50:10Z fomels $	 */

