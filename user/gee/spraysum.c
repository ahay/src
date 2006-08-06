/* Spray or sum over 1, 2, and/or 3-axis. */
/*
  Copyright (C) 2006 University of Texas at Austin
  
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

static int n1,n2,n3, m1,m2,m3;

void spraysum_init(int n1_in, int n2_in, int n3_in, 
		   int m1_in, int m2_in, int m3_in)
/*< initialize >*/
{
    n1 = n1_in; n2 = n2_in; n3 = n3_in;
    m1 = m1_in; m2 = m2_in; m3 = m3_in;

    if( n1 != 1  &&  m1 != 1  &&  n1 != m1) sf_error("%s: n1,m1",__FILE__);
    if( n2 != 1  &&  m2 != 1  &&  n2 != m2) sf_error("%s: n2,m2",__FILE__);
    if( n3 != 1  &&  m3 != 1  &&  n3 != m3) sf_error("%s: n3,m3",__FILE__);
}

void spraysum_lop(bool adj, bool add, int nx, int ny, float *xx, float *yy)
/*< linear operator >*/
{ 
    int i1,i2,i3,  x,y;

    sf_adjnull(adj,add,nx,ny,xx,yy);
    
    for (i3=0; i3 < SF_MAX(n3,m3); i3++) {  
	x = SF_MIN(i3,n3-1); 
	y = SF_MIN(i3,m3-1);
	for (i2=0; i2 < SF_MAX(n2,m2); i2++) {   
	    x = x*n2 + SF_MIN(i2,n2-1);   
	    y = y*m2 + SF_MIN(i2,m2-1); 
	    for (i1=0; i1 < SF_MAX(n1,m1); i1++) {   
		x = x*n1 + SF_MIN(i1,n1-1);   
		y = y*m1 + SF_MIN(i1,m1-1); 

		if( adj)  xx[x]  +=  yy[y];
		else      yy[y]  +=  xx[x];
	    }
	}
    }
}
