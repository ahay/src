/* Simple data extrapolation. */
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

#include "extend.h"

static const int nw = 3;
static const float a[] = {7./3., -5./3., 1./3.};

void extend (int ne     /* padding */, 
	     int nd     /* data length */, 
	     float *dat /* data [nd] */, 
	     float *ext /* extension [nd+2*ne] */)
/*< 1-D extension >*/
{
    int i, j;
    float s;

    for (i=0; i < nd; i++) {
	ext[ne+i] = dat[i];
    }
    for (i=ne-1; i >= 0; i--) {
	for (s=0., j=0; j < nw; j++) {
	    s += a[j]*ext[i+j+1];
	}
	ext[i] = s;
    }
    for (i=nd+ne; i < nd+2*ne; i++) {
	for (s=0., j=0; j < nw; j++) {
	    s += a[j]*ext[i-j-1];
	}
	ext[i] = s;
    }
}

void extend2 (int ne         /* padding */, 
	      int n1, int n2 /* data size */, 
	      float** dat    /* data [n2][n1] */, 
	      float** ext    /* extension [n2+2*ne][n1+2*ne] */, 
	      float* tmp1    /* temporary storage [n2] */, 
	      float* tmp2    /* temporary storage [n2+2*ne] */)
/*< 2-D extension >*/
{
    int i1, i2;
    for (i2=0; i2 < n2; i2++) {
	extend (ne,n1,dat[i2],ext[i2+ne]);
    }
    for (i1=0; i1 < n1+2*ne; i1++) {
	for (i2=0; i2 < n2; i2++) {
	    tmp1[i2] = ext[i2+ne][i1];
	}
	extend (ne,n2,tmp1,tmp2);
	for (i2=0; i2 < n2+2*ne; i2++) {
	    ext[i2][i1] = tmp2[i2];
	} 
    }
}

/* 	$Id: extend.c 7107 2011-04-10 02:04:14Z ivlad $	 */
