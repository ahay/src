/* Normal operator for CR decomposition */
/*
  Copyright (C) 2015 University of Texas at Austin
  
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

static int nr, nc;

void cr_init(int nr1 /* number of rows */, 
	     int nc1 /* number of columns */)
/*< initialize >*/
{
    nr = nr1;
    nc = nc1;
}

void cr_prec(int nx, const float *x, float *y)
/*< apply preconditioning >*/
{
    int i;
    for (i=0; i < nr; i++) {
	y[i] = x[i]/nc;
    }
    for (i=nr; i < nr+nc; i++) {
	y[i] = x[i]/nr;
    }
}

void cr_apply(int nx, const float *x, float *y)
/*< apply normal operator >*/
{
    int i;
    float sr, sc;

    if (nx != nr+nc) sf_error("%s: wrong size nx != %d",__FILE__,nr+nc);

    /* stack rows and columns */
    sr = 0.0f;
    for (i=0; i < nr; i++) {
	sr += x[i];
    }
    sc = 0.0f;
    for (i=nr; i < nx; i++) {
	sc += x[i];
    }
    
    for (i=0; i < nr; i++) {
	y[i] = nc*x[i] + sc;
    }
    for (i=nr; i < nx; i++) {
	y[i] = nr*x[i] + sr;
    }
}

