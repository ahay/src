/* Division by downgoing wave operators */
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

static int nn3, nn1, nn2;
float* dowav;


void divk_init (float* ff       /* wavelet [na] */,
                 int nk1         /* number of components */,
                 int n2          /* size of one component 1D*/,
                 int n1          /* size of one component 2D*/) 
/*< initialize with a wavelet>*/
{
    
    dowav  = ff;
    nn3   = nk1;
    nn1  = n1;
    nn2  = n2;
    
    
}

void divk_lop (bool adj, bool add, int nx, int ny, float* xx, float* yy) 
/*< 
linear operator
convolution with known wavelet:
yy = wav*xx
xx = wav'*yy
xx and yy are 2D or more
wav is 1D
>*/
{
    int i, n;
    n = nn1*nn2*nn3;

    if (nx != n || ny != n) sf_error("%s: wrong size",__FILE__);

    sf_adjnull (adj, add, nx, ny, xx, yy);

    for (i=0; i < n; i++) {
         if (adj) {
		xx[i] += dowav[i]*yy[i];
	    } else {
		yy[i] += dowav[i]*xx[i];
	    }

       
        
    } 
}












