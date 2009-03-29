/* Row of identity operators */
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

static int nn3, nn1, nn2, nwav, ori_wav;
float* wav;


void convk_init (int na          /* wavelet size */, 
		 float* ff       /* wavelet [na] */,
                 int ori_wavelet     /*origin coordinate of wavelet*/,
                 int nk1         /* number of components */,
                 int n2          /* size of one component 1D*/,
                 int n1          /* size of one component 2D*/) 
/*< initialize with a wavelet>*/
{
    nwav = na;
    wav  = ff;
    nn3   = nk1;
    nn1  = n1;
    nn2  = n2;
    ori_wav = ori_wavelet;
    
}

void convk_lop (bool adj, bool add, int nx, int ny, float* xx, float* yy) 
/*< 
linear operator
convolution with known wavelet:
yy = wav*xx
xx = wav'*yy
xx and yy are 2D or more
wav is 1D
>*/
{
    int ik, in1, in2, iwav, n, y;
    n = nn1*nn2;

    if (nx != n*nn3 || ny != n*nn3) sf_error("%s: wrong size",__FILE__);

    sf_adjnull (adj, add, nx, ny, xx, yy);

    for (ik=0; ik < nn3; ik++) {
        for (in2=0; in2 < nn2; in2++) {
            for (iwav=0; iwav < nwav; iwav++){
	    for (in1=0; in1 < nn1; in1++) {
                y = iwav+in1+ori_wav;
                if (y < 0 || y >=nn1) continue; /*zero value boundary condition*/

	        if (adj) {
		    xx[ik*n+in2*nn1+in1] += yy[ik*n+in2*nn1+y] * wav[iwav];
	        } else {
		    yy[ik*n+in2*nn1+y] += xx[ik*n+in2*nn1+in1] * wav[iwav];
	        }


	    }
            } /*1d*/
        } /*2d*/
    } /*3d*/
}











