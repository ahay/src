/* Smooth gradient operations */
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

#include "edge.h"

void sf_grad2 (int n          /* data size */, 
	       const float *x /* input trace [n] */, 
	       float *w       /* output gradient squared [n] */)
/*< centered finite-difference gradient >*/
{
    int i;
    float ww;

    w[0] = 0.;
    for (i=1; i < n-1; i++) {
	ww = 0.5*(x[i+1]-x[i-1]);
	w[i] = ww*ww;
    }
    w[n-1] = 0.;
}

void sf_sobel (int n1, int n2         /* data size */, 
	       float **x              /* input data [n2][n1] */, 
	       float **w1, float **w2 /* output gradient [n2][n1] */)
/*< Sobel's 9-point gradient >*/
{
    int i1, i2;

    for (i2=0; i2 < n2; i2++) {
	for (i1=0; i1 < n1; i1++) {
	    if (i2 == 0 || i2 == n2-1 || i1 == 0 || i1 == n1-1) {
		w1[i2][i1] = 0.;
		w2[i2][i1] = 0.;
	    } else {
		w1[i2][i1] =
		    x[i2+1][i1+1] - x[i2+1][i1-1] +
		    2*(x[i2][i1+1] - x[i2][i1-1]) +
		    x[i2-1][i1+1] - x[i2-1][i1-1];
		w2[i2][i1] =
		    x[i2+1][i1+1] - x[i2-1][i1+1] +
		    2*(x[i2+1][i1] - x[i2-1][i1]) +
		    x[i2+1][i1-1] - x[i2-1][i1-1]; 
	    }
	}
    }
}

void sf_sobel2 (int n1, int n2  /* data size */, 
		float **x        /* input data [n2][n1] */, 
		float **w        /* output gradient squared [n2][n1] */)
/*< Sobel's gradient squared >*/
{
    int i1, i2;
    float w1, w2;

    for (i2=0; i2 < n2; i2++) {
	for (i1=0; i1 < n1; i1++) {
	    if (i2 == 0 || i2 == n2-1 || i1 == 0 || i1 == n1-1) {
		w[i2][i1] = 0.;
	    } else {
		w1 =
		    x[i2+1][i1+1] - x[i2+1][i1-1] +
		    2*(x[i2][i1+1] - x[i2][i1-1]) +
		    x[i2-1][i1+1] - x[i2-1][i1-1];
		w2 =
		    x[i2+1][i1+1] - x[i2-1][i1+1] +
		    2*(x[i2+1][i1] - x[i2-1][i1]) +
		    x[i2+1][i1-1] - x[i2-1][i1-1];
		w[i2][i1] = w1*w1 + w2*w2;
	    }
	}
    }
}

void sf_sobel32 (int n1, int n2, int n3  /* data size */, 
		 float ***x              /* input data [n3][n2][n1] */, 
		 float ***w              /* output gradient squared */)
/*< Sobel's gradient squared in 3-D>*/
{
    int i1, i2, i3;
    float w1, w2, w3;

    for (i3=0; i3 < n2; i3++) {
	for (i2=0; i2 < n2; i2++) {
	    for (i1=0; i1 < n1; i1++) {
		if (i3 == 0 || i3 == n3-1 || 
		    i2 == 0 || i2 == n2-1 || 
		    i1 == 0 || i1 == n1-1) {
		    w[i3][i2][i1] = 0.;
		} else {
		    w1 =
			x[i3-1][i2-1][i1+1] - x[i3-1][i2-1][i1-1] +
			x[i3+1][i2-1][i1+1] - x[i3+1][i2-1][i1-1] +
			x[i3-1][i2+1][i1+1] - x[i3-1][i2+1][i1-1] +
			x[i3+1][i2+1][i1+1] - x[i3+1][i2+1][i1-1] +
			4*(x[i3-1][i2][i1+1] - x[i3-1][i2][i1-1] +
			   x[i3+1][i2][i1+1] - x[i3+1][i2][i1-1] +
			   x[i3][i2-1][i1+1] - x[i3][i2-1][i1-1] +
			   x[i3][i2+1][i1+1] - x[i3][i2+1][i1-1] +
			   x[i3][i2][i1+1] - x[i3][i2][i1-1]);
		    w2 =
			x[i3-1][i2+1][i1-1] - x[i3-1][i2-1][i1-1] +
			x[i3+1][i2+1][i1-1] - x[i3+1][i2-1][i1-1] +
			x[i3-1][i2+1][i1+1] - x[i3-1][i2-1][i1+1] +
			x[i3+1][i2+1][i1+1] - x[i3+1][i2-1][i1+1] +
			4*(x[i3-1][i2+1][i1] - x[i3-1][i2-1][i1] +
			   x[i3+1][i2+1][i1] - x[i3+1][i2-1][i1] +
			   x[i3][i2+1][i1-1] - x[i3][i2-1][i1-1] +
			   x[i3][i2+1][i1+1] - x[i3][i2-1][i1+1] +
			   x[i3][i2+1][i1] - x[i3][i2-1][i1]);
		    w3 =
			x[i3+1][i2-1][i1-1] - x[i3-1][i2-1][i1-1] +
			x[i3+1][i2+1][i1-1] - x[i3-1][i2+1][i1-1] +
			x[i3+1][i2-1][i1+1] - x[i3-1][i2-1][i1+1] +
			x[i3+1][i2+1][i1+1] - x[i3-1][i2+1][i1+1] +
			4*(x[i3+1][i2-1][i1] - x[i3-1][i2-1][i1] +
			   x[i3+1][i2+1][i1] - x[i3-1][i2+1][i1] +
			   x[i3+1][i2][i1-1] - x[i3-1][i2][i1-1] +
			   x[i3+1][i2][i1+1] - x[i3-1][i2][i1+1] +
			   x[i3+1][i2][i1] - x[i3-1][i2][i1]);		    
		    w[i3][i2][i1] = (w1*w1 + w2*w2 + w3*w3)/36.;
		}
	    }
	}
    }
}
