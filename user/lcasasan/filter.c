
/*
  Copyright (C) 2009 Politecnico di Milano
  Author: Emmanuel Spadavecchia     
  <emmanuel.spadavecchia@gmail.com> 

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

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <rsf.h>

#include "filter.h"

#define MY_MIN(X,Y) ((X) < (Y) ? (X) : (Y))
#define MY_MAX(X,Y) ((X) > (Y) ? (X) : (Y))
void hanning_win (int n, int flag, float* w) 
/*<Hanning window >*/{

	int i, semiLen;
	float a;
	
	semiLen = (int) (n/2);

	for (i = 0, a = 0; i < semiLen; i++) {
		w[i] = .5*(1 - cos(2*SF_PI*i/(n-1)));
		w[n-1-i] = w[i];
		a += 2*w[i];
	}

    if (n%2){
		w[semiLen]=1.0;
		a=+1.0;
	}
	
	if (flag == 0) {
		for (i = 0; i < n; i++)
			w[i] = w[i] / a;
	}

	// *** for (i = 0; i < n; i++) sf_warning("\nw(%d)=%f;", i+1, w[i]);
}

void gaussian_win (int n, int flag, float* w) 
/*<Gaussian window>*/{
	int i;
	float a,b,c;
	float sigmasq, sigma;

	sigma = n;

	b = (n-1)/2;
	b = roundf(b);
	sigmasq=sigma*sigma;

	for (i = 0, a = 0; i < n; i++) {
		c = (i-b)*(i-b)/(2*sigmasq);
		w[i] = expf(-c);
		a += w[i];
	}

	if (flag == 0) {
		for (i = 0; i < n; i++)
			w[i] = w[i] / a;
	}

	// *** for (i = 0; i < n; i++) printf("\nw(%d)=%.100f;", i+1, w[i]);
}

void boxcar_win (int n, int flag, float* w) {
/*<rectangular window>*/
	int i;
	float a;
	
	a = 1.0;
	if (flag == 0)
		a = 1 / (float)n;

	for (i = 0; i < n; i++) {
		w[i] = a;
	}

	// *** for (i = 0; i < n; i++) printf("\nw(%d)=%.100f;", i+1, w[i]);
}

void tri_win (int n, int flag, float* w) {
/*<triangular window>*/
	int i, semiLen;
	float a;
	
	semiLen = (int) (n/2);
    
	for (i = 0, a = 0; i < semiLen; i++) {
		w[i] = 2*i/(float)(n-1);
		
		w[n-1-i] = w[i];
			
		a += 2*w[i];
	}

	if (n%2){
		w[semiLen]=1.0;
		a=+1.0;
	}

	if (flag == 0) {
		for (i = 0; i < n; i++)
			w[i] = w[i] / a;
	}

	// *** for (i = 0; i < n; i++) printf("\nw(%d)=%.100f;", i+1, w[i]);
}

void conv (const float* x, int xLen, const float* h, int hLen, float* y) {
/*<1D convolution>*/
	int i,j,hSemiLen,ind;
	
	if (hLen == 1) {
		for (i = 0; i < xLen; i++)
			y[i] = x[i];
	} else {

		hSemiLen=(int) ((hLen-1)/2);
		for (i = 0; i < xLen; i++) {
			y[i] = 0.0;
			for(j = -1*hSemiLen; j <= hSemiLen; j++) {

				if (i+j <0 || i+j >= xLen)
					ind = i-j;	
				else
					ind = i+j;

				y[i] += x[ind] * h[j + hSemiLen];
			}
		}

	}
	// *** for (i = 0; i < xLen; i++) printf("\ny(%d)=%.100f;", i+1, y[i]);
}

/* 		for (i = 0; i < xLen; i++) { */
/* 			y[i] = 0.0; */
/* 			for(j = -1*hSemiLen; j <= hSemiLen; j++) { */
/* 				if (i+j >=0 && i+j < xLen) */
/* 					y[i] += x[i+j] * h[j + hSemiLen]; */
/* 			} */
/* 		} */

void transpose (float* A, int N, int M, float* At) {
/*<matrix transpose>*/	
	int i, j;
	
	for ( i = 0; i < N; i++) {
		for ( j = 0; j < M; j++) {
			At[(i*M)+j] = A[(j*N)+i]; 
		}
	}
}

void conv_2D (const float* x, int N, int M, const float* hCol, int hColLen, const float* hRow, int hRowLen, float* y) {
/*<2D convolution>*/	
	int i, j;
	
	// alloca la memoria temporanea
	float *tmp  = (float*) malloc (N*M*sizeof(float));    // N x M
	float *tmp2 = (float*) malloc (N*M*sizeof(float));    // M x N

	// filtra le colonne
	for (j = 0; j < M; j++)
		conv (&x[j*N], N, hCol, hColLen, &tmp[j*N]);

	// trasponi (inverti ordine righe/colonne)
	transpose (tmp, N, M, tmp2);

	// filtra per righe
	for (i = 0; i < N; i++)
		conv (&tmp2[i*M], M, hRow, hRowLen, &tmp[i*M]);

	// trasponi
	transpose (tmp, M, N, y);

	// libera la memoria
	free(tmp);
	free(tmp2);
}

void median_filter (const float* x, int xLen, int N, float* y) {
/*<median filter>*/
       int i,j,k,ind;
       float tmp, *TMP;

       TMP  = (float*) malloc ((2*N+1)*sizeof(float));

       for (k = 0; k < xLen; k++)
       {
               for (i = -N; i <= N; i++)
               {
                       ind = k+i;
                       ind = MY_MAX (ind, 0);
                       ind = MY_MIN (ind,xLen-1);
                       TMP[i+N] = x[ind];
                       // *** printf("\n\t%d) x(%d)=%f ", i, ind, TMP[i+N]);
               }

               for (i = 0; i <= 2*N; i++)
               {
                       for (j = i+1; j <= 2*N; j++)
                       {
                               if (TMP[i] > TMP[j])
                               {
                                       tmp = TMP[i];
                                       TMP[i] = TMP[j];
                                       TMP[j] = tmp;
                               }
                       }
                       // *** printf("\n\t%d) x_ord=%f;", i, TMP[i]);
               }
               y[k] = TMP[N];
       }
       free(TMP);
}
