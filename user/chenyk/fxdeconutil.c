/* f-x deconolution utilities. */
/*
  Copyright (C) 2013 University of Texas at Austin
  
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

  References:      							
  Canales(1984):'Random noise reduction' 54th. SEGM	
  Gulunay(1986):'FXDECON and complex Wiener Predicition filter' 56th. SEGM	                
  Galbraith(1991):'Random noise attenuation by F-X prediction: a tutorial' 61th. SEGM
  Seismic Unix (SU) software package: http://www.cwp.mines.edu/cwpcodes/.
*/

#include<rsf.h>
#include"fxdeconutil.h"

kiss_fft_cpx cmplx(float re, float im)
/*< complex number >*/
{
    kiss_fft_cpx c;
    c.r = re;
    c.i = im;
    return c;
}

void backward_substitution (int nrows, /* number of rows (and columns) of input matrix */
			    float **matrix, /* matrix of coefficients (after LU decomposition) */ 
			    int *idx, /* permutation vector obtained from routine LU decomposition */
			    float *b /* right hand side vector in equation Ax=b */)
/*< backward_substitution >*/
{
	int i,ii = -1,j;	/* loop counters */
	int ip;			/* index of first nonvanishing element of b */
	float sum;		/* auxiliary variable for partial sums */

	for (i=0; i<nrows; i++) {
		ip = idx[i];

		/* do forward substitution */
		sum = b[ip];
		b[ip] = b[i];
		if (ii !=-1)
			for (j=ii; j<i; j++)
				sum -= matrix[i][j]*b[j];
		else if (sum !=0.0)
			ii = i;
		b[i] = sum;
	}

	/* now, do the backward substitution */
	for (i=nrows-1; i>=0; i--) {
		sum = b[i];
		for (j=i+1; j<nrows; j++)
			sum -= matrix[i][j]*b[j];

		/* store results in output vector */
		b[i] = sum/matrix[i][i];
	}
}

void lu_decomposition (int nrows, /* number of rows of matrix to invert */
		       float **matrix, /* matrix to invert */
		       int *idx,  /* vector recording the row permutations effected by partial pivoting */
		       float *d /*+/- 1 depending on whether the number of row interchanges was even or odd*/)
/*< LU decomposition >*/
{
	int i,j,k;	/* loop counters for rows, columns and summs */
	int imax=0;	/* index of maximum pivot */
	float big;	/* largest number in input matrix */
	float dum;	/* pivot scale factor */
	float sum;	/* auxiliary variable for summ  */
	float temp;	/* auxiliary variable */
	float *vv;	/* vector to store implicit scaling for each row */

	/* allocate working space */
	vv = sf_floatalloc(nrows);

	/* initialize interchanges counter */
	*d = 1.0;

	/* loop over rows to get implicit scaling information */
	for (i=0; i<nrows; i++) {
		big = 0.0;
		for (j=0; j<nrows; j++) {
		    temp=SF_ABS(matrix[i][j]);
		    if(temp>big) big=temp;
		}
		if (big == 0.0) 
			fprintf(stderr,"error, singular matrix in LU decomposition\n");
		
		/* save the scaling */
		vv[i] = 1.0/big;
	}

	/* loop over columns (Crout's method) */
	for (j=0; j<nrows; j++) {
		for (i=0; i<j; i++) {
			sum = matrix[i][j];
			for (k=0; k<i; k++)
				sum -= matrix[i][k]*matrix[k][j];
			matrix[i][j] = sum;
		}
	
		/* initialize for the search for largest pivot element */
		big = 0.0;
		for (i=j; i<nrows; i++) {
			sum = matrix[i][j];
			for (k=0; k<j; k++)
				sum -= matrix[i][k]*matrix[k][j];
			matrix[i][j] = sum;

			/*  Is new pivot better than best so far? */
			if ((dum=vv[i]*SF_ABS(sum)) >= big) {
 				big = dum;
				imax = i;
			}
		}
	
		/* Do we need to interchange rows */
		if (j != imax) {
			for (k=0; k<nrows; k++) {
				dum = matrix[imax][k];
				matrix[imax][k] = matrix[j][k];
				matrix[j][k] = dum;
			}

			/* change the parity of d */
			*d = -(*d);

			/* interchange the scale factor */
			vv[imax] = vv[j];
		}
		idx[j] = imax;

		/* if matrix becomes singular don't use pivot=0 */
		if (matrix[j][j] == 0.0f) matrix[j][j]= SF_EPS;
	
		if (j !=nrows) {

			/* divide by the pivot element */
			dum = 1.0/matrix[j][j];
			for (i=j+1; i<nrows; i++)
				matrix[i][j] *= dum;
		}
	}

	/* free workspace */
	free(vv);
}


void cxcor (int lx, int ifx, kiss_fft_cpx *x,
            int ly, int ify, kiss_fft_cpx *y, 
            int lz, int ifz, kiss_fft_cpx *z)
/*< complex correlation >*/
{
        int i,j;
	kiss_fft_cpx *xr;
	xr = (kiss_fft_cpx*) sf_complexalloc(lx);
	for (i=0,j=lx-1; i<lx; ++i,--j)
	        xr[i] = sf_conjf (x[j]);
	cconv(lx,1-ifx-lx,xr,ly,ify,y,lz,ifz,z);
	free(xr);
}


void cconv (int lx, int ifx, kiss_fft_cpx *x,
       	    int ly, int ify, kiss_fft_cpx *y,
	    int lz, int ifz, kiss_fft_cpx *z)
/*< complex convolution >*/
{
        int ilx=ifx+lx-1,ily=ify+ly-1,ilz=ifz+lz-1,i,j,jlow,jhigh;
	kiss_fft_cpx sum;
  
	x -= ifx;  y -= ify;  z -= ifz; 
	for (i=ifz; i<=ilz; ++i) {
	        jlow = i-ily;  if (jlow<ifx) jlow = ifx;
		jhigh = i-ify;  if (jhigh>ilx) jhigh = ilx;
		for (j=jlow,sum=cmplx(0.,0.); j<=jhigh; ++j){
		        sum = sf_cadd(sum,sf_cmul(x[j],y[i-j]));
		}
		z[i] = sum;
	}
}
