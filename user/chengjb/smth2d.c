/*    Smoothing 2D array */
/*
  Copyright (C) 2012 Tongji University (Jiubing Cheng)
  
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
#include "_cjb.h"


/* Prototype for function used internally */
/*****************************************************************************
*****************************************************************************/
void tripd(float *d, float *e, float *b, int n)
/*<****************************************************************************
Given an n-by-n symmetric, tridiagonal, positive definite matrix A and
 n-vector b, following algorithm overwrites b with the solution to Ax = b
Input:
d[]     the diagonal of A 
e[]     the superdiagonal of A
b[]     the rhs of Ax=b

Output:
b[]     b[] is overwritten with the solution to Ax=b

*****************************************************************************
Notes:

Given an n-by-n symmetric, tridiagonal, positive definite matrix A and
 n-vector b, following algorithm overwrites b with the solution to Ax = b

*****************************************************************************
Author: Zhenyue Liu, Colorado School of Mines, 1993.
****************************************************************************>*/
{
        int k;
        float temp;

        /* decomposition */
        for(k=1; k<n; ++k){
           temp = e[k-1];
           e[k-1] = temp/d[k-1];
           d[k] -= temp*e[k-1];
        }

        /* substitution */
        for(k=1; k<n; ++k)  b[k] -= e[k-1]*b[k-1];

        b[n-1] /=d[n-1];
        for(k=n-1; k>0; --k)  b[k-1] = b[k-1]/d[k-1] - e[k-1]*b[k];

 }

/***************************************************************************
***************************************************************************/
void smooth2d(float **v, int n1, int n2, float r1, float r2, float rw)
/*<**************************************************************************
	int n1;		 number of points in x1 (fast) dimension
	int n2;		 number of points in x1 (fast) dimension 
        float **v0;       array of input velocities 
        float r1;        smoothing parameter for x1 direction
        float r2;        smoothing parameter for x2 direction
        float rw;        smoothing parameter for window
" Notes:								",
" Larger r1 and r2 result in a smoother data. Recommended ranges of r1 	", 
" and r2 are from 1 to 20.						",
"									",
" Smoothing can be implemented in a selected window. The range of 1st   ",
" dimension for window is from win[0] to win[1]; the range of 2nd   	",
" dimension is from win[2] to win[3]. 					",
"									",
" Smoothing the window function (i.e. blurring the edges of the window)	",
" may be done by setting a nonzero value for rw, otherwise the edges	",
" of the window will be sharp.						",
" 									",
**************************************************************************>*/
{
	int nmax;	/* max of n1 and n2 */
	int ix, iz;	/* counters */
	int win[4];	/* 1d array defining the corners of smoothing window */
	float **v0;	/* array of output velocities */
	float **w;	/* intermediate array */
	float *d, *e;	/* input arrays for subroutine tripd */
	float *f;	/* intermediate array */

	/* scale the smoothing parameter */
	r1 = r1*r1*0.25;
	r2 = r2*r2*0.25;

	/* allocate space */
	nmax = (n1<n2)?n2:n1;

	v0 = sf_floatalloc2(n1,n2);
	w = sf_floatalloc2(n1,n2);

	d = sf_floatalloc(nmax);
	e = sf_floatalloc(nmax);
	f = sf_floatalloc(nmax);

	for(ix=0; ix<nmax; ++ix)
	    d[ix] = e[ix] = f[ix] = 0.0f;

	/* save the original velocity */
        for(ix=0; ix<n2; ++ix)
	 	for(iz=0; iz<n1; ++iz)
			v0[ix][iz]=v[ix][iz];

	/* get parameters for window function */
	win[0] = 0;
	win[1] = n1;
	win[2] = 0;
	win[3] = n2;
	rw = rw*rw*0.25;
 
	/* define the window function */
	for(ix=0; ix<n2; ++ix)
	 	for(iz=0; iz<n1; ++iz)
			w[ix][iz] = 0;	
	for(ix=win[2]; ix<win[3]; ++ix)
	 	for(iz=win[0]; iz<win[1]; ++iz)
			w[ix][iz] = 1;	

	if(win[0]>0 || win[1]<n1 || win[2]>0 || win[3]<n2){
	/*	smooth the window function */
         	for(iz=0; iz<n1; ++iz){
	 		for(ix=0; ix<n2; ++ix){
				d[ix] = 1.0+2.0*rw;
				e[ix] = -rw;
				f[ix] = w[ix][iz];
			}
        		d[0] -= rw;
         		d[n2-1] -= rw;
         		tripd(d,e,f,n2);
	 		for(ix=0; ix<n2; ++ix)
				w[ix][iz] = f[ix];
		}
         	for(ix=0; ix<n2; ++ix){
	 		for(iz=0; iz<n1; ++iz){
				d[iz] = 1.0+2.0*rw;
				e[iz] = -rw;
				f[iz] = w[ix][iz];
		}
        		d[0] -= rw;
         		d[n1-1] -= rw;
         		tripd(d,e,f,n1);
	 		for(iz=0; iz<n1; ++iz)
				w[ix][iz] = f[iz];
		}
	}

	/*      solving for the smoothing velocity */
        for(iz=0; iz<n1; ++iz){
	 	for(ix=0; ix<n2-1; ++ix){
			d[ix] = 1.0+r2*(w[ix][iz]+w[ix+1][iz]);
			e[ix] = -r2*w[ix+1][iz];
			f[ix] = v[ix][iz];
		}
        	d[0] -= r2*w[0][iz];
         	d[n2-1] = 1.0+r2*w[n2-1][iz];
		f[n2-1] = v[n2-1][iz];
         	tripd(d,e,f,n2);
	 	for(ix=0; ix<n2; ++ix)
			v[ix][iz] = f[ix];
	}
         for(ix=0; ix<n2; ++ix){
	 	for(iz=0; iz<n1-2; ++iz){
			d[iz] = 1.0+r1*(w[ix][iz+1]+w[ix][iz+2]);
			e[iz] = -r1*w[ix][iz+2];
			f[iz] = v[ix][iz+1];
		}
		f[0] += r1*w[ix][1]*v[ix][0];
         	d[n1-2] = 1.0+r1*w[ix][n1-1];
		f[n1-2] = v[ix][n1-1];
         	tripd(d,e,f,n1-1);
	 	for(iz=0; iz<n1-1; ++iz)
			v[ix][iz+1] = f[iz];
	}

	free(d);
	free(e);
	free(f);
	free(*v0);
	free(*w);
}
