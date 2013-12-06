/* Copyright (c) Colorado School of Mines, 2011.*/
/* All rights reserved.                       */

/*********************** self documentation **********************/
/*****************************************************************************
YXTOXY - Compute a regularly-sampled, monotonically increasing function x(y)
	from a regularly-sampled, monotonically increasing function y(x) by
	inverse linear interpolation.

yxtoxy		compute a regularly sampled function x(y) from a regularly
		sampled, monotonically increasing function y(x)

******************************************************************************
Function Prototype:
void yxtoxy (int nx, float dx, float fx, float y[], 
	int ny, float dy, float fy, float xylo, float xyhi, float x[]);

******************************************************************************
Input:
nx		number of samples of y(x)
dx		x sampling interval; dx>0.0 is required
fx		first x
y		array[nx] of y(x) values; y[0] < y[1] < ... < y[nx-1] required
ny		number of samples of x(y)
dy		y sampling interval; dy>0.0 is required
fy		first y
xylo		x value assigned to x(y) when y is less than smallest y(x)
xyhi		x value assigned to x(y) when y is greater than largest y(x)

Output:
x		array[ny] of x(y) values

******************************************************************************
 Notes:
 User must ensure that:
 (1) dx>0.0 && dy>0.0
 (2) y[0] < y[1] < ... < y[nx-1] */

/*******************************************************************************
Author:  Dave Hale, Colorado School of Mines, 06/02/89
*****************************************************************************/
/**************** end self doc ********************************/
#include "yxtoxy.h"
void yxtoxy (int nx, float dx, float fx, float y[], 
	int ny, float dy, float fy, float xylo, float xyhi, float x[])
/*< compute regularly sampled x(y) from  regularly sampled y(x)>*/
/*****************************************************************************
Compute a regularly-sampled, monotonically increasing function x(y) from a 
regularly-sampled, monotonically increasing function y(x) by inverse linear 
interpolation.
******************************************************************************
Input:
nx		number of samples of y(x)
dx		x sampling interval; dx>0.0 is required
fx		first x
y		array[nx] of y(x) values; y[0] < y[1] < ... < y[nx-1] required
ny		number of samples of x(y)
dy		y sampling interval; dy>0.0 is required
fy		first y
xylo		x value assigned to x(y) when y is less than smallest y(x)
xyhi		x value assigned to x(y) when y is greater than largest y(x)

Output:
x		array[ny] of x(y) values
******************************************************************************
Notes:
User must ensure that:
(1) dx>0.0 && dy>0.0
(2) y[0] < y[1] < ... < y[nx-1]
******************************************************************************
Author:  Dave Hale, Colorado School of Mines, 06/02/89
*****************************************************************************/
{
	int nxi,nyo,jxi1,jxi2,jyo;
	float dxi,fxi,dyo,fyo,fyi,yo,xi1,yi1,yi2; 

	nxi = nx; dxi = dx; fxi = fx;
	nyo = ny; dyo = dy; fyo = fy;
	fyi = y[0];

	/* loop over output y less than smallest input y */
	for (jyo=0,yo=fyo; jyo<nyo; jyo++,yo+=dyo) {
		if (yo>=fyi) break;
		x[jyo] = xylo;
	}

	/* loop over output y between smallest and largest input y */
	if (jyo==nyo-1 && yo==fyi) {
		x[jyo++] = fxi;
		yo += dyo;
	}
	jxi1 = 0;
	jxi2 = 1;
	xi1 = fxi;
	while (jxi2<nxi && jyo<nyo) {
		yi1 = y[jxi1];
		yi2 = y[jxi2];
		if (yi1<=yo && yo<=yi2) {
			x[jyo++] = xi1+dxi*(yo-yi1)/(yi2-yi1);
			yo += dyo;
		} else {
			jxi1++;
			jxi2++;
			xi1 += dxi;
		}
	}

	/* loop over output y greater than largest input y */
	while (jyo<nyo)
		x[jyo++] = xyhi;
}
