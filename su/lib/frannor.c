/* Copyright (c) Colorado School of Mines, 2003.*/
/* All rights reserved.                       */

/*********************** self documentation **********************/
/*****************************************************************************
FRANNOR - functions to generate a pseudo-random float normally distributed
		with N(0,1); i.e., with zero mean and unit variance.

frannor		return a normally distributed random float
srannor		seed random number generator for normal distribution

******************************************************************************
Function Prototypes:
float frannor (void);
void srannor (int seed);

******************************************************************************
frannor:
Input:		(none)
Returned:	normally distributed random float

srannor:
Input:
seed		different seeds yield different sequences of random numbers.

******************************************************************************
Notes:
Adapted from subroutine rnor in Kahaner,  Moler and Nash (1988)
which in turn was based on  an algorithm by 
Marsaglia and Tsang (1984).

******************************************************************************
References:
"Numerical Methods and Software", D.
Kahaner, C. Moler, S. Nash, Prentice Hall, 1988.

Marsaglia G. and Tsang, W. W., 1984,
A fast, easily implemented method for sampling from decreasing or symmetric
unimodal density functions:  SIAM J. Sci. Stat. Comput., v. 5, no. 2,
p. 349-359.

*****************************************************************************
Author:  Dave Hale, Colorado School of Mines, 01/21/89
*****************************************************************************/
/**************** end self doc ********************************/

#include "frannor.h"
#include <math.h>

/* constants for normal random number generator */
#define AA 12.37586
#define B 0.4878992
#define C 12.67706
#define C1 0.9689279
#define C2 1.301198
#define PC 0.01958303
#define XN 2.776994
#define OXN 0.3601016
#define NBITS 24
#define ABS(x) ((x) < 0 ? -(x) : (x))
#define SGN(x) ((x) < 0 ? -1.0 : 1.0)
#define MAX(x,y) ((x) > (y) ? (x) : (y))
#define MIN(x,y) ((x) < (y) ? (x) : (y))
static float v[]={
	0.3409450, 0.4573146, 0.5397793, 0.6062427, 0.6631691,
	0.7136975, 0.7596125, 0.8020356, 0.8417227, 0.8792102, 0.9148948,
	0.9490791, 0.9820005, 1.0138492, 1.0447810, 1.0749254, 1.1043917,
	1.1332738, 1.1616530, 1.1896010, 1.2171815, 1.2444516, 1.2714635,
	1.2982650, 1.3249008, 1.3514125, 1.3778399, 1.4042211, 1.4305929,
	1.4569915, 1.4834526, 1.5100121, 1.5367061, 1.5635712, 1.5906454,
	1.6179680, 1.6455802, 1.6735255, 1.7018503, 1.7306045, 1.7598422,
	1.7896223, 1.8200099, 1.8510770, 1.8829044, 1.9155830, 1.9492166,
	1.9839239, 2.0198430, 2.0571356, 2.0959930, 2.1366450, 2.1793713,
	2.2245175, 2.2725185, 2.3239338, 2.3795007, 2.4402218, 2.5075117,
	2.5834658, 2.6713916, 2.7769943, 2.7769943, 2.7769943, 2.7769943
};

/* internal state variables for uniform random number generator */
static int i=16,j=4;
static float u[]={
	0.8668672834288,  0.3697986366357,  0.8008968294805,
	0.4173889774680,  0.8254561579836,  0.9640965269077,
	0.4508667414265,  0.6451309529668,  0.1645456024730,
	0.2787901807898,  0.06761531340295, 0.9663226330820,
	0.01963343943798, 0.02947398211399, 0.1636231515294,
	0.3976343250467,  0.2631008574685
};

/* macro to generate a random number uni uniform on [0,1) */
#define UNI(uni) \
	uni = u[i]-u[j]; if (uni<0.0) uni += 1.0; u[i] = uni; \
	if (--i<0) i = 16;  if (--j<0) j = 16

float
frannor(void)
/**************************************************************************
return a normally distributed random float
***************************************************************************
Returned:	normally distributed random float
**************************************************************************/
{
	int k;
	float uni,vni,rnor,x,y,s,bmbx,xnmx;
	
	/* uni is uniform on [0,1) */
	UNI(uni);
		
	/* vni is uniform on [-1,1) */
	vni = uni+uni-1.0;
	
	/* k is in range [0,63] */
	k = ((int)(u[i]*128))%64;
	
	/* fast part */
	rnor = vni*v[k+1];
	if (ABS(rnor)<=v[k]) return rnor;
	
	/* slow part */
	x = (ABS(rnor)-v[k])/(v[k+1]-v[k]);
	UNI(y);
	s = x+y;
	if (s<=C2) {
		if (s<=C1) return rnor;
		bmbx = B-B*x;
		if (y<=C-AA*exp(-0.5*bmbx*bmbx)) {
			if (exp(-0.5*v[k+1]*v[k+1])+y*PC/v[k+1] <= 
				exp(-0.5*rnor*rnor)) return rnor;
			do {
				UNI(y);
				x = OXN*log(y);
				UNI(y);
			} while (-2.0*log(y)<=x*x);
			xnmx = XN-x;
			return (rnor>=0.0 ? ABS(xnmx) : -ABS(xnmx));
		}
	}
	bmbx = B-B*x;	
	return (rnor>=0.0 ? ABS(bmbx) : -ABS(bmbx));
}			

void
srannor (int seed)
/*****************************************************************************
seed random number generator
******************************************************************************
Input:
seed		different seeds yield different sequences of random numbers.
*****************************************************************************/
{
	int ii,jj,ia,ib,ic,id;
	float s,t;
	
	i = 16;
	j = 4;
	ia=ABS(seed)%32707;
	ib=1111;
	ic=1947;
	for (ii=0; ii<17; ii++) {
		s = 0.0;
		t = 0.5;
		for (jj=0; jj<64; jj++) {
			id = ic-ia;
			if (id<0) {
				id += 32707;
				s += t;
			}
			ia = ib;
			ib = ic;
			ic = id;
			t *= 0.5;
		}
		u[ii] = s;
	}
}
