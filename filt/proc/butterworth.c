#include <math.h>

#include <rsf.h>

#include "butterworth.h"

static const float pi=3.141592653589793;

/* Copyright (c) Colorado School of Mines, 2000.*/
/* All rights reserved.                       */

/*********************** self documentation **********************/
/*****************************************************************************
BUTTERWORTH - Functions to design and apply Butterworth filters:

bfdesign	design a Butterworth filter
bfhighpass	apply a high-pass Butterworth filter 
bflowpass	apply a low-pass Butterworth filter 

******************************************************************************
Function Prototypes:
void bfhighpass (int npoles, float f3db, int n, float p[], float q[]);
void bflowpass (int npoles, float f3db, int n, float p[], float q[]);
void bfdesign (float fpass, float apass, float fstop, float astop,
	int *npoles, float *f3db);

	******************************************************************************
bfdesign:
Input:
fpass		frequency in pass band at which amplitude is >= apass
apass		amplitude in pass band corresponding to frequency fpass
fstop 		frequency in stop band at which amplitude is <= astop
astop		amplitude in stop band corresponding to frequency fstop

Output:
npoles		number of poles
f3db		frequency at which amplitude is sqrt(0.5) (-3 db)

bfhighpass and bflowpass:
Input:
npoles		number of poles (and zeros); npoles>=0 is required
f3db		3 db frequency; nyquist = 0.5; 0.0<=f3db<=0.5 is required
n		length of p and q
p		array[n] to be filtered

Output:
q		filtered array[n] (may be equivalent to p)

******************************************************************************
Notes:
(1) Nyquist frequency equals 0.5

(2) The following conditions must be true:
	(0.0<fpass && fpass<0.5) &&
	(0.0<fstop && fstop<0.5) &&
	(fpass!=fstop) &&
	(0.0<astop && astop<apass && apass<1.0)

(3) if (fpass<fstop)

bfdesign:
Butterworth filter:  compute number of poles and -3 db frequency
for a low-pass or high-pass filter, given a frequency response
constrained at two frequencies.

******************************************************************************
Author:  Dave Hale, Colorado School of Mines, 06/02/89
*****************************************************************************/
/**************** end self doc ********************************/

    void
	bfdesign (float fpass, float apass, float fstop, float astop,
		  int *npoles, float *f3db)
/*****************************************************************************
Butterworth filter:  compute number of poles and -3 db frequency
for a low-pass or high-pass filter, given a frequency response
constrained at two frequencies.
******************************************************************************
Input:
fpass		frequency in pass band at which amplitude is >= apass
apass		amplitude in pass band corresponding to frequency fpass
fstop 		frequency in stop band at which amplitude is <= astop
astop		amplitude in stop band corresponding to frequency fstop

Output:
npoles		number of poles
f3db		frequency at which amplitude is sqrt(0.5) (-3 db)
******************************************************************************
Notes:
(1) Nyquist frequency equals 0.5

(2) The following conditions must be true:
	(0.0<fpass && fpass<0.5) &&
	(0.0<fstop && fstop<0.5) &&
	(fpass!=fstop) &&
	(0.0<astop && astop<apass && apass<1.0)

(3) if (fpass<fstop)
		a low-pass filter is assumed
	else
		a high-pass filter is assume
******************************************************************************
Author:  Dave Hale, Colorado School of Mines, 06/02/89
*****************************************************************************/
{
    float wpass,wstop,fnpoles,w3db;

    /* warp frequencies according to bilinear transform */
    wpass = 2.0*tanf(pi*fpass);
    wstop = 2.0*tanf(pi*fstop);

    /* if lowpass filter, then */
    if (fstop>fpass) {
	fnpoles = 0.5*logf((1.0/(apass*apass)-1.0)/(1.0/(astop*astop)-1.0))
	    / logf(fabsf(wpass/wstop));
	w3db = wpass/powf((1.0/(apass*apass)-1.0),0.5/fnpoles);

	/* else, if highpass filter, then */
    } else {
	fnpoles = 0.5*logf((1.0/(apass*apass)-1.0)/(1.0/(astop*astop)-1.0))
	/ logf(fabsf(wstop/wpass));
	w3db = wpass*powf((1.0/(apass*apass)-1.0),0.5/fnpoles);
    }

    /* determine integer number of poles */
    *npoles = 1+(int)fnpoles;

    /* determine (unwarped) -3 db frequency */
    *f3db = atanf(0.5*w3db)/pi;
}

void
bfhighpass (int npoles, float f3db, int n, float p[], float q[])
/*****************************************************************************
Butterworth filter:  high-pass
******************************************************************************
Input:
npoles		number of poles (and zeros); npoles>=0 is required
f3db		3 db frequency; nyquist = 0.5; 0.0<=f3db<=0.5 is required
n		length of p and q
p		array[n] to be filtered

Output:
q		filtered array[n] (may be equivalent to p)
******************************************************************************
Author:  Dave Hale, Colorado School of Mines, 06/02/89
*****************************************************************************/
{
    int jpair,j;
    float r,scale,theta,a,b1,b2,pj,pjm1,pjm2,qjm1,qjm2;

    r = 2.0*tan(pi*fabs(f3db));
    if (npoles%2!=0) {
	scale = r+2.0;
	a = 2.0/scale;
	b1 = (r-2.0)/scale;
	pj = 0.0;
	qjm1 = 0.0;
	for (j=0; j<n; j++) {
	    pjm1 = pj;
	    pj = p[j];
	    q[j] = a*(pj-pjm1)-b1*qjm1;
	    qjm1 = q[j];
	}
    } else {
	for (j=0; j<n; j++)
	    q[j] = p[j];
    }
    for (jpair=0; jpair<npoles/2; jpair++) {
	theta = pi*(2*jpair+1)/(2*npoles);
	scale = 4.0+4.0*r*sinf(theta)+r*r;
	a = 4.0/scale;
	b1 = (2.0*r*r-8.0)/scale;
	b2 = (4.0-4.0*r*sinf(theta)+r*r)/scale;
	pjm1 = 0.0;
	pj = 0.0;
	qjm2 = 0.0;
	qjm1 = 0.0;
	for (j=0; j<n; j++) {
	    pjm2 = pjm1;
	    pjm1 = pj;
	    pj = q[j];
	    q[j] = a*(pj-2.0*pjm1+pjm2)-b1*qjm1-b2*qjm2;
	    qjm2 = qjm1;
	    qjm1 = q[j];
	}
    }
}

void
bflowpass (int npoles, float f3db, int n, float p[], float q[])
/*****************************************************************************
Butterworth filter:  low-pass
******************************************************************************
Input:
npoles		number of poles (and zeros); npoles>=0 is required
f3db		3 db frequency; nyquist = 0.5; 0.0<=f3db<=0.5 is required
n		length of p and q
p		array[n] to be filtered

Output:
q		filtered array[n] (may be equivalent to p)
******************************************************************************
Author:  Dave Hale, Colorado School of Mines, 06/02/89
*****************************************************************************/
{
    int jpair,j;
    float r,scale,theta,a,b1,b2,pj,pjm1,pjm2,qjm1,qjm2;

    r = 2.0*tan(pi*fabs(f3db));
    if (npoles%2!=0) {
	scale = r+2.0;
	a = r/scale;
	b1 = (r-2.0)/scale;
	pj = 0.0;
	qjm1 = 0.0;
	for (j=0; j<n; j++) {
	    pjm1 = pj;
	    pj = p[j];
	    q[j] = a*(pj+pjm1)-b1*qjm1;
	    qjm1 = q[j];
	}
    } else {
	for (j=0; j<n; j++)
	    q[j] = p[j];
    }
    for (jpair=0; jpair<npoles/2; jpair++) {
	theta = pi*(2*jpair+1)/(2*npoles);
	scale = 4.0+4.0*r*sinf(theta)+r*r;
	a = r*r/scale;
	b1 = (2.0*r*r-8.0)/scale;
	b2 = (4.0-4.0*r*sinf(theta)+r*r)/scale;
	pjm1 = 0.0;
	pj = 0.0;
	qjm2 = 0.0;
	qjm1 = 0.0;
	for (j=0; j<n; j++) {
	    pjm2 = pjm1;
	    pjm1 = pj;
	    pj = q[j];
	    q[j] = a*(pj+2.0*pjm1+pjm2)-b1*qjm1-b2*qjm2;
	    qjm2 = qjm1;
	    qjm1 = q[j];
	}
    }
}

/* 	$Id: butterworth.c,v 1.4 2003/10/01 22:45:56 fomels Exp $	 */

