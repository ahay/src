/* Copyright (c) Colorado School of Mines, 2007.*/
/* All rights reserved.                       */

#include <rsf.h>

float fsinc (float x)
/*< Return sinc(x) = sin(PI*x)/(PI*x) (float version) >*/
/******************************************************************************
Author:  Dave Hale, Colorado School of Mines, 06/02/89
*****************************************************************************/
{
	float pix;

	if (x==0.0) {
		return 1.0;
	} else {
	    pix = SF_PI*x;
	    return sinf(pix)/pix;
	}
}

double dsinc (double x)
/*< Return sinc(x) = sin(PI*x)/(PI*x) (double version) >*/
/******************************************************************************
Author:  Dave Hale, Colorado School of Mines, 06/02/89
*****************************************************************************/
{
	double pix;

	if (x==0.0) {
		return 1.0;
	} else {
	    pix = SF_PI*x;
	    return sin(pix)/pix;
	}
}
