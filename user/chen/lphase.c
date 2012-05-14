/* Linear phase filters */

#include <rsf.h>
#include "lphlag.h"
#include "lphmf.h"
#include "lphbsp.h"

static	int itp;	// 0: maxflat 	1: Lagrange 2: bspline
static	void *h;	// handle


void lphase_init(int n, bool causal, int interp)
/*< initiliaze >*/
{
	
	itp = interp;
	switch(interp)
	{
	case 1:
		h = lphlag_init(n, causal);
		break;
	case 2:
		h = lphbsp_init(n, causal);
		break;
	default:
		h = lphmf_init(n, causal);
		break;
	}
}

void lphase_filt(float delay, float *out)
/*< filter >*/
{
	switch(itp)
	{
	case 1:
		lphlag_filt(h, delay, out);
		break;
	case 2:
		lphbsp_filt(h, delay, out);
		break;
	default:
		lphmf_filt(h, delay, out);
		break;
	}
}

void lphase_dfilt(float delay, float *out)
/*< derivative filter >*/
{
	switch(itp)
	{
	case 1:
		lphlag_dfilt(h, delay, out);
		break;
	case 2:
		lphbsp_dfilt(h, delay, out);
		break;
	default:
		lphmf_dfilt(h, delay, out);
		break;
	}
}


void lphase_freq(float delay, int nk, sf_complex *out)
/*< frequency response >*/
{
	switch(itp)
	{
	case 1:
		lphlag_freq(h, delay, nk, out);
		break;
	case 2:
		lphbsp_freq(h, delay, nk, out);
		break;
	default:
		lphmf_freq(h, delay, nk, out);
		break;
	}
}

void lphase_close()
/*< release memory >*/
{
	switch(itp)
	{
	case 1:
		lphlag_close(h);
		break;
	case 2:
		lphbsp_close(h);
		break;
	default:
		lphmf_close(h);
		break;
	}
}

