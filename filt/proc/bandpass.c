#include <math.h>
#include <stdio.h>

#include "bandpass.h"
#include "butter.h"

static butter blo, bhi;

void bandpass_init (void)
{
    const float fhi=0.275314, flo=0.0141445;
    const int nphi=10, nplo=5;

    blo = butter_init(false,flo,nplo);
    bhi = butter_init(true, fhi,nphi);
}

void bandpass_close (void)
{
    butter_close(blo);
    butter_close(bhi);
}

void bandpass (int n1, float* trace)
{
    int i1;
    float t;

    butter_apply(blo,n1,trace);
    for (i1=0; i1 < n1/2; i1++) { 
	t=trace[i1];
	trace[i1]=trace[n1-1-i1];
	trace[n1-1-i1]=t;
    }
    butter_apply(blo,n1,trace);
    for (i1=0; i1 < n1/2; i1++) { 
	t=trace[i1];
	trace[i1]=trace[n1-1-i1];
	trace[n1-1-i1]=t;
    }
    butter_apply(bhi,n1,trace);
    for (i1=0; i1 < n1/2; i1++) { 
	t=trace[i1];
	trace[i1]=trace[n1-1-i1];
	trace[n1-1-i1]=t;
    }
    butter_apply(bhi,n1,trace);
    for (i1=0; i1 < n1/2; i1++) { 
	t=trace[i1];
	trace[i1]=trace[n1-1-i1];
	trace[n1-1-i1]=t;
    }
}

/* 	$Id$	 */
