#include <math.h>
#include <stdio.h>

#include "bandpass.h"
#include "butterworth.h"

static int hipoles, lopoles;
static float hidb, lodb;

void bandpass_init (int verb)
{
    const float lopass=0.25, apass=0.95, lostop=0.3, astop=0.05;
    const float hipass=0.02, histop=0.01;

    bfdesign (lopass, sqrt(apass), lostop, sqrt(astop), &lopoles, &lodb);
    bfdesign (hipass, sqrt(apass), histop, sqrt(astop), &hipoles, &hidb);

    if (verb) fprintf(stderr,"lopoles=%d hipoles=%d lodb=%g hidb=%g\n",
		      lopoles,hipoles,lodb,hidb);
}

void bandpass (int n1, float* trace)
{
    int i1;
    float t;

    bflowpass (lopoles, lodb, n1, trace, trace);
    for (i1=0; i1 < n1/2; i1++) { 
	t=trace[i1];
	trace[i1]=trace[n1-1-i1];
	trace[n1-1-i1]=t;
    }
    bflowpass (lopoles, lodb, n1, trace, trace);
    for (i1=0; i1 < n1/2; i1++) { 
	t=trace[i1];
	trace[i1]=trace[n1-1-i1];
	trace[n1-1-i1]=t;
    }
    
    bfhighpass (hipoles, hidb, n1, trace, trace);
    for (i1=0; i1 < n1/2; i1++) { 
	t=trace[i1];
	trace[i1]=trace[n1-1-i1];
	trace[n1-1-i1]=t;
    }
    bfhighpass (hipoles, hidb, n1, trace, trace);
    for (i1=0; i1 < n1/2; i1++) { 
	t=trace[i1];
	trace[i1]=trace[n1-1-i1];
	trace[n1-1-i1]=t;
    }
}



 
