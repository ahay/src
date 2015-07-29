#ifndef _aal_h
#define _aal_h

void aal_doubint(int nt, float *trace);
/* causal and anticausal integration */

float aal_pick(float ti, float deltat, 
	       const float *trace,
	       int nt, float dt, float t0);
/* pick a traveltime sample from a trace */

#endif
