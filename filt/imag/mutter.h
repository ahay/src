#ifndef _mutter_h
#define _mutter_h

void mutter_init (int n1, float o1, float d1, 
		  bool abs, bool inner, bool hyper);
void mutter (float tp, float slope0, float slopep, float x, float *data);

#endif

/* 	$Id$	 */
