#ifndef _fastmarch_h
#define _fastmarch_h

void fastmarch_init (int n3,int n2,int n1);
void fastmarch (float* time, float* v, int* in,
		int n3,int n2,int n1,
		float o3,float o2,float o1,
		float d3,float d2,float d1,
		float s3,float s2,float s1,
		int b3, int b2, int b1,
		int order);
void fastmarch_close (void);

#endif

/* 	$Id: fastmarch.h,v 1.2 2003/09/30 14:30:52 fomels Exp $	 */
