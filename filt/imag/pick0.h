#ifndef _pick0_h
#define _pick0_h

void pick0_init (int n1_in, int n2_in, int order);
void pick0_set (int i2, float* dip);
void pick0_close (void);
void pick0_step (float t0, float* t);
void pick0_step0 (float t0, float* t);
void pick0_delta (int k2, float* t);

#endif

/* 	$Id: pick0.h,v 1.2 2003/09/30 14:30:53 fomels Exp $	 */
