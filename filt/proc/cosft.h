#ifndef _cosft_h
#define _cosft_h

void cosft_init(int n1 /*, float o1, float d1 */);
void cosft_close(void);
void cosft_frw (float *q, int o1, int d1);
void cosft_inv (float *q, int o1, int d1);

#endif

/* 	$Id: cosft.h,v 1.3 2003/12/04 05:13:21 fomels Exp $	 */
