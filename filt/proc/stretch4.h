#ifndef _stretch4_h
#define _stretch4_h

typedef struct Map4 *map4;

map4 stretch4_init (int n1, float o1, float d1, int nd, float eps);
void stretch4_define (map4 str, float* coord);
void stretch4_apply (map4 str, float* ord, float* mod);
void stretch4_invert (map4 str, float* ord, float* mod);
void stretch4_close (map4 str);

#endif

/* 	$Id: stretch4.h,v 1.2 2004/06/03 05:35:51 fomels Exp $	 */
