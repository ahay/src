#ifndef _stretch_h
#define _stretch_h

typedef struct Map *map;

map stretch_init (int n1, float o1, float d1, int nd, float eps, bool narrow);
void stretch_define (map str, float* coord);
void stretch_apply (map str, float* ord, float* mod);
void stretch_invert (map str, float* ord, float* mod);
void stretch_close (map str);

#endif

/* 	$Id: stretch.h,v 1.3 2004/03/18 03:23:49 fomels Exp $	 */
