#ifndef _allp2_h
#define _allp2_h

typedef struct Allpass2 *allpass2;

allpass2 allpass2_init(int nw, int nj, int nx, int ny, float **pp);
void allpass21 (bool der, const allpass2 ap, float** xx, float** yy);

#endif

/* 	$Id: allp2.h,v 1.1 2004/02/14 06:59:24 fomels Exp $	 */
