#ifndef _callp2_h
#define _callp2_h

typedef struct Callpass2 *callpass2;

callpass2 callpass2_init(int nw, int nj, int nx, int ny);
void callpass21 (bool der, const callpass2 ap, float** xx, float** yy);
void callpass21_set (callpass2 ap, float p);

#endif

/* 	$Id: callp2.h 704 2004-07-13 18:22:06Z fomels $	 */
