#ifndef _allp2_h
#define _allp2_h

typedef struct Allpass2 *allpass2;

allpass2 allpass2_init(int nw, int nj, int nx, int ny, float **pp);
void allpass21 (bool der, const allpass2 ap, float** xx, float** yy);

void allpass22_init (allpass2 ap1);
void allpass21_lop (bool adj, bool add, int n1, int n2, float* xx, float* yy);

#endif

/* 	$Id$	 */
