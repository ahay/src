#ifndef _allp3_h
#define _allp3_h

typedef struct Allpass *allpass;

allpass allpass_init(int nw, int nj, int nx, int ny, int nz, float ***pp);
void allpass1 (bool der, const allpass ap, float*** xx, float*** yy);
void allpass2 (bool der, const allpass ap, float*** xx, float*** yy);

#endif
