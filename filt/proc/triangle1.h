#ifndef _triangle1_h
#define _triangle1_h

void triangle1_init (int nbox, int ndat);
void triangle1_lop (bool adj, bool add, int nx, int ny, float* x, float* y);
void triangle1_close(void);

#endif
