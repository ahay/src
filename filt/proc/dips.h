#ifndef _dips_h
#define _dips_h

void dips_init(int nd1, int nw, int nj, int nx, int ny, float** x1);
void dips_close(void);
void dips(const float *d, float *b, float **aa);

#endif
