#ifndef _mask_h
#define _mask_h

void mask_init(const bool *m_in);
void mask_lop(bool adj, bool add, int nx, int ny, float *x, float *y);

#endif
