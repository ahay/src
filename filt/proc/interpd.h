#ifndef _interpd_h
#define _interpd_h

void interp_init (int n, float e, int verb);
void interp_close (void);
void interp2(int n2, float** in, float** out, float** pp);

#endif
