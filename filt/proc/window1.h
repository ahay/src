#ifndef _window1_h
#define _window1_h

void window1_init (int w_in, int nw_in, int n_in, float h_in);
int window1_apply (int iw, float* dat, bool left, bool right, float *win);

#endif
