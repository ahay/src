#ifndef _predict_h
#define _predict_h

void predict_init (int nx, int ny, float eps, int verb);
void predict_close (void);
void predict_flat (float** d, float** m, float** pp);

#endif
