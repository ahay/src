#ifndef _predict_h
#define _predict_h

void predict_init (int nx, int ny, float eps);
void predict_close (void);
void predict_flat (int i0, float** d, float** m, float** pp);

#endif
