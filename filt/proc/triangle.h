#ifndef _triangle_h
#define _triangle_h

#include <rsf.h>

void triangle_init (int nbox, int ndat);
void triangle (int o, int d, bool der, float *x);
void  triangle_close(void);

#endif
