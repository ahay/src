#ifndef _pwd_h
#define _pwd_h

#include <rsf.h>

typedef struct Pwd *pwd;

pwd pwd_init(int n1, int nw);
void pwd_close (pwd w);
void pwd_define (bool adj, pwd w, const float* pp, float* diag, float** offd);
void pwd_set (pwd w, float* inp, float* out, float* tmp);

#endif
