#include "extend.h"

static const int nw = 3;
static const float a[] = {7./3., -5./3., 1./3.};

void extend (int ne, int nd, float *dat, float *ext)
{
    int i, j;
    float s;

    for (i=0; i < nd; i++) {
	ext[ne+i] = dat[i];
    }
    for (i=ne-1; i >= 0; i--) {
	for (s=0., j=0; j < nw; j++) {
	    s += a[j]*ext[i+j+1];
	}
	ext[i] = s;
    }
    for (i=nd+ne; i < nd+2*ne; i++) {
	for (s=0., j=0; j < nw; j++) {
	    s += a[j]*ext[i-j-1];
	}
	ext[i] = s;
    }
}

void extend2 (int ne, int n1, int n2, float** dat, float** ext, 
	      float* tmp1, float* tmp2)
{
    int i1, i2;
    for (i2=0; i2 < n2; i2++) {
	extend (ne,n1,dat[i2],ext[i2+ne]);
    }
    for (i1=0; i1 < n1+2*ne; i1++) {
	for (i2=0; i2 < n2; i2++) {
	    tmp1[i2] = ext[i2+ne][i1];
	}
	extend (ne,n2,tmp1,tmp2);
	for (i2=0; i2 < n2+2*ne; i2++) {
	    ext[i2][i1] = tmp2[i2];
	} 
    }
}
