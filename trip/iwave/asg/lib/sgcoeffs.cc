#include "utils.h"


float * sgcoeffs(int k) {
  int n,i;
  float * c = (float *)usermalloc_(k*sizeof(float));
  for (n=1;n<=k;n++) {
    c[n-1]=1.0f;
    for (i=1;i<=k;i++) {
      if (i != n) {
	c[n-1] *= ((2.0f*i-1.0f)*(2.0f*i-1.0f))
	  /(((2.0f*n-1.0f)*(2.0f*n-1.0f)) - ((2.0f*i-1.0f)*(2.0f*i-1.0f)));
      }
    }
    c[n-1]=fabs(c[n-1])/(2.0f*n-1.0f);
    if (((n-(n/2)*2) == 0)) c[n-1]=-c[n-1];
  }
  return c;
}
