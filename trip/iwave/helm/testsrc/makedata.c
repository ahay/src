#include "cstd.h"

int main() {

  /* 2D 10 x 10 suitable for 2x2 decomp */
  FILE * fph = NULL;
  FILE * fpd = NULL;
  float * a  = NULL;
  int i;

  fph=fopen("pt2d.rsf","w");
  fprintf(fph,"n1=100\nn2=100\nd1=1\nd2=1\no1=0\no2=0\ndata_format=native_float\nscale=0\nin=pt2d.rsf@");
  fclose(fph);
  fpd=fopen("pt2d.rsf@","w");
  a = (float *) malloc(10000*sizeof(float));
  for (i=0;i<10000;i++) a[i]=0.0f;
  a[5050]=1.0f;
  fwrite(a,sizeof(float),10000,fpd);
  fclose(fpd);

}
