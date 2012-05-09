#include "cstd.h"

int main() {

  /* 2D 10 x 10 suitable for 2x2 decomp */
  FILE * fph = NULL;
  FILE * fpd = NULL;
  float * a  = NULL;
  int i;

  fph=fopen("test2d.rsf","w");
  fprintf(fph,"n1=10\nn2=10\nd1=1\nd2=1\no1=0\no2=10\ndata_format=native_float\nscale=3\nin=test2d.rsf@");
  fclose(fph);
  fpd=fopen("test2d.rsf@","w");
  a = (float *) malloc(100*sizeof(float));
  for (i=0;i<100;i++) a[i]=(float)i;
  fwrite(a,sizeof(float),100,fpd);
  fclose(fpd);

  fph=fopen("proc0test2d.rsf","w");
  fprintf(fph,"n1=10\nn2=10\nd1=1\nd2=1\no1=0\no2=10\ndata_format=native_float\nscale=3\nin=proc0test2d.rsf@");
  fclose(fph);
  fpd=fopen("proc0test2d.rsf@","w");
  for (i=0;i<100;i++) a[i]=0.0f;
  fwrite(a,sizeof(float),100,fpd);
  fclose(fpd);

  /* par file for this test */
  fph=fopen("test2d.par","w");
  fprintf(fph,"Data=test2d.rsf\nmpi_np1=2\nmpi_np2=2\nmpi_np3=1\n");
  fclose(fph);

}
