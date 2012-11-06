#include "rarray.h"

int main(int argc, char ** argv) {
  
  int ndim=2;
  IPNT gs0, ge0, gs, ge, i, n0, n;
  RARR arr;
  int err=0;

  IASN(gs0,IPNT_0);
  ge0[0]=10;
  ge0[1]=7;
  gs0[2]=1;
  ge0[2]=0;
  gs[0]=2; 
  ge[0]=8;
  gs[1]=2;
  ge[1]=6;
  gs[2]=1;
  ge[2]=0;

  ra_setnull(&arr);
  fprintf(stderr,"main->ra_create\n");
  err=ra_create(&arr,ndim,gs0,ge0);
  fprintf(stderr,"main<-ra_create err=%d\n",err);

  fprintf(stderr,"main->ra_a_size\n");
  ra_a_size(&arr,n0);
  fprintf(stderr,"main<-ra_a_size\n");
  
  fprintf(stderr,"set array\n");
  for (i[1]=gs0[1]; i[1]<=ge0[1]; i[1]++) {
    for (i[0]=gs0[0]; i[0]<=ge0[0]; i[0]++) {
      arr._s0[i[0]+i[1]*n0[0]] = 
	i[0]+10*i[1];
    }
  }
  fprintf(stderr,"main->ra_greset\n");  
  ra_greset(&arr,gs,ge);
  fprintf(stderr,"main<-ra_greset\n");  

  fprintf(stderr,"************\n");
  ra_dump(&arr,stderr);
  fprintf(stderr,"************\n");
  ra_print(&arr,stderr);
  fprintf(stderr,"************\n");
  fprintf(stderr,"main->ra_size\n");  
  ra_size(&arr,n);
  fprintf(stderr,"main<-ra_size\n");  
  fprintf(stderr,"************\n");
  for (i[1]=0;i[1]<n[1];i[1]++) {
    for (i[0]=0;i[0]<n[0];i[0]++) {
      fprintf(stderr,"%+11.3e ",arr._s2[i[1]][i[0]]);
    }
    fprintf(stderr,"\n");
  }

  fprintf(stderr,"************\n");
  fprintf(stderr,"s2[0][-1]=%+11.3e\n",arr._s2[0][-1]);
  fprintf(stderr,"s2[-1][0]=%+11.3e\n",arr._s2[-1][0]);
  fprintf(stderr,"s2[-1][-1]=%+11.3e\n",arr._s2[-1][-1]);
  fprintf(stderr,"s02[2][2]=%+11.3e\n",arr._s02[2][2]);
  ra_destroy(&arr);

}
