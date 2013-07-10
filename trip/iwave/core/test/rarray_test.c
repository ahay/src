#include "rarray.h"

int main(int argc, char ** argv) {
  
  int ndim=3;
  IPNT gs0, ge0, gs, ge, s, e, i, n0, n;
  RARR arr;

  IASN(gs0,IPNT_0);
  ge0[0]=10;
  ge0[1]=7;
  ge0[2]=4;
  gs[0]=2; 
  ge[0]=8;
  gs[1]=2;
  ge[1]=6;
  gs[2]=1;
  ge[2]=3;

  fprintf(stderr,"main->ra_create\n");
  ra_create(&arr,ndim,gs0,ge0);
  fprintf(stderr,"main<-ra_create\n");
  fprintf(stderr,"main->ra_a_size\n");
  ra_a_size(&arr,n0);
  fprintf(stderr,"main<-ra_a_size\n");
  
  fprintf(stderr,"set array\n");
  for (i[2]=gs0[2]; i[2]<=ge0[2]; i[2]++) {
    for (i[1]=gs0[1]; i[1]<=ge0[1]; i[1]++) {
      for (i[0]=gs0[0]; i[0]<=ge0[0]; i[0]++) {
	arr._s0[i[0]+i[1]*n0[0]+i[2]*n0[0]*n0[1]] = 
	  i[0]+10*i[1]+100*i[2];
      }
    }
  }
  fprintf(stderr,"main->ra_greset\n");  
  ra_greset(&arr,gs,ge);
  fprintf(stderr,"main<-ra_greset\n");  
  ra_se(&arr,s,e);

  fprintf(stderr,"************\n");
  ra_dump(&arr,stderr);
  fprintf(stderr,"********** alloc\n");
  for (i[2]=0;i[2]<n0[2];i[2]++) {
    fprintf(stderr,"slice %d\n",i[2]);
    for (i[1]=0;i[1]<n0[1];i[1]++) {
      for (i[0]=0;i[0]<n0[0];i[0]++) {
	fprintf(stderr,"%+11.3e ",
		arr._s0[i[2]*n0[0]*n0[1] +i[1]*n0[0] +i[0]]);
      }
      fprintf(stderr,"\n");
    }
  }
  fprintf(stderr,"************ comp\n");
  ra_print(&arr,stderr);
  fprintf(stderr,"************\n");
  fprintf(stderr,"main->ra_size\n");  
  ra_size(&arr,n);
  fprintf(stderr,"main<-ra_size\n");  
  fprintf(stderr,"********** alloc\n");
  for (i[2]=0;i[2]<n0[2];i[2]++) {
    fprintf(stderr,"slice %d\n",i[2]);
    for (i[1]=0;i[1]<n0[1];i[1]++) {
      for (i[0]=0;i[0]<n0[0];i[0]++) {
	fprintf(stderr,"%+11.3e ",arr._s3[i[2]][i[1]][i[0]]);
      }
      fprintf(stderr,"\n");
    }
  }
  fprintf(stderr,"********** comp\n");
  for (i[2]=s[2];i[2]<=e[2];i[2]++) {
    fprintf(stderr,"slice %d\n",i[2]);
    for (i[1]=s[1];i[1]<=e[1];i[1]++) {
      for (i[0]=s[0];i[0]<=e[0];i[0]++) {
	fprintf(stderr,"%+11.3e ",arr._s3[i[2]][i[1]][i[0]]);
      }
      fprintf(stderr,"\n");
    }
  }
  fprintf(stderr,"********** test out-of-bounds indices\n");
  fprintf(stderr,"s3[%d][%d][%d]=%+11.3e\n",s[2]-1,s[1],s[0],arr._s3[s[2]-1][s[1]][s[0]]);
  fprintf(stderr,"s3[%d][%d][%d]=%+11.3e\n",s[2],s[1]-1,s[0],arr._s3[s[2]][s[1]-1][s[0]]);
  fprintf(stderr,"s3[%d][%d][%d]=%+11.3e\n",s[2],s[1],s[0]-1,arr._s3[s[2]][s[1]][s[0]-1]);
  fprintf(stderr,"s3[%d][%d][%d]=%+11.3e\n",s[2]-1,s[1]-1,s[0]-1,arr._s3[s[2]-1][s[1]-1][s[0]-1]);
  fprintf(stderr,"s3[%d][%d][%d]=%+11.3e\n",e[2]+1,e[1],e[0],arr._s3[e[2]+1][e[1]][e[0]]);
  fprintf(stderr,"s3[%d][%d][%d]=%+11.3e\n",e[2],e[1]+1,e[0],arr._s3[e[2]][e[1]+1][e[0]]);
  fprintf(stderr,"s3[%d][%d][%d]=%+11.3e\n",e[2],e[1],e[0]+1,arr._s3[e[2]][e[1]][e[0]+1]);
  fprintf(stderr,"s3[%d][%d][%d]=%+11.3e\n",e[2]+1,e[1]+1,e[0]+1,arr._s3[e[2]+1][e[1]+1][e[0]+1]);

  ra_destroy(&arr);

  fprintf(stderr,"exit\n");
}
