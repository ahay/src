#include <traceio.h>

#define HDRFILE "test/hdr.su"

char ** xargv;

int main(int argc, char ** argv) {

  int err=0;
  tracegeom tg;
  RPNT og;
  IPNT n;
  IPNT axord;
  RPNT o;
  RPNT d;
  int i;
  int ndim=2;
  ireal dt=1.4;
  int usernt=0;
  ireal usert0=0.0;
  int initbuf=0;
  int order=1;

  for (i=0;i<RARR_MAX_NDIM;i++) {
    og[i]=REAL_ZERO;
    o[i]=REAL_ZERO;
    d[i]=20;
    axord[i]=i;
  }
  n[0]=91;
  n[1]=391;
  n[2]=1;

  /* construct trace geometry object */
  if (err=construct_tracegeom(&tg,HDRFILE,HDRFILE,SRC_TOL,stderr)) {
    fprintf(stderr,"Error from construct_tracegeom, err=%d\n",err);
    exit(1);
  }

  if (err=init_tracegeom(&tg,og,n,d,o,axord,order,dt,ndim,usernt,usert0,initbuf,stderr)) {
    fprintf(stderr,"Error from init_tracegeom, err=%d\n",err);
    exit(1);
  }

  print_tracegeom(&tg);

  destroy_tracegeom(&tg);

  iwave_fdestroy();
}
