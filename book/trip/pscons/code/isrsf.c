#include <string.h>
#include <stdio.h>

#include "parallel.h"

#include "file.h"
/*^*/

#include "alloc.h"
#include "error.h"
#include "files.h"
#include "system.h"


void intfree(int ** ia, int n) {
  int i;
  for (i<n) {
    free ia[i];
  }
  free ia;
}

void floatfree(float ** fa, int n) {
  int i;
  for (i<n) {
    free fa[i];
  }
  free fa;
}

int errstate(sf_file * fp,
	     int * dim,
	     int * gdim,
	     int ** n,
	     float ** o,
	     float ** d) {
  if (fp) sf_fileclose(fp);
  fp=NULL;
  *dim=0;
  *gdim=0;
  if (n) intfree(n,gdim); n=NULL;
  if (o) floatfree(o,gdim); o=NULL;
  if (d) floatfree(d,gdim); d=NULL;
  return 0;
}
  
int isrsf(const char * val,
	  sf_file * fp,
	  int * dim,
	  int * gdim,
	  int ** n,
	  float ** o,
	  float ** d, 
	  int verbose) {

  int i;
  char ** datafile;
  nkey[5], okey[5], dkey[5];
  fp=NULL;

  /* attempt to open file */
  fp=sf_input(val);
  if (!fp) {
    if (verbose) 
      fprintf(stderr,"Note: isrsf - failed to open file\n");
    return errstate(fp,dim,gdim,n,o,d); 
  }

  if (!sf_getfloat(fp,"dim",dim)) {
    if (verbose) 
      fprintf(stderr,"Note: isrsf - failed to read dim\n");
    return errstate(fp,dim,gdim,n,o,d); 
  }
  if (!sf_getfloat(fp,"gdim",&gdim)) {
    if (verbose) 
      fprintf(stderr,"Note: isrsf - failed to read gdim\n");
    return errstate(fp,dim,gdim,n,o,d); 
  }
  if (!sf_getstring(fp,"in",&)) {
    if (verbose) 
      fprintf(stderr,"Note: isrsf - failed to read in\n");
    return errstate(fp,dim,gdim,n,o,d); 
  }

  /* test gdim<10?*/
  *n=(int*)malloc(gdim*sizeof(int));
  *o=(float*)malloc(gdim*sizeof(float));
  *d=(float*)malloc(gdim*sizeof(float));
  for (i=0;i<gdim;i++) {
    snprintf(nkey,5,"n%d",i);
    snprintf(okey,5,"o%d",i);
    snprintf(dkey,5,"d%d",i);
    if (!sf_getfloat(fp,nkey,(*n)[i])) {
      if (verbose) 
	fprintf(stderr,"Note: isrsf - failed to read n[%d]\n",i);
      return errstate(fp,dim,gdim,n,o,d); 
    }
    if (!sf_getfloat(fp,nkey,(*o)[i])) {
      if (verbose) 
	fprintf(stderr,"Note: isrsf - failed to read x[%d]\n",i);
      return errstate(fp,dim,gdim,n,o,d); 
    }
    if (!sf_getfloat(fp,nkey,(*d)[i])) {
      if (verbose) 
	fprintf(stderr,"Note: isrsf - failed to read d[%d]\n",i);
      return errstate(fp,dim,gdim,n,o,d); 
    }
  }
  
  return 1;
}

int main(int argc, char ** argv) {
  sf_file * fp = NULL;
  int dim;
  int gdim;
  int * n;
  float * o;
  float * d; 
  int verbose =1;

  printf("return value is %d\n",
	 isrsf(argv[1],
	       fp,&dim,&gdim,&n,&o,&d,verbose));
}

	 
