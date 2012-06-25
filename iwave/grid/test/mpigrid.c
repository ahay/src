#ifdef IWAVE_USE_MPI
#include "mpigridio.h"
#else
#include "gridio.h"
#endif
#include "parser.h"
#include <time.h>

#define N1 10
#define N2 10
#define N3 1
#define MPI_N1 2
#define MPI_N2 2
#define MPI_N3 1
#define DIM 2

#define VERBOSE

int main(int argc, char **argv) {
  
  FILE * fp = NULL;
  FILE * fpd = NULL;
  char * fname;
  char * dname;
  char * new_fname;
  char * new_dname;
  IPNT gs,n,ran,rags;
  int dim;
  int ndata;
  int ii;
  ireal * a;
  ireal * b;
  ireal ferr = REAL_ZERO;
#ifdef IWAVE_USE_MPI
  int mpi_n1=MPI_N1;
  int mpi_n2=MPI_N2; 
  int mpi_n3=MPI_N3;
#else
  // serial defaults
  int mpi_n1=1;
  int mpi_n2=1; 
  int mpi_n3=1;
#endif
  int wrank=0;

#ifdef IWAVE_USE_MPI
  double mpitime0;
  double mpitime1;
#endif

  PARARRAY * hdrpar;

#ifdef IWAVE_USE_MPI
  MPI_Init(&argc, &argv);
  MPI_Comm_rank(MPI_COMM_WORLD, &wrank);
  storeComm(MPI_COMM_WORLD);
  storeRank(wrank);
  mpitime0=MPI_Wtime();
#endif

  fname = (char *)usermalloc_(sizeof(char)*(strlen("test2d.rsf")+1));
  strcpy(fname,"test2d.rsf");
  dname = (char *)usermalloc_(sizeof(char)*(strlen("test2d.rsf@")+1));
  strcpy(dname,"test2d.rsf@");

  if (wrank==0) {
    fp=iwave_fopen(&fname,"w",NULL,stderr);
    fprintf(fp,"n1=%d\nn2=%d\nd1=1\nd2=1\no1=0\no2=10\ndata_format=native_float\nscale=3\nin=test2d.rsf@",N1,N2);
    fflush(fp);
    iwave_fclose(fp);
    fpd=iwave_fopen(&dname,"w",NULL,stderr);
    a = (ireal *) usermalloc_(N1*N2*N3*sizeof(ireal));
    
    for (ii=0;ii<N1*N2*N3;ii++) {
      a[ii]=(ireal)ii;
    }
    fwrite(a,sizeof(ireal),N1*N2*N3,fpd);
    fflush(fpd);
    iwave_fclose(fpd);
  }

  IASN(gs,IPNT_0);
  n[0]=N1;
  n[1]=N2;
  n[2]=N3;
  dim=DIM;

  if (dim == 2) {
    for (ii=0;ii<RARR_MAX_NDIM;ii++) {
      rags[ii]=gs[ii];
      ran[ii]=n[ii];
    }
    ran[dim-2]=n[dim-2]/mpi_n1;
    ran[dim-1]=n[dim-1]/mpi_n2;
    rags[dim-2]=wrank%mpi_n1*(n[dim-2]/mpi_n1);
    rags[dim-1]=wrank/mpi_n1*(n[dim-1]/mpi_n2);
    if (wrank%mpi_n1 == mpi_n1-1)  ran[dim-2]=n[dim-2]-(mpi_n1-1)*(n[dim-2]/mpi_n1);
    if (wrank/mpi_n1 == mpi_n2-1)  ran[dim-1]=n[dim-1]-(mpi_n2-1)*(n[dim-1]/mpi_n2);
    /*
    fprintf(stderr,"rank: %d\n", wrank);
    for (ii=0;ii<RARR_MAX_NDIM;ii++) {
      fprintf(stderr,"rags[%d]=%d\n",ii,rags[ii]);
    }
    for (ii=0;ii<RARR_MAX_NDIM;ii++) {
      fprintf(stderr,"ran[%d]=%d\n",ii,ran[ii]);
    }
    */
  }
  else if (dim == 3) {
    ran[0]=n[0]/mpi_n1;
    ran[1]=n[1]/mpi_n2;
    ran[2]=n[2]/mpi_n3;
    rags[0]=wrank%mpi_n1*(n[0]/mpi_n1)+gs[0];
    rags[1]=(wrank/mpi_n1)%mpi_n2*(n[1]/mpi_n2)+gs[1];
    rags[2]=(wrank/mpi_n1/mpi_n2)*(n[2]/mpi_n3)+gs[2];
    if (wrank%mpi_n1 == mpi_n1-1)  
      ran[0]=n[0]-(mpi_n1-1)*(n[0]/mpi_n1);
    if ((wrank/mpi_n1)%mpi_n2 == mpi_n2-1)  
      ran[1]=n[1]-(mpi_n2-1)*(n[1]/mpi_n2);
    if (wrank/mpi_n1/mpi_n2 == mpi_n3-1) 
      ran[2]=n[2]-(mpi_n3-1)*(n[2]/mpi_n3);
    /*
    fprintf(stderr,"rank: %d\n", wrank);
    for (ii=0;ii<RARR_MAX_NDIM;ii++) {
      fprintf(stderr,"rags[%d]=%d\n",ii,rags[ii]);
    }
    for (ii=0;ii<RARR_MAX_NDIM;ii++) {
      fprintf(stderr,"ran[%d]=%d\n",ii,ran[ii]);
    }
    */
  }
  else {
#ifdef VERBOSE
    fprintf(stderr,"Error: dim must be 2 or 3\n");
#endif
    userfree_(a);
    userfree_(fname);
    userfree_(dname);
    exit(1);
  }

  if (rsfread(a,rags,ran,"test2d.rsf",0,stderr,0)) {
#ifdef VERBOSE
    fprintf(stderr,"Error: from rsfread\n");
#endif
    userfree_(a);
    userfree_(fname);
    userfree_(dname);
    return 1;
  }

  // write data to new file 
  new_fname = (char *)usermalloc_(sizeof(char)*(strlen("test2d.rsf")+10));
  new_dname = (char *)usermalloc_(sizeof(char)*(strlen("test2d.rsf")+11));
  sprintf(new_fname, "proc%d", wrank);
  strcat(new_fname,"test2d.rsf");
  strcpy(new_dname,new_fname);
  strcat(new_dname,"@");
  
  if (!(hdrpar=ps_new())) {
#ifdef VERBOSE
    fprintf(stderr,"Error: failed to create new par ptr\n");
#endif
    userfree_(a);
    userfree_(fname);
    userfree_(dname);
    userfree_(new_fname);
    userfree_(new_dname);
    exit(1);
  }

  ps_slint(*hdrpar,"n1",N1);
  ps_slint(*hdrpar,"n2",N2);
  ps_slint(*hdrpar,"n3",N3);
  ps_slfloat(*hdrpar,"d1",1.25);
  ps_slfloat(*hdrpar,"d2",1.25);
  ps_slfloat(*hdrpar,"d3",1.25);
  ps_slfloat(*hdrpar,"o1",0.0);
  ps_slfloat(*hdrpar,"o2",0.0);
  ps_slfloat(*hdrpar,"o3",0.0);
  ps_slcstring(*hdrpar,"data_format","native_float");
  ps_slint(*hdrpar,"scale",3);
  ps_slcstring(*hdrpar,"in",new_dname);

  if (!(fp = iwave_fopen(&new_fname,"w",NULL,stderr))) {
#ifdef VERBOSE
    fprintf(stderr,"Error: from iwave_fopen\n");
#endif
    ps_delete(&hdrpar);
    userfree_(a);
    userfree_(fname);
    userfree_(dname);
    userfree_(new_fname);
    userfree_(new_dname);
    exit(1);
  }
  ps_printall(*hdrpar,fp);
  fflush(fp);
  iwave_fclose(fp);

  // open data file so that write works
  if (!(fpd = iwave_fopen(&new_dname,"w",dname,stderr))) {
#ifdef VERBOSE
    fprintf(stderr,"Error: from iwave_fopen\n");
#endif
    ps_delete(&hdrpar);
    userfree_(a);
    userfree_(fname);
    userfree_(dname);
    userfree_(new_fname);
    userfree_(new_dname);
    exit(1);
  }
  iwave_fclose(fpd);

  if (rsfwrite(a,rags,ran,new_fname,0,stderr,0)) {
#ifdef VERBOSE
    fprintf(stderr,"Error: from rsfwrite\n");
#endif
    ps_delete(&hdrpar);
    userfree_(a);
    userfree_(fname);
    userfree_(dname);
    userfree_(new_fname);
    userfree_(new_dname);
    exit(1);
  }

  // now read it back in again
  b = (ireal *) usermalloc_(N1*N2*N3*sizeof(ireal));
  if (rsfread(b,rags,ran,new_fname,0,stderr,0)) {
#ifdef VERBOSE
    fprintf(stderr,"Error: from rsfread, 2nd time\n");
#endif
    ps_delete(&hdrpar);
    userfree_(a);
    userfree_(b);
    userfree_(fname);
    userfree_(dname);
    userfree_(new_fname);
    userfree_(new_dname);
    exit(1);
  }

  // check results
  ferr=REAL_ZERO;
  ndata=ran[0]*ran[1]*ran[2];
  for (ii=0;ii<ndata;ii++) 
    ferr = iwave_max(ferr,iwave_abs(a[ii]-b[ii]));
  if (ferr<10*REAL_EPS) printf("Gridio Test: data same after write, re-read - error less than 10 * macheps\n");
  else printf("Gridio Test: data differ after write, re-read - error greater than 10 * macheps\n");

#ifdef IWAVE_USE_MPI
  mpitime1=MPI_Wtime();
  fprintf(stderr,"#%d time consumption: %7.5f sec\n",wrank,(mpitime1-mpitime0));
#endif
  userfree_(a);
  userfree_(b);
  userfree_(fname);
  userfree_(dname);
  userfree_(new_fname);
  userfree_(new_dname);
  ps_delete(&hdrpar);

#ifdef IWAVE_USE_MPI
  MPI_Finalize();
#endif
  exit(0);
}
