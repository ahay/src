#include "mpigridio.h"
#include "parser.h"
#include <time.h>

int main(int argc, char **argv) {
  
  int wrank, wsize;
  PARARRAY par;
  FILE * fp;
  char * parfname;
  char * fname;
  char * data_format;
  char new_fname[64];
  IPNT gs,n,ran,rags;
  grid g, l_g;
  int dim;
  int ii;
  ireal * a;
  int err;
  int mpi_n1, mpi_n2, mpi_n3;
  /*int version=0;*/

  /*xargc=argc; xargv=argv;*/
  /*sleep(2);*/
  MPI_Init(&argc, &argv);
  MPI_Comm_rank(MPI_COMM_WORLD, &wrank);
  MPI_Comm_size(MPI_COMM_WORLD, &wsize);

  double mpitime0=MPI_Wtime();
  double mpitime1;

  storeComm(MPI_COMM_WORLD);
  storeRank(wrank);
  storeSize(wsize);

  fprintf(stderr,"N_SLOWEST_AXIS_SLICE:%d\n",N_SLOWEST_AXIS_SLICE); 
  ps_createargs(&par, argc, argv);
  if (ps_ffcstring(par,"par",&parfname)) {
    fprintf(stderr, "Error: par file is not provided\n");
    exit(1);
  }
  ps_setnull(&par);
  /* read parameter table from file */
  if (ps_createfile(&par,parfname)) {
    fprintf(stderr,"mpigrid: failed to parse file = %s\n",fname);
    return E_FILE;
  }
  ps_ffint(par,"mpi_np1",&mpi_n1);
  ps_ffint(par,"mpi_np2",&mpi_n2);
  ps_ffint(par,"mpi_np3",&mpi_n3);
  
  if (!ps_ffcstring(par,"Data",&fname)) {
    if (!(fp=fopen(fname,"r"))) {
      fprintf(stderr,"Error: mpigrid\n");
      fprintf(stderr,"failed to open existing header file %s\n",fname);
      exit(1);
    }
  }
  else {
    fprintf(stderr, "Error: Data is not provided\n");
    exit(1);
  }

  ps_setnull(&par);
  if (ps_createfile(&par,fname)) {
    fprintf(stderr,"mpigrid: failed to parse file = %s\n",fname);
    return E_FILE;
  }
  if (ps_ffcstring(par,"data_format",&data_format)) {
    fprintf(stderr, "Error: data format is not provided. Default: native_float\n");
    data_format = "native_float";
  }
  if (err = par_grid(&g, par, stderr)) {
    fprintf(stderr,"Error: read from read_grid\n");
    ps_destroy(&par);
    return err;
  }
  IASN(gs,IPNT_0);
  /*get_gs(gs,g);*/
  get_n(n,g);
  fprintf(stderr,"gs[0]:%d  gs[1]:%d  gs[2]:%d\n",
          gs[0],gs[1],gs[2]);
  dim=g.dim;
  fprintf(stderr,"dim: %d\n",dim);
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
    fprintf(stderr,"rank: %d\n", wrank);
    for (ii=0;ii<RARR_MAX_NDIM;ii++) {
      fprintf(stderr,"rags[%d]=%d\n",ii,rags[ii]);
    }
    for (ii=0;ii<RARR_MAX_NDIM;ii++) {
      fprintf(stderr,"ran[%d]=%d\n",ii,ran[ii]);
    }
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
    fprintf(stderr,"rank: %d\n", wrank);
    for (ii=0;ii<RARR_MAX_NDIM;ii++) {
      fprintf(stderr,"rags[%d]=%d\n",ii,rags[ii]);
    }
    for (ii=0;ii<RARR_MAX_NDIM;ii++) {
      fprintf(stderr,"ran[%d]=%d\n",ii,ran[ii]);
    }
  }
  a=(ireal *)malloc(ran[0]*ran[1]*ran[2]*sizeof(ireal));
  /*
  if (!version) {
    mpirsfread(a,rags,ran,fname,1,stderr);
    if (wrank == 0)
      fprintf(stderr,"using mpirfsread\n");
  }
  else {
    mpirsfread_v1(a,rags,ran,fname,1,stderr);
    if (wrank == 0)
      fprintf(stderr,"using mpirfsread_v1\n");
  }*/
  rsfread(a,rags,ran,fname,1,stderr,0);
  
  sprintf(new_fname, "proc%d", wrank);
  strcat(new_fname,fname);
  l_g.dim=dim;
  for(ii=0;ii<RARR_MAX_NDIM;ii++) {
    l_g.axes[ii].n=ran[ii];
    l_g.axes[ii].d=g.axes[ii].d;
    l_g.axes[ii].o=g.axes[ii].o+rags[ii]*g.axes[ii].d;
    l_g.axes[ii].id=g.axes[ii].id;
  }
  rsfwrite(a,rags,ran,new_fname,stderr,0);
  mpitime1=MPI_Wtime();
  fprintf(stderr,"#%d time consumption: %7.5f sec\n",wrank,(mpitime1-mpitime0));
  free(a);
  fclose(fp);
  fflush(fp);
  MPI_Finalize();
  exit(1);
}
