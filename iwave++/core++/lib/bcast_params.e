#include "bcast_params.h"

int bcast_params(PARARRAY * par, FILE * stream) {

#ifdef IWAVE_USE_MPI
  int rk;
  int i;
  int ioff;
  PARARRAY newpar;
  MPI_Comm cml = retrieveComm();    /* local communicator */
  ps_setnull(&newpar);

  if (par->buffersize==0) return 0;
  /* get rank */
  MPI_Comm_rank(cml,&rk);
  //  MPI_Comm_rank(MPI_COMM_WORLD,&rk);
  /* set buffersize */
  if (rk==0) newpar.buffersize=par->buffersize;
  /* broadcast buffersize */
  //  if (!(MPI_SUCCESS==MPI_Bcast(&(newpar.buffersize),1,MPI_INT,0,MPI_COMM_WORLD))) 
  if (!(MPI_SUCCESS==MPI_Bcast(&(newpar.buffersize),1,MPI_INT,0,cml))) 
    return E_OTHER;
  /* allocated buffer */
  if (!(newpar.buffer=(char *)malloc(newpar.buffersize * sizeof(char)))) 
    return E_ALLOC;
  /* copy buffer on root */
  if (rk==0)
    memcpy(newpar.buffer,par->buffer,par->buffersize*sizeof(char));
  /* broadcast buffer */
  //  if (!(MPI_SUCCESS==MPI_Bcast(newpar.buffer,newpar.buffersize,MPI_CHAR,0,MPI_COMM_WORLD))) {
  if (!(MPI_SUCCESS==MPI_Bcast(newpar.buffer,newpar.buffersize,MPI_CHAR,0,cml))) {
    free(newpar.buffer);
    return E_ALLOC;
  }
  fprintf(stream,"buffer bcast rk=%d\n",rk);

  /* assign number of pairs on root */
  if (rk==0) newpar.npar=par->npar;
  /* broadcast number of pairs */
  // if (!(MPI_SUCCESS==MPI_Bcast(&(newpar.npar),1,MPI_INT,0,MPI_COMM_WORLD))) {
  if (!(MPI_SUCCESS==MPI_Bcast(&(newpar.npar),1,MPI_INT,0,cml))) {
    free(newpar.buffer);
    return E_OTHER;
  }
  /* allocated buffer */
  /* allocate sized string arrays */
  if (!(newpar.ns=(SIZEDSTRING *)malloc(2L* newpar.npar * sizeof(SIZEDSTRING)))) {
    free(newpar.buffer);
    return E_ALLOC;
  }
  newpar.vs=newpar.ns+newpar.npar;

  /* foreach sized string */
  for (i=0;i<2*newpar.npar;++i) {
    /* key: copy length, offset on root */
    if (rk==0) {
      newpar.ns[i].n=par->ns[i].n;
      ioff=(intptr_t)(par->ns[i].s-par->buffer);
    }
    /* key: broadcast length */
    //    if (!(MPI_SUCCESS==MPI_Bcast(&(newpar.ns[i].n),1,MPI_INT,0,MPI_COMM_WORLD))) {
    if (!(MPI_SUCCESS==MPI_Bcast(&(newpar.ns[i].n),1,MPI_INT,0,cml))) {
      ps_destroy(&newpar);
      return E_ALLOC;
    }
    /* key: broadcast offset */
    //  if (!(MPI_SUCCESS==MPI_Bcast(&ioff,1,MPI_INT,0,MPI_COMM_WORLD))) {
    if (!(MPI_SUCCESS==MPI_Bcast(&ioff,1,MPI_INT,0,cml))) {
      ps_destroy(&newpar);
      return E_ALLOC;
    }
    /* key: set pointer into buffer */
    newpar.ns[i].s=newpar.buffer+ioff;
  }
  /*
  fprintf(stream,"CREATED IN BCAST rk=%d\n",rk);
  ps_printall(newpar,stream);
  */

  /* copy over input */
  ps_copy(newpar,par);

  /* trash workspace */
  ps_destroy(&newpar);
#endif
  return 0;
}
