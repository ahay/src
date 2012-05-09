#include "pio.h"

FILE * pio_fopen(const char * restrict path, const char * restrict mode) {
  FILE * pfp = NULL;
#ifdef IWAVE_USE_MPI
  int rk=0;
  rk = retrieveRank();
  if (rk==PIO_IOPROC) pfp=fopen(path,mode);
#else
  pfp=fopen(path,mode);
#endif
  return pfp;
}

int pio_fclose(FILE * pfp) {
#ifdef IWAVE_USE_MPI
  int rk=0;
#endif
  int err=0;
  if (pio_isopen(pfp)) {
#ifdef IWAVE_USE_MPI
    rk = retrieveRank();
    if (rk==PIO_IOPROC) err=fclose(pfp);
#else
    err=fclose(pfp);
#endif
  }
  return err;
}

int pio_isopen(FILE * pfp) {
  int ver=0;
#ifdef IWAVE_USE_MPI
  MPI_Comm cm = retrieveComm();
  int rk = retrieveRank();
  if (rk==PIO_IOPROC && pfp) ver=1;
  /* if Bcast fails, returns 0 */
  MPI_Bcast(&ver,1,MPI_INT,PIO_IOPROC,cm);
#else
  if (pfp) ver=1;
#endif
  return ver;
}

int pio_fread(void * restrict ptr, 
	      size_t size, 
	      size_t nmemb, 
	      FILE * restrict pstream) {
  int err=0;
#ifdef IWAVE_USE_MPI
  int rk = retrieveRank();
  MPI_Comm cm = retrieveComm();
  
  if (rk==PIO_IOPROC) {
    if (nmemb != fread(ptr,size,nmemb,pstream)) err=E_FILE;
  }
  if (!err) {
    if (MPI_SUCCESS != MPI_Bcast((char*)ptr,nmemb*size,MPI_CHAR,PIO_IOPROC,cm)) 
      err=E_OTHER;
  }
#else
  if (nmemb != fread(ptr,size,nmemb,pstream)) err=E_FILE;
#endif
  if (err) {
    return -1;
  }
  return nmemb;
}  

int pio_fwrite(void * restrict ptr, 
	       size_t size, 
	       size_t nmemb, 
	       FILE * restrict pstream,
	       int src) {
  int err=0;
#ifdef IWAVE_USE_MPI
  int rk = retrieveRank();
  MPI_Comm cm = retrieveComm();
  MPI_Status status;
  int sendtag=src;
  int recvtag=PIO_IOPROC;
  if (MPI_SUCCESS != MPI_Sendrecv((char *)ptr,size*nmemb,MPI_CHAR,PIO_IOPROC,sendtag,
				  (char *)ptr,size*nmemb,MPI_CHAR,src,recvtag,
				  cm,&status)) err=E_OTHER;
  if (!err && (rk==PIO_IOPROC)) {
#endif
    if (nmemb != fwrite(ptr,size,nmemb,pstream)) err=E_FILE;
#ifdef IWAVE_USE_MPI
  }
#endif
  if (err) return -1;
  return nmemb;
}

int pio_fseeko(FILE *stream, off_t offset, int whence) {
  int err=0;
#ifdef IWAVE_USE_MPI
  int rk = retrieveRank();
  if (rk==PIO_IOPROC) {
#endif
    err=fseeko(stream,offset,whence);
#ifdef IWAVE_USE_MPI
  }
#endif
  return err;
}
