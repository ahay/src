/* mpi rsf i/o 
   X.W., 09/12/09

   WWS 09/09: make compilation conditional

   WWS 09.11.09: correct assignment of gs throughout
*/

#include "mpigridio.h"

/* all i/o through world rk0
   #define WC_ONLY
*/

//#define VERBOSE

// UPDATE - if set, rsfwrite scales output, rsfread performs saxpy
// if not set, rsfwrite does not scale, rsfread overwrites
#define UPDATE

#ifdef IWAVE_USE_MPI

/* common interface for serial and parallel rsfread */

int rsfread(ireal * a, 
	    const IPNT rags, 
	    const IPNT ran, 
	    const char * fname, 
	    float scale,
	    FILE * stream,
	    int	panelindex) {
  return mpirsfread(a,rags,ran,fname,scale,stream, panelindex);
}

int mpirsfread(ireal * a, 
               const IPNT rags, 
               const IPNT ran, 
               const char * fname, 
	       float scale,
               FILE * stream,
	       int panelindex) {
  
  /**************************
   * BEGIN DECLARATIONS     *
   **************************/

  /* rags = start indices of target array */
  /* ran  = lengths of target array axes  */

  /* strings */
  char * type;
  char * dname;

  /* workspace */  
  int ii;            /* axis counter */
  off_t j, i=0;      /* offset counters */
  int err=0;         /* error flag */
  grid g;            /* grid workspace, init from file */
  int dim;           /* grid dimension */
  int gdim;          /* grid global dimension */
  IPNT gs;           /* start indices of grid (file) */
  IPNT gsa;          /* start indices of grid intersection */
  IPNT g_gsa;        /* start indices, global */
  IPNT l_gsa;        /* start indices, local  */
  IPNT gea;          /* end   indices of grid intersection */
  IPNT n;            /* lengths of grid (file) axes */
  IPNT na;           /* lengths of grid intersection axes */
  IPNT gl_na;        /* lengths of grid intersection axes, local or global */
  size_t recsize_b;  /* length of 1D chunk to be read (bytes) */
  size_t recsize;    /* length of 1D chunk to be read (ireal) */
  size_t nbr;        /* length of 1D chunk actually read */

  IPNT read_gs;       /* start indices of grid (memory) */
  IPNT read_gn;       /* lengths of grid (memory) axes */ 

  FILE * fp = NULL;
  PARARRAY * par;
  float * fbuf;      /* input buffer for read */

  MPI_Comm wcomm;    /* global communicator */
  MPI_Comm lcomm;    /* local communicator */
  int lsize;         /* size of the local communicator */
  int lrank;         /* rank in the local communicator  */

  /* XDR buffers etc. */
#ifdef SUXDR
  XDR  xdrs;         /* xdr structure */
  int xflag;         /* i/o success flag */
  xflag = 0;         /* To avoid "uninitialized variable" warning */
  char * buf;        /* input buffer for XDR stream */
  buf = NULL;        /* To avoid "uninitialized variable" warning */
#endif

  /* grid offsets - allocated dynamically */
  /*off_t * goffs;*/ /* global offsets - into file data */
  off_t * read_goffs; /* global offsets - into memory data */
  off_t * loffs;     /* local offsets - into RARR data */
  size_t noffs;      /* number of offsets */ 

  off_t file_goffs;    /* global offsets - into file data */

  /* vars used to read extended model*/  
  off_t cur_pos = 0;

  /**************************
   * END DECLARATIONS       *
   **************************/

  fflush(stream);
  /* local variables for mpi info */
  wcomm=retrieveGlobalComm();

#ifdef WC_ONLY
  lrank=retrieveGlobalRank();
  lsize=retrieveGlobalSize();
  lcomm=wcomm;
#else
  lrank=retrieveRank();
  lsize=retrieveSize();
  lcomm=retrieveComm();
#endif

  /*  MPI_Errhandler_set(wcomm,MPI_ERRORS_RETURN);*/

  /* rank 0 - read parameter table, create grid, open files 
   *
   * local (rank 0) variables initialized:
   * 
   * par         = parameter table
   * g           = grid
   * gs          = global grid start = ALWAYS = 0 (IPNT)
   * dname       = data filename
   * data_format = native_float or xdr_float
   * fp          = data FILE*
   *
   * global variables initialied, to be broadcast or computed
   * on exit from rank 0 scope:
   *
   * dim         = grid dimension
   * n           = axis lengths (IPNT) 
   * read_gs     = grid start for current block
   * read_gn     = axis lengths for read block
   * recsize     = size of read block (words)
   * recsize_b   = size of read block (bytes)
   * scale       = power of 10 used in scaling during read
   */
  if (lrank==0) {

    /* read parameter table from file */
    par = ps_new();
    if (ps_createfile(par,fname)) {
      fprintf(stream,"read_grid: failed to parse file = %s\n",fname);
      fflush(stream);
      err=E_FILE;
    }
    
    /* create grid from parameter table */
    if ( err || (err = par_grid(&g, *par, stream)) ) {
      fprintf(stream,"Error [mpirsfread]: read from read_grid\n");
      ps_delete(&par);
    }

    if (!err)  {
      /* out 02.14 WWS
	 scale=0;
	 ps_flint(*par,"scale",&scale);
      */
      /* get global array params from grid
       */
      
      get_n(n, g);
      dim = get_dimension_grid(g); // WWS 23.03.15: changed from g.dim
      gdim = g.gdim;
      
    }

    /* an error on rk 0 is fatal here */
    if (err) {
      fflush(stream);
      MPI_Abort(wcomm,err);
    }

  }

  /* next: broadcast grid info, allocated buffers */
  // scale is now a float - WWS 02.14
  if (lsize>1) {
    err = (MPI_Bcast(&dim,1,MPI_INT,0,lcomm)) ||
      (MPI_Bcast(&gdim,1,MPI_INT,0,lcomm)) ||
      (MPI_Bcast(n,RARR_MAX_NDIM,MPI_INT,0,lcomm)) ||
      //      (MPI_Bcast(&scale,1,MPI_INT,0,lcomm));
      (MPI_Bcast(&scale,1,MPI_FLOAT,0,lcomm));
    /*  mpi_err(err,wcomm,stream);*/
    if (err) {
      fflush(stream);
      MPI_Abort(wcomm,err);
    }
  }

  /* initialize global gs, gs, gn for read block. NOTE: global gs is
     ALWAYS IPNT_0 - it is the index array of the grid origin
  */
  IASN(gs,IPNT_0);
  IASN(read_gn,n);
  for (ii=dim;ii<gdim;ii++) read_gn[ii]=1;
  // 09.01.14: only reduce last axis if dim > 2
  if (dim>2) 
    read_gn[dim-1]=iwave_min(N_SLOWEST_AXIS_SLICE, n[dim-1]);
  
  /* start indices begin at 0 */
  IASN(read_gs, IPNT_0);
  
  /* initial record size in words & bytes */
  recsize = read_gn[0]*read_gn[1]*read_gn[2];
  /* set up for float types only - any extension to 
     other types will need to move this assignment inside
     the rank=0 branch and add a bcast */
  recsize_b = recsize*sizeof(float);

#ifdef VERBOSE
  fprintf(stream,"mpirsfread: initialization\n");
  fprintf(stream,"n[0]=%d n[1]=%d n[2]=%d\n",n[0],n[1],n[2]);
  fprintf(stream,"rags[0]=%d rags[1]=%d rags[2]=%d\n",rags[0],rags[1],rags[2]);
  fprintf(stream,"ran[0]=%d ran[1]=%d ran[2]=%d\n",ran[0],ran[1],ran[2]);
  fprintf(stream,"read_gn[0]=%d read_gn[1]=%d read_gn[2]=%d\n", 
	  read_gn[0],read_gn[1],read_gn[2]);
  fprintf(stream,"recsize=%ld recsize_b=%ld\n",recsize,recsize_b);
#endif

  /* read buffer - will always be adequate, as record 
     size will never increase 
  */
  fbuf=(float *)usermalloc_(recsize_b);
  if (!fbuf) {
    fprintf(stream,"Fatal Error [mpirsfread]:\n");
    fprintf(stream,"  failed to allocate %ld bytes on process %d\n",recsize_b,lrank);
    fflush(stream);
    MPI_Abort(wcomm,E_ALLOC);
  }

  /* prepare binary input on rk 0 */
  if (lrank==0) {

    /* get filename */
    if ( (err=ps_flcstring(*par,"in",&dname)) ) {
      fprintf(stream,"Error [mpirsfread]: read from ps_flcstring\n");
      fprintf(stream,"failed to extract in param\n");
    }

    if ( err || (err=ps_flcstring(*par,"data_format",&type)) ) {
      fprintf(stream,"Error [mpirsfread]: read from ps_flcstring\n");
      fprintf(stream,"failed to extract type param\n");
    }

    /* sanity test - float data only for this implementation */
    if (strcmp(type,"xdr_float") &&
	strcmp(type,"native_float")) {
      fprintf(stream,"Fatal Error [mpirsfread]:\n");
      fprintf(stream,"  type=%s not float type - this version only accepts float data from file\n",type);
      fflush(stream);
      MPI_Abort(wcomm,E_BADINPUT);
    }
    
    /* now finished with PARARRAY, trash it */
    ps_delete(&par);

    /* open file */
    if (!err) {
      fp=NULL;
#ifdef IWAVE_USE_FMGR
      fp=iwave_const_fopen(dname,"r",NULL,stream);
#else
      fp=fopen(dname,"r");
#endif
      /* done with filename */
      userfree_(dname);
      if (!fp) {
	fprintf(stream,"Error [mpirsfread]: read from fopen\n");
	err=E_FILE;
      }
    }

    if (!err) {
      /* initial file offset */
      file_goffs = 0;
           
      /* first sanity check */
      if (panelindex < 0) {
	fprintf(stream,"Error: rsfread\n");
	fprintf(stream,"panelindex = %d < 0\n",panelindex);
	return E_OTHER;
      }
      
      /* seek to panel at input panelindex modulo number of panels */
      panelindex = panelindex % get_panelnum_grid(g);
      /* compute current starting position */
      cur_pos = panelindex * get_datasize_grid(g) * sizeof(float);
      
#ifdef VERBOSE
      fprintf(stream, "--- mpirsfread: current file cursor : %lld , moving to cur_pos: %lld \n", ftello(fp), cur_pos);
#endif
      fseeko(fp,cur_pos,SEEK_SET); 
      /*<-- D.S. 01.01.11: extended-model related */

#ifdef SUXDR
      if (!strcmp(type,"xdr_float")) {
	/* allocate char buffer - again, initial value of recsize_b is always enough */
	buf=usermalloc_(recsize_b);
	if (buf) {
	  /* (re)initialize xdr stream */
	  xdrmem_create(&xdrs,buf,recsize_b,XDR_DECODE);
	}
	else {
	  fprintf(stream,"Error [mpirsfread]: failed to allocate %ld bytes to XDR buffer\n",recsize_b);
	  err=E_ALLOC;
	}
      }
#endif

    }

    /* an error on rk 0 is fatal here */
    if (err) {
      fflush(stream);
      MPI_Abort(wcomm,err);
    }
  }

  /**********************************************************
   * READ LOOP - until updated start index exceeds axis len
   **********************************************************/

  while (read_gs[dim-1] < n[dim-1] ) {


    /*
     * read data on rk 0
     * broadcast data
     * compute intersection of grids
     * insert data
     * update read_gs, read_gn
     */

    if (lrank==0) {

      file_goffs=ftello(fp);

      if (!strcmp(type,"native_float")) {
	
	nbr = fread(fbuf,sizeof(float), recsize, fp);
	if (nbr != recsize) {
	  fprintf(stream,"Error [mpirsfread]: read from fread at file offset %lld\n", file_goffs);
	  fprintf(stream,"       attempted to read %ld floats\n",recsize);
	  fprintf(stream,"       succeeded in reading %ld floats\n",nbr);
	  fprintf(stream,"       index of %d-diml slice = %d\n",dim-1,read_gs[dim-1]);
	  fprintf(stream,"       length of axis %d = %d\n",dim-1,n[dim-1]);
	  err=E_FILE;
	}
#ifdef VERBOSE
	fprintf(stream,"Note [mpirsfread]: read from fread at file offset %lld\n", file_goffs);
	fprintf(stream,"       attempted to read %ld floats\n",recsize);
	fprintf(stream,"       succeeded in reading %ld floats\n",nbr);
	fprintf(stream,"       index of %d-diml slice = %d\n",dim-1,read_gs[dim-1]);
	fprintf(stream,"       length of axis %d = %d\n",dim-1,n[dim-1]);
	for (int iii=0;iii<recsize;iii++) fprintf(stream,"i=%d f=%e\n",iii,fbuf[iii]);
#endif
      
      }

#ifdef SUXDR
      
      else if (!strcmp(type,"xdr_float")) {
	
	/* reset xdr stream */
	xdr_setpos(&xdrs,0);
	/* read in byte array */
	nbr = fread(buf,sizeof(char),recsize_b,fp);
	if (nbr != recsize_b) {
	  fprintf(stream,"Error [mpirsfread]: read from fread at file offset %lld\n", file_goffs);
	  fprintf(stream,"       attempted to read %ld bytes\n",recsize_b);
	  fprintf(stream,"       succeeded in reading %ld bytes\n",nbr);
	  fprintf(stream,"       index of %d-diml slice = %d\n",dim-1,read_gs[dim-1]);
	  fprintf(stream,"       length of axis %d = %d\n",dim-1,n[dim-1]);
	  err=E_FILE;
	}        
	/* convert to ireal - need only one word of float storage */ 
	if (!err) {
	  xflag=1;
	  for (j=0;j<recsize;j++) {
	    xflag=xflag && xdr_float(&xdrs,fbuf+j);
	  }
	}
	if (!err && !xflag) {
	  fprintf(stream,"Error [mpirsfread]: read - failed xdr conversion\n");
	  err=E_OTHER;
	}
      }
      else {
	fprintf(stream,"Error [mpirsfread]: read - data format must be either\n");
	fprintf(stream,"native_float or xdr_float\n");
	err = E_OTHER;
      }
    
#else
      else {
	fprintf(stream,"Error [mpirsfread]: read - data format must be native_float\n");
	err= E_OTHER;
      }
#endif

      /* an error on rk 0 is fatal here */
      if (err) {
	fflush(stream);
	MPI_Abort(wcomm,err);
      }
    }
  
    /* broadcast data */
    if (lsize > 1) {
      err=MPI_Bcast(fbuf,recsize,IWAVE_MPI_REAL,0,lcomm);
      if (err != MPI_SUCCESS) {
	fprintf(stream,"Error [mpirsfread]: MPI_Bcast error #%d. ABORT\n", lrank);
	fflush(stream);
	MPI_Abort(wcomm,err);
      }
    }

    IASN(g_gsa,IPNT_0);
    IASN(l_gsa,IPNT_0);
    IASN(gl_na,IPNT_0);

    /* compute grid intersections */
    for (ii=0; ii < dim; ii++) {
      gsa[ii] = iwave_max(read_gs[ii], rags[ii]);
      gea[ii] = iwave_min(read_gs[ii] + read_gn[ii] - 1, rags[ii] + ran[ii] - 1);
      na[ii]  = iwave_max(gea[ii] - gsa[ii] + 1, 0);
      /* intersection nonempty */
      if ((rags[ii] + ran[ii] - 1 >= read_gs[ii]) &&
	  (rags[ii] <= read_gs[ii]+read_gn[ii]-1)) {
	g_gsa[ii] = gsa[ii];
	l_gsa[ii] = gsa[ii];
	gl_na[ii] = na[ii];
      }
    }
#ifdef VERBOSE
    for (ii=0;ii<dim;ii++) {
      fprintf(stream,"#%d: gsa[%d]=%d, gea[%d]=%d, na[%d]=%d, rags[%d]=%d, ran[%d]=%d \n",
	      lrank,ii,gsa[ii],ii,gea[ii],ii,na[ii],ii,rags[ii],ii,ran[ii]);
      fprintf(stream,"#%d: read_gs[%d]=%d, read_gn[%d]=%d, g_gsa[%d]=%d, l_gsa[%d]=%d, gl_na[%d]=%d\n", 
	      lrank,ii,read_gs[ii],ii,read_gn[ii],ii,g_gsa[ii],ii,l_gsa[ii],ii,gl_na[ii]);
    }
#endif
    
    if ( (err = get_array_offsets(&read_goffs, &noffs, dim, read_gs, read_gn, g_gsa, gl_na)) ||
	 (err = get_array_offsets(&loffs, &noffs, dim, rags, ran, l_gsa, gl_na)) ) {
      if (err == E_OUTOFBOUNDS) 
	fprintf(stream,"NOTE [mpirsfread]: grid intersection empty\n");
      else 
	fprintf(stream,"Error [mpirsfread]: read from get_array_offsets, err=%d\n",err);
      fflush(stream);
      MPI_Abort(wcomm,err);
    }

    /* load data into subgrid */
    if (noffs > 0) {
      for (i=0;i<noffs;i++) {
#ifdef UPDATE
	for (j=0;j<gl_na[0];j++) { *(a+loffs[i]+j)+=scale*(*(fbuf+read_goffs[i]+j)); }
#else
	for (j=0;j<gl_na[0];j++) { *(a+loffs[i]+j)=*(fbuf+read_goffs[i]+j); }
#endif
      }

      if (read_goffs) userfree_(read_goffs);
      if (loffs) userfree_(loffs);     
    }

    /* update start, size of read block */
    read_gs[dim-1] += read_gn[dim-1];
    read_gn[dim-1] = iwave_min(N_SLOWEST_AXIS_SLICE, n[dim-1]-read_gs[dim-1]);

    /* update record size in words & bytes */
    recsize = read_gn[0]*read_gn[1]*read_gn[2];
    recsize_b = recsize*sizeof(float);
#ifdef VERBOSE
    if (lrank==0) {
      file_goffs=ftello(fp);
      fprintf(stream,"read_gs[%d] updated to %d\n",dim-1,read_gs[dim-1]);
      fprintf(stream,"read_gn[%d] updated to %d\n",dim-1,read_gn[dim-1]);
      fprintf(stream,"recsize updated to %ld\n",recsize);
      fprintf(stream,"recsize_b updated to %ld\n",recsize_b);
      fprintf(stream,"file_goffs updated to %lld\n",file_goffs);
    }
#endif
  }  

  /* cleanup - close file, delete read buffer */

  userfree_(fbuf);
  if (lrank==0) {
#ifdef IWAVE_USE_FMGR
    iwave_fclose(fp);
#else
    fclose(fp);
#endif
#ifdef SUXDR
    if (!strcmp(type,"xdr_float")) {
      /* disengage xdr stream */
      userfree_(buf);
      xdr_destroy(&xdrs); 
    }
#endif
    /* done with type designator */
    userfree_(type);
  }

  return err;
}

int rsfwrite(ireal * a, 
	     const IPNT rags, 
	     const IPNT ran, 
	     const char * fname, 
	     float scale,
	     FILE * stream,
	     int panelindex) {
  return mpirsfwrite(a,rags,ran,fname,scale,stream,panelindex);
}

int mpirsfwrite(ireal * a, 
		const IPNT rags, 
		const IPNT ran, 
		const char * fname, 
		float scale,
		FILE * stream,
		int panelindex) {
  
  /**************************
   * BEGIN DECLARATIONS     *
   **************************/

  /* rags = start indices of target array */
  /* ran  = lengths of target array axes  */

  /* strings */
  char * type;
  char * dname;

  /* workspace */  
  int ii;            /* axis counter */
  off_t j, i=0;      /* offset counters */
  int err=0;         /* error flag */
  grid g;            /* grid workspace, init from file */
  int dim;           /* grid dimension */
  int gdim;          /* grid global dimension */
  IPNT gs;           /* start indices of grid (file) */
  IPNT gsa;          /* start indices of grid intersection */
  IPNT g_gsa;        /* start indices, global */
  IPNT l_gsa;        /* start indices, local  */
  IPNT gea;          /* end   indices of grid intersection */
  IPNT n;            /* lengths of grid (file) axes */
  IPNT na;           /* lengths of grid intersection axes */
  IPNT gl_na;        /* lengths of grid intersection axes, local or global */
  size_t recsize_b;  /* length of 1D chunk to be read (bytes) */
  size_t recsize;    /* length of 1D chunk to be read (ireal) */
  size_t nbr;        /* length of 1D chunk actually read */

  IPNT read_gs;       /* start indices of grid (memory) */
  IPNT read_gn;       /* lengths of grid (memory) axes */ 

  FILE * fp = NULL;
  PARARRAY * par;
  float * fbuf;      /* local output buffer */
  float * rbuf;      /* rank 0 reduction buffer */

  MPI_Comm wcomm;    /* local communicator */
  int wrank;

  /* XDR buffers etc. */
#ifdef SUXDR
  XDR  xdrs;         /* xdr structure */
  int xflag;         /* i/o success flag */
  xflag = 0;         /* To avoid "uninitialized variable" warning */
  char * buf;        /* input buffer for XDR stream */
  buf = NULL;        /* To avoid "uninitialized variable" warning */
#endif

  /* grid offsets - allocated dynamically */
  /*off_t * goffs;*/ /* global offsets - into file data */
  off_t * read_goffs; /* global offsets - into memory data */
  off_t * loffs;     /* local offsets - into RARR data */
  size_t noffs;      /* number of offsets */ 

  off_t file_goffs;    /* global offsets - into file data */

  /* vars used to read extended model*/  
  off_t cur_pos = 0;

  /**************************
   * END DECLARATIONS       *
   **************************/

  fflush(stream);
  /* local variables for mpi info */
  wrank=retrieveRank();
  wcomm=retrieveComm();

  fflush(stream);

  /*  MPI_Errhandler_set(wcomm,MPI_ERRORS_RETURN);*/

  /* rank 0 - read parameter table, create grid, open files 
   *
   * local (rank 0) variables initialized:
   * 
   * par         = parameter table
   * g           = grid
   * gs          = global grid start = ALWAYS = 0 (IPNT)
   * dname       = data filename
   * data_format = native_float or xdr_float
   * fp          = data FILE*
   *
   * global variables initialied, to be broadcast or computed
   * on exit from rank 0 scope:
   *
   * dim         = grid dimension
   * n           = axis lengths (IPNT) 
   * read_gs     = grid start for current block
   * read_gn     = axis lengths for read block
   * recsize     = size of read block (words)
   * recsize_b   = size of read block (bytes)
   */
  if (wrank==0) {

    /* read parameter table from file */
    par=ps_new();
    if (ps_createfile(par,fname)) {
      fprintf(stream,"read_grid: failed to parse file = %s\n",fname);
      fflush(stream);
      err=E_FILE;
    }
    /*    fprintf(stream,"NOTE: par file %s opened in mpirsfwrite = \n",fname); */
    /*    ps_printall(par,stream); */
    
    /* create grid from parameter table */
    if ( err || (err = par_grid(&g, *par, stream)) ) {
      fprintf(stream,"Error [mpirsfwrite]: read from read_grid\n");
    }

    if (!err)  {
      /*
	scale=0;
	ps_flint(*par,"scale",&scale);
      */
      /* get global array params from grid
       */
      
      get_n(n, g);
      dim = get_dimension_grid(g); // WWS 23.03.15: changed from g.dim
      gdim = g.gdim;
      
    }

    /* an error on rk 0 is fatal here */
    if (err) {
      fprintf(stream,"xsfwrite: Error in rank 0 setup\n");
      fflush(stream);
      MPI_Abort(wcomm,err);
    }
  }

  /* next: broadcast grid info, allocated buffers, and use this
     information to compute various global dependent quantities.
  */
  err = (MPI_Bcast(&dim,1,MPI_INT,0,wcomm)) ||
    (MPI_Bcast(&gdim,1,MPI_INT,0,wcomm)) ||
    (MPI_Bcast(n,RARR_MAX_NDIM,MPI_INT,0,wcomm)) ||
    //    (MPI_Bcast(&scale,1,MPI_INT,0,wcomm));
    (MPI_Bcast(&scale,1,MPI_FLOAT,0,wcomm));

  if (err) {
    fprintf(stream,"xsfwrite: Bcast error\n");
    fflush(stream);
    MPI_Abort(wcomm,err);
  }

  /* initialize global gs, gs, gn for read block. NOTE: global gs is
     ALWAYS IPNT_0 - it is the index array of the grid origin
  */
  IASN(gs,IPNT_0);
  IASN(read_gn,n);
  for (ii=dim;ii<gdim;ii++) read_gn[ii]=1;
  // 09.01.14 only subslice for dim > 2
  if (dim>2) read_gn[dim-1]=iwave_min(N_SLOWEST_AXIS_SLICE, n[dim-1]);
  
  /* start indices begin at 0 */
  IASN(read_gs, IPNT_0);
  
  /* initial record size in words & bytes - note that record buffer always consists of floats */
  recsize = read_gn[0]*read_gn[1]*read_gn[2];
  recsize_b = recsize*sizeof(float);
  
  /* read buffer - will always be adequate, as record 
     size will never increase 
  */
  fbuf=(float *)usermalloc_(recsize_b);
  if (!fbuf) {
    fprintf(stream,"Fatal Error [mpirsfwrite]:\n");
    fprintf(stream,"  failed to allocate %ld bytes on process %d\n",recsize_b,wrank);
    fflush(stream);
    MPI_Abort(wcomm,E_ALLOC);
  }
  /* zero out buffer - only part for local process will be updated, then MPI_Reduce
     produces global data on root */
  for (i=0;i<recsize;i++) fbuf[i]=0.0f;

  /* reduction buffer on rank 0 */
  if (wrank==0) {
    rbuf=(float *)usermalloc_(recsize_b);
  }
  else {
    rbuf=NULL;
  }

  /* compute scale factor - reciprocal */
  //ntot=1;
  //for (ii=0;ii<dim;ii++) ntot*=ran[ii];
  //if (scale>0) for (ii=0;ii< scale;ii++) scfac *= 0.10; /* 10.0; */
  //if (scale<0) for (ii=0;ii<-scale;ii++) scfac *= 10.0; /* 0.10; */
  
  /* prepare binary output on rk 0 */
  if (wrank==0) {

    /* get filename */
    if ( (err=ps_flcstring(*par,"in",&dname)) ) {
      fprintf(stream,"Error [mpirsfwrite]: from ps_flcstring\n");
      fprintf(stream,"failed to extract in param\n");
    }

    if ( err || (err=ps_flcstring(*par,"data_format",&type)) ) {
      fprintf(stream,"Error [mpirsfwrite]: from ps_flcstring\n");
      fprintf(stream,"failed to extract type param\n");
    }

    /* now finished with PARARRAY, trash it */
    ps_delete(&par);

    /* open file - since file must exist but needs to be written, use "r+" rather than "w" */
    if (!err) {
      fp=NULL;
#ifdef IWAVE_USE_FMGR
      fp=iwave_const_fopen(dname,"r+",NULL,stream);
#else
      fp=fopen(dname,"r+");
#endif
      if (!fp) {
	fprintf(stream,"Error [mpirsfwrite]: from fopen\n");
	err=E_FILE;
      }
      else {
	/* done with filename */
	userfree_(dname);
      }
    }

    if (!err) {

      /* initial file offset */
      file_goffs = 0;
           
      /* first sanity check */
      if (panelindex < 0) {
	fprintf(stream,"Error: rsfread\n");
	fprintf(stream,"panelindex = %d < 0\n",panelindex);
	return E_OTHER;
      }
      
      /* seek to panel at input panelindex modulo number of panels */
      panelindex = panelindex % get_panelnum_grid(g);
      /* compute current starting position */
      cur_pos = panelindex * get_datasize_grid(g) * sizeof(float);
	
      //      fprintf(stream, "--- mpirsfwrite: rank=%d current file cursor : %lld, moving to cus_pos: %lld\n",wrank, ftello(fp), cur_pos);
      fseeko(fp,cur_pos,SEEK_SET); 

#ifdef SUXDR
      if (!strcmp(type,"xdr_float")) {
	/* allocate char buffer - again, initial value of recsize_b is
	   always enough */
	buf=usermalloc_(recsize_b);
	if (buf) {
	  /* (re)initialize xdr stream */
	  xdrmem_create(&xdrs,buf,recsize_b,XDR_ENCODE);
	}
	else {
	  fprintf(stream,"Error [mpirsfwrite]: failed to allocate %ld bytes to XDR buffer\n",recsize_b);
	  err=E_ALLOC;
	}
      }
#endif

    }

    /* an error on rk 0 is fatal here */
    if (err) { 
      fprintf(stream,"rank 0 error\n");
      fflush(stream);
      MPI_Abort(wcomm,err);
    }
  }

  /**********************************************************
   * WRITE LOOP - until updated start index exceeds axis len
   **********************************************************/

  while (read_gs[dim-1] < n[dim-1]) {

    /*
     * compute intersection of grids
     * insert data
     * reduce data
     * update read_gs, read_gn
     * write data on rk 0
     */

    /* compute grid intersections */
    for (ii=0; ii < dim; ii++) {
      gsa[ii] = iwave_max(read_gs[ii], rags[ii]);
      gea[ii] = iwave_min(read_gs[ii] + read_gn[ii] - 1, rags[ii] + ran[ii] - 1);
      na[ii]  = iwave_max(gea[ii] - gsa[ii] + 1, 0);
      /*      if (extend) { */
      /* rarray to left of garray */
      /*
	if (rags[ii] + ran[ii] - 1 < read_gs[ii]) {
	g_gsa[ii] = read_gs[ii];
	l_gsa[ii] = rags[ii] + ran[ii] - 1;
	if (read_gs[ii] == gs[ii]) { 
	gl_na[ii] = 1;
	}
	else {
	gl_na[ii] = 0;
	}
	}
      */
      /* rarray to right of garray */
      /*
	else if (rags[ii] > read_gs[ii] + read_gn[ii] - 1) {
	g_gsa[ii] = read_gs[ii] + read_gn[ii] - 1;
	l_gsa[ii] = rags[ii];
	if (read_gs[ii]+read_gn[ii] == gs[ii]+n[ii]) {
	gl_na[ii] = 1;
	}
	else{
	gl_na[ii] = 0;
	}
      */      
      /*      } */
      /* intersection nonempty */
      if ((rags[ii] + ran[ii] - 1 >= read_gs[ii]) &&
	  (rags[ii] <= read_gs[ii]+read_gn[ii]-1)) {
	g_gsa[ii] = gsa[ii];
	l_gsa[ii] = gsa[ii];
	gl_na[ii] = na[ii];
      }
    }
#ifdef VERBOSE
    for (ii=0;ii<dim;ii++) {
      fprintf(stream,"#%d: read_gs[%d]=%d, read_gn[%d]=%d, g_gsa[%d]=%d, l_gsa[%d]=%d, gl_na[%d]=%d\n", 
	      wrank,ii,read_gs[ii],ii,read_gn[ii],ii,g_gsa[ii],ii,l_gsa[ii],ii,gl_na[ii]);
    }
#endif
    
    if ( (err = get_array_offsets(&read_goffs, &noffs, dim, read_gs, read_gn, g_gsa, gl_na)) ||
	 (err = get_array_offsets(&loffs, &noffs, dim, rags, ran, l_gsa, gl_na)) ) {
      fprintf(stream,"Error [mpirsfwrite]: from get_array_offsets, err=%d\n",err);
      fflush(stream);
      MPI_Abort(wcomm,err);
    }

    /* adjoint extension loop - adj_extend by const along axes 0,...,dim-1.
       added 09.03.12 to make rsfwrite l2-adjoint to rsfread WWS
    */
    /*
      if (extend) {
      
      for (ii=0;ii<g.dim;ii++) {
      if (adj_extend_array(a,rags,ran,l_gsa,gl_na,g.dim,ii)) {
      fprintf(stream,"Error: rsfwrite - adjoint extension failed at axis %d\n",ii);
      return E_OTHER;
      }
      }
      } 
    */
    /* load data into write block, with scale */
    if (noffs > 0) {
      for (i=0;i<noffs;i++) {
#ifdef UPDATE
	for (j=0;j<gl_na[0];j++) { *(fbuf+read_goffs[i]+j) = scale * (*(a+loffs[i]+j)); }
#else
	for (j=0;j<gl_na[0];j++) { *(fbuf+read_goffs[i]+j) = *(a+loffs[i]+j); }
#endif
      }
      if (read_goffs) userfree_(read_goffs);
      if (loffs) userfree_(loffs);     
    }
    else{
      for (i=0;i<recsize;i++) fbuf[i]=0.0f;
    }

    if(rbuf){
      for (i=0;i<recsize;i++) rbuf[i]=0.0f;
    }
    /* reduce data - in MPI2 could use MPI_IN_PLACE */
    err=MPI_Reduce(fbuf,rbuf,recsize,MPI_FLOAT,MPI_SUM,0,wcomm);
    if (err != MPI_SUCCESS) {
      fprintf(stream,"Error [mpirsfwrite]: MPI_Reduce error #%d. ABORT\n", wrank);
      /*      mpi_err(err,wcomm,stream);*/
      fflush(stream);
      MPI_Abort(wcomm,err);
    }
    
    /* write data on rk 0 */
    if (wrank==0) {

      file_goffs=ftello(fp);

      if (!strcmp(type,"native_float")) {	
	nbr = fwrite(rbuf,sizeof(float), recsize, fp);
	if (nbr != recsize) {
	  fprintf(stream,"Error [mpirsfwrite]: from fwrite at file offset %lld\n", file_goffs);
	  fprintf(stream,"       attempted to write %ld bytes\n",recsize);
	  fprintf(stream,"       succeeded in writing %ld bytes\n",nbr);
	  err=E_FILE;
	}        
	// added 21.03.14 WWS
	fflush(fp);
      }
      
#ifdef SUXDR
      
      else if (!strcmp(type,"xdr_float")) {

	/* reset xdr stream */
	xdr_setpos(&xdrs,0);

	/* convert from ireal - need only one word of float storage */ 
	if (!err) {
	  xflag=1;
	  for (j=0;j<recsize;j++) {
	    xflag=xflag && xdr_float(&xdrs,rbuf+j);
	  }
	}
	if (!err && !xflag) {
	  fprintf(stream,"Error [mpirsfwrite]: read - failed xdr conversion\n");
	  err=E_OTHER;
	}
	
	/* write out byte array */
	nbr = fwrite(buf,sizeof(char),recsize_b,fp);
	if (nbr != recsize_b) {
	  fprintf(stream,"Error [mpirsfwrite]: from fread at file offset %lld\n", file_goffs);
	  fprintf(stream,"       attempted to write %ld bytes\n",recsize_b);
	  fprintf(stream,"       succeeded in writing %ld bytes\n",nbr);
	  err=E_FILE;
	}        
	// added 21.03.14 WWS
	fflush(fp);
      }
      else {
	fprintf(stream,"Error [mpirsfwrite]: data format must be either\n");
	fprintf(stream,"native_float or xdr_float\n");
	err = E_OTHER;
      }
#else
      else {
	fprintf(stream,"Error [mpirsfwrite] - data format must be native_float\n");
	err= E_OTHER;
      }
#endif

      /* an error on rk 0 is fatal here */
      if (err) {
	fflush(stream);
	MPI_Abort(wcomm,err);
      }
    }

    /* update start, size of read block */
    read_gs[dim-1] += read_gn[dim-1];
    read_gn[dim-1] = iwave_min(N_SLOWEST_AXIS_SLICE, n[dim-1]-read_gs[dim-1]);

    /* update record size in words & bytes */
    recsize = read_gn[0]*read_gn[1]*read_gn[2];
    recsize_b = recsize*sizeof(float);

  }    

  /* cleanup - close file, delete read buffer */

  userfree_(fbuf);
  userfree_(rbuf);
  if (wrank==0) {
#ifdef IWAVE_USE_FMGR
    iwave_fclose(fp);
#else
    fclose(fp);
#endif
#ifdef SUXDR
    if (!strcmp(type,"xdr_float")) {
      /* disengage xdr stream */
      userfree_(buf);
      xdr_destroy(&xdrs); 
    }
    userfree_(type);
#endif
  }

  return err;
}
    
#else
  
/* inserted just to avoid compiler warnings  */
void frufru() { fprintf(stderr,"ha ha ha\n"); }

#endif
