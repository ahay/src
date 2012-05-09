/* mpi rsf i/o 
X.W., 09/12/09

WWS 09/09: make compilation conditional

WWS 09.11.09: correct assignment of gs throughout
*/

#include "mpigridio.h"

/* all i/o through world rk0
#define WC_ONLY
*/

#ifdef IWAVE_USE_MPI

/*#define VERBOSE*/

/* an MPI error handler from 
 * 
 * http://beige.ucs.indiana.edu/I590/node85.html
 *
 * should probably be in utils. Relies on a prior
 * call to MPI_Errhandler_set - which should probably 
 * be in iwave_construct
 */
/*
void mpi_err(int err, MPI_Comm comm,FILE * stream) {
  char err_str[MPI_MAX_ERROR_STRING];
  int err_len;
  int err_cls;

  if (err != MPI_SUCCESS) {
    MPI_Error_class(err,&err_cls);
    MPI_Error_string(err_cls, err_str, &err_len);
    fprintf(stream,"%s\n",err_str);
    MPI_Error_string(err, err_str, &err_len);
    fprintf(stream,"%s\n",err_str);
    fflush(stream);
    MPI_Abort(comm,err);
  }
}
*/
/* common interface for serial and parallel rsfread */

int rsfread(ireal * a, 
	    IPNT rags, 
	    IPNT ran, 
	    char * fname, 
	    int extend, 
	    FILE * stream,
	    int	panelindex   /* D.S. 01.01.11: extended-model related */
	    ) {
  return mpirsfread(a,rags,ran,fname,extend,stream
		    , panelindex /* D.S. 01.01.11: extended-model related */
		    );
}

int mpirsfread(ireal * a, 
               IPNT rags, 
               IPNT ran, 
               char * fname, 
               int extend, 
               FILE * stream,
	       int panelindex  /* D.S. 01.01.11: extended-model related */
	       ) {
  
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
  IPNT gs;           /* start indices of grid (file) */
  IPNT gsa;          /* start indices of grid intersection */
  IPNT g_gsa;        /* start indices, global */
  IPNT l_gsa;        /* start indices, local  */
  IPNT gea;          /* end   indices of grid intersection */
  IPNT g_gea;        /* end   indices, global */
  IPNT l_gea;        /* end   indices, local */
  IPNT n;            /* lengths of grid (file) axes */
  IPNT na;           /* lengths of grid intersection axes */
  IPNT gl_na;        /* lengths of grid intersection axes, local or global */
  size_t recsize_b;  /* length of 1D chunk to be read (bytes) */
  size_t recsize;    /* length of 1D chunk to be read (ireal) */
  size_t nbr;        /* length of 1D chunk actually read */

  IPNT read_gs;       /* start indices of grid (memory) */
  IPNT read_gn;       /* lengths of grid (memory) axes */ 

  FILE * fp = NULL;
  PARARRAY par;
  float * fbuf;      /* input buffer for read */

  MPI_Comm wcomm;    /* global communicator */
  int wsize;         /* size of the global communicator */
  int wrank;         /* rank in the global communicator  */
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

  /* scale flag - like scalco, scalel in SEGY; factor */
  int scale=0;       /* flag - read from parameters */
  float scfac=1.0;   /* factor workspace */
  size_t ntot;       /* total length of constructed local array */
  /*
  float a_max;
  float a_min;
  */

  /* vars used to read extended model*/  
  int panel_size = 1; 
  off_t cur_pos = 0;

  /**************************
   * END DECLARATIONS       *
   **************************/

  fflush(stream);
  /* local variables for mpi info */
  wrank=retrieveGlobalRank();
  wsize=retrieveGlobalSize();
  wcomm=retrieveGlobalComm();

#ifdef WC_ONLY
  lrank=wrank;
  lsize=wsize;
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

    ps_setnull(&par);
    /* read parameter table from file */
    if (ps_createfile(&par,fname)) {
      fprintf(stream,"read_grid: failed to parse file = %s\n",fname);
      fflush(stream);
      err=E_FILE;
    }
    
    /* create grid from parameter table */
    if ( err || (err = par_grid(&g, par, stream)) ) {
      fprintf(stream,"Error [mpirsfread]: read from read_grid\n");
    }

    if (!err)  {

      scale=0;
      ps_ffint(par,"scale",&scale);
  
      /* get global array params from grid
      */
      
      get_n(n, g);
      dim = g.dim;
      
    }

    /* an error on rk 0 is fatal here */
    if (err) {
      fflush(stream);
      MPI_Abort(wcomm,err);
    }

  }

  /* next: broadcast grid info, allocated buffers */
    
  if (lsize>1) {
    err = (MPI_Bcast(&dim,1,MPI_INT,0,lcomm)) ||
      (MPI_Bcast(n,RARR_MAX_NDIM,MPI_INT,0,lcomm)) ||
      (MPI_Bcast(&scale,1,MPI_INT,0,lcomm));
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
  fprintf(stream,"read_gn[0]=%d read_gn[1]=%d read_gn[2]=%d\n", 
	  read_gn[0],read_gn[1],read_gn[2]);
  fprintf(stream,"recsize=%ld recsize_b=%ld\n",recsize,recsize_b);
#endif

  /* read buffer - will always be adequate, as record 
     size will never increase 
  */
  fbuf=(float *)malloc(recsize_b);
  if (!fbuf) {
    fprintf(stream,"Fatal Error [mpirsfread]:\n");
    fprintf(stream,"  failed to allocate %ld bytes on process %d\n",recsize_b,lrank);
    fflush(stream);
    MPI_Abort(wcomm,E_ALLOC);
  }

  /* prepare binary input on rk 0 */
  if (lrank==0) {

    /* get filename */
    if ( (err=ps_ffcstring(par,"in",&dname)) ) {
      fprintf(stream,"Error [mpirsfread]: read from ps_ffcstring\n");
      fprintf(stream,"failed to extract in param\n");
    }

    if ( err || (err=ps_ffcstring(par,"data_format",&type)) ) {
      fprintf(stream,"Error [mpirsfread]: read from ps_ffcstring\n");
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
    ps_destroy(&par);

    /* open file */
    if (!err) {
      fp=NULL;
#ifdef IWAVE_USE_FMGR
      fp=iwave_const_fopen(dname,"r",NULL,stream);
#else
      fp=fopen(dname,"r");
#endif
      if (!fp) {
	fprintf(stream,"Error [mpirsfread]: read from fopen\n");
	err=E_FILE;
      }
      else {
	/* done with filename */
	free(dname);
      }
    }

    if (!err) {
      /* initial file offset */
      file_goffs = 0;
           
      /* D.S. 01.01.11: extended-model related --> */
      /* get size of every model panel */
      for (ii=0; ii< g.dim; ii++){
	panel_size *= (n[ii]);   
      }

      /* compute current starting position
	 Note: for extended model, don't move cursor to the begining
	       when cur_pos is not less than file size (should 
	       throw an error in this case)*/
      cur_pos = panelindex * panel_size * sizeof(float);
      
#ifdef VERBOSE
      fprintf(stream, "--- mpirsfread: current file cursor : %lld , moving to cus_pos: %lld \n", ftello(fp), cur_pos);
#endif
      fseeko(fp,cur_pos,SEEK_SET); 
      /*<-- D.S. 01.01.11: extended-model related */

#ifdef SUXDR
      if (!strcmp(type,"xdr_float")) {
	/* allocate char buffer - again, initial value of recsize_b is always enough */
	buf=malloc(recsize_b);
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
      else {
	/* if you get this far, you've read a block -
	   update file offset (only on rk 0 of course) */
	/*	file_goffs += recsize_b;*/
      }
    }

    /* broadcast data */
    if (lsize > 1) {
      err=MPI_Bcast(fbuf,recsize,IWAVE_MPI_REAL,0,lcomm);
      if (err != MPI_SUCCESS) {
	fprintf(stream,"Error [mpirsfread]: MPI_Bcast error #%d. ABORT\n", lrank);
	/*      mpi_err(err,wcomm,stream);*/
	fflush(stream);
	MPI_Abort(wcomm,err);
      }
    }
    
    /* compute grid intersections */
    for (ii=0; ii < dim; ii++) {
      gsa[ii] = iwave_max(read_gs[ii], rags[ii]);
      gea[ii] = iwave_min(read_gs[ii] + read_gn[ii] - 1, rags[ii] + ran[ii] - 1);
      na[ii]  = iwave_max(gea[ii] - gsa[ii] + 1, 0);
      if (extend) {
	/* rarray to left of garray */
	if (rags[ii] + ran[ii] - 1 < read_gs[ii]) {
	  g_gsa[ii] = read_gs[ii];
	  g_gea[ii] = read_gs[ii];
	  l_gsa[ii] = rags[ii] + ran[ii] - 1;
	  l_gea[ii] = rags[ii] + ran[ii] - 1;
	  if (read_gs[ii] == gs[ii]) { 
	    gl_na[ii] = 1;
	  }
	  else {
	    gl_na[ii] = 0;
	  }
	}
	/* rarray to right of garray */
	else if (rags[ii] > read_gs[ii] + read_gn[ii] - 1) {
	  g_gsa[ii] = read_gs[ii] + read_gn[ii] - 1;
	  g_gea[ii] = read_gs[ii] + read_gn[ii] - 1;
	  l_gsa[ii] = rags[ii];
	  l_gea[ii] = rags[ii];
	  if (read_gs[ii]+read_gn[ii] == gs[ii]+n[ii]) {
	    gl_na[ii] = 1;
	  }
	  else{
	    gl_na[ii] = 0;
	  }
	}      
      }
      /* intersection nonempty */
      if ((rags[ii] + ran[ii] - 1 >= read_gs[ii]) &&
	  (rags[ii] <= read_gs[ii]+read_gn[ii]-1)) {
	g_gsa[ii] = gsa[ii];
	g_gea[ii] = gea[ii];
	l_gsa[ii] = gsa[ii];
	l_gea[ii] = gea[ii];
	gl_na[ii] = na[ii];
      }
    }
#ifdef VERBOSE
    for (ii=0;ii<dim;ii++) {
      fprintf(stream,"#%d: read_gs[%d]=%d, read_gn[%d]=%d, g_gsa[%d]=%d, l_gsa[%d]=%d, gl_na[%d]=%d\n", 
	      lrank,ii,read_gs[ii],ii,read_gn[ii],ii,g_gsa[ii],ii,l_gsa[ii],ii,gl_na[ii]);
    }
#endif
    
    if ( (err = get_array_offsets(&read_goffs, &noffs, dim, read_gs, read_gn, g_gsa, gl_na)) ||
	 (err = get_array_offsets(&loffs, &noffs, dim, rags, ran, l_gsa, gl_na)) ) {
      fprintf(stream,"Error [mpirsfread]: read from get_array_offsets, err=%d\n",err);
      fflush(stream);
      MPI_Abort(wcomm,err);
    }

    /* load data into subgrid */
    if (noffs > 0) {
      for (i=0;i<noffs;i++) {
	for (j=0;j<gl_na[0];j++) { *(a+loffs[i]+j)=*(fbuf+read_goffs[i]+j); }
      }
      if (read_goffs) free(read_goffs);
      if (loffs) free(loffs);     
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
      fprintf(stream,"file_goffs updated to %ld\n",file_goffs);
    }
      
#endif


  }  

  /* cleanup - close file, delete read buffer */

  free(fbuf);
  if (lrank==0) {
#ifdef IWAVE_USE_FMGR
    iwave_fclose(fp);
#else
    fclose(fp);
#endif
#ifdef SUXDR
    if (!strcmp(type,"xdr_float")) {
      /* disengage xdr stream */
      free(buf);
      xdr_destroy(&xdrs); 
    }
    free(type);
#endif
  }

  /* extension outside of file-defined grid */

  for (ii=0; ii < dim; ii++) {
    gsa[ii] = iwave_max(gs[ii], rags[ii]);
    gea[ii] = iwave_min(gs[ii] + n[ii] - 1, rags[ii] + ran[ii] - 1);
    na[ii]  = iwave_max(gea[ii] - gsa[ii] + 1, 0);
    
    /* rarray to left of garray */
    if ( (rags[ii] + ran[ii] - 1 < gs[ii]) && extend ) { 
      g_gsa[ii] = gs[ii];
      g_gea[ii] = gs[ii];
      l_gsa[ii] = rags[ii] + ran[ii] - 1;
      l_gea[ii] = rags[ii] + ran[ii] - 1;
      gl_na[ii] = 1;
    }
    /* rarray to right of garray */
    else if ( (rags[ii] > gs[ii] + n[ii] -1) && extend ) {
      g_gsa[ii] = gs[ii] + n[ii] - 1;
      g_gea[ii] = gs[ii] + n[ii] - 1;
      l_gsa[ii] = rags[ii];
      l_gea[ii] = rags[ii];
      gl_na[ii] = 1;
    }
    /* intersection nonempty */
    else {
      g_gsa[ii] = gsa[ii];
      g_gea[ii] = gea[ii];
      l_gsa[ii] = gsa[ii];
      l_gea[ii] = gea[ii];
      gl_na[ii] = na[ii];
    }
  } 
  /* extension loop - extend by const along axes dim-1,..,0.
   */
  if (extend) {
    for (ii=dim-1;ii>-1;ii--) {
      if (extend_array(a,rags,ran,l_gsa,gl_na,dim,ii)) {
        fprintf(stream,"Error [mpirsfread]: rsfread: extension failed at axis %d\n",ii);
        return E_OTHER;
      }
    }
  } 
  
  /* scaling loop */
  /*fprintf(stream,"gridio: SCALE=%d\n",scale);*/
  ntot=1;
  for (ii=0;ii<dim;ii++) ntot*=ran[ii];
  if (scale>0) for (ii=0;ii<scale;ii++) scfac *= 10.0;
  if (scale<0) for (ii=0;ii<-scale;ii++) scfac *= 0.10;
  /*
  a_max=scfac*a[0];
  a_min=scfac*a[0];
  */
  if (scale) {
    for (i=0;i<ntot;i++) {
      a[i]*=scfac;
      /*
      a_max=iwave_max(a[i],a_max);
      a_min=iwave_min(a[i],a_min);
      */
    }
  }
  /*  fprintf(stream,"mpigridio: fac=%e max=%e min=%e\n",scfac,a_max,a_min);*/
  
  return err;
}

int rsfwrite(ireal * a, 
	     IPNT rags, 
	     IPNT ran, 
	     char * fname, 
	     int extend,
	     FILE * stream,
	     int panelindex  /* D.S. 01.01.11: extended-model related */
	     ) {
  return mpirsfwrite(a,rags,ran,fname,extend,stream,panelindex);
}

int mpirsfwrite(ireal * a, 
		IPNT rags, 
		IPNT ran, 
		char * fname, 
		int extend,
		FILE * stream,
		int panelindex  /* D.S. 01.01.11: extended-model related */
		) {
  
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
  IPNT gs;           /* start indices of grid (file) */
  IPNT gsa;          /* start indices of grid intersection */
  IPNT g_gsa;        /* start indices, global */
  IPNT l_gsa;        /* start indices, local  */
  IPNT gea;          /* end   indices of grid intersection */
  IPNT g_gea;        /* end   indices, global */
  IPNT l_gea;        /* end   indices, local */
  IPNT n;            /* lengths of grid (file) axes */
  IPNT na;           /* lengths of grid intersection axes */
  IPNT gl_na;        /* lengths of grid intersection axes, local or global */
  size_t recsize_b;  /* length of 1D chunk to be read (bytes) */
  size_t recsize;    /* length of 1D chunk to be read (ireal) */
  size_t nbr;        /* length of 1D chunk actually read */

  IPNT read_gs;       /* start indices of grid (memory) */
  IPNT read_gn;       /* lengths of grid (memory) axes */ 

  FILE * fp = NULL;
  PARARRAY par;
  float * fbuf;      /* local output buffer */
  float * rbuf;      /* rank 0 reduction buffer */

  MPI_Comm wcomm;    /* local communicator */
  int wsize;         /* size of the local communicator */
  int wrank;         /* rank in the local communicator  */

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

  /* scale flag - like scalco, scalel in SEGY; factor */
  int scale=0;       /* flag - read from parameters */
  float scfac=1.0;   /* factor workspace */
  size_t ntot;       /* total length of constructed local array */

  /* vars used to read extended model*/  
  int panel_size = 1; 
  off_t cur_pos = 0;

  /**************************
   * END DECLARATIONS       *
   **************************/

  fflush(stream);
  /* local variables for mpi info */
  wrank=retrieveRank();
  wsize=retrieveSize();
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
   * scale       = power of 10 used in scaling during read
   */
  if (wrank==0) {

    ps_setnull(&par);
    /* read parameter table from file */
    if (ps_createfile(&par,fname)) {
      fprintf(stream,"read_grid: failed to parse file = %s\n",fname);
      fflush(stream);
      err=E_FILE;
    }
    //    fprintf(stream,"NOTE: par file %s opened in mpirsfwrite = \n",fname);
    //    ps_printall(par,stream);
    
    /* create grid from parameter table */
    if ( err || (err = par_grid(&g, par, stream)) ) {
      fprintf(stream,"Error [mpirsfwrite]: read from read_grid\n");
    }

    if (!err)  {

      scale=0;
      ps_ffint(par,"scale",&scale);
  
      /* get global array params from grid
      */
      
      get_n(n, g);
      dim = g.dim;
      
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
    (MPI_Bcast(n,RARR_MAX_NDIM,MPI_INT,0,wcomm)) ||
    (MPI_Bcast(&scale,1,MPI_INT,0,wcomm));

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
  read_gn[dim-1]=iwave_min(N_SLOWEST_AXIS_SLICE, n[dim-1]);
  
  /* start indices begin at 0 */
  IASN(read_gs, IPNT_0);
  
  /* initial record size in words & bytes - note that record buffer always consists of floats */
  recsize = read_gn[0]*read_gn[1]*read_gn[2];
  recsize_b = recsize*sizeof(float);
  
  /* read buffer - will always be adequate, as record 
     size will never increase 
  */
  fbuf=(float *)malloc(recsize_b);
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
    rbuf=(float *)malloc(recsize_b);
  }
  else {
    rbuf=NULL;
  }

  /* compute scale factor - reciprocal */
  ntot=1;
  for (ii=0;ii<dim;ii++) ntot*=ran[ii];
  if (scale>0) for (ii=0;ii< scale;ii++) scfac *= 0.10; // 10.0;
  if (scale<0) for (ii=0;ii<-scale;ii++) scfac *= 10.0; // 0.10;

  /* prepare binary output on rk 0 */
  if (wrank==0) {

    /* get filename */
    if ( (err=ps_ffcstring(par,"in",&dname)) ) {
      fprintf(stream,"Error [mpirsfwrite]: from ps_ffcstring\n");
      fprintf(stream,"failed to extract in param\n");
    }

    if ( err || (err=ps_ffcstring(par,"data_format",&type)) ) {
      fprintf(stream,"Error [mpirsfwrite]: from ps_ffcstring\n");
      fprintf(stream,"failed to extract type param\n");
    }

    /* now finished with PARARRAY, trash it */
    ps_destroy(&par);

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
	free(dname);
      }
    }

    if (!err) {

      /* initial file offset */
      file_goffs = 0;
           
      /* WWS 20.01.11 -- turn off ext model feature by setting panelindex=-1 */

      if (panelindex >-1) {

	/* D.S. 01.01.11: extended-model related --> */
	/* get size of every model panel */
	for (ii=0; ii< g.dim; ii++){
	  panel_size *= (n[ii]);   
	}
	/* compute current starting position
	   Note: for extended model, don't move cursor to the begining
	   when cur_pos is not less than file size (should 
	   throw an error in this case)*/
	cur_pos = panelindex * panel_size * sizeof(float);
	
	fprintf(stream, "--- mpirsfwrite: rank=%d current file cursor : %lld, moving to cus_pos: %lld\n",wrank, ftello(fp), cur_pos);
	fseeko(fp,cur_pos,SEEK_SET); 
	/*<-- D.S. 01.01.11: extended-model related */
	
      }
#ifdef SUXDR
      if (!strcmp(type,"xdr_float")) {
	/* allocate char buffer - again, initial value of recsize_b is always enough */
	buf=malloc(recsize_b);
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
   * READ LOOP - until updated start index exceeds axis len
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
      //      if (extend) {
	/* rarray to left of garray */
	if (rags[ii] + ran[ii] - 1 < read_gs[ii]) {
	  g_gsa[ii] = read_gs[ii];
	  g_gea[ii] = read_gs[ii];
	  l_gsa[ii] = rags[ii] + ran[ii] - 1;
	  l_gea[ii] = rags[ii] + ran[ii] - 1;
	  if (read_gs[ii] == gs[ii]) { 
	    gl_na[ii] = 1;
	  }
	  else {
	    gl_na[ii] = 0;
	  }
	}
	/* rarray to right of garray */
	else if (rags[ii] > read_gs[ii] + read_gn[ii] - 1) {
	  g_gsa[ii] = read_gs[ii] + read_gn[ii] - 1;
	  g_gea[ii] = read_gs[ii] + read_gn[ii] - 1;
	  l_gsa[ii] = rags[ii];
	  l_gea[ii] = rags[ii];
	  if (read_gs[ii]+read_gn[ii] == gs[ii]+n[ii]) {
	    gl_na[ii] = 1;
	  }
	  else{
	    gl_na[ii] = 0;
	  }
	}      
	//      }
      /* intersection nonempty */
      if ((rags[ii] + ran[ii] - 1 >= read_gs[ii]) &&
	  (rags[ii] <= read_gs[ii]+read_gn[ii]-1)) {
	g_gsa[ii] = gsa[ii];
	g_gea[ii] = gea[ii];
	l_gsa[ii] = gsa[ii];
	l_gea[ii] = gea[ii];
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
    if (extend) {
      
      for (ii=0;ii<g.dim;ii++) {
	if (adj_extend_array(a,rags,ran,l_gsa,gl_na,g.dim,ii)) {
	  fprintf(stream,"Error: rsfwrite - adjoint extension failed at axis %d\n",ii);
	  return E_OTHER;
	}
      }
    } 

    /* load data into write block, with scale */
    if (noffs > 0) {
      for (i=0;i<noffs;i++) {
	for (j=0;j<gl_na[0];j++) { *(fbuf+read_goffs[i]+j) = scfac * (*(a+loffs[i]+j)); }
      }
      if (read_goffs) free(read_goffs);
      if (loffs) free(loffs);     
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

  free(fbuf);
  free(rbuf);
  if (wrank==0) {
#ifdef IWAVE_USE_FMGR
    iwave_fclose(fp);
#else
    fclose(fp);
#endif
#ifdef SUXDR
    if (!strcmp(type,"xdr_float")) {
      /* disengage xdr stream */
      free(buf);
      xdr_destroy(&xdrs); 
    }
    free(type);
#endif
  }

  return err;
}
    
#else
  
/* inserted just to avoid compiler warnings  */
void frufru() { fprintf(stderr,"ha ha ha\n"); }

#endif
