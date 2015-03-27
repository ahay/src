#include "gridio.h"

// UPDATE - if set, rsfwrite scales output, rsfread performs saxpy
// if not set, rsfwrite does not scale, rsfread overwrites
#define UPDATE
// if you want to see what's happening...
#define VERBOSE

int read_grid(grid * g, const char * fname, FILE * fp) {

  int err=0;
  PARARRAY * par = ps_new();

  /* parser file for key=value pairs */
  if ((err=ps_createfile(par,fname))) {
    fprintf(fp,"ERROR: read_grid from ps_createfile, err=%d\n",err); 
    fprintf(fp,"failed to parse file = %s\n",fname);
    return E_FILE;
  }

  if ( (err=par_grid(g,*par,fp)) ) {
    fprintf(fp,"ERROR: read_grid from par_grid, err=%d\n",err); 
    return err;
  }

  ps_delete(&par);

  return err;
}

int par_grid(grid * g, PARARRAY par, FILE * fp) {
  
  int i;
  int tmp;
  char key[4];
  size_t kl=4;
  /*
  fprintf(stderr,"in par_grid:\n");
  ps_printall(par,stderr);
  */
  if (IARR_MAX_NDIM > 9) {
    fprintf(fp,"ERROR: par_grid\n");
    fprintf(fp,"  Oooops...grid parser assumes max dim < 10\n");
    fprintf(fp,"  check macro defn in include/utils.h\n");
    return E_OTHER;
  }

  /*  fprintf(stderr,"initialize grid\n");*/
  init_default_grid(g);

  /*  fprintf(stderr,"look for n's, d's, and o's, also order\n");*/
  for (i=0;i<IARR_MAX_NDIM;i++) {
    snprintf(key,kl,"n%d",i+1);
    tmp=1;
    ps_flint(par,key,&tmp);
    //////
    /*
    fprintf(stderr,"i=%d tmp=%d key=%s\n",i,tmp,key);
    */
    g->axes[i].n=tmp;
    snprintf(key,kl,"d%d",i+1);
    if (ps_flreal(par,key,&(g->axes[i].d))) g->axes[i].d=1.0;
    snprintf(key,kl,"o%d",i+1);
    if (ps_flreal(par,key,&(g->axes[i].o))) g->axes[i].o=0.0;
    // added 30.11.13 WWS
    snprintf(key,kl,"id%d",i+1);
    if (ps_flint(par,key,&(g->axes[i].id))) g->axes[i].id=-1;
    /* determine global dim by finding least axis index with n>1 */
    // added 15.12.13 WWS - regard axes with assigned id's as precious
    if (g->axes[i].n>1 || g->axes[i].id > -1) g->gdim=iwave_max(g->gdim,i+1);
  }

  /* mod of 10.10.12: keyword dim = physical dimension = record
     dimension - default = global dim */
  g->dim = g->gdim;
  ps_flint(par,"dim",&(g->dim));
  if (g->dim < 0 || g->dim > RARR_MAX_NDIM) {
    fprintf(fp,"ERROR: par_grid\n");
    fprintf(fp,"keyword = dim value = %d out of bounds [0, %d]\n",g->dim,RARR_MAX_NDIM);
    return E_OTHER;
  }
  /* allow some top dims to be = 1 */	    
  if (g->dim > g->gdim) g->gdim = g->dim;

  // next block of code gives primacy to already 
  // assigned axis ids, but also allows deprecated 
  // assignment by name for spatial axes.
  if (g->dim > 0) {
    // search for already allocated id=0
    tmp=1;
    for (i=0;i<g->gdim;i++) {
      if (g->axes[i].id==0) tmp=-1;
    }
    if (tmp==1) {
      ps_flint(par,"z_axis",&tmp);
      tmp--;
      if (tmp<0 || tmp>g->dim-1) {
	fprintf(fp,"ERROR: par_grid\n");
	fprintf(fp,"  z_axis index = %d out of range for dim = %d\n",tmp,g->dim);
	return E_OTHER;
      }
      g->axes[tmp].id=0;
    }
  }
  if (g->dim > 1) { 
    tmp=2;
    for (i=0;i<g->gdim;i++) {
      if (g->axes[i].id==1) tmp=-1;
    }
    if (tmp==2) {
      ps_flint(par,"x_axis",&tmp);
      tmp--;
      if (tmp<0 || tmp>g->dim-1) {
	fprintf(fp,"ERROR: par_grid\n");
	fprintf(fp,"  x_axis index = %d out of range for dim = %d\n",tmp,g->dim);
	return E_OTHER;
      }
      g->axes[tmp].id=1;
    }
  }
  if (g->dim > 2) {
    tmp=3;
    for (i=0;i<g->gdim;i++) {
      if (g->axes[i].id==2) tmp=-1;
    }
    if (tmp==3) {
      ps_flint(par,"y_axis",&tmp);
      tmp--;
      if (tmp<0 || tmp>g->dim-1) {
	fprintf(fp,"ERROR: par_grid\n");
	fprintf(fp,"  y_axis index = %d out of range for dim = %d\n",tmp,g->dim);
	return E_OTHER;
      }
      g->axes[tmp].id=2;
    }
  }
  /*
  fprintf(stderr,"in par_grid:\n");
  fprint_grid(stderr,*g);
  */
  return 0;
}

int extend_loop(ireal * a,
		const IPNT ran,
		int dim,
		int ax,
		int kxmin,
		int kxmax,
		int kxlim) {

  int kx, k0, k1, k2;  /* counters for dim'l loops */

  for (kx=kxmin;kx<kxmax;kx++) {

    if (dim==1) {
      a[kx]=a[kxlim];
    }    
    
    if (dim==2) {
      if (ax==0) {
	for (k1=0;k1<ran[1];k1++) {
	  a[kx+k1*ran[0]]=a[kxlim+k1*ran[0]];	    
	}
      }
      else {
	for (k0=0;k0<ran[0];k0++) {
	  a[k0+kx*ran[0]]=a[k0+kxlim*ran[0]];
	}
      }
    }
    
    if (dim==3) {
      if (ax==0) {
	for (k2=0;k2<ran[2];k2++) {
	  for (k1=0;k1<ran[1];k1++) {
	    a[kx+(k1+k2*ran[1])*ran[0]]
	      =a[kxlim+(k1+k2*ran[1])*ran[0]];	    
	  }
	}
      }
      else if (ax==1) {
	for (k2=0;k2<ran[2];k2++) {
	  for (k0=0;k0<ran[0];k0++) {
	    a[k0+(kx+k2*ran[1])*ran[0]]
	      =a[k0+(kxlim+k2*ran[1])*ran[0]];
	  }
	}
      }
      else {
	for (k1=0;k1<ran[1];k1++) {
	  for (k0=0;k0<ran[0];k0++) {
	    a[k0+(k1+kx*ran[1])*ran[0]]
	      =a[k0+(k1+kxlim*ran[1])*ran[0]];
	  }
	}
      }      
    }
  }

  return 0;
}

int adj_extend_loop(ireal * a,
		    const IPNT ran,
		    int dim,
		    int ax,
		    int kxmin,
		    int kxmax,
		    int kxlim) {

  int kx, k0, k1, k2;  /* counters for dim'l loops */

  for (kx=kxmin;kx<kxmax;kx++) {

    if (dim==1) {
      a[kxlim]+=a[kx];
      if (kx != kxlim) a[kx]=REAL_ZERO;
    }    
    
    if (dim==2) {
      if (ax==0) {
	for (k1=0;k1<ran[1];k1++) {
	  a[kxlim+k1*ran[0]] += a[kx+k1*ran[0]];
	  if (kxlim != kx) a[kx+k1*ran[0]]=REAL_ZERO;
	}
      }
      else {
	for (k0=0;k0<ran[0];k0++) {
	  a[k0+kxlim*ran[0]]+=a[k0+kx*ran[0]];
	  if (kxlim != kx) a[k0+kx*ran[0]]=REAL_ZERO;
	}
      }
    }
    
    if (dim==3) {
      if (ax==0) {
	for (k2=0;k2<ran[2];k2++) {
	  for (k1=0;k1<ran[1];k1++) {
	    a[kxlim+(k1+k2*ran[1])*ran[0]]
	      +=a[kx+(k1+k2*ran[1])*ran[0]];
	    if (kxlim != kx) a[kx+(k1+k2*ran[1])*ran[0]] = REAL_ZERO;
	  }
	}
      }
      else if (ax==1) {
	for (k2=0;k2<ran[2];k2++) {
	  for (k0=0;k0<ran[0];k0++) {
	    a[k0+(kxlim+k2*ran[1])*ran[0]]
	      +=a[k0+(kx+k2*ran[1])*ran[0]];
	    if (kxlim != kx) a[k0+(kx+k2*ran[1])*ran[0]] = REAL_ZERO;
	  }
	}
      }
      else {
	for (k1=0;k1<ran[1];k1++) {
	  for (k0=0;k0<ran[0];k0++) {
	    a[k0+(k1+kxlim*ran[1])*ran[0]]
	      +=a[k0+(k1+kx*ran[1])*ran[0]];
	    if (kxlim != kx) a[k0+(k1+kx*ran[1])*ran[0]] = REAL_ZERO;
	  }
	}
      }      
    }
  }

  return 0;
}
      
/* array extension happens only over physical dimensions, assumed to
   be at most 3 */
int extend_array(ireal * a, 
		 const IPNT rags, 
		 const IPNT ran, 
		 const IPNT gs, 
		 const IPNT n, 
		 int dim, 
		 int ax) {

  int kxlim;        /* source index for ax loop */
  int kxmin, kxmax; /* limits for ax loop */
  int err=0;        /* return value */

  /* sanity */
  /*  if (RARR_MAX_NDIM > 3) err=E_OUTOFBOUNDS; */
  if (dim > 3 || dim < 1) err=E_OUTOFBOUNDS;
  if (ax < 0 || ax > dim-1) err=E_OUTOFBOUNDS;
  
  if (err) return err;

  if (!err && (rags[ax]<gs[ax])) {

    kxlim=iwave_min(ran[ax],gs[ax]-rags[ax]);
    kxmin=0;
    kxmax=kxlim;

    err=extend_loop(a,ran,dim,ax,kxmin,kxmax,kxlim);

  }

  if (!err && (rags[ax]+ran[ax]>gs[ax]+n[ax])) {

    kxlim=iwave_max(0,gs[ax]-rags[ax]+n[ax]-1);
    kxmin=kxlim+1;
    kxmax=ran[ax];

    err=extend_loop(a,ran,dim,ax,kxmin,kxmax,kxlim);
 
  }
   
  return err;
}

int adj_extend_array(ireal * a, 
		     const IPNT rags, 
		     const IPNT ran, 
		     const IPNT gs, 
		     const IPNT n, 
		     int dim, 
		     int ax) {

  int kxlim;        /* source index for ax loop */
  int kxmin, kxmax; /* limits for ax loop */
  int err=0;        /* return value */

  /* sanity */
  /*  if (RARR_MAX_NDIM > 3) err=E_OUTOFBOUNDS; */
  if (dim > 3 || dim < 1) err=E_OUTOFBOUNDS;
  if (ax < 0 || ax > dim-1) err=E_OUTOFBOUNDS;
  
  if (err) return err;

  if (!err && (rags[ax]<gs[ax])) {

    kxlim=iwave_min(ran[ax],gs[ax]-rags[ax]);
    kxmin=0;
    kxmax=kxlim;

    err=adj_extend_loop(a,ran,dim,ax,kxmin,kxmax,kxlim);

  }

  if (!err && (rags[ax]+ran[ax]>gs[ax]+n[ax])) {

    kxlim=iwave_max(0,gs[ax]-rags[ax]+n[ax]-1);
    kxmin=kxlim+1;
    kxmax=ran[ax];

    err=adj_extend_loop(a,ran,dim,ax,kxmin,kxmax,kxlim);
 
  }
   
  return err;
}

/* serial implementation of rsfread */

#ifndef IWAVE_USE_MPI

/* private - not prototyped in .h */
int rsfread_base(ireal * a, 
		 const IPNT rags, 
		 const IPNT ran,
		 const char * dname,
		 const char * type,
		 const char * protodata,
		 int dim,
		 int panelnum,
		 int datasize,
		 IPNT n,
		 ireal scale,
		 int extend, 
		 FILE * stream,
		 int panelindex   /* D.S. 01.01.11: extended-model related */
		 ) {

  /**************************
   * BEGIN DECLARATIONS     *
   **************************/

  /* rags = start indices of target array */
  /* ran  = lengths of target array axes  */

  /* workspace */  
  int ii;          /* axis counter */
  off_t i,j;       /* offset counters */
  int err=0;       /* error flag */
  IPNT gs;         /* start indices of grid (file) - always 0! */
  IPNT gsa;        /* start indices of grid intersection */
  IPNT g_gsa;      /* start indices, global */
  IPNT l_gsa;      /* start indices, local  */
  IPNT gea;        /* end   indices of grid intersection */
  IPNT na;         /* lengths of grid intersection axes */
  IPNT gl_na;      /* lengths of grid intersection axes, local or global */
  size_t recsize_b;/* length of 1D chunk to be read (bytes) */

  FILE * fp = NULL;
  float * fbuf;    /* input buffer for read */

  /* XDR buffers etc. */
#ifdef SUXDR
  XDR  xdrs;       /* xdr structure */
  int xflag;       /* i/o success flag */
  xflag = 0;       /* To avoid "uninitialized variable" warning */
  char * buf;      /* input buffer for XDR stream */
#endif

  /* grid offsets - allocated dynamically */
  off_t * goffs;   /* global offsets - into file data */
  off_t * loffs;   /* local offsets - into RARR data */
  size_t noffs;    /* number of offsets */

  // scale outsource to calling function - 02.14
  /* scale flag - like scalco, scalel in SEGY; factor */
  //  float scfac=1.0; /* factor workspace */
  //  size_t ntot;     /* total length of constructed local array */
  //  float a_max;
  //  float a_min;
  
  /* to write movie panel - may need to seek */
  off_t cur_pos = 0;

  /**************************
   * END DECLARATIONS       *
   **************************/

  /* open data file */
#ifdef IWAVE_USE_FMGR
  fp=iwave_const_fopen(dname,"r",protodata,stream);
  if (!fp) {
    fprintf(stream,"Error: rsfread - read from iwave_fopen, file = %s\n",dname);
    return E_FILE;
  }
#else
  fp=fopen(dname,"r");
  if (!fp) {
    fprintf(stream,"Error: rsfread - read from fopen, file = %s\n",dname);
    return E_FILE;
  }
#endif

  /* NOTE: the global index origin is always IPNT_0, 
     regardless of what the physical grid coordinates are
  */
  IASN(gs,IPNT_0);
  /*  get_gs(gs, g);*/

  /* figger out array intersection. special case: if there is no
     intersection, read nearest part of data */
  for (ii=0; ii < dim; ii++) {
    gsa[ii] = iwave_max(gs[ii], rags[ii]);
    gea[ii] = iwave_min(gs[ii] + n[ii] - 1, rags[ii] + ran[ii] - 1);
    na[ii]  = iwave_max(gea[ii] - gsa[ii] + 1, 0);
    
    /* rarray to left of garray */
    if ( (rags[ii] + ran[ii] - 1 < gs[ii]) && extend ) { 
      g_gsa[ii] = gs[ii];
      l_gsa[ii] = rags[ii] + ran[ii] - 1;
      gl_na[ii] = 1;
    }
    /* rarray to right of garray */
    else if ( (rags[ii] > gs[ii] + n[ii] -1) && extend ) {
      g_gsa[ii] = gs[ii] + n[ii] - 1;
      l_gsa[ii] = rags[ii];
      gl_na[ii] = 1;
    }
    /* intersection nonempty */
    else {
      g_gsa[ii] = gsa[ii];
      l_gsa[ii] = gsa[ii];
      gl_na[ii] = na[ii];
    }
  }

  /* use array data to work out offsets of 1D segments in both
     grid array (global) and rarray (local) */
  if ( (err = get_array_offsets(&goffs, &noffs, dim, gs, n, g_gsa, gl_na)) ||
       (err = get_array_offsets(&loffs, &noffs, dim, rags, ran, l_gsa, gl_na)) ) {
    fprintf(stream,"Error: rsfread - read from get_array_offsets, err=%d\n",err);
    return err;
  }  

  /* compute current starting position */
  /* first sanity check */
  if (panelindex < 0) {
    fprintf(stream,"Error: rsfread\n");
    fprintf(stream,"panelindex = %d < 0\n",panelindex);
    return E_OTHER;
  }
      
  /* seek to panel at input panelindex modulo number of panels */
  panelindex = panelindex % panelnum;
  cur_pos = panelindex * datasize * sizeof(float);

  /*** from here on, memory is allocated, so save return until
       deallocation */
  
  recsize_b = gl_na[0]*sizeof(float);

  /* native read loop */
  
  if (!strcmp(type,"native_float")) {

    /* allocate read buffer */
    fbuf=(float *)usermalloc_(recsize_b);
    for (i=0;i<(int)noffs;i++) {
      /* seek to read segment */
      if (!err && fseeko(fp,goffs[i]*sizeof(float) + cur_pos,SEEK_SET)) {
	fprintf(stream,"Error: rsfread from fseeko at file offset %ld\n",(intmax_t)(goffs[i]*sizeof(float)+cur_pos)); 
	err=E_FILE;
      }
      /* read in byte string */
      if (!err && (gl_na[0] != (int)fread(fbuf,sizeof(float),gl_na[0],fp))) {
	fprintf(stream,"Error: rsfread from fread for array offset %ld\n",(intmax_t)loffs[i]);
	fprintf(stream,"-- failed to read %d words at file offset %ld\n",
		gl_na[0],(intmax_t)(goffs[i]*sizeof(float)+cur_pos)); 
	err=E_FILE;
      }
      /* convert to ireal, scale, add */
      if (!err) {
#ifdef UPDATE
	for (j=0;j<gl_na[0];j++) { *(a+loffs[i]+j)+=scale*(*(fbuf+j)); }
#else
	for (j=0;j<gl_na[0];j++) { *(a+loffs[i]+j)=*(fbuf+j); }
#endif
      }
    }
    userfree_(fbuf);
  }
  
#ifdef SUXDR
  /* xdr read loop */
  else if (!strcmp(type,"xdr_float")) {
    fprintf(stream,"rsfread: XDR option\n");

    /* allocate char buffer */
    buf=usermalloc_(recsize_b);
    /* allocate one float */
    fbuf=usermalloc_(sizeof(float));

    /* (re)initialize xdr stream */
    xdrmem_create(&xdrs,buf,recsize_b,XDR_DECODE);
    
    for (i=0; i<noffs; i++) {
      
      /* reset xdr stream */
      if (!err) xdr_setpos(&xdrs,0);

      /* seek to read segment */
      if (!err && fseeko(fp,goffs[i]*sizeof(float) + cur_pos,SEEK_SET)) {
	fprintf(stream,"Error: rsfread from fseeko at file offset %ld\n",(intmax_t)goffs[i]);
	err=E_FILE;
      }
      /* read in byte array */
      if (!err && (recsize_b != fread(buf,sizeof(char),recsize_b,fp))) {
	fprintf(stream,"Error: rsfread from fread at array offset %ld\n",(intmax_t)loffs[i]);
	err=E_FILE;
      }
      
      /* convert to ireal - need only one word of float storage */ 
      if (!err) {
	xflag=1;
	for (j=0;j<gl_na[0];j++) {
	  xflag=xflag && xdr_float(&xdrs,fbuf);
#ifdef UPDATE
	  *(a+loffs[i]+j)+=scale * (*fbuf); 
#else
	  *(a+loffs[i]+j)=*fbuf; 
#endif
	}
      }
      if (!err && !xflag) {
	fprintf(stream,"Error: rsfread - failed xdr conversion\n");
	err=E_OTHER;
      }
    }
    /* disengage xdr stream */
    userfree_(buf);
    userfree_(fbuf);
    xdr_destroy(&xdrs);

  }
  else {
    fprintf(stream,"Error: rsfread - data format must be either\n");
    fprintf(stream,"native_float or xdr_float\n");

    return E_OTHER;
  }
#else
  else {
    fprintf(stream,"Error: rsfread - data format must be native_float\n");

    return E_OTHER;
  }
#endif

  fflush(fp);

#ifdef IWAVE_USE_FMGR
  iwave_fclose(fp);
#else
  fclose(fp);
#endif

  /* extension loop - extend by const along axes dim-1,..,0.
   */
  if (extend) {

    for (ii=dim-1;ii>-1;ii--) {
      if (extend_array(a,rags,ran,l_gsa,gl_na,dim,ii)) {
	fprintf(stream,"Error: rsfread - extension failed at axis %d\n",ii);
	return E_OTHER;
      }
    }
  } 

  /* scaling loop */
  /* REMOVED 3 Feb 2014 WWS - outsource scaling to calling routine! */
  /*
  ntot=1;
  for (ii=0;ii<dim;ii++) ntot*=ran[ii];
  if (scale>0) for (ii=0;ii<scale;ii++) scfac *= 10.0;
  if (scale<0) for (ii=0;ii<-scale;ii++) scfac *= 0.10;
  a_max=scfac*a[0];
  a_min=scfac*a[0];
  if (scale) { 
    for (i=0;i<ntot;i++) {
      a[i]*=scfac;
      a_max=iwave_max(a[i],a_max);
      a_min=iwave_min(a[i],a_min);
    }
  }
  */
  if (goffs) userfree_(goffs);
  if (loffs) userfree_(loffs);

  return err;
}

/* PRIVATE - docs removed from header file 04.02.14 WWS
 * read data from subcube of SEP77/RSF data cube, scale, add to data buffer
 *
 *  This version of rsfread has two legal use cases: 
 *  <ol>
 *  <li>prototype files, grid not provided (corresponding pointer args = NULL); in this
 *  case metadata file (fname) is source of all info about RSF data structure - data file,
 *  grid, scale exponent, data type are extracted from fname, so consistency is guaranteed
 *  if fname is the path to a legitimate metadata file (which does not need to have been
 *  opened in the current process). <p>This option is safe for read-only apps like IWAVE.</li>
 *  <li>prototype files and grid are provided - args point to existing objects. Then 
 *  consistency is responsibility of calling function. Assuming that metadata files are
 *  already opened using iwave_fopen with appropriate (r, r+, or w+) permissions on appropriate
 *  prototypes, existing file pointers will be used. Prototype grid, scale, and type arguments
 *  provide metadata, and must also be consistent - which they are if calling function has
 *  generated them from a legit metadata file.<p> This option provides a means to guarantee
 *  data integrity of temporary data in read/write apps such as IWAVE++: 
 *  GridDC opens temporary (and archival) files and keeps them open for
 *  the life of the GridDC object. The reference count feature of the IWAVE file manager 
 *  prevents the corresponding FILE*s from being re-used for writes of other DC data, despite
 *  calls to iwave_fclose within rsfread.
 *  <p>
 *  Preconditions:
 *  <ol>
 *  <li>pathname fname describes existing SEP77/RSF metadata file, including (n,d,o) grid
 *  parameters, "in" key with data file name as value, and
 *  "data_format" key with either "native_float" or "xdr_float" as
 *  parameters.</li>
 *  <li>pathname dname describes existing SEP77/RSF data file, compatible with fname</li>
 *  <li>pointers to prototype metadata file, data file, and grid either (a) all exist 
 *  and are compatible with the specified RSF files, or (b) are all NULL. In first
 *  case, compatibility is NOT checked and is responsibility of calling function.</li>
 *  <li> fname may be opened for read access, and either (a) if already
 *  open, was opened by iwave_fopen with prototype specified by proto (or
 *  any prototype if proto=NULL), or (b) has not been previously opened.
 *  <li>file pointed to by "in=" contains float data, either native or
 *  xdr. No other data types are currently admitted.</li>
 * </ol>
 *  Postconditions:
 *  <ol>
 *  <li>intersection of target and file cubes computed; data from
 *  intersection read, scaled, and added to target cube
 *  corresponding to intersection</li>
 *  </ol>
 *  @param[out]  a             - pointer to first word of target cube data array
 *  @param[in]   ags           - global indices of axis starts - target cube
 *  @param[in]   an            - global axis lengths - target cube
 *  @param[in]   fname         - pathname of RSF metadata source file
 *  @param[in]   dname         - pathname of RSF data source file
 *  @param[in]   type          - data type string, for insertion in RSF metadata
 *  @param[in]   scale         - input data scale 
 *  @param[in]   protohdr      - pathname of RSF metadata prototype, or NULL
 *  @param[in]   protodata     - pathname of RSF data prototype, or NULL
 *  @param[in]   protog        - grid, presumed to be defined by protohdr
 *  @param[in]   extens        - extension flag - extend by const along all axes in 
 *                               decreasing axis order if set (DEPRECATED)
 *  @param[in]   fp            - verbose output parameter
 *  @param[in]   panelindex    - panel index of extended model (always 
 *                               0 for non-extended model) 
 *  @return 0 on success, else nonzero error code (see \ref utils.h)
 */

int rsfread_proto(ireal * a, 
		  const IPNT rags, 
		  const IPNT ran, 
		  const char * fname,
		  const char * dname,
		  const char * type,
		  float scale,
		  const char * protohdr,
		  const char * protodata,
		  const grid * protogrid,
		  int extend, 
		  FILE * stream,
		  int panelindex   /* D.S. 01.01.11: extended-model related */
		  ) {

  /* data to be fished out of metadata file */
  FILE * fph = NULL;
  PARARRAY * par = NULL;
  grid g;
  char * ltype = NULL;
  char * ldname = NULL;
  //  int lscale=0;     /* local scaling index */
  IPNT n;          /* lengths of grid (file) axes */
  int err=0;

  if (!(protohdr) && !(protodata) && !(protogrid)) {
    par = ps_new();
    
    /* read parameter table from file */
    if (ps_createfile_fproto(par,&fph,protohdr,fname)) {
      fprintf(stream,"read_grid from ps_creatfile_fproto: failed to parse file = %s\n",fname);
      ps_delete(&par);
      return E_FILE;
    }

    /* having created par, can close stream */
#ifdef IWAVE_USE_FMGR
    iwave_fclose(fph);
#else
    fclose(fph);
#endif

    /* create grid from parameter table */
    if ( (err = par_grid(&g, *par, stream)) ) {
      fprintf(stream,"Error: read from read_grid\n");
      ps_delete(&par);
      return err;
    }
  
    /* get data filename */
    if ( (err=ps_flcstring(*par,"in",&ldname)) ) {
      fprintf(stream,"Error: rsfread - read from ps_flcstring\n");
      fprintf(stream,"failed to extract in param\n");
      ps_delete(&par);
      return err;
    }    

    if ( (err=ps_flcstring(*par,"data_format",&ltype)) ) {
      fprintf(stream,"Error: rsfread - read from ps_flcstring\n");
      fprintf(stream,"failed to extract type param\n");
      ps_delete(&par);
      return err;
    }
    //    ps_flint(*par,"scale",&lscale);
    ps_delete(&par);
  }

  else if ((protohdr) && (protodata) && (protogrid)) { 
    ldname = (char *)usermalloc_((strlen(dname)+10)*sizeof(char));
    strcpy(ldname,dname);
    ltype = (char *)usermalloc_((strlen(type)+10)*sizeof(char));
    strcpy(ltype,type);
    //    lscale=scale;
    copy_grid(&g,protogrid);
  }
  else {
    fprintf(stream,"Error: rsfread_proto - inconsistent values of\n");
    fprintf(stream,"protohdr, protodata, protog args\n");
    return E_BADINPUT;
  }
    
  /* get global array lengths from grid 
     */
  get_n(n, g);

  err=rsfread_base(a,rags,ran,
		   ldname,ltype,protodata,
		   get_dimension_grid(g),get_panelnum_grid(g),get_extended_datasize_grid(g),n,
		   scale,extend,stream,panelindex);
  if (err) 
    fprintf(stream,"Error: rsfread_proto from rsfread_base err=%d\n",err);
  
  userfree_(ldname);
  userfree_(ltype);

  return err;

}

int rsfread(ireal * a, 
	    const IPNT rags, 
	    const IPNT ran, 
	    const char * fname,
	    float scale, 
	    //	    int extend,
	    FILE * stream,
	    int panelindex   /* D.S. 01.01.11: extended-model related */
	    ) {
  // note 04.02.14: extend option now deprecated, const = false
  // NOTE scale is now float WWS 03.02.14
  return rsfread_proto(a,rags,ran,fname,NULL,NULL,scale,NULL,NULL,NULL,0,stream,panelindex);
}

// NOTE scale is now float WWS 03.02.14
int rsfwrite_proto(ireal * a, 
		   const IPNT rags, 
		   const IPNT ran, 
		   const char * fname,
		   const char * dname,
		   const char * type,
		   float scale,
		   const char * protohdr,
		   const char * protodata,
		   const grid * protogrid,
		   int extend, 
		   FILE * stream,
		   int panelindex   /* D.S. 01.01.11: extended-model related */
		   ) {

  /**************************
   * BEGIN DECLARATIONS     *
   **************************/

  /* strings */
  char * ltype;
  char * ldname;

  /* workspace */  
  grid g;
  int ii;
  off_t i,j;       /* counters */
  int err=0;       /* error flag */
  IPNT g_gsa;      /* start indices, global */
  IPNT l_gsa;      /* start indices, local  */
  IPNT gea;        /* end   indices of grid intersection */
  IPNT gs;         /* start indices of grid (file) */
  IPNT gsa;        /* start indices of grid intersection */
  IPNT n;          /* lengths of grid (file) axes */
  IPNT na;         /* lengths of grid intersection axes */
  IPNT gl_na;      /* lengths of grid intersection axes, local or global */
  FILE * fp = NULL;
  FILE * fph = NULL;
  PARARRAY * par = NULL;

  /*  char * schar;*/
  float * fbuf;    /* input buffer for read */
  size_t recsize_b;/* length of 1D chunk to be read (bytes) */

  /* XDR buffers etc. */
#ifdef SUXDR
  XDR  xdrs;       /* xdr structure */
  int xflag;       /* i/o success flag */
  char * buf;      /* input buffer for xdr stream */
#endif

  // scale outsourced to calling function - 02.14
  /* scale flag - like scalco, scalel in SEGY; factor */
  //  int ntot;
  //  int lscale=0;     /* flag - read from parameters */
  //  float scfac=1.0; /* factor workspace */

  /* grid offsets - allocated dynamically */
  off_t * goffs=NULL;   /* global offsets - into file data */
  off_t * loffs=NULL;   /* local offsets - into RARR data */
  size_t noffs;    /* number of offsets */

  /* possible nonzero offset for movie writes */
  off_t cur_pos = 0;

  /**************************
   * END DECLARATIONS       *
   **************************/

  if (!(protohdr) && !(protodata) && !(protogrid)) {

    par = ps_new();
    
    /* read parameter table from file */
    if (ps_createfile_fproto(par,&fph,protohdr,fname)) {
      fprintf(stream,"read_grid from ps_creatfile_fproto: failed to parse file = %s\n",fname);
      ps_delete(&par);
      return E_FILE;
    }

    /* having created par, can close stream */
#ifdef IWAVE_USE_FMGR
    iwave_fclose(fph);
#else
    fclose(fph);
#endif

    /* create grid from parameter table */
    if ( (err = par_grid(&g, *par, stream)) ) {
      fprintf(stream,"Error: read from read_grid\n");
      ps_delete(&par);
      return err;
    }
  
    /* get data filename */
    if ( (err=ps_flcstring(*par,"in",&ldname)) ) {
      fprintf(stream,"Error: rsfwrite - read from ps_flcstring\n");
      fprintf(stream,"failed to extract in param\n");
      ps_delete(&par);
      return err;
    }    

    if ( (err=ps_flcstring(*par,"data_format",&ltype)) ) {
      fprintf(stream,"Error: rsfwrite - read from ps_flcstring\n");
      fprintf(stream,"failed to extract type param\n");
      ps_delete(&par);
      return err;
    }
    //    ps_flint(*par,"scale",&lscale);
    ps_delete(&par);
  }

  else if ((protohdr) && (protodata) && (protogrid)) { 
    ldname = (char *)usermalloc_((strlen(dname)+10)*sizeof(char));
    strcpy(ldname,dname);
    ltype = (char *)usermalloc_((strlen(type)+10)*sizeof(char));
    strcpy(ltype,type);
    //    lscale=scale;
    copy_grid(&g,protogrid);
    fprintf(stderr,"IN RSFWRITE: hdr=%s data=%s\n",fname,dname);
    iwave_fprintall(stderr);
  }
  else {
    fprintf(stream,"Error: rsfwrite_proto - inconsistent values of\n");
    fprintf(stream,"protohdr, protodata, protog args\n");
    return E_BADINPUT;
  }

  /* DATA FILE SECTION */
  
  // WWS 03.02.14 - outsource!      
  /* scaling loop note that scale should be inverse of read */
  /*  ntot=1;
  for (ii=0;ii<get_dimension_grid(g);ii++) ntot*=ran[ii];
  lscale=-lscale;
  if (lscale>0) for (ii=0;ii<lscale;ii++) scfac *= 10.0;
  if (lscale<0) for (ii=0;ii<-lscale;ii++) scfac *= 0.10;
  if (lscale) 
    for (i=0;i<ntot;i++) a[i]*=scfac;
  */

  /* get global array data from grid 
     WWS 09.11.09 - gs is ALWAYS IPNT_0 */
  IASN(gs,IPNT_0);
  get_n(n,g);

  /* figger out array intersection. special case: if there is no
     intersection, read nearest part of data */
  for (ii=0; ii < get_dimension_grid(g); ii++) {
    gsa[ii] = iwave_max(gs[ii], rags[ii]);
    gea[ii] = iwave_min(gs[ii] + n[ii] - 1, rags[ii] + ran[ii] - 1);
    na[ii]  = iwave_max(gea[ii] - gsa[ii] + 1, 0);
    
    /* rarray to left of garray */
    if ( (rags[ii] + ran[ii] - 1 < gs[ii]) && extend ) { 
      g_gsa[ii] = gs[ii];
      l_gsa[ii] = rags[ii] + ran[ii] - 1;
      gl_na[ii] = 1;
    }
    /* rarray to right of garray */
    else if ( (rags[ii] > gs[ii] + n[ii] -1) && extend ) {
      g_gsa[ii] = gs[ii] + n[ii] - 1;
      l_gsa[ii] = rags[ii];
      gl_na[ii] = 1;
    }
    /* intersection nonempty */
    else {
      g_gsa[ii] = gsa[ii];
      l_gsa[ii] = gsa[ii];
      gl_na[ii] = na[ii];
    }
  }

  /* use array data to work out offsets of 1D segments in both
     grid array (global) and rarray (local) */
  if ( (err = get_array_offsets(&goffs, &noffs, get_dimension_grid(g), gs, n, g_gsa, gl_na)) ||
       (err = get_array_offsets(&loffs, &noffs, get_dimension_grid(g), rags, ran, l_gsa, gl_na)) ) {
    fprintf(stream,"ERROR: rsfwrite from get_array_offsets, err=%d\n",err);
    ps_delete(&par);
    return err;
  }
  if (!goffs || !loffs) {
    fprintf(stream,"Error: rsfwrite - failed to allocate goffs or loffs\n");
    return E_ALLOC;
  }

  /* adjoint extension loop - adj_extend by const along axes 0,...,dim-1.
     added 06.03.12 to make rsfwrite l2-adjoint to rsfwrite WWS
  */
  if (extend) {
    for (ii=0;ii<get_dimension_grid(g);ii++) {
      if (adj_extend_array(a,rags,ran,l_gsa,gl_na,get_dimension_grid(g),ii)) {
	fprintf(stream,"Error: rsfwrite - adjoint extension failed at axis %d\n",ii);
	ps_delete(&par);
	return E_OTHER;
      }
    }
  } 

  recsize_b=na[0]*sizeof(float);

  /* open data file - "w" option works if header file already exists but
     data file does not
     WWS 19.01.11: data file must exist prior to call, so use r+ exclusively
  */
  fp=NULL;
#ifdef IWAVE_USE_FMGR
  fp=iwave_const_fopen(ldname,"r+",NULL,stream);
#else
  p=fopen(ldname,"r+")
#endif
    if (!fp) {
      fprintf(stream,"Error: rsfwrite from fopen: cannot open file=%s\n",ldname);
      return E_FILE;
    }

  fflush(stream);

  /* compute current starting position */
  /* first sanity check */
  if (panelindex < 0) {
    fprintf(stream,"Error: rsfwrite\n");
    fprintf(stream,"panelindex = %d < 0\n",panelindex);
    return E_OTHER;
  }
      
  /* seek to panel at input panelindex modulo number of panels */
  panelindex = panelindex % get_panelnum_grid(g);
  cur_pos = panelindex * get_extended_datasize_grid(g) * sizeof(float);

  if (!strcmp(ltype,"native_float")) {

    /* allocate float buffer */
    fbuf=(float *)usermalloc_(recsize_b);

    for (i=0;i<(int)noffs;i++) {
      /* seek to write segment */
      if (!err && fseeko(fp,goffs[i]*sizeof(float) + cur_pos,SEEK_SET)) {
	fprintf(stream,"Error: rsfwrite from fseeko at file offset %ld\n",(intmax_t)goffs[i]); 
	fprintf(stream,"possible cause: attempt to write off end of file\n");
	fseeko(fp,0L,SEEK_END);
	fprintf(stream,"file length = %ld:\n",(intmax_t)ftello(fp));
	fprintf(stream,"note that new file can only be written in contiguous,\n");
	fprintf(stream,"consecutive blocks\n");
	err=E_FILE;
      }
      /* convert from ireal */
      if (!err) {
#ifdef UPDATE
	for (j=0;j<na[0];j++) { *(fbuf+j)=scale * (*(a+loffs[i]+j)); }
#else
	for (j=0;j<na[0];j++) { *(fbuf+j)=*(a+loffs[i]+j); }
#endif
      }
      /* write out float buffer */
      if (!err && (na[0] != (int)fwrite(fbuf,sizeof(float),na[0],fp))) {
	fprintf(stream,"Error: rsfwrite from fwrite at array offset %ld\n",(intmax_t)loffs[i]);
	fprintf(stream,"failed to write %d words\n",na[0]);
	err=E_FILE;
      }
      // added 21.03.14 WWS
      fflush(fp);

      //      fprintf(stderr,"rsfwrite: wrote %d floats at offset %jd\n",na[0],goffs[i]*sizeof(float) + cur_pos);
    }
    userfree_(fbuf);
    
  }

#ifdef SUXDR
  /* xdr write loop */
  else if (!strcmp(ltype,"xdr_float")) {

    /* allocate char buffer */
    buf=usermalloc_(recsize_b);
    /* allocate one float */
    fbuf=usermalloc_(sizeof(float));

    /* (re)initialize xdr stream */
    xdrmem_create(&xdrs,buf,recsize_b,XDR_ENCODE);
    
    for (i=0;i<noffs;i++) {

      /* reset xdr stream */
      if (!err) xdr_setpos(&xdrs,0);

      /* seek to write segment */
      if (!err && fseeko(fp,goffs[i]*sizeof(float) + cur_pos, SEEK_SET)) {
	fprintf(stream,"Error: rsfwrite from fseeko at file offset %ld\n",(intmax_t)goffs[i]);
	fprintf(stream,"possible cause: attempt to write off end of file\n");
	j=0;
	fseeko(fp,j,SEEK_END);
	fprintf(stream,"file length = %ld:\n",(intmax_t)ftello(fp));
	fprintf(stream,"note that new file can only be written in contiguous,\n");
	fprintf(stream,"consecutive blocks\n");
	err=E_FILE;
      }

      /* convert from ireal - need only one word of float storage */ 
      xflag = 0;
      if (!err) {
	xflag=1;
	for (j=0;j<na[0];j++) {
#ifdef UPDATE
	  *fbuf=scale * (*(a+loffs[i]+j));
#else
	  *fbuf=*(a+loffs[i]+j);
#endif
	  xflag=xflag && xdr_float(&xdrs,fbuf);
	}
      }
      /* write out byte array */
      if (!err && (recsize_b != fwrite(buf,sizeof(char),recsize_b,fp))) {
	fprintf(stream,"Error: rsfwrite from fwrite at array offset %ld\n",(intmax_t)loffs[i]);
	err=E_FILE;
      }
      // added 21.03.14 WWS
      fflush(fp);
      if (!err && !xflag) {
	fprintf(stream,"Error: rsfwrite - failed xdr conversion\n");
	err=E_OTHER;
      }
    }
    /* disengage xdr stream */
    userfree_(buf);
    userfree_(fbuf);
    xdr_destroy(&xdrs);

  }
  else {
    fprintf(stream,"Error: rsfwrite - data format must be either\n");
    fprintf(stream,"native_float or xdr_float\n");
    return E_OTHER;
  }

#else
  else {
    fprintf(stream,"Error: rsfwrite - data format must be native_float\n");
    return E_OTHER;
  }
#endif

  /* free workspace */
  if (goffs) userfree_(goffs);
  if (loffs) userfree_(loffs);
  
  userfree_(ldname);
  userfree_(ltype);
  ps_delete(&par);
  fflush(fp);
#ifdef IWAVE_USE_FMGR
  iwave_fclose(fp);
#else
  fclose(fp);
#endif

  return err;
}

int rsfwrite(ireal * a, 
	     const IPNT rags, 
	     const IPNT ran, 
	     const char * fname, 
	     float scale, // added 04.02.14
	     //int extend,       /* WWS 06.03.11 */ /* masked 04.02.14 */
	     FILE * stream,
	     int panelindex  /* D.S. 01.01.11: extended-model related */
	     ) {
  // note 04.02.14: extend option now deprecated, const = false
  // NOTE scale is now float WWS 03.02.14
  return rsfwrite_proto(a,rags,ran,fname,NULL,NULL,scale,NULL,NULL,NULL,0,stream,panelindex);
}

#endif


  
