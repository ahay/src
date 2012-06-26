/*
WWS 09.09: make compilation conditional on MPI
*/
#include <stdio.h>

#include "gridio.h"

int read_grid(grid * g, char * fname, FILE * fp) {

  int err=0;
  PARARRAY * par = ps_new();

  /* parser file for key=value pairs */
  if (ps_createfile(par,fname)) {
    //    fprintf(fp,"read_grid: failed to parse file %s\n",fname);
    return E_FILE;
  }

  if ( (err=par_grid(g,*par,fp)) ) {
    //    fprintf(fp,"read_grid: failed to parse table from file %s\n",fname);
    return err;
  }

  ps_delete(&par);

  return err;
}

int par_grid(grid * g, PARARRAY par, FILE * fp) {
  
  int i;
  int tmp;
  char key[3];
  size_t kl=3;

  /*  fprintf(stderr,"in par_grid\n");*/
  if (RARR_MAX_NDIM > 9) {
    fprintf(fp,"ERROR: par_grid\n");
    fprintf(fp,"  Oooops...grid parser assumes max dim < 10\n");
    fprintf(fp,"  check macro defn in include/utils.h\n");
    return E_OTHER;
  }

  /*  fprintf(stderr,"initialize grid\n");*/
  init_default_grid(g);

  /*  fprintf(stderr,"look for n's, d's, and o's, also order\n");*/
  for (i=0;i<RARR_MAX_NDIM;i++) {
    snprintf(key,kl,"n%d",i+1);
    tmp=1;
    ps_flint(par,key,&tmp);
    g->axes[i].n=tmp;
    snprintf(key,kl,"d%d",i+1);
    if (ps_flreal(par,key,&(g->axes[i].d))) g->axes[i].d=1.0;
    snprintf(key,kl,"o%d",i+1);
    if (ps_flreal(par,key,&(g->axes[i].o))) g->axes[i].o=0.0;
    /* determine dim by finding least axis index with n>1 */
    if (g->axes[i].n>1) g->dim=iwave_max(g->dim,i+1);
  }
  /*  fprintf(stderr,"order params - dim=%d\n",g->dim);*/
  if (g->dim > 0) { 
    tmp=1;
    ps_flint(par,"z_axis",&tmp);
    tmp--;
    /*    fprintf(stderr,"z axis tmp=%d\n",tmp);*/
    if (tmp<0 || tmp>g->dim-1) {
      fprintf(fp,"ERROR: par_grid\n");
      fprintf(fp,"  z_axis index = %d out of range for dim = %zu\n",tmp,g->dim);
      return E_OTHER;
    }
    g->axes[tmp].id=0;
  }
  if (g->dim > 1) { 
    tmp=2;
    ps_flint(par,"x_axis",&tmp);
    tmp--;
    /*    fprintf(stderr,"x axis tmp=%d\n",tmp);*/
    if (tmp<0 || tmp>g->dim-1) {
      fprintf(fp,"ERROR: par_grid\n");
      fprintf(fp,"  x_axis index = %d out of range for dim = %zu\n",tmp,g->dim);
      return E_OTHER;
    }
    g->axes[tmp].id=1;
  }
  if (g->dim > 2) {
    tmp=3;
    ps_flint(par,"y_axis",&tmp);
    tmp--;
    /*    fprintf(stderr,"y axis tmp=%d\n",tmp);*/
    if (tmp<0 || tmp>g->dim-1) {
      fprintf(fp,"ERROR: par_grid\n");
      fprintf(fp,"  y_axis index = %d out of range for dim = %zu\n",tmp,g->dim);
      return E_OTHER;
    }
    g->axes[tmp].id=2;
  }
  
  /*  fprintf(stderr,"exit par_grid\n");*/
  return 0;
}

/*
int ra_create_g(RARR * _ra, grid g, FILE * fp) {

  IPNT n;
  IPNT gs;
  int err=0;

  if (g.dim > RARR_MAX_NDIM) {
    fprintf(fp,"Error: init_rarr\n");
    fprintf(fp,"input grid dim out of bounds\n");
    return E_OUTOFBOUNDS;
  }
  err=err||get_n(n,g);
  err=err||get_gs(gs,g);
  if (err) {
    fprintf(fp,"Error: init_rarr\n");
    fprintf(fp,"failed to get params from grid input\n");
    print_grid(g);
    
    return err;
  }
  return ra_create_s(_ra,g.dim,gs,n);
}
*/

int extend_loop(ireal * a,
		IPNT ran,
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
		    IPNT ran,
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
      
int extend_array(ireal * a, 
		 IPNT rags, 
		 IPNT ran, 
		 IPNT gs, 
		 IPNT n, 
		 int dim, 
		 int ax) {

  int kxlim;        /* source index for ax loop */
  int kxmin, kxmax; /* limits for ax loop */
  int err=0;        /* return value */

  /* sanity */
  if (RARR_MAX_NDIM > 3) err=E_OUTOFBOUNDS;
  if (dim > RARR_MAX_NDIM || dim < 1) err=E_OUTOFBOUNDS;
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
		     IPNT rags, 
		     IPNT ran, 
		     IPNT gs, 
		     IPNT n, 
		     int dim, 
		     int ax) {

  int kxlim;        /* source index for ax loop */
  int kxmin, kxmax; /* limits for ax loop */
  int err=0;        /* return value */

  /* sanity */
  if (RARR_MAX_NDIM > 3) err=E_OUTOFBOUNDS;
  if (dim > RARR_MAX_NDIM || dim < 1) err=E_OUTOFBOUNDS;
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

int rsfread(ireal * a, 
	    IPNT rags, 
	    IPNT ran, 
	    char * fname, 
	    int extend, 
	    FILE * stream,
	    int panelindex   /* D.S. 01.01.11: extended-model related */
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
  int ii;          /* axis counter */
  off_t i,j;       /* offset counters */
  int err=0;       /* error flag */
  grid g;          /* grid workspace, init from file */
  IPNT gs;         /* start indices of grid (file) */
  IPNT gsa;        /* start indices of grid intersection */
  IPNT g_gsa;      /* start indices, global */
  IPNT l_gsa;      /* start indices, local  */
  IPNT gea;        /* end   indices of grid intersection */
  IPNT g_gea;      /* end   indices, global */
  IPNT l_gea;      /* end   indices, local */
  IPNT n;          /* lengths of grid (file) axes */
  IPNT na;         /* lengths of grid intersection axes */
  IPNT gl_na;      /* lengths of grid intersection axes, local or global */
  size_t recsize_b;/* length of 1D chunk to be read (bytes) */

  FILE * fp = NULL;
  PARARRAY * par;
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

  /* scale flag - like scalco, scalel in SEGY; factor */
  int scale=0;     /* flag - read from parameters */
  float scfac=1.0; /* factor workspace */
  size_t ntot;     /* total length of constructed local array */
  float a_max;
  float a_min;
  
  /* to write movie panel - may need to seek */
  off_t cur_pos = 0;

  /**************************
   * END DECLARATIONS       *
   **************************/

  par = ps_new();

  /* read parameter table from file */
  if (ps_createfile(par,fname)) {
    fprintf(stream,"read_grid: failed to parse file = %s\n",fname);
    ps_delete(&par);
    return E_FILE;
  }

  /* create grid from parameter table */
  if ( (err = par_grid(&g, *par, stream)) ) {
    fprintf(stream,"Error: read from read_grid\n");
    //    ps_destroy(&par);
    ps_delete(&par);
    return err;
  }
  
  /* get global array params from grid 
     NOTE: the global index origin is always IPNT_0, 
     regardless of what the physical grid coordinates are
   */
  IASN(gs,IPNT_0);
  /*  get_gs(gs, g);*/
  get_n(n, g);

 /*
  fprintf(stream,"RSFREAD:\n");
  for (ii=0;ii<RARR_MAX_NDIM;ii++) {
    fprintf(stream,"ran[%d]=%d rags[%d]=%d n[%d]=%d gs[%d]=%d\n ",ii,ran[ii],ii,rags[ii],ii,n[ii],ii,gs[ii]);
  }
  */

  /* figger out array intersection. special case: if there is no
     intersection, read nearest part of data */
  for (ii=0; ii < g.dim; ii++) {
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

  /* use array data to work out offsets of 1D segments in both
     grid array (global) and rarray (local) */
  if ( (err = get_array_offsets(&goffs, &noffs, g.dim, gs, n, g_gsa, gl_na)) ||
       (err = get_array_offsets(&loffs, &noffs, g.dim, rags, ran, l_gsa, gl_na)) ) {
    fprintf(stream,"Error: rsfread - read from get_array_offsets, err=%d\n",err);
    //    ps_destroy(&par);
    ps_delete(&par);
    return err;
  }
  
  /* get filename */
  if ( (err=ps_flcstring(*par,"in",&dname)) ) {
    fprintf(stream,"Error: rsfread - read from ps_flcstring\n");
    fprintf(stream,"failed to extract in param\n");
    ps_delete(&par);
    return err;
  }
  
  /* open file */
  /* file mgr version - no prototype file is possible */
  fp=NULL;
#ifdef IWAVE_USE_FMGR
  //  fprintf(stream,"RSFREAD->IWAVE_FOPEN, file=%s\n",dname);
  fp=iwave_const_fopen(dname,"r",NULL,stream);
  if (!fp) {
    fprintf(stream,"Error: rsfread - read from iwave_fopen, file = %s\n",dname);
    //    ps_destroy(&par);
    ps_delete(&par);
    return E_FILE;
  }
#else
  fp=fopen(dname,"r");
  if (!fp) {
    fprintf(stream,"Error: rsfread - read from fopen, file = %s\n",dname);
    //    ps_destroy(&par);
    ps_delete(&par);
    return E_FILE;
  }
#endif
 
  if ( (err=ps_flcstring(*par,"data_format",&type)) ) {
    fprintf(stream,"Error: rsfread - read from ps_flcstring\n");
    fprintf(stream,"failed to extract type param\n");
    ps_delete(&par);
    return err;
  }
  ps_flint(*par,"scale",&scale);

  /* compute current starting position */
  cur_pos = panelindex * get_datasize_grid(g) * sizeof(float);

  /*** from here on, memory is allocated, so save return until
       deallocation */
  
  recsize_b = gl_na[0]*sizeof(float);

  /* native read loop */
  
  if (!strcmp(type,"native_float")) {

    /* allocate read buffer */
    fbuf=(float *)usermalloc_(recsize_b);
    for (i=0;i<noffs;i++) {
      /* seek to read segment */
      if (!err && fseeko(fp,goffs[i]*sizeof(float) + cur_pos,SEEK_SET)) {
	fprintf(stream,"Error: rsfread from fseeko at file offset %jd\n",(intmax_t)(goffs[i]*sizeof(float)+cur_pos)); 
	err=E_FILE;
      }
      /* read in byte string */
      if (!err && (gl_na[0] != fread(fbuf,sizeof(float),gl_na[0],fp))) {
	fprintf(stream,"Error: rsfread from fread for array offset %jd\n",(intmax_t)loffs[i]);
	fprintf(stream,"-- failed to read %d words at file offset %jd\n",
		gl_na[0],(intmax_t)(goffs[i]*sizeof(float)+cur_pos)); 
	err=E_FILE;
      }
      /* convert to ireal */
      if (!err) {
	for (j=0;j<gl_na[0];j++) { *(a+loffs[i]+j)=*(fbuf+j); }
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
	fprintf(stream,"Error: rsfread from fseeko at file offset %jd\n",(intmax_t)goffs[i]);
	err=E_FILE;
      }
      /* read in byte array */
      if (!err && (recsize_b != fread(buf,sizeof(char),recsize_b,fp))) {
	fprintf(stream,"Error: rsfread from fread at array offset %jd\n",(intmax_t)loffs[i]);
	err=E_FILE;
      }
      
      /* convert to ireal - need only one word of float storage */ 
      if (!err) {
	xflag=1;
	for (j=0;j<gl_na[0];j++) {
	  xflag=xflag && xdr_float(&xdrs,fbuf);
	  *(a+loffs[i]+j)=*fbuf; 
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
    //    ps_destroy(&par);
    ps_delete(&par);

    return E_OTHER;
  }
#else
  else {
    fprintf(stream,"Error: rsfread - data format must be native_float\n");
    //    ps_destroy(&par);
    ps_delete(&par);

    return E_OTHER;
  }
#endif

  fflush(fp);
#ifdef IWAVE_USE_FMGR
  //fprintf(stream,"RSFREAD -> IWAVE_FCLOSE, file=%s\n",dname);
  iwave_fclose(fp);
#else
  fclose(fp);
#endif

  
  /* extension loop - extend by const along axes dim-1,..,0.
  */
  if (extend) {

    for (ii=g.dim-1;ii>-1;ii--) {
      if (extend_array(a,rags,ran,l_gsa,gl_na,g.dim,ii)) {
	fprintf(stream,"Error: rsfread - extension failed at axis %d\n",ii);
	return E_OTHER;
      }
    }
  } 

  /* scaling loop */
  //  fprintf(stream,"gridio: SCALE=%d\n",scale);
  ntot=1;
  for (ii=0;ii<g.dim;ii++) ntot*=ran[ii];
  if (scale>0) for (ii=0;ii<scale;ii++) scfac *= 10.0;
  if (scale<0) for (ii=0;ii<-scale;ii++) scfac *= 0.10;
  a_max=scfac*a[0];
  a_min=scfac*a[0];
  if (scale) for (i=0;i<ntot;i++) {
    a[i]*=scfac;
    a_max=iwave_max(a[i],a_max);
    a_min=iwave_min(a[i],a_min);
  }
  /*fprintf(stream,"gridio: fac=%e max=%e min=%e\n",scfac,a_max,a_min);*/

  if (goffs) userfree_(goffs);
  if (loffs) userfree_(loffs);
  userfree_(dname);
  userfree_(type);

  ps_delete(&par);

  return err;
}

int rsfwrite(ireal * a, IPNT rags, IPNT ran, char * fname, 
	     int extend,       /* WWS 06.03.11 */
	     FILE * stream
	     , int panelindex  /* D.S. 01.01.11: extended-model related */
	     ) {

  /**************************
   * BEGIN DECLARATIONS     *
   **************************/

  /* strings */
  char * dname;

  /* workspace */  
  grid g;
  char * type; 
  int ii;
  off_t i,j;       /* counters */
  int err=0;       /* error flag */
  IPNT g_gsa;      /* start indices, global */
  IPNT l_gsa;      /* start indices, local  */
  IPNT gea;        /* end   indices of grid intersection */
  IPNT g_gea;      /* end   indices, global */
  IPNT l_gea;      /* end   indices, local */
  IPNT gs;         /* start indices of grid (file) */
  IPNT gsa;        /* start indices of grid intersection */
  IPNT n;          /* lengths of grid (file) axes */
  IPNT na;         /* lengths of grid intersection axes */
  IPNT gl_na;      /* lengths of grid intersection axes, local or global */
  FILE * fp = NULL;
  PARARRAY * par;

  /*  char * schar;*/
  float * fbuf;    /* input buffer for read */
  size_t recsize_b;/* length of 1D chunk to be read (bytes) */

  /* XDR buffers etc. */
#ifdef SUXDR
  XDR  xdrs;       /* xdr structure */
  int xflag;       /* i/o success flag */
  char * buf;      /* input buffer for xdr stream */
#endif

  /* scale flag - like scalco, scalel in SEGY; factor */
  int ntot;
  int scale=0;     /* flag - read from parameters */
  float scfac=1.0; /* factor workspace */

  /* grid offsets - allocated dynamically */
  off_t * goffs;   /* global offsets - into file data */
  off_t * loffs;   /* local offsets - into RARR data */
  size_t noffs;    /* number of offsets */

  /* possible nonzero offset for movie writes */
  off_t cur_pos = 0;

  /**************************
   * END DECLARATIONS       *
   **************************/

  /*  fprintf(stderr,"enter rsfwrite\n");*/

  /* HEADER FILE SECTION */

  par = ps_new();

  if (!ps_createfile(par,fname)) {
    //    fprintf(stderr,"in rsfwrite: par file = \n");
    //    ps_printall(*par,stderr);
    //    fprintf(stderr,"next\n");
    if ((ii=par_grid(&g,*par,stream))) {
      fprintf(stream,"Error rsfwrite err=%d\n",ii);
      fprintf(stream,"file %s exists but does not define RSF/SEP data structure\n",fname);
      return E_FILE;
    }

    if (ps_flcstring(*par,"in",&dname)) {
      fprintf(stream,"Error: rsfwrite from ps_flcstring\n");
      fprintf(stream,"failed to read data filename\n");
      return E_FILE;
    }
    if (ps_flcstring(*par,"data_format",&type)) {
      fprintf(stream,"Error: rsfwrite from ps_flcstring\n");
      fprintf(stream,"failed to read data_format from header file %s\n",fname);
      return E_FILE;
    }
    ps_flint(*par,"scale",&scale);
    //    fprintf(stderr,"rsfwrite: header=%s data=%s data_format=%s scale=%d parfile=\n",fname,dname,type,scale);
    //      ps_printall(*par,stderr);
  }
  else {
    fprintf(stream,"Error: rsfwrite - header file %s not opened\n",fname);
    return E_FILE;
  }

  /* DATA FILE SECTION */
        
  /* scaling loop note that scale should be inverse of read */
  /*  fprintf(stream,"gridio: SCALE=%d\n",scale);*/
  ntot=1;
  for (ii=0;ii<g.dim;ii++) ntot*=ran[ii];
  /* invert scale */
  scale=-scale;
  if (scale>0) for (ii=0;ii<scale;ii++) scfac *= 10.0;
  if (scale<0) for (ii=0;ii<-scale;ii++) scfac *= 0.10;
  if (scale) for (i=0;i<ntot;i++) {
    a[i]*=scfac;
  }

  /*  fprintf(stderr,"get global array data from grid \n");*/
  /* get global array data from grid 
     WWS 09.11.09 - gs is ALWAYS IPNT_0 */
  IASN(gs,IPNT_0);
  get_n(n,g);

  /* figger out array intersection. special case: if there is no
     intersection, read nearest part of data */
  for (ii=0; ii < g.dim; ii++) {
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

  /* use array data to work out offsets of 1D segments in both
     grid array (global) and rarray (local) */
  if ( (err = get_array_offsets(&goffs, &noffs, g.dim, gs, n, g_gsa, gl_na)) ||
       (err = get_array_offsets(&loffs, &noffs, g.dim, rags, ran, l_gsa, gl_na)) ) {
    fprintf(stream,"Error: rsfread - read from get_array_offsets, err=%d\n",err);
    //    ps_destroy(&par);
    ps_delete(&par);
    return err;
  }
  /*
  fprintf(stream,"noffs=%zd\n",noffs);
  for (i=0;i<noffs;i++) 
    fprintf(stream,"i=%zd goffs=%zd loffs=%zd\n",i,goffs[i],loffs[i]);
  */
  /*  fprintf(stderr,"from here on, memory is allocated, so save return until deallocation\n");*/

  /* adjoint extension loop - adj_extend by const along axes 0,...,dim-1.
     added 06.03.12 to make rsfwrite l2-adjoint to rsfread WWS
  */
  if (extend) {
    for (ii=0;ii<g.dim;ii++) {
      if (adj_extend_array(a,rags,ran,l_gsa,gl_na,g.dim,ii)) {
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
  fp=iwave_const_fopen(dname,"r+",NULL,stream);
#else
  p=fopen(dname,"r+")
#endif
  if (!fp) {
    fprintf(stream,"Error: rsfwrite from fopen: cannot open file=%s\n",dname);
    return E_FILE;
  }

  fflush(stream);

  /* compute current starting position */
  cur_pos = panelindex * get_datasize_grid(g) * sizeof(float);
  /* <-- D.S. 01.01.11: extended-model related */

  /*  fprintf(stderr,"native write loop\n");*/

  if (!strcmp(type,"native_float")) {

    /* allocate float buffer */
    fbuf=(float *)usermalloc_(recsize_b);

    for (i=0;i<noffs;i++) {
      /* seek to write segment */
      /*      fprintf(stderr,"seek trace %d\n",i);*/
      if (!err && fseeko(fp,goffs[i]*sizeof(float) + cur_pos,SEEK_SET)) {
	fprintf(stream,"Error: rsfwrite from fseeko at file offset %jd\n",(intmax_t)goffs[i]); 
	fprintf(stream,"possible cause: attempt to write off end of file\n");
	fseeko(fp,0L,SEEK_END);
	fprintf(stream,"file length = %jd:\n",(intmax_t)ftello(fp));
	fprintf(stream,"note that new file can only be written in contiguous,\n");
	fprintf(stream,"consecutive blocks\n");
	err=E_FILE;
      }
      /* convert from ireal */
      if (!err) {
	for (j=0;j<na[0];j++) { *(fbuf+j)=*(a+loffs[i]+j); }
      }
      /* write out float buffer */
      if (!err && (na[0] != fwrite(fbuf,sizeof(float),na[0],fp))) {
	fprintf(stream,"Error: rsfwrite from fwrite at array offset %jd\n",(intmax_t)loffs[i]);
	fprintf(stream,"failed to write %d words\n",na[0]);
	err=E_FILE;
      }
    }
    userfree_(fbuf);
  }

#ifdef SUXDR
  /* xdr write loop */
  else if (!strcmp(type,"xdr_float")) {

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
	fprintf(stream,"Error: rsfwrite from fseeko at file offset %jd\n",(intmax_t)goffs[i]);
	fprintf(stream,"possible cause: attempt to write off end of file\n");
	j=0;
	fseeko(fp,j,SEEK_END);
	fprintf(stream,"file length = %jd:\n",(intmax_t)ftello(fp));
	fprintf(stream,"note that new file can only be written in contiguous,\n");
	fprintf(stream,"consecutive blocks\n");
	err=E_FILE;
      }

      /* convert from ireal - need only one word of float storage */ 
      xflag = 0;
      if (!err) {
	xflag=1;
	for (j=0;j<na[0];j++) {
	  *fbuf=*(a+loffs[i]+j);
	  xflag=xflag && xdr_float(&xdrs,fbuf);
	}
      }
      /* write out byte array */
      if (!err && (recsize_b != fwrite(buf,sizeof(char),recsize_b,fp))) {
	fprintf(stream,"Error: rsfwrite from fwrite at array offset %jd\n",(intmax_t)loffs[i]);
	err=E_FILE;
      }
      
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
  userfree_(dname);
  userfree_(type);
  ps_delete(&par);
  fflush(fp);
#ifdef IWAVE_USE_FMGR
  iwave_fclose(fp);
#else
  fclose(fp);
#endif
  return err;
}

#endif

  
