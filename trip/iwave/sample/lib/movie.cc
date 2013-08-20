#include "movie.h"

int movie_init(MOVIE * mt, 
	       IMODEL *m, 
	       PARARRAY *par, 
	       tracegeom *tg,  
	       FILE * stream) {
  ireal dt;
  int i,j,k,ntot;
  int kl=KEYLEN;
  char key[KEYLEN];
  ireal slice;        /* actual slice coord, 3d case */
  RPNT d;             /* workspace for grid build-up */
  RPNT o;             /* workspace for grid build-up */
  IPNT n;             /* workspace for grid build-up */
  char * fname;       /* workspace for data file name */
  FILE * fp;          /* for output initializiation */
  float * buf;        /* null buffer for data file init */

  /* remind user that there are no 1D movies */
  if ((m->g).dim<2) {
    fprintf(stream,"NOTE: no movie created for 1D\n");
    return 0;
  }

  /* read movie choices */
  mt->maxslice=0;
  mt->nmovie=0;
  for (i=0;i<RDOM_MAX_NARR;i++) {
    j=snprintf(key,kl,"movie%d",i+1);
    if (j>kl || j<0) {
      fprintf(stream,"ERROR: movie_init from snprintf\n");
      return E_ALLOC;
    }
    if (!ps_flcstring(*par,key,&(fname))) {
      mt->imovie[mt->nmovie]=mt->iselect(fname);
      if (mt->imovie[mt->nmovie] < 0) {
	fprintf(stream,"NOTE: movie_init: key %s does not correspond to legit field, so no movie\n",fname);
      }
      else {
	mt->smovie[mt->nmovie] = (char *)usermalloc_(sizeof(char)*(KEYLEN+strlen(fname)));
	strcpy(mt->smovie[mt->nmovie],"movie_");
	strcat(mt->smovie[mt->nmovie],fname);
	strcat(mt->smovie[mt->nmovie],".rsf");
	rd_size(&(m->ld_a),mt->imovie[mt->nmovie],n);
	/* compute size of biggest 2d slice */
	if (m->g.dim > 1) {
	  for (k=0;k<RARR_MAX_NDIM;k++) {
	    ntot=1;
	    for (j=0;j<RARR_MAX_NDIM;j++) { if (k!=j) ntot*=n[j]; }
	    mt->maxslice=iwave_max(mt->maxslice,ntot);
	  }
	}
	mt->nmovie++;
      }
    }
  }  

  /* no movies - return */
  if (!mt->nmovie) return 0;

  /* defaults for 3D slice - last axis, middle slice */
  mt->dim3d=2;
  mt->slice3d = 0;
  /* default step = 0.0 - no movie! */
  dt=0.0;

  ps_flreal(*par, "moviestep", &dt);
  if (m->g.dim>2) {
    ps_flint(*par, "movieaxis3D",&(mt->dim3d));
    ps_flreal(*par, "movieslice3D",&slice);
  }

  /* sanity check - note that if */
  if (mt->dim3d < 0 || ((m->g.dim==3)&&(mt->dim3d > 2))) {
    fprintf(stream,"Error: movie_init - invalid dim3d=%d\n",mt->dim3d);
    return E_OTHER;
  }

  /* integerize step */  
  mt->framestep = (int)(0.1+(dt/tg->dt)); 

  /* if step<=0 then no movie - clean up */
  if (mt->framestep<1) {
    for (i=0;i<mt->nmovie;i++) {
      userfree_(mt->smovie[i]);
    }
    mt->nmovie=0;
    return 0;
  }

  /* integerize slice */
  if (m->g.dim>2) {
    get_d(d,m->g);
    get_o(o,m->g);
    mt->slice3d=(int)(0.1+((slice-o[mt->dim3d])/d[mt->dim3d]));
  }
  
  /* integerize start, stop */
  mt->framestart = (int)(((tg->t0)/((m->tsind).dt))+0.1);
  mt->framestop = mt->framestart + tg->nt-1;
  if (mt->framestep==0)
    mt->framestart=mt->framestop;
  else
    mt->framestart = (int)(((tg->t0)/m->tsind.dt)+0.1);
  
  mt->hname=(char *)usermalloc_(10*sizeof(char));
  strcpy(mt->hname,"./");

  /* compute next step */
  mt->it = mt->framestart;
  mt->frameindex=0;

  /* create movie grid */
  /* correction 03.12.12: actually movie frames are extended axis! */
  /* also movie dim is ALWAYS = 2 and movie gdim ALWAYS = 3 ! */
  /*
    if (m->g.dim > 2) (mt->mg).dim=m->g.dim;
    else (mt->mg).dim=m->g.dim+1;
  */
  (mt->mg).dim=2;
  (mt->mg).gdim=3;
  for (i=0;i<mt->dim3d; i++) {
    (mt->mg).axes[i].n  = (m->g).axes[i].n;
    (mt->mg).axes[i].d  = (m->g).axes[i].d;
    (mt->mg).axes[i].o  = (m->g).axes[i].o;    
    (mt->mg).axes[i].id = (m->g).axes[i].id;
  }
  for (i=mt->dim3d;i<(mt->mg).dim-1; i++) {
    (mt->mg).axes[i].n  = (m->g).axes[i+1].n;
    (mt->mg).axes[i].d  = (m->g).axes[i+1].d;
    (mt->mg).axes[i].o  = (m->g).axes[i+1].o;    
    (mt->mg).axes[i].id = (m->g).axes[i+1].id;
  }
  i=(mt->mg).gdim-1;
  (mt->mg).axes[i].n  = 1+(mt->framestop-mt->framestart)/mt->framestep;
  (mt->mg).axes[i].d  = dt;
  (mt->mg).axes[i].o  = tg->t0;
  /*  if ((mt->mg).dim>2) (mt->mg).axes[i].id = (m->g.axes[mt->dim3d].id);  */
  (mt->mg).axes[i].id = i+1;

  /* compute length of each frame, create null buffer */
  ntot=1;
  for (j=0;j<(mt->mg).dim-1;j++) ntot *= (mt->mg).axes[j].n;
  buf = (float *)usermalloc_(ntot*sizeof(float));
  memset(buf,0,ntot*sizeof(float));

  /* write out the header, data files (latter initialized to zero) */
  for (i=0;i<mt->nmovie;i++) {

    fp=NULL;
#ifdef IWAVE_USE_FMGR
    fp=iwave_const_fopen(mt->smovie[i],"w",NULL,stream);
#else
    fp=fopen(mt->smovie[i],"w");
#endif

    if (!fp) {
      fprintf(stream,"Error: movie_init\n");
      fprintf(stream,"failed to create movie header file %s\n",mt->smovie[i]);
      return E_FILE;
    }

    for (j=0;j<(mt->mg).dim;j++) {
      fprintf(fp,"n%d=%lu d%d=%g o%d=%g\n",
	      j+1,(mt->mg).axes[j].n,
	      j+1,(mt->mg).axes[j].d,
	      j+1,(mt->mg).axes[j].o);
    }
    fprintf(fp,"n%d=%lu d%d=%g o%d=%g\n",
	    (mt->mg).gdim,(mt->mg).axes[(mt->mg).gdim-1].n,
	    (mt->mg).gdim,(mt->mg).axes[(mt->mg).gdim-1].d,
	    (mt->mg).gdim,(mt->mg).axes[(mt->mg).gdim-1].o);

    fprintf(fp,"data_format=native_float\n");
    fprintf(fp,"scale=0\n");
    fprintf(fp,"dim=2\n");
    fname=(char *)usermalloc_(strlen(mt->smovie[i])+10);
    strcpy(fname,mt->smovie[i]);
    strcat(fname,"@");
    fprintf(fp,"in=%s\n",fname);

    fseek(fp,0L,SEEK_SET);

#ifdef IWAVE_USE_FMGR
    iwave_fclose(fp);
#else
    fclose(fp);
#endif

#ifdef IWAVE_USE_FMGR
    fp=iwave_const_fopen(fname,"w",NULL,stream);
#else
    fp=fopen(fname,"w");
#endif

    if (!fp) {
      fprintf(stream,"Error: movie_init\n");
      fprintf(stream,"failed to create movie data file %s\n",fname);
      return E_FILE;
    }
    
    /* fill data file with zeros */
    for (j=0;j<(mt->mg).axes[(mt->mg).gdim-1].n;j++) {
      if (ntot != fwrite(buf,sizeof(float),ntot,fp)) {
	fprintf(stream,"Error: movie_init\n");
	fprintf(stream,"failed to write %d bytes to data file %s\n",ntot,fname);
      }
    }

    fseek(fp,0L,SEEK_SET);

#ifdef IWAVE_USE_FMGR
    iwave_fclose(fp);
#else
    fclose(fp);
#endif
    userfree_(fname);
  }

  userfree_(buf);

  return 0;
}

int movie_run(MOVIE * mt, 
	      IMODEL * m, 
	      FILE * stream,
	      int _irec          /* D.S. 01.01.11: extended-model related */
	      ) {

  int i, k, ii, jj, kk;
  /* int j; */
  IPNT gs, ge, n, np;      /* workspace for size info */
  ireal * aptr, * aptr0;            /* target data pointer */
  off_t aoff;              /* offset into target data */
  IPNT gs_a, ge_a, gs_c, ge_c;
  int ndim;
  int err=0;
  
  /* time test - not there yet */
  if (m->tsind.it < mt->it) return 0;
  /* time test - beyond end */
  if (m->tsind.it > mt->framestop) return 0;
  /* time test - not at integral time */
  if (m->tsind.iv) return 0;

  /* at next step */

  for (i=0;i<mt->nmovie;i++) {

    /* gussy up gs and ge starting with domain */
    rd_size(&(m->ld_a),mt->imovie[i],n);
    rd_size(&(m->ld_a),mt->imovie[i],np);
    rd_gse(&(m->ld_a),mt->imovie[i],gs,ge);
    gs[mt->dim3d]=(mt->it-mt->framestart)/(mt->framestep);
    np[mt->dim3d]=1;

    /* if mt->dim3d==2 then life is easy */
    aptr=m->ld_a._s[mt->imovie[i]]._s0;

    /* XW 4-23-11: set value in non-computational domain zero */
    rd_gse(&(m->ld_c),mt->imovie[i],gs_c,ge_c);
    rd_gse(&(m->ld_a),mt->imovie[i],gs_a,ge_a);
    aptr0 = aptr;
    ndim = (m->g).dim;

    /* 1D */
    if (ndim == 1) {
      for (ii = gs_a[0];ii <= ge_a[0];ii ++) {
        if (ii < gs_c[0] || ii > ge_c[0])
          *aptr0 = REAL_ZERO;
        aptr0 ++;
      }
    }
    /* 2D */
    if (ndim == 2) {
      for (jj = gs_a[1];jj <= ge_a[1];jj ++) {
        for (ii = gs_a[0];ii <= ge_a[0];ii ++) {
          if (ii < gs_c[0] || ii > ge_c[0] || jj < gs_c[1] || jj > ge_c[1])
            *aptr0 = REAL_ZERO;
          aptr0 ++;
        }
      }
    }
    /* 3D */
    if (ndim == 3) {
      for (kk = gs[2];kk <= ge[2];kk ++) {
        for (jj = gs[1];jj <= ge[1];jj ++) {
          for (ii = gs[0];ii <= ge[0];ii ++) {
            if (ii < gs_c[0] || ii > ge_c[0] || jj < gs_c[1] || jj > ge_c[1] ||
                kk < gs_c[2] || kk > ge_c[2])
              *aptr0 = REAL_ZERO;
            aptr0 ++;
          }
        }
      }
    }

    if (ndim==2) {
	if ((err=rsfwrite(aptr, gs, n, mt->smovie[i], 0, stream, mt->frameindex))) {
	fprintf(stream,"ERROR: movie_run from rsfwrite, err=%d\n",err);
	fflush(stream);
	return err;
      }
    }
    else if (mt->dim3d==2) {
      aoff=mt->slice3d*n[0]*n[1];
      if ((err=rsfwrite(&(aptr[aoff]), gs, np, mt->smovie[i], 0, stream, mt->frameindex))) {
	fprintf(stream,"ERROR: movie_run from rsfwrite, err=%d\n",err);
	fflush(stream);
	return err;
      }
    }
    /* otherwise must copy data to buffer first */
    else if (mt->dim3d==1) {
      if (!mt->buf) mt->buf=(ireal *)usermalloc_(mt->maxslice*sizeof(ireal));
      for (i=0;i<n[2];i++) {
	for (k=0;i<n[0];i++) {
	  (mt->buf)[k+i*n[0]]=aptr[k+mt->slice3d*n[0]+i*n[0]*n[1]];
	}
      }
    }
    /* mt->dim3d=0 */
    else { 
      if (!mt->buf) mt->buf=(ireal *)usermalloc_(mt->maxslice*sizeof(ireal));      
      for (i=0;i<n[2];i++) {
	for (k=0;k<n[1];k++) {
	  (mt->buf)[k+i*n[1]]=aptr[mt->slice3d+k*n[0]+i*n[0]*n[1]];
	}
      }
      /* D.S. 01.01.11: extended-model related --> */
      rsfwrite(mt->buf, gs, np, mt->smovie[i], 0, stream, mt->frameindex); 
      /* <-- extended-model related */
    }
	
  }

  /* advance next time */
  mt->it += mt->framestep;
  mt->frameindex++;

  return 0;
}	 

int movie_setnull(MOVIE * mt) {
  int i;
  mt->hname = NULL;
  for (i=0;i<RDOM_MAX_NARR;i++) {
    mt->imovie[i]     = 0;
    mt->smovie[i]     = NULL;
  }
  mt->nmovie        = 0;
  mt->framestep     = 0;
  mt->framestart    = 0;
  mt->framestop     = 0;
  mt->it            = 0;
  mt->dim3d         = 0;
  mt->slice3d       = 0;			
  mt->maxslice      = 0;
  mt->buf=NULL;
  
  return 0;
}

int movie_destroy(MOVIE * mt) {
  int i;
  if (mt->hname) userfree_(mt->hname);
  for (i=0;i<mt->nmovie;i++) {
    if (mt->smovie[i]) userfree_(mt->smovie[i]);
  }
  if (mt->buf) userfree_(mt->buf);
  movie_setnull(mt);
  return 0;
}

void movie_fprint(MOVIE * mt, FILE * fp) {

  int i;
  if (mt->nmovie>0) {
    fprintf(fp,"/*---------------------------------------------------------*/\n");
    fprintf(fp,  "IWAVE MOVIE WRITER\n");
    fprintf(fp,  "number of movies = %d\n",mt->nmovie);
    for (i=0;i<mt->nmovie;i++) 
      fprintf(fp,"%d file = %s field index = %d\n",
	      i,mt->smovie[i],mt->imovie[i]);
    fprintf(fp,  "framestart       = %d\n", mt->framestart);
    fprintf(fp,  "framestop        = %d\n", mt->framestop);
    fprintf(fp,  "framestep        = %d\n", mt->framestep);
    fprintf(fp,  "next step        = %d\n", mt->it);
    fprintf(fp,  "normal axis      = %d\n", mt->dim3d);
    fprintf(fp,  "slice index      = %d\n", mt->slice3d);
    fprintf(fp,  "hname            = %s\n", mt->hname);
    fflush(fp);
  }

} 

