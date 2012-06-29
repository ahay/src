#include "pointsrc.h"

#define DT_DBG

/* point sampling */
int pointsource(IPNT is,
		RPNT rs,
		int order,
		ireal scramp,
		ireal src,
		RARR arr) {

  IPNT gs, ge;
  ireal fac;
  int ndim;
  ireal p;
  
  FILE *stream;
  
  stream = retrieveOutstream();

  ra_ndim(&arr, &ndim);
  ra_gse(&arr, gs, ge);
  
  /* order 0 option */
  if ( order==0 ) {
    if ((ndim == 1 &&
	 is[0]>gs[0]-1 &&
	 is[0]<ge[0]+1) ||
	(ndim == 2 && 
	 is[0]>gs[0]-1 &&
	 is[0]<ge[0]+1 && 
	 is[1]>gs[1]-1 &&
	 is[1]<ge[1]+1) ||
	(ndim == 3 && 
	 is[0]>gs[0]-1 &&
	 is[0]<ge[0]+1 &&
	 is[1]>gs[1]-1 &&
	 is[1]<ge[1]+1 &&
	 is[2]>gs[2]-1 &&
	 is[2]<ge[2]+1)) {
      p = ra_gget(&arr, is);
      ra_gset(&arr, is, p + src * scramp);
    }
  }
  else if ( order==1 ) {
    if (ndim==1) {
      /* 0,0,0 */
      if (is[0]>gs[0]-1 &&
	  is[0]<ge[0]+1 ) {
	p = ra_gget(&arr, is);
	fac=1.0-rs[0];
	ra_gset(&arr, is, p + fac * src * scramp); 
      }
      /* 1,0,0 */
      is[0]++;
      if (is[0]>gs[0]-1 &&
	  is[0]<ge[0]+1 ) {
	p = ra_gget(&arr, is);
	fac=rs[0];
	ra_gset(&arr, is, p + fac * src * scramp); 
      }
      is[0]--;
    }
    if (ndim==2) {
      /* 0,0,0 */
      if (is[0]>gs[0]-1 &&
	  is[1]>gs[1]-1 &&
	  is[0]<ge[0]+1 &&
	  is[1]<ge[1]+1) {
	p = ra_gget(&arr, is);
	fac=(1.0-rs[0])*(1.0-rs[1]);
	ra_gset(&arr, is, p + fac * src * scramp); 
      }
      /* 1,0,0 */
      is[0]++;
      if (is[0]>gs[0]-1 &&
	  is[1]>gs[1]-1 &&
	  is[0]<ge[0]+1 &&
	  is[1]<ge[1]+1) {;
	p = ra_gget(&arr, is);
	fac=rs[0]*(1.0-rs[1]);
	ra_gset(&arr, is, p + fac * src * scramp);
      }
      is[0]--;
      /* 0,1,0 */
      is[1]++;
      if (is[0]>gs[0]-1 &&
	  is[1]>gs[1]-1 &&
	  is[0]<ge[0]+1 &&
	  is[1]<ge[1]+1) {
	p = ra_gget(&arr, is);
	fac=(1.0-rs[0])*rs[1];
	ra_gset(&arr, is, p + fac * src * scramp);
      }
      is[1]--;
      /* 1,1,0 */
      is[0]++;is[1]++;
      if (is[0]>gs[0]-1 &&
	  is[1]>gs[1]-1 &&
	  is[0]<ge[0]+1 &&
	  is[1]<ge[1]+1) {
	p = ra_gget(&arr, is);
	fac=rs[0]*rs[1];
	ra_gset(&arr, is, p + fac * src * scramp);
      }
      is[0]--;is[1]--;
    }
    if (ndim==3) {

      /* dy=0 plane */
      /* 0,0,0 */
      if (is[0]>gs[0]-1 &&
	  is[1]>gs[1]-1 &&
	  is[2]>gs[2]-1 &&
	  is[0]<ge[0]+1 &&
	  is[1]<ge[1]+1 &&
	  is[2]<ge[2]+1) {
	p = ra_gget(&arr, is);
	fac=(1.0-rs[0])*(1.0-rs[1])*(1.0-rs[2]);
	ra_gset(&arr, is, p + fac * src * scramp); 
      }
      /* 1,0,0 */
      is[0]++;
      if (is[0]>gs[0]-1 &&
	  is[1]>gs[1]-1 &&
	  is[2]>gs[2]-1 &&
	  is[0]<ge[0]+1 &&
	  is[1]<ge[1]+1 &&
	  is[2]<ge[2]+1) {
	p = ra_gget(&arr, is);
	fac=rs[0]*(1.0-rs[1])*(1.0-rs[2]);
	ra_gset(&arr, is, p + fac * src * scramp); 
      }
      is[0]--;
      /* 0,1,0 */
      is[1]++;
      if (is[0]>gs[0]-1 &&
	  is[1]>gs[1]-1 &&
	  is[2]>gs[2]-1 &&
	  is[0]<ge[0]+1 &&
	  is[1]<ge[1]+1 &&
	  is[2]<ge[2]+1) {
	p = ra_gget(&arr, is);
	fac=(1.0-rs[0])*rs[1]*(1.0-rs[2]);
	ra_gset(&arr, is, p + fac * src * scramp); 
      }
      is[1]--;
      /* 1,1,0 */
      is[0]++;is[1]++;
      if (is[0]>gs[0]-1 &&
	  is[1]>gs[1]-1 &&
	  is[2]>gs[2]-1 &&
	  is[0]<ge[0]+1 &&
	  is[1]<ge[1]+1 &&
	  is[2]<ge[2]+1) {
	p = ra_gget(&arr, is);
	fac=rs[0]*rs[1]*(1.0-rs[2]);
	ra_gset(&arr, is, p + fac * src * scramp); 
      }
      is[0]--;is[1]--;

      /* dy=1 plane */
      is[2]++;
      /* 0,0,1 */
      if (is[0]>gs[0]-1 &&
	  is[1]>gs[1]-1 &&
	  is[2]>gs[2]-1 &&
	  is[0]<ge[0]+1 &&
	  is[1]<ge[1]+1 &&
	  is[2]<ge[2]+1) {
	p = ra_gget(&arr, is);
	fac=(1.0-rs[0])*(1.0-rs[1])*rs[2];
	ra_gset(&arr, is, p + fac * src * scramp); 
      }
      /* 1,0,1 */
      is[0]++;
      if (is[0]>gs[0]-1 &&
	  is[1]>gs[1]-1 &&
	  is[2]>gs[2]-1 &&
	  is[0]<ge[0]+1 &&
	  is[1]<ge[1]+1 &&
	  is[2]<ge[2]+1) {
	p = ra_gget(&arr, is);
	fac=rs[0]*(1.0-rs[1])*rs[2];
	ra_gset(&arr, is, p + fac * src * scramp); 
      }
      is[0]--;
      /* 0,1,1 */
      is[1]++;
      if (is[0]>gs[0]-1 &&
	  is[1]>gs[1]-1 &&
	  is[2]>gs[2]-1 &&
	  is[0]<ge[0]+1 &&
	  is[1]<ge[1]+1 &&
	  is[2]<ge[2]+1) {
	p = ra_gget(&arr, is);
	fac=(1.0-rs[0])*rs[1]*rs[2];
	ra_gset(&arr, is, p + fac * src * scramp); 
      }
      is[1]--;
      /* 1,1,1 */
      is[0]++;is[1]++;
      if (is[0]>gs[0]-1 &&
	  is[1]>gs[1]-1 &&
	  is[2]>gs[2]-1 &&
	  is[0]<ge[0]+1 &&
	  is[1]<ge[1]+1 &&
	  is[2]<ge[2]+1) {
	p = ra_gget(&arr, is);
	fac=rs[0]*rs[1]*rs[2];
	ra_gset(&arr, is, p + fac * src * scramp); 
      }
      is[0]--;is[1]--;
      is[2]--;
    }
  }
 
  return 0;
}

int pointsrc_init(POINTSRC * tr, IMODEL * m, PARARRAY * par, tracegeom *tg, FILE * stream) {

  int err = 0;                /* error flag */
  char * wp;                /* wavelet phase */
  int iw;                   /* half-width */
  int i;                    /* counter */
  int ndim;                 /* dimension of grids */
  RPNT d;                   /* steps */
  ireal prod_d;             /* cell volume scale factor for delta */
  char *srcfile;            /* workspace for file name */
  segy trsrc;               /* segy workspace for reading file */
  Value val;                /* workspace for reading segy headers */
  int tmpnt;                /* length of time series read from file */
  int lnt;                  /* length of extended time series, for integration */
  ireal tmpdt;              /* time step for time series read from file */
  ireal tmpt0;              /* time origin for time series read from file */
  ireal tmax;               /* max time for either trace or source time series */
  int wl;                   /* length of workspace for cubic spline interp */
  ireal * wk;               /* workspace for cubic spline interp */
  int iend = 1;             /* end condition for cubic spline interp */
  ireal t0;                 /* time origin on model dt grid */
  ireal tmp0, tmp1, q;      /* workspace for in-place trapezoidal rule */
  IPNT tis;                 /* buffer for coefficient sampling near source */
  int iflag;                /* test flag */
  ireal tdt;                /* ireal buffer for dt */

  ireal refvel;             /* reference velocity for target prop wavelet */
  ireal refbou;		    /* reference buoyancy for target prop wavelet */
  ireal refkappa;           /* near-source bulk modulus */
  ireal refdist;            /* reference distance for target prop wavelet */
  ireal refamp;             /* reference amplitude for target prop wavelet */
  ireal fpeak;              /* peak frequency for Ricker computation */

  segy trdbg;               /* workspace for building output segy */
  IPNT gs, ge;              /* workspace for buoyancy exchange comps */
  ireal *resc;

  /* MPI workspace */
#ifdef IWAVE_USE_MPI
  int rk, sz;
  MPI_Comm cm;
  ireal *procbuf;
	
#endif
  
  /* end declarations */

  stream = retrieveOutstream();

  /* assign default reference values */
  refvel   = CREF_DEF;
  refbou   = REAL_ZERO;
  refkappa = REAL_ZERO;
  refdist  = RREF_DEF;
  refamp   = REAL_ONE;
  fpeak    = FPEAK_DEF;        

  /* get array of grid steps */
  get_d(d, m->gl);

  /* extract dimension */
  rd_ndim(&m->ld_a, D_MP0, &ndim);

  /* TV */
  IASN(tr->is,IPNT_0);
  RASN(tr->rs,RPNT_0);
  tr->is[0]=tg->is[0]; tr->rs[0]=tg->rs[0]; tis[0]=tr->is[0];
  if (ndim > 1) { tr->is[1]=tg->is[1]; tr->rs[1]=tg->rs[1]; tis[1]=tr->is[1]; }
  if (ndim > 2) { tr->is[2]=tg->is[2]; tr->rs[2]=tg->rs[2]; tis[2]=tr->is[2]; } 

  /* read sampling order */
  tr->order = 0;
  ps_flint(*par, "sampord", &(tr->order));

  /* either read from file, or it's a gaussian */
  tr->fpsrc = NULL;
  tr->fpdbg = NULL;

  /* extract near-source bulk modulus from grid - necessary for
     several (but not all) cases, and not enough of an expense to 
     avoid.

     WWS, 29.11.08: to avoid failure if integer part of source point is in 
     physical grid (including physical boundary points) but not in comp.
     grid, check neighboring grid points - in that case, at least one
     of these should be in part of comp grid, where kappa is initialized.
  */

  rd_gse(&(m->ld_a), D_MP0, gs, ge);
  for ( i = 0; i < ndim; ++i ) {
    if ( (tr->is[i] < gs[i]) && (tr->is[i]+1 > gs[i]-1) ) tis[i]++;
    if ( (tr->is[i] > ge[i]) && (tr->is[i]-1 < ge[i]+1) ) tis[i]--;
  } 
  iflag=1;
  for (i = 0; i < ndim; ++i) 
    if ( (tis[i] < gs[i]) || (tis[i] > ge[i]) ) iflag=0;

  if (iflag) refkappa = rd_gget(&(m->ld_a), D_MP0, tis);

  fprintf(stream,"NOTE: proc=%d sample kappa at ",retrieveRank());
  for (i=0;i<ndim;++i) 
    fprintf(stream,"index[%d]=%d ",i,tis[i]);
  fprintf(stream,"sample flag=%d\n",iflag);

  /* extract near-source buoyancy from grid - necessary for
     several (but not all) cases, and not enough of an expense to 
     avoid.

     WWS, 04.03.09: need to do this here too - to avoid failure if
     integer part of source point is in physical grid (including
     physical boundary points) but not in comp.  grid, check
     neighboring grid points - in that case, at least one of these
     should be in part of comp grid, where kappa is initialized.
  */
  
  rd_gse(&(m->ld_a), D_MV0, gs, ge);
  for ( i = 0; i < ndim; ++i ) {
    tis[i] = tr->is[i];
    if ( (tr->is[i] < gs[i]) && (tr->is[i]+1 > gs[i]-1) ) tis[i]++;
    if ( (tr->is[i] > ge[i]) && (tr->is[i]-1 < ge[i]+1) ) tis[i]--;
  } 
  iflag=1;
  for (i = 0; i < ndim; ++i) 
    if ( (tis[i] < gs[i]) || (tis[i] > ge[i]) ) iflag=0;

  if (iflag) refbou = rd_gget(&(m->ld_a), D_MV0, tis);
  fprintf(stream,"NOTE: proc=%d sample buoyancy at ",retrieveRank());
  for (i=0;i<ndim;++i) 
    fprintf(stream,"index[%d]=%d ",i,tis[i]);
  fprintf(stream,"sample flag=%d\n",iflag);

#ifdef IWAVE_USE_MPI
  rk = retrieveRank();
  sz = retrieveSize();
  cm = retrieveComm();
  if ( rk == 0 ) {
    procbuf = (ireal*)usermalloc_(sz * sizeof(ireal));
    if ( procbuf == NULL ) return E_ALLOC;
  }
  MPI_Gather(&refkappa, 1, IWAVE_MPI_REAL, procbuf, 1, IWAVE_MPI_REAL, 0, cm);
  if(rk ==0 ) {
    refkappa = 0.0;
    for ( i = 0; i < sz; ++i ) refkappa=iwave_max(refkappa,procbuf[i]);
  }
  MPI_Bcast(&refkappa, 1, IWAVE_MPI_REAL, 0, cm);
  MPI_Gather(&refbou, 1, IWAVE_MPI_REAL, procbuf, 1, IWAVE_MPI_REAL, 0, cm);
  if(rk ==0 ) {
    refbou = 0.0;
    for ( i = 0; i < sz; ++i ) refbou=iwave_max(refbou,procbuf[i]);
    userfree_(procbuf);
  }
  MPI_Bcast(&refbou, 1, IWAVE_MPI_REAL, 0, cm);
#endif	

  if (refbou > REAL_ZERO) {
    fprintf(stream,"NOTE: in pointsrc, using  buoyancy at source location = %e\n", refbou);
  }
  else {
    fprintf(stream,"ERROR: in pointsrc, ref buoyancy nonpositive, = %e\n",refbou);
    return E_OUTOFBOUNDS;
  }

  if (refkappa > REAL_ZERO) {
    fprintf(stream,"NOTE: in pointsrc, using  bulk mod at source location = %e\n", refkappa);
  }
  else {
    fprintf(stream,"ERROR: in pointsrc, ref bulk mod nonpositive, = %e\n",refkappa);
    return E_OUTOFBOUNDS;
  }

  /* get reference velocity, either from parameters or from buoyancy
     and bulk modulus near source location. Since velocity is not
     stored in RDOM, must be computed from buoyancy and bulk mod at
     slightly different points - this of course does not matter if
     source is located in homogeneous region.
  */
  
  if (ps_flreal(*par,"refvel",&refvel)) {  

    /* compute velocity from buoyancy and bulk modulus */

    refvel = sqrt(refkappa * refbou);
    
    fprintf(stream,"NOTE: in pointsrc, using velocity computed from \n");
    fprintf(stream,"      bulk mod and buoyancy near source location; \n");
    fprintf(stream,"      computed value = %e\n", refvel);

  }
  else {
    fprintf(stream,"NOTE: in pointsrc, using velocity from param table = %e\n", refvel);
  }

  /* Either read reference distance from parameters, or use default. */
  ps_flreal(*par,"refdist",&refdist);
  if (refdist>0) {
    fprintf(stream,"NOTE: in pointsrc, using reference distance = %e\n", refdist);
  }
  else {
    fprintf(stream,"NOTE: in pointsrc, read nonpos. reference distance = %e\n", refdist);
    fprintf(stream,"      this implies that wavelet will be read from file and used \n");
    fprintf(stream,"      directly on RHS as multiplier of spatial delta, rather than\n");
    fprintf(stream,"      to produce target propagating pulse.\n");
  }

  /* Either read reference amplitude from parameters, or use default. */
  
  ps_flreal(*par,"refamp",&refamp);
  fprintf(stream,"NOTE: in pointsrc, using reference amplitude = %e\n", refamp);

  /* read peak frequency from parameters, or use default (only used in
     Gaussian option II) */
  if (ps_flreal(*par,"fpeak", &fpeak)) {
    fprintf(stream,"NOTE: pointsrc_init - using default ");
    fprintf(stream,"peak frequency (fpeak) = %e\n",fpeak);
  }
    
  /* Overall scale factor for source insertion to produce target pulse, per 
     paper by WWS&TV: contains factors of
     - 4 pi c^2 dt (from rhs of equation 9, using kappa/rho=c^2 - RHS of 
       difference scheme is multiplied by kappa * dt);
     - r (reference distance for normalization, per eqn 13);
     - reference amplitude;
     - reciprocal of cell volume, for delta function.
     Note: for Option Ib (direct insertion of source wavelet, signalled by 
     refdist<=0) this factor accounts only for time step, bulk modulus, and 
     cell volume, as this case defines no target scale.
  */
  prod_d = 1.0;
  for (i = 0; i <  ndim; i++) prod_d *= d[i];
  if (refdist>0) {
    tr->scramp =  4.0 * 3.1415927 * refvel * refvel * 
      refdist * refamp * ((m->tsind).dt) / prod_d;
  }
  else {
    tr->scramp =  refkappa * refamp * ((m->tsind).dt) / prod_d;
  }

  /* Option I: read source from file */

  if (!ps_flcstring(*par,"source", &srcfile))  {
    if (!(tr->fpsrc = iwave_const_fopen(srcfile, "r",NULL,stream))) {
      fprintf(stream, "Error: pointsrc_init - failed to open source file\n");
      return E_FILE;
    }
    if (fseek(tr->fpsrc,0L,SEEK_SET)) {
      fprintf(stream,"Error: pointsrc_init - failed to seek to start of file\n");
      return E_FILE;
    }
    if (!fgettr(tr->fpsrc, &trsrc)) {
      fprintf(stream,"Error: pointsrc_init - failed to read source file\n");
      return E_FILE;
    }
    iwave_fclose(tr->fpsrc);

    /* at this point, the data member of trsrc contains the source wavelet
       at an external sample rate - read headers relevant for single dilat
       point source.
    */
    gethdval(&trsrc, "ns", &val);
    tmpnt = vtoi(hdtype("ns"), val);
    gethdval(&trsrc, "dt", &val);
    tmpdt = 0.001 * vtof(hdtype("dt"), val);
    gethdval(&trsrc, "delrt", &val);
    tmpt0 = vtof(hdtype("delrt"), val);		
        
    /* calculate istart, length of resamp wavelet,
       allocate buffer to hold it.
    */
    tr->istart = (int)(tmpt0/((m->tsind).dt));
    t0 = (m->tsind).dt * tr->istart;
    /*    tr->n = (int)(tmpnt*tmpdt/((m->tsind).dt))+1;*/
    /* final time for ext source set to max (final time for
       trace, final time for input source) */
    tmax = iwave_max(tmpt0 + tmpnt * tmpdt, tg->t0 + tg->nt * ((m->tsind).dt));
    /* ext src array length set to number of samples from istart
       to max time */
    tr->n = (int)(tmax/((m->tsind).dt)) + 1 - tr->istart;
    tr->w = (ireal *)usermalloc_(sizeof(ireal)*(tr->n));
    
    /* allocate buffer for integrated wavelet 
       at input sample rate */
    lnt  = (int)(tr->n * ((m->tsind).dt) / tmpdt) + 1;
    resc = (ireal *)usermalloc_(sizeof(ireal) * lnt);
    for (i = 0; i < tmpnt; i++) resc[i] = trsrc.data[i];
    for (i = tmpnt; i < lnt; i++) resc[i] = trsrc.data[tmpnt-1];

    /* interpolation workspace */
    wl = cubic_getworksize(lnt);
    wk=(ireal *)usermalloc_(sizeof(ireal) * wl);
    if (!wk) return E_ALLOC;
    
    /* Option Ia: normalize to produce propagating wavelet at distance
       refdist, speed refvel, if refdist > 0. */
    
    if (refdist > 0) {

      fprintf(stream,"NOTE: in pointsrc_init, compute wavelet for target pulse (Option Ia)\n");

      /* integrate once - trapezoidal rule in-place 
	 w[j] <- sum_{i=1}^{i=j} 0.5*dt*(w[i-1]+w[i])
	 leave multiplication by dt for scaling step
      */
      
      tmp0 = resc[0];
      q = 0.0;
      resc[0] = 0.0;
      for (i = 1; i < lnt; i++) {
	tmp1 = resc[i];
	q += 0.5 * tmpdt * (tmp0 + tmp1);
	resc[i] = q;
	tmp0 = tmp1;
      }
      
      tdt = (ireal)((m->tsind).dt);
      /* interpolate */
      if ((err=cubic_(&tmpt0, &tmpdt, resc, &lnt, &t0, &(tdt), tr->w, &(tr->n), &iend,wk,&wl))) {
	fprintf(stream,"Error: pointsrc_init - from cubic\n");
	return err;
      }
    }

    /* Option Ib: use wavelet as is, if refdist <= 0 */

    else {

      fprintf(stream,"NOTE: in pointsrc_init, using file wavelet directly on RHS (Option Ib)\n");
      tdt = (ireal)((m->tsind).dt);

      /* interpolate */
      if ((err=cubic_(&tmpt0, &tmpdt, resc, &lnt, &t0, &(tdt), tr->w,&(tr->n),&iend,wk,&wl))) {
	fprintf(stream,"Error: pointsrc_init - from cubic\n");
	return err;
      }
    }
    
    /* clean up */
    userfree_(wk);
    userfree_(srcfile);
    userfree_(resc);
  }

  /* Option II: create Gaussian derivative wavelet to produce Ricker
     propagating pulse with specified peak frequency, amplitude at
     spec'd distance. */

  else {

    /* check that reference distance is positive - only legal option */
    if (!(refdist>0)) {
      fprintf(stream,"Error: pointsrc_init, gaussian case:\n");
      fprintf(stream,"negative reference distance = %e no legal in this case\n",refdist);
      return E_OTHER;
    }
     
    fprintf(stream,"NOTE: in pointsrc_init, computing Gaussian derivative RHS for target Ricker pulse at f=%e, r=%e\n",fpeak,refdist);
 
    /* RHS wavelet is derivative of gaussian */
    tr->w = igetdgauss(&iw, (m->tsind).dt, fpeak);
    tr->n = 2 * iw+1;
    
    /* source phase - default is zero-phase */
    tr->istart = -iw;
    if (!ps_flcstring(*par,"waveletphase",&wp)) {
      if (!strcmp(wp,"zerophase")) tr->istart=-iw;
      else if (!strcmp(wp,"causal")) tr->istart=0;
      else if (!strcmp(wp,"anticausal")) tr->istart=-2*iw;
      else {
	fprintf(stream,"Error: pointsrc_init, gaussian case:\n");
	fprintf(stream,"only legit values of waveletphase are \n");
	fprintf(stream,"zerophase, causal, anticausal\n");
	return E_OTHER;
      }
      userfree_(wp);
    }
  }

  /* optionally write out wavelet appearing on RHS */
  tr->idbg = 0;
  ps_flint(*par, "dump_wavelet", &(tr->idbg));

  if (tr->idbg) {
    memcpy(trdbg.data,tr->w,tr->n*sizeof(ireal));
    val.u=1000.0*((m->tsind).dt);
    puthdval(&trdbg,"dt",&val);
    val.h=tr->n;
    puthdval(&trdbg,"ns",&val);
    val.h=((m->tsind).dt)*tr->istart;
    puthdval(&trdbg,"delrt",&val);
    
    if (!(tr->fpdbg=iwave_const_fopen("wavelet.debug","w",NULL,stream))) {
      fprintf(stream,"Error: init_point: failed to open test wavelet file\n");
      return E_FILE;
    }
    else {
      fprintf(stream,"write wavelet trace\n");
    }
    fputtr(tr->fpdbg,&trdbg);
    fflush(tr->fpdbg);
    iwave_fclose(tr->fpdbg);
  }

  return 0;
  
}

int pointsrc_destroy(POINTSRC * tr) {
  if (tr->w)       userfree_(tr->w); 
  if ( tr->fpsrc ) iwave_fclose(tr->fpsrc);
  if ( tr->fpdbg ) iwave_fclose(tr->fpdbg);
  return 0;
}

int pointsrc_run(POINTSRC * tr, IMODEL * m) {

  int i;

  /* key dimn off pressure field - NO-OP if iv!=0 */
  if ( ((m->tsind).it >= tr->istart) && 
       ((m->tsind).it < tr->istart + tr->n) &&
       ((m->tsind).iv == 0) ) {
    for (i=0; i < (m->ld_a)._s[D_P0].ndim; i++) {
      pointsource(tr->is,
		  tr->rs,
		  tr->order,
		  tr->scramp,
		  (tr->w)[(m->tsind).it-tr->istart],
       		  (m->ld_c)._s[D_P[i]]);
    }
  }
  return 0;
}

void pointsrc_fprint(POINTSRC const * tr, FILE * fp) {
  int i;
  fprintf(fp,"/*---------------------------------------------------------*/\n");
  fprintf(fp,"POINT SOURCE\n");
  fprintf(fp,"pulse length = %d\n",tr->n);
  fprintf(fp,"istart       = %d\n",tr->istart);
  fprintf(fp,"order        = %d\n",tr->order);
  for (i=0;i<RARR_MAX_NDIM;i++)
      fprintf(fp,"is[%d]=%d rs[%d]=%e\n", i, tr->is[i], i, tr->rs[i]);
    
  fflush(fp);
}
