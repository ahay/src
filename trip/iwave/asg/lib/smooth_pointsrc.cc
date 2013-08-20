#include "smooth_pointsrc.h"

#define DT_DBG

static inline ireal cubic_interp(ireal t, ireal dt, ireal *y, int n, ireal t0)
{
	register int i, i0, i1, i2, i3;
	register ireal dt3, ti0, ti1, ti2, ti3, l0, l1, l2, l3;
	
  t -= t0;
  
	dt3 = dt * dt * dt;
  i = (int)(floor( t / dt) );
  
  if ( (i < 0) || (i >= n-1) ) return 0.0;
  
  if (i <= 0) 
  {
    i0 = 0;
    i1 = 1;
    i2 = 2;
    i3 = 3;
  }
  else if (i >= n-2)
  {
    i0 = n - 4;
    i1 = n - 3;
    i2 = n - 2;
    i3 = n - 1;
  }
  else
  {
    i0 = i - 1;
    i1 = i;
    i2 = i + 1;
    i3 = i + 2;
  }
  
  ti0 = i0 * dt;
  ti1 = i1 * dt;
  ti2 = i2 * dt;
  ti3 = i3 * dt;
			
  l0 = -						(t - ti1 )* (t - ti2) * (t - ti3) / 6.0 / dt3;
  l1 =  (t - ti0)							*	(t - ti2) * (t - ti3) / 2.0 / dt3;
  l2 = -(t - ti0) * (t - ti1 )						* (t - ti3) / 2.0 / dt3;
  l3 =  (t - ti0) * (t - ti1 )* (t - ti2)             / 6.0 / dt3;

  return y[i0] * l0 + y[i1] * l1 + y[i2] * l2 + y[i3] * l3;
}

static inline ireal pfun(IPNT ind, RPNT d,  int ndim,  SPOINTSRC *tr, ireal t, int iv, ireal dt)
{
	register int i;
	register ireal r, x, xd=0.0, p=0.0, gradphi;
  register ireal pi = 3.1415927;
  
  r = 0.0;
	for ( i = 0; i < ndim; ++i )
	{
    x = (ireal)(ind[i]);
		if ( i == iv ) x += 0.5;
		x = x * d[i] - tr->xs[i];
    if (i == iv)  xd = x;
		r += x * x;
	}  
  
  r = sqrt(r);
  
  if (r < tr->rad)  return 0.0;
	
  gradphi = -compdgauss(r - tr->rad,  tr->phipeak )
          * (2.0 * pi * pi * tr->phipeak * tr->phipeak)  * xd / r;
  
	if (tr->tmpflag == 3) /* Option Ia */
	{
		p = cubic_interp( t - r / tr->c, tr->dt, tr->w, tr->n, tr->t0)/ r;
	}
  else if (tr->tmpflag == 2) /* Option Ib */
  {
    p = cubic_interp(t - r / tr->c, tr->dt, tr->w1, tr->n, tr->t0)/ r;

  }
	else if (tr->tmpflag == 1)  /* Ricker */
	{
		p = comprick(t + tr->istart * dt - r / tr->c, tr->fpeak) / r;
	}
	return p * gradphi;
}

static inline ireal vfun(IPNT ind, RPNT d, int ndim, SPOINTSRC *tr, ireal t, ireal dt)
{
	register int i;
	register ireal r, x, gradphi, f0=0.0, f1=0.0;
  ireal pi = 3.1415927;
    	
	r = 0.0;
	for ( i = 0; i < ndim; ++i )
	{
		x = ((ireal)(ind[i])) * d[i] - tr->xs[i];
		r += x * x;
	}
  r = sqrt(r);
  
  if (r < tr->rad)  return 0.0;
  
  gradphi = -compdgauss(r - tr->rad,  tr->phipeak)
          * (2.0 * pi * pi * tr->phipeak * tr->phipeak) / r;
                  
	if (tr->tmpflag == 3)  /* Option Ia */
	{
		f0 = cubic_interp( t - r / tr->c, tr->dt, tr->w, tr->n, tr->t0) / tr->c;
		f1 = cubic_interp( t - r / tr->c, tr->dt, tr->w1, tr->n, tr->t0) / r;
	}
  else if (tr->tmpflag == 2)  /* Option Ib */
  {
    f0 = cubic_interp( t - r / tr->c, tr->dt, tr->w1, tr->n, tr->t0) / tr->c;
		f1 = cubic_interp( t - r / tr->c, tr->dt, tr->w, tr->n, tr->t0) / r;
  }
	else if (tr->tmpflag == 1) /* Ricker */
	{
		f0 = comprick( t + tr->istart * dt - r / tr->c, tr->fpeak) / tr->c;
		f1 = compdgauss( t + tr->istart * dt - r / tr->c, tr->fpeak) / r;
	}
	
	return tr->bou * ( f0 + f1) * gradphi;;
}

int spointsrc_init(SPOINTSRC * tr, IMODEL * m, PARARRAY * par, tracegeom *tg, FILE * stream) 
{

  char * wp;                /* wavelet phase */
  int iw;                   /* half-width */
  int i;                    /* counter */
  int ndim;                 /* dimension of grids */
  RPNT d;                   /* steps */
	IPNT tis;                 /* buffer for coefficient sampling near source */
  
  IPNT is;                  /* source indices */
  RPNT rs;                  /* source in-cell offsets */

  ireal refvel;             /* reference velocity for target prop wavelet */
  ireal refbou;							/* reference buoyancy for target prop wavelet */
  ireal refkappa;           /* near-source bulk modulus */
  ireal refdist;            /* reference distance for target prop wavelet */
  ireal refamp;             /* reference amplitude for target prop wavelet */

  IPNT gs, ge;              /* workspace for buoyancy exchange comps */
	IPNT phiw;
	int iflag; 
	
	char *srcfile;						/* workspace for file name */
	segy trsrc;               /* segy workspace for reading file */
  Value val;                /* workspace for reading segy headers */
  ireal tmpt0;              /* time origin for time series read from file */
  segy trdbg;               /* workspace for building output segy */
	double q;

  int rk;
	
  /* MPI workspace */
#ifdef IWAVE_USE_MPI
  int sz = retrieveSize();
  MPI_Comm cm;
  ireal *procbuf;
	
#endif
  
  /* end declarations */

  stream = retrieveOutstream();
  rk = retrieveRank();

  /* assign default reference values */
  refvel   = CREF_DEF;
  refbou   = REAL_ZERO;
  refkappa = REAL_ZERO;
  refdist  = RREF_DEF;
  refamp   = REAL_ONE;
  tr->fpeak = FPEAK_DEF;  
  tr->phipeak = FPEAK_DEF;        
      

  /* get array of grid steps */
  get_d(d, m->gl);

  /* extract dimension */
  rd_ndim(&m->ld_a, D_MP0, &ndim);

  /* TV */
  IASN(is, IPNT_0);
  RASN(rs, RPNT_0);
  is[0] = tg->is[0]; rs[0]=tg->rs[0]; tis[0] = is[0];
  if (ndim > 1) { is[1]=tg->is[1]; rs[1]=tg->rs[1]; tis[1]=is[1]; }
  if (ndim > 2) { is[2]=tg->is[2]; rs[2]=tg->rs[2]; tis[2]=is[2]; } 
  
  /* get physical coordinates of the source */
  RASN(tr->xs, RPNT_0);
  for (i = 0; i < ndim; i++) tr->xs[i] = ( (float)(tg->is[i]) + tg->rs[i] ) * d[i];


	/* extract near-source bulk modulus from grid - necessary for
     several (but not all) cases, and not enough of an expense to 
     avoid.

     WWS, 29.11.08: to avoid failure if integer part of source point is in 
     physical grid (including physical boundary points) but not in comp.
     grid, check neighboring grid points - in that case, at least one
     of these should be in part of comp grid, where kappa is initialized.
 */

  rd_gse(&(m->ld_a), D_MP0, gs, ge);
  for ( i = 0; i < ndim; ++i ) 
	{
    if ( (is[i] < gs[i]) && (is[i] + 1 > gs[i]-1) ) tis[i]++;
    if ( (is[i] > ge[i]) && (is[i] - 1 < ge[i]+1) ) tis[i]--;
  } 
  iflag=1;
  for (i = 0; i < ndim; ++i) 
    if ( (tis[i] < gs[i]) || (tis[i] > ge[i]) ) iflag=0;

  if (iflag) refkappa = rd_gget(&(m->ld_a), D_MP0, tis);
  fprintf(stream,"NOTE: proc=%d sample kappa at ", retrieveRank());
  for (i=0;i<ndim;++i) fprintf(stream, "index[%d]=%d ", i, tis[i]);
  fprintf(stream, "sample flag=%d\n", iflag);

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
  for ( i = 0; i < ndim; ++i ) 
	{
    tis[i] = is[i];
    if ( (is[i] < gs[i]) && (is[i]+1 > gs[i]-1) ) tis[i]++;
    if ( (is[i] > ge[i]) && (is[i]-1 < ge[i]+1) ) tis[i]--;
  } 
  iflag=1;
  for (i = 0; i < ndim; ++i) 
    if ( (tis[i] < gs[i]) || (tis[i] > ge[i]) ) iflag=0;

  if (iflag) refbou = rd_gget(&(m->ld_a), D_MV0, tis);
  fprintf(stream,"NOTE: proc=%d sample buoyancy at ",retrieveRank());
  for (i=0;i<ndim;++i) fprintf(stream,"index[%d]=%d ",i,tis[i]);
  fprintf(stream, "sample flag=%d\n", iflag);

#ifdef IWAVE_USE_MPI
  cm = retrieveComm();
  if ( rk == 0 ) 
	{
    procbuf = (ireal*)malloc(sz * sizeof(ireal));
    if ( procbuf == NULL ) return E_ALLOC;
  }
  MPI_Gather(&refkappa, 1, IWAVE_MPI_REAL, procbuf, 1, IWAVE_MPI_REAL, 0, cm);
  if(rk == 0 ) 
	{
    refkappa = 0.0;
    for ( i = 0; i < sz; ++i ) refkappa=iwave_max(refkappa,procbuf[i]);
  }
  MPI_Bcast(&refkappa, 1, IWAVE_MPI_REAL, 0, cm);
 
	MPI_Gather(&refbou, 1, IWAVE_MPI_REAL, procbuf, 1, IWAVE_MPI_REAL, 0, cm);
  if(rk ==0 ) 
	{
    refbou = 0.0;
    for ( i = 0; i < sz; ++i ) refbou=iwave_max(refbou,procbuf[i]);
    free(procbuf);
  }
  MPI_Bcast(&refbou, 1, IWAVE_MPI_REAL, 0, cm);
#endif	

  if (refbou > REAL_ZERO) 
	{
    fprintf(stream,"NOTE: in pointsrc, using  buoyancy at source location = %e\n", refbou);
  }
  else 
	{
    fprintf(stream,"ERROR: in pointsrc, ref buoyancy nonpositive, = %e\n",refbou);
    return E_OUTOFBOUNDS;
  }

  if (refkappa > REAL_ZERO) 
	{
    fprintf(stream,"NOTE: in pointsrc, using  bulk mod at source location = %e\n", refkappa);
  }
  else 
	{
    fprintf(stream,"ERROR: in pointsrc, ref bulk mod nonpositive, = %e\n",refkappa);
    return E_OUTOFBOUNDS;
  }

  /* get reference velocity, either from parameters or from buoyancy
     and bulk modulus near source location. Since velocity is not
     stored in RDOM, must be computed from buoyancy and bulk mod at
     slightly different points - this of course does not matter if
     source is located in homogeneous region.
  */
  
  if (ps_ffreal(*par,"refvel", &refvel)) 
	{  

    /* compute velocity from buoyancy and bulk modulus */

    refvel = sqrt(refkappa * refbou);
    
    fprintf(stream,"NOTE: in pointsrc, using velocity computed from \n");
    fprintf(stream,"      bulk mod and buoyancy near source location; \n");
    fprintf(stream,"      computed value = %e\n", refvel);

  }
  else 
	{
    fprintf(stream,"NOTE: in pointsrc, using velocity from param table = %e\n", refvel);
  }
	
  /* Either read reference distance from parameters, or use default. */
  ps_ffreal(*par,"refdist",&refdist);
  if (refdist>0) 
	{
    fprintf(stream,"NOTE: in pointsrc, using reference distance = %e\n", refdist);
  }
  else 
	{
    fprintf(stream,"NOTE: in pointsrc, read nonpos. reference distance = %e\n", refdist);
    fprintf(stream,"      this implies that wavelet will be read from file and used \n");
    fprintf(stream,"      directly on RHS as multiplier of spatial delta, rather than\n");
    fprintf(stream,"      to produce target propagating pulse.\n");
  }
	
	tr->c = refvel;
	tr->bou = refbou;
	
  /* Either read reference amplitude from parameters, or use default. */
  ps_ffreal(*par,"refamp", &refamp);
  fprintf(stream,"NOTE: in pointsrc, using reference amplitude = %e\n", refamp);

  /* Either read  peak for the cutoff function, or use default. */
	if (ps_ffreal(*par, "cutoff_peak", &(tr->phipeak)))
	{
    fprintf(stream,"NOTE: pointsrc_init - using default ");
    fprintf(stream,"parameter for the cutoff function (cutoff_width) = %e\n", tr->phipeak);
  }
	
  /* Either read the radius of the flat area of the cutoff function, or use default */
  tr->rad = 1.4 / tr->phipeak; 
  if (ps_ffreal(*par, "cutoff_rad", &(tr->rad)))
	{
    fprintf(stream,"NOTE: pointsrc_init - using default ");
    fprintf(stream, "rad of the flat area in the cutoff function = %e\n", tr->rad);
  }

  /* compute start and end indecies for the cutoff funciton */
  IASN(tr->ixs, IPNT_0);
  IASN(tr->ixe, IPNT_0);
	for (i = 0; i < ndim; ++i)
	{
		phiw[i] = 1 + (int)(floor( (tr->rad +  1.4 / tr->phipeak) / d[i] + 0.1));  
		tr->ixs[i] = is[i] - phiw[i];
		tr->ixe[i] = is[i] + phiw[i];    
	}
  
  /* Overall scale factor for source insertion to produce target pulse, per 
     paper by WWS&TV: contains factors of
     - dt - RHS of the difference scheme is multiplied by the time step
     - r (reference distance for normalization, per eqn 13);
     - reference amplitude;
     Note: for Option Ib (direct insertion of source wavelet, signalled by 
     refdist<=0) this factor accounts only for time step as this case defines no target scale.
  */

  if (refdist>0) 
  {
    tr->scramp = refamp * refdist * m->tsind.dt;
  }
  else 
  {
    tr->scramp = m->tsind.dt;
  }
	/* Option I: read source from file */
	
  if (!ps_ffcstring(*par, "source", &srcfile))  
	{
    if (!(tr->fpsrc = fopen(srcfile, "r"))) {
      fprintf(stream, "Error: pointsrc_init - failed to open source file\n");
      return E_FILE;
    }
    if (!fgettr(tr->fpsrc, &trsrc)) {
      fprintf(stream,"Error: pointsrc_init - failed to read source file\n");
      return E_FILE;
    }
		
		
    /* at this point, the data member of trsrc contains the source wavelet
		at an external sample rate - read headers relevant for single dilat
		point source.
    */
    mygethdval(&trsrc, "ns", &val);
    tr->n = vtoi(hdtype("ns"), val);
    mygethdval(&trsrc, "dt", &val);
    tr->dt = 0.001 * vtof(hdtype("dt"), val);
    mygethdval(&trsrc, "delrt", &val);
    tmpt0 = vtof(hdtype("delrt"), val);		
		
    tr->istart = (int)(floor(tmpt0 / m->tsind.dt));
    tr->t0 = (m->tsind).dt * tr->istart;

    tr->w = (ireal *)malloc(sizeof(ireal)*(tr->n));
		tr->w1 = (ireal *)malloc(sizeof(ireal)*(tr->n));
    
    /* allocate buffer for integrated wavelet 
			at input sample rate */
    
    for (i = 0; i < tr->n; i++)  tr->w[i] = trsrc.data[i];
    
    /* Option Ia: normalize to produce propagating wavelet at distance
			refdist, speed refvel, if refdist > 0. */
    
    if (refdist > 0) 
		{
    
			tr->tmpflag = 3;
      fprintf(stream,"NOTE: in pointsrc_init, compute wavelet for target pulse (Option Ia)\n");

      /* integrate once - trapezoidal rule in-place 
          w[j] <- sum_{i=1}^{i=j} 0.5*dt*(w[i-1]+w[i])
      */

      q = 0.0;
      tr->w1[0] = 0.0;
      for (i = 1; i < tr->n; i++) 
			{
				q += 0.5 * tr->dt * (tr->w[i] + tr->w[i - 1]);
				tr->w1[i] = q;
			}
  
    }
    /* Option Ib: use wavelet as is, if refdist <= 0 */
    else
    {
      tr->tmpflag = 2;
      fprintf(stream,"NOTE: in pointsrc_init, using file wavelet directly on RHS (Option Ib)\n");
           
      /* differentiate once - trapezoidal rule in-place 
          w[j] <- sum_{i=1}^{i=j} 0.5 * (w[i+1]+w[i - 1])/dt
      */
       
      tr->w1[0] = 0.0;

      for (i = 1; i < tr->n - 1; i++) 
				tr->w1[i] = (tr->w[i + 1] - tr->w[i - 1]) / (2.0 *  tr->dt);
      
      tr->w1[tr->n- 1] = 0.0;
      
    }
    free(srcfile);
	}
	else
	{
  
   /* check that reference distance is positive - only legal option */
    if ( !(refdist>0) ) 
    {
      fprintf(stream,"Error: pointsrc_init, gaussian case:\n");
      fprintf(stream,"negative reference distance = %e no legal in this case\n", refdist);
      return E_OTHER;
    }
     
		tr->tmpflag = 1;
		/* read peak frequency from parameters, or use default (only used in
		Gaussian option II) */
		if (ps_ffreal(*par,"fpeak", &(tr->fpeak))) 
		{
			fprintf(stream,"NOTE: pointsrc_init - using default ");
			fprintf(stream,"peak frequency (fpeak) = %e\n", tr->fpeak);
		}
		
		/* source phase - default is zero-phase */
    tr->w = getrick(&iw, (m->tsind).dt, tr->fpeak);
    tr->w1 = getdgauss(&iw, (m->tsind).dt, tr->fpeak);

		tr->istart = -iw;
		tr->n = 2 * iw + 1;
    tr->dt = m->tsind.dt;


		if (!ps_ffcstring(*par, "waveletphase", &wp)) 
		{
			if (!strcmp(wp,"zerophase")) tr->istart = -iw;
			else if (!strcmp(wp,"causal")) tr->istart = 0;
			else if (!strcmp(wp,"anticausal")) tr->istart = -2 * iw;
			else 
			{
				fprintf(stream,"Error: pointsrc_init, gaussian case:\n");
				fprintf(stream,"only legit values of waveletphase are \n");
				fprintf(stream,"zerophase, causal, anticausal\n");
				return E_OTHER;
			}
			free(wp);
		}			
	}
  
  tr->sn = (int)((tr->rad + 1.4 / tr->phipeak) / tr->c / m->tsind.dt) + 1 + tr->n;
  
  /* optionally write out wavelet appearing on RHS */
  tr->idbg = 0;
  ps_ffint(*par, "dump_wavelet", &(tr->idbg));
  if ( (tr->idbg) && (rk==0) )
  {
    if (tr->idbg == 1)  memcpy(trdbg.data, tr->w, tr->n * sizeof(ireal));
		if (tr->idbg == 2)  memcpy(trdbg.data, tr->w1, tr->n * sizeof(ireal));
    val.u=1000.0*(tr->dt);
    myputhdval(&trdbg,"dt",&val);
    val.h=tr->n;
    myputhdval(&trdbg,"ns",&val);
    val.h=tr->dt * tr->istart;
    myputhdval(&trdbg,"delrt",&val);
    
    if (!(tr->fpdbg=fopen("wavelet.debug","w"))) {
      fprintf(stream,"Error: init_point: failed to open test wavelet file\n");
      return E_FILE;
    }
    else {
      fprintf(stream,"write wavelet trace\n");
    }
    fputtr(tr->fpdbg,&trdbg);
    fflush(tr->fpdbg);
  }
	
  return 0;
  
}

int spointsrc_destroy(SPOINTSRC * tr) 
{
   return 0;
}

int spointsrc_run(SPOINTSRC * tr, IMODEL * m) 
{

  RPNT d;  
  int ndim, i;
  IPNT ix, gs, ge;
  ireal bm, bou, p0, p1=0.0, p2=0.0, ps, t; 
/*  FILE * stream;
    stream = retrieveOutstream(); */


  
  if ( (m->tsind.it >= tr->istart) &&  (m->tsind.it <= tr->istart + tr->sn) )
  {
    
		t = (ireal)(m->tsind.it) * m->tsind.dt;
    /* extract dimension */
    rd_ndim(&m->ld_a, D_MP0, &ndim);
    get_d(d,m->gl);
    
    if ( (m->tsind).iv == 0 )
    {
      
      IASN(gs, IPNT_0);
      IASN(ge, IPNT_0);
			rd_gse(&(m->ld_a), D_MP0, gs, ge);
						
      #if RARR_MAX_NDIM > 2
      for (ix[2] = iwave_max(tr->ixs[2], gs[2]); ix[2] < iwave_min(tr->ixe[2] + 1, ge[2] + 1); ix[2]++) 
      {
      #endif
          #if RARR_MAX_NDIM > 1
          for (ix[1] = iwave_max(tr->ixs[1], gs[1]); ix[1] < iwave_min(tr->ixe[1] + 1, ge[1] + 1); ix[1]++) 
          {
          #endif
              for (ix[0] = iwave_max(tr->ixs[0], gs[0]); ix[0] < iwave_min(tr->ixe[0] + 1, ge[0] + 1); ix[0]++) 
              {
                  p0 = rd_gget(&(m->ld_a), D_P[0], ix);
                  if (ndim > 1) p1 = rd_gget(&(m->ld_a), D_P[1], ix);
                  if (ndim > 2) p2 = rd_gget(&(m->ld_a), D_P[2], ix);


                  bm  = rd_gget(&(m->ld_a), D_MP0, ix);
                  
                  ps = bm * tr->scramp * vfun(ix, d, ndim, tr, t + 0.5 * m->tsind.dt, m->tsind.dt);
                  
                  rd_gset(&(m->ld_a), D_P[0], ix, p0 - ps);
                  if (ndim  > 1) rd_gset(&(m->ld_a), D_P[1], ix, p1 - ps);
                  if (ndim  > 2) rd_gset(&(m->ld_a), D_P[2], ix, p2 - ps);

              }
          #if RARR_MAX_NDIM > 1
          }
          #endif
      #if RARR_MAX_NDIM > 2
      }
      #endif    
    }	

		for (i=0; i < ndim; i++) 
    {
      
			if ((m->tsind).iv == 3 * ndim + 1 + 4 * i + 1) 

      {
        IASN(gs, IPNT_0);
        IASN(ge, IPNT_0);
        rd_gse(&(m->ld_a), D_MV[i], gs, ge);
       
        #if RARR_MAX_NDIM > 2
        for (ix[2] = iwave_max(tr->ixs[2], gs[2]); ix[2] < iwave_min(tr->ixe[2] + ((i==2) ? 0 : 1), ge[2] + 1); ix[2]++) 
        {
        #endif
            #if RARR_MAX_NDIM > 1
            for (ix[1] = iwave_max(tr->ixs[1], gs[1]); ix[1] < iwave_min(tr->ixe[1] + ((i==1) ? 0 : 1), ge[1] + 1); ix[1]++) 
            {
            #endif
                for (ix[0] = iwave_max(tr->ixs[0], gs[0]); ix[0] < iwave_min(tr->ixe[0] + ((i==0) ? 0 : 1), ge[0] + 1); ix[0]++) 
                {
                      bou = rd_gget(&(m->ld_a), D_MV[i], ix);
                      p0 = rd_gget(&(m->ld_a), D_V[i], ix);
                      ps = bou * tr->scramp * pfun(ix, d, ndim, tr, t + m->tsind.dt, i, m->tsind.dt);
                      rd_gset(&(m->ld_a), D_V[i], ix, p0-ps);
                  }
            #if RARR_MAX_NDIM > 1
            }
            #endif
          #if RARR_MAX_NDIM > 2
        }
        #endif
      }
    }
  }
  
	return 0;
}

void spointsrc_fprint(SPOINTSRC const * tr, FILE * fp) 
{
  int i;
  fprintf(fp,"/*---------------------------------------------------------*/\n");
  fprintf(fp,"SMOOTH POINT SOURCE\n");
  fprintf(fp,"source sound vel = %f\n",tr->c);
	fprintf(fp,"source bou       = %f\n", tr->bou);
  fprintf(fp,"istart           = %d\n", tr->istart);
	fprintf(fp,"n                = %d\n", tr->n);
  fprintf(fp,"dt                = %f\n", tr->dt);
  fprintf(fp,"sn               = %d\n", tr->sn);
  fprintf(fp,"fpeak            = %f\n", tr->fpeak);
	fprintf(fp,"phipeak          = %f\n", tr->phipeak);
  for (i=0;i<RARR_MAX_NDIM;i++)
    fprintf(fp,"xs[%d]=%f\n", i, tr->xs[i]);
	for (i=0;i<RARR_MAX_NDIM;i++)
    fprintf(fp,"ixs[%d]=%d ixe[%d]=%d\n", i, tr->ixs[i], i, tr->ixe[i]);
  fprintf(fp,"/*---------------------------------------------------------*/\n");
}
