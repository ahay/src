#include "acd.hh"

using RVL::parse;

//#define IWAVE_VERBOSE

/*--- time step functions ---------------------------------------------------*/

extern void acd_2d_2(float **, 
		     float **, 
		     float **, 
		     int *, 
		     int *, 
		     float, 
		     float *);

extern void acd_2d_4(float **, 
		     float **, 
		     float **, 
		     int *, 
		     int *, 
		     float, 
		     float *,
		     float *,
		     int *,
		     int *);

extern void acd_2d_8(float **, 
		     float **, 
		     float **, 
		     int *, 
		     int *, 
		     float, 
		     float *,
		     float *,
		     float *,
		     float *,
		     int *,
		     int *);

extern void acd_3d_2(float ***, 
		     float ***, 
		     float ***, 
		     int *, 
		     int *, 
		     float, 
		     float *);

extern void acd_3d_4(float ***, 
		     float ***, 
		     float ***, 
		     int *, 
		     int *, 
		     float, 
		     float *,
		     float *,
		     int *,
		     int *);

extern void acd_3d_8(float ***, 
		     float ***, 
		     float ***, 
		     int *, 
		     int *, 
		     float, 
		     float *,
		     float *,
		     float *,
		     float *,
		     int *,
		     int *);

extern int acd_step(RDOM*, int, void *);

/*----------------------------------------------------------------------------*/
/* no-ops for this implementation                                             */
/*----------------------------------------------------------------------------*/
int acd_build_sten_dep(FILE * stream, 
		       int ndim, 
		       int stendep[RDOM_MAX_NARR][RDOM_MAX_NARR]) {
    return 0;
}

/*----------------------------------------------------------------------------*/
/* working functions                                                          */
/*----------------------------------------------------------------------------*/

/*int acd_modelinit(PARARRAY *pars, 
		  FILE *stream, 
		  grid const & g,
		  ireal dt,
		  std::vector<std::string> & active,
		  void ** fdpars) {
*/
int acd_modelinit(PARARRAY pars, 
		  FILE *stream,
		  IMODEL & model) {

    int err=0;           /* return value */
    ACD_TS_PARS *acdpars;   /* model pars */

    IPNT cdims;          /* workspace for cartesian grid dim info */
    IPNT crank;          /* workspace for cartesian grid rank */
#ifdef IWAVE_USE_MPI
    IPNT cpers;          /* workspace for periodic wrap info  - currently not used */
#endif

    int i;       /* counter */
    RPNT dxs;    /* grid steps */
    ireal lam;   /* slownesses dt/dx */
    int idim;    /* counter */

    /* allocate sgn model ----------------------------------------------------*/   
    acdpars = (ACD_TS_PARS*)usermalloc_(sizeof(ACD_TS_PARS));
    if ( acdpars == NULL ) { 
      err=E_ALLOC;
      fprintf(stream,"ERROR: acd_modelinit\n");
      fprintf(stream,"failed to allocate SGN_TS_PARS object\n");
      return err;
    }
    
    /* assign heritage pars
    // 11.04.12: remove these
    //  sgnm->psingle = 0; // TODO: single pressures in Step 3 
    //  IASN(sgnm->eflags,IPNT_1); // this is related - with code from old sg added, could try using this device */
    
    /* decode dimensions, parallel rank - read grid dimn on rank 0, broadcast */
    
    IASN(cdims, IPNT_1); /* default grid size */ 
    IASN(crank, IPNT_0); /* default cartisian ranks */ 
    
#ifdef IWAVE_USE_MPI
    MPI_Comm cm=retrieveComm();
    
    if ( MPI_Cart_get(cm, RARR_MAX_NDIM, cdims, cpers, crank) != MPI_SUCCESS )  {
      fprintf(stream, "ERROR. Internal: cannot get Cartesian coordinates.\n");
      return E_INTERNAL;
    }
    
    MPI_Bcast((void*)(&(model.g.dim)),1,MPI_INT,0,cm);
#endif
    
    /* set boundary flags */
    IASN(acdpars->lbc, IPNT_0); /* default left bc flag */ 
    IASN(acdpars->rbc, IPNT_0); /* default right bc flag */ 

    for (i=0;i<model.g.dim;i++) {
      if (crank[i]==0) acdpars->lbc[i]=1;
      if (crank[i]==cdims[i]-1) acdpars->rbc[i]=1;
    }

    /* decode order - with version 2.0, deprecated syntax "scheme_phys" etc. is dropped */
    acdpars->k=1;
    parse(pars,"order",acdpars->k);
#ifdef IWAVE_VERBOSE
    fprintf(stream,"NOTE: initializing ACD with half-order = %d\n",acdpars->k);
#endif

    /* initialize scaled Courant arrays */
    acdpars->c0=REAL_ONE;
    RASN(acdpars->c1,RPNT_0);
    RASN(acdpars->c2,RPNT_0);
    RASN(acdpars->c3,RPNT_0);
    RASN(acdpars->c4,RPNT_0);

    /* initialize bound check params */
    if (!parse(pars,"cmax",acdpars->cmax)) {
      RVLException e;
      e<<"Error: acd_modelinit\n";
      e<<"  failed to find cmax in param table\n";
      throw e;
    }
    if (!parse(pars,"cmin",acdpars->cmin)) {
      RVLException e;
      e<<"Error: acd_modelinit\n";
      e<<"  failed to find cmin in param table\n";
      throw e;
    }

    /* extract grid steps from grid */
    get_d(dxs, model.g);

    /* set model dimn par */
    acdpars->ndim = model.g.dim;

    acdpars->c0 = 0.0;
    for (idim = 0;idim < acdpars->ndim;idim ++) {

      if (dxs[idim] <= 0.0) {
	fprintf(stream, "Error: bad input: wrong grid space step, dim=%d, step=%g\n",
		idim, dxs[idim]);
	return E_BADINPUT;
      }
      lam = model.tsind.dt / dxs[idim];
      
      /* assign scaled Courant numbers for orders 2, 4, and 8 - these are the only */
      /* choices implemented */
      if (acdpars->k==1) {
	acdpars->c1[idim]   = lam*lam;
	acdpars->c0        += lam*lam*(-2.0);
      }
      else if (acdpars->k==2) {
	acdpars->c1[idim]   = lam*lam*(4.0/3.0);
	acdpars->c2[idim]   = lam*lam*(-1.0/12.0);
	acdpars->c0        += lam*lam*(-5.0/2.0);
      }
      else if (acdpars->k==4) {
	acdpars->c1[idim]   = lam*lam*(8.0/5.0);
	acdpars->c2[idim]   = lam*lam*(-1.0/5.0);
	acdpars->c3[idim]   = lam*lam*(8.0/315.0);
	acdpars->c4[idim]   = lam*lam*(-1.0/560.0);
	acdpars->c0        += lam*lam*(-205.0/72.0);
      }
      else {
	fprintf(stream,"ERROR: acd_readschemeinfo\n");
	fprintf(stream,"assigned scheme half-order = %d not defined\n",acdpars->k);
	fprintf(stream,"currently defined schemes: half-orders 1, 2, and 4\n");
	fflush(stream);
	return E_BADINPUT;
      }
#ifdef IWAVE_VERBOSE
      fprintf(stderr, "k=%d lam[%d] = %g\n", acdpars->k, idim, lam);
      if (acdpars->k==1)
	fprintf(stderr,"c1[%d]=%g\n",
		idim,acdpars->c1[idim]);
      if (acdpars->k==2) 
	fprintf(stderr,"c1[%d]=%g c2[%d]=%g\n",
		idim,acdpars->c1[idim],idim,acdpars->c2[idim]);
      if (acdpars->k==4) 
	fprintf(stderr,"c1[%d]=%g c2[%d]=%g c3[%d]=%g c4[%d]=%g\n",
		idim,acdpars->c1[idim],idim,acdpars->c2[idim],
		idim,acdpars->c3[idim],idim,acdpars->c4[idim]);
#endif
    }
#ifdef IWAVE_VERBOSE
    fprintf(stderr,"c0=%g\n",acdpars->c0);
#endif
    /* reserve a copy of dt for use in source scaling */
    acdpars->dt = model.tsind.dt;
    
    /* identify active fields */
    model.active.resize(3);
    model.active[0]="csq";
    model.active[1]="uc";
    model.active[2]="up";
    
    /* assign param object pointer */
    model.specs = (void*)acdpars;
    return 0;
}

/*----------------------------------------------------------------------------*/
void acd_modeldest(void ** fdpars) {

    /* destroy acdpars - all data members allocated on stack */
    userfree_(*fdpars);
}

/*----------------------------------------------------------------------------*/
/* this variant so simple that might as well write it from
// scratch, so no need for sten_dep_mat - note that all arrays in this
// app are primal, so no
// need for gtype. Clearly this interface should be refactored so that
// these things can be hidden.
// ndim (from fdpars) gtype and sten_dep_mat should be internal details */

int acd_create_sten(void * fdm,
		    FILE * stream,
		    //      		    IWaveInfo const & ic,
		    int ndim,
		    IPNT gtype[RDOM_MAX_NARR], 
		    STENCIL * sten) {
    ACD_TS_PARS * acdpars = (ACD_TS_PARS *)(fdm);
    STENCIL_MASK mask;/* workspace */
    int nmask;        /* number of masks - dependent pairs of dynamic arrays */
    int ipair[2][2];  /* index pairs for use in defining mask */
    int imask;        /* mask counter */
    int idim;         /* dim counter */
    int iv;           /* mask entry counter */
    int len;          /* length of mask - number of entries */
    int j;            /* counter */
    int k;            /* scheme order */
    IPNT ind;         /* workspace for mask entry */
    int err = 0;

    /* set order variable */
    k = acdpars->k;

    /* initialize index pairs */
    /* first - uc->up */
    ipair[0][0] = D_UC;
    ipair[0][1] = D_UP;
    /* second - up->uc */
    ipair[1][0] = D_UP;
    ipair[1][1] = D_UC;

    /* initialize STENCIL to null stencil */
    sten_setnull(sten);
  
    /* sanity check */
    if (k < 1) {
	fprintf(stream,"ERROR: acd_create_sten - illegal value of k = %d\n",k);
	return E_BADINPUT;
    }

    /* declare number of masks: one for shape of uc->up stencil, one for
    // shape of up->uc stencil */
    nmask=2;

    /* nontrivial STENCIL initialization */
    if ((err = sten_create(sten,nmask))) {
	fprintf(stream,"ERROR: acd_create_sten - failed to create stencil\n");
	return err;
    }

    /* length of stencil is 2k+1 in each direction, but origin is common to all
    // directions, so */
    len = 2*k*ndim+1;
    for (imask=0;imask<nmask;imask++) {
	if ((err = mask_create(&mask, ipair[imask][0], ipair[imask][1], len))) {
	    fprintf(stream,"ERROR: acd_create_sten from mask_create\n");
	    sten_destroy(sten);
	    return err;
	}
	/* "cross stencil" - same in every dimension
	// 2d 4th order - k=0
	// idim = 0
	//   j = 0
	//     iv=0: ind[0]=-1
	//     iv=2: ind[0]= 1
	//     iv=1: ind[0]=-2
	//     iv=3  ind[0]= 2
	// etc
	// eventually iv = ndim*4, ind=IPNT_0 */

	for (idim=0;idim<ndim;idim++) {
	    IASN(ind,IPNT_0);
	    for (j=0;j<k;j++) {
		/* left half of mask on axis idim */
		ind[idim]=-j-1;
		iv = idim*2*k+j;
		if ((err = mask_set(&mask,iv,ind))) {
		    fprintf(stream,"ERROR: acd_create_sten from mask_set\n");	
		    sten_destroy(sten);
		    return err;
		}	
		/* right half of mask on axis idim */
		ind[idim]=j+1;
		iv = idim*2*k+j+k;
		if ((err = mask_set(&mask,iv,ind))) {
		    fprintf(stream,"ERROR: acd_create_sten from mask_set\n");	
		    sten_destroy(sten);
		    return err;
		}	
	    }
    
	    IASN(ind,IPNT_0);
	    iv=2*k*ndim;
	    if ((err = mask_set(&mask,iv,ind))) {
		fprintf(stream,"ERROR: acd_create_sten from mask_set\n");	
		sten_destroy(sten);
		return err;
	    }	
	}	
	/*
	  fprintf(stderr,"setting mask %d\n",imask);
	  fprintf(stderr,"ip = %d\n",mask.ip);
	  fprintf(stderr,"ir = %d\n",mask.ir);
	  fprintf(stderr,"n  = %d\n",mask.n);
	  for (j=0;j<mask.n;j++) 
	  fprintf(stderr,"s[%d] = [%d,%d,%d]\n",j,(mask.s[j])[0],(mask.s[j])[1],(mask.s[j])[2]);
	*/
	if ((err=sten_set(sten,imask, &mask))) {
	    fprintf(stream,"ERROR: acd_create_sten from sten_set\n");	
	    sten_destroy(sten);
	    return err;
	}
    }
    /*
      if (err=sten_out(sten,stderr,acd_ind2str)) {
      fprintf(stream,"ERROR: acd_create_sten from sten_out\n");	
      sten_destroy(sten);
      return err;
      }
    */
    return 0;
}

/*----------------------------------------------------------------------------*/
/* implements new time grid logic: choose stable time step (max_step set), 
   optionally with reduced cfl - to be called in iwave_construct BEFORE
   any data are read, so must depend on max velo, cannot check here.
*/

//int acd_readtimegrid(PARARRAY *pars, FILE * stream, IMODEL * model) {
//int acd_readtimegrid(PARARRAY *pars, FILE * stream, grid const & g, ireal & dt) {
int acd_timegrid(PARARRAY *pars, 
		 FILE * stream, 
		 grid const & g, 
		 ireal & dt,
		 ireal & rhs) {

    ireal cmax;                   /* max velo, computed or read from params */
    ireal cfl = CFL_DEF;          /* default cfl fraction */
    ireal a;                      /* accumulator for max space step */
    int i;                        /* counter */

    /* branch on presence of parameter dt - if dt set in param table,
       use it */
    if (parse(*pars,"dt", dt) ){
#ifdef IWAVE_VERBOSE
	fprintf(stream,"NOTE: sg_readtimegrid - dt=%12.4e read from param table\n", dt);	
	fprintf(stream,"NOTE: NOT CHECKED FOR STABILITY!\n");
#endif
	rhs=dt*dt;
	return 0;
    }
    	
    /* for either method of computing dt, need cfl = proportion
       of max time step to use */
    cfl=CFL_DEF;
    if (!parse(*pars,"cfl",cfl)) {
#ifdef IWAVE_VERBOSE
	fprintf(stream,"NOTE: acd_readtimegrid\n");
	fprintf(stream,"  using default cfl fraction %g\n",cfl);;
#endif
    }

    if (!parse(*pars,"cmax",cmax)) { 
      fprintf(stream,"ERROR: acd_readtimegrid - failed to read cmax from parameters \n");
      return E_BADINPUT;
    }

    a = g.axes[0].d;
  
    for ( i = 1; i < g.dim; ++i ) a = iwave_min(a, g.axes[i].d);
#ifdef IWAVE_VERBOSE
    fprintf(stream,"NOTE: acd_readtimegrid: ndim=%d min dx=%e cfl=%e\n",g.dim,a,cfl);
    fflush(stream);
#endif

    if ( (a < REAL_EPS) || (cfl < REAL_EPS) ) {
	fprintf(stream,"ERROR: acd_readtimegrid - either min dx=%e or cfl=%e "
		" too small\n", a, cfl);
	return E_BADINPUT;
    }

    dt = a*cfl/(cmax*sqrt((float)(g.dim)));
    //cerr << "dt="<< dt << endl;
#ifdef IWAVE_VERBOSE
    fprintf(stream,"NOTE: acd_readtimegrid: on return, dt=%e\n",dt);
    fflush(stream);
#endif
    rhs=dt*dt;
    return 0;
}

int acd_step(RDOM* dom, int iv, void * tspars) {

    /* pointers for 2D case */
    ireal ** restrict uc2;
    ireal ** restrict up2;
    ireal ** restrict csq2;
    /* pointers for 3D case */
    ireal *** restrict uc3;
    ireal *** restrict up3;
    ireal *** restrict csq3;
    int ndim;                       /* problem dmn */
    IPNT s, s0;                     /* loop starts  */
    IPNT e, e0;                     /* loop ends */

    ireal tmp;
    IPNT i;

    /* acd struct */
    ACD_TS_PARS * acdpars = (ACD_TS_PARS *)tspars;

    /* extract dimn info */
    ra_ndim(&(dom->_s[D_UC]),&ndim);
    ra_gse(&(dom->_s[D_UC]),s,e);
    ra_a_gse(&(dom->_s[D_UC]),s0,e0);
    if (ndim == 2) {


	/* 2D computational arrays */
	uc2   = (dom->_s)[D_UC ]._s2;
	up2   = (dom->_s)[D_UP ]._s2;
	csq2  = (dom->_s)[D_CSQ]._s2;

	/* 2nd order case */
	if (acdpars->k == 1) {
	  //	    acd_2d_2(uc2, up2, csq2, 
	  acd_2d_2((dom->_s)[D_UC ]._s2,
		   (dom->_s)[D_UP ]._s2,
		   (dom->_s)[D_CSQ]._s2,
		     s, e,
		     acdpars->c0, 
		     acdpars->c1);
	}
	/* 4th order case */
	else if (acdpars->k == 2) {
	    acd_2d_4(uc2, up2, csq2,
		     s, e, 
		     acdpars->c0, 
		     acdpars->c1, acdpars->c2,
		     acdpars->lbc, acdpars->rbc);
	}
	/* 8th order case */
	else if (acdpars->k == 4) {
	    acd_2d_8(uc2, up2, csq2,
		     s, e, 
		     acdpars->c0, 
		     acdpars->c1, acdpars->c2,
		     acdpars->c3, acdpars->c4,
		     acdpars->lbc, acdpars->rbc);
	}
	else {
	    fprintf(stderr,"ERROR: acd_step\n");
	    fprintf(stderr,"called with half-order != 1, 2, 4\n");
	    return E_BADINPUT;
	}

	int s00=s0[0]; int e00=e0[0];
	for (i[1]=s0[1];i[1]<=e0[1];i[1]++) {
#if defined(__INTEL_COMPILER)
#pragma ivdep
#elif defined(__GNUC__) && !defined(__clang__)
#pragma GCC ivdep
#endif
	    for (i[0]=s00;i[0]<=e00;i[0]++) {
		tmp=((dom->_s)[D_UC]._s2)[i[1]][i[0]];
		((dom->_s)[D_UC]._s2)[i[1]][i[0]]=((dom->_s)[D_UP]._s2)[i[1]][i[0]];
		((dom->_s)[D_UP]._s2)[i[1]][i[0]]=tmp;
	    }
	}

    }
    else if (ndim == 3) {
    
	uc3   = (dom->_s)[D_UC ]._s3;
	up3   = (dom->_s)[D_UP ]._s3;
	csq3  = (dom->_s)[D_CSQ]._s3;

	/* 2nd order case */
	if (acdpars->k == 1) {
	    acd_3d_2(uc3, up3, csq3, 
		     s, e, 
		     acdpars->c0, 
		     acdpars->c1);
	}
	/* 4th order case */
	else if (acdpars->k == 2) {
	    acd_3d_4(uc3, up3, csq3,
		     s, e, 
		     acdpars->c0, 
		     acdpars->c1, acdpars->c2,
		     acdpars->lbc, acdpars->rbc);
	}
	/* 8th order case */
	else if (acdpars->k == 4) {
	    acd_3d_8(uc3, up3, csq3,
		     s, e, 
		     acdpars->c0, 
		     acdpars->c1, acdpars->c2,
		     acdpars->c3, acdpars->c4,
		     acdpars->lbc, acdpars->rbc);
	}
	else {
	    fprintf(stderr,"ERROR: acd_step\n");
	    fprintf(stderr,"called with half-order != 1, 2, 4\n");
	    return E_BADINPUT;
	}

	for (i[2]=s0[2];i[2]<=e0[2];i[2]++) {
	    for (i[1]=s0[1];i[1]<=e0[1];i[1]++) {
		for (i[0]=s0[0];i[0]<=e0[0];i[0]++) {
		    tmp=((dom->_s)[D_UC]._s3)[i[2]][i[1]][i[0]];
		    ((dom->_s)[D_UC]._s3)[i[2]][i[1]][i[0]]=((dom->_s)[D_UP]._s3)[i[2]][i[1]][i[0]];
		    ((dom->_s)[D_UP]._s3)[i[2]][i[1]][i[0]]=tmp;
		}
	    }
	}
    }
    else {
	fprintf(stderr,"ERROR: acd_step\n");
	fprintf(stderr,"called with space dim != 2 or 3\n");
	return E_BADINPUT;
    }
  
    return 0;
}

/*----------------------------------------------------------------------------*/
void acd_check(RDOM * dom, void * specs, FILE * stream) {
  // arrays to test: only csq = rarr[0]
  size_t n = 0; // length of data array
  ireal vmax = -FLT_MAX;
  ireal vmin = FLT_MAX;

  /* acd struct */
  ACD_TS_PARS * acdpars = (ACD_TS_PARS *)specs;
  
  /* extract dimn info */
  ra_a_datasize(&(dom->_s[0]),&n);

  /* max & min */
  for (int i=0;i<(int)n;i++) {
    vmax=iwave_max((dom->_s[0]._s0)[i],vmax);
    vmin=iwave_min((dom->_s[0]._s0)[i],vmin);
    if (vmax<0.0f || sqrt(vmax)>acdpars->cmax || 
	vmin<0.0f || sqrt(vmin)<acdpars->cmin) {
      RVLException e;
      e<<"i = " << i << "\n";
      e<<"sqrt(vmax) = " << sqrt(vmax) << "\n";
      e<<"sqrt(vmin) = " << sqrt(vmin) << "\n";
      e<<"Error: input csq at index "<<i<<" = "<<(dom->_s[0]._s0)[i]<<"\n";
      e<<"  out of bounds ["<<acdpars->cmin<<", "<<acdpars->cmax<<"]\n";
      throw e;
    }
  }
}
