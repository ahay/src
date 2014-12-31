#include "asg.hh"

using RVL::valparse;
using RVL::RVLException;

extern "C" void asg_pstep2d(float ** restrict bulk, 
			    float ** restrict p0, float ** restrict p1,
			    float ** restrict v0, float ** restrict v1,
			    float ** restrict ep, float ** restrict epp,
			    float * restrict sdiv,
			    int * gsc, int * gec, 
			    int maxoff, float * restrict c[RARR_MAX_NDIM]);

extern "C" void asg_pstep3d(float *** restrict bulk, 
			    float *** restrict p0, float *** restrict p1, float *** restrict p2,
			    float *** restrict v0, float *** restrict v1, float *** restrict v2,
			    float ** restrict ep, float ** restrict epp,
			    float * restrict sdiv,
			    int * gsc, int * gec, 
			    int maxoff, float ** restrict c);

extern "C" void asg_vstep2d(float ** restrict buoy,
			    float ** restrict p0,float ** restrict p1,
			    float ** restrict v0,float ** restrict v1,
			    float ** restrict ev, float ** restrict evp,
			    float ** restrict gradp,
			    int * gsc_v0, int * gec_v0,
			    int * gsc_v1, int * gec_v1,
			    int maxoff,float ** restrict c);

extern "C" void asg_vstep3d(float *** restrict buoy,
			    float *** restrict p0,float *** restrict p1, float *** restrict p2,
			    float *** restrict v0,float *** restrict v1, float *** restrict v2,
			    float ** restrict ev, float ** restrict evp,
			    float ** restrict gradp,
			    int * gsc_v0, int * gec_v0,
			    int * gsc_v1, int * gec_v1,
			    int * gsc_v2, int * gec_v2,
			    int maxoff,float ** restrict c);

/*
int asg_modelinit(PARARRAY *pars,
		  FILE *stream,
		  grid const & g,
		  ireal dt,
		  std::vector<std::string> & active,
		  void ** specs) {
*/
int asg_modelinit(PARARRAY pars,
		  FILE *stream,
		  IMODEL & model) {
  try {

    int err=0;           /* return value */
    ASG_TS_PARS *asgpars;   /* model pars */

    IPNT cdims;          /* workspace for cartesian grid dim info */
    IPNT crank;          /* workspace for cartesian grid rank */
#ifdef IWAVE_USE_MPI
    IPNT cpers;          /* workspace for periodic wrap info  - currently not used */
#endif

    RPNT dxs;    /* grid steps */
    ireal lam;   /* slownesses dt/dx */

    /* allocate sgn model ----------------------------------------------------*/   
    asgpars = (ASG_TS_PARS*)usermalloc_(sizeof(ASG_TS_PARS));
    if ( asgpars == NULL ) { 
      err=E_ALLOC;
      fprintf(stream,"ERROR: asg_modelinit\n");
      fprintf(stream,"failed to allocate SGN_TS_PARS object\n");
      return err;
    }

    /************ The next block of code determines which boundary faces of the
     * domain assigned to the current process are external, and which internal.
     * this code should be moved to the core, along with the boundary flags 
     * lbc and rbc.
     ************                         BEGIN                     *************/	    
    /* decode dimensions, parallel rank - read grid dimn on rank 0, broadcast */

    IASN(cdims, IPNT_1); /* default grid size */ 
    IASN(crank, IPNT_0); /* default cartisian ranks */ 

#ifdef IWAVE_USE_MPI
    MPI_Comm cm=retrieveComm();

    if ( MPI_Cart_get(cm, RARR_MAX_NDIM, cdims, cpers, crank) != MPI_SUCCESS )  {
      RVLException e;
      e<<"ERROR: asg_modelinit\n";
      e<<"  cannot get Cartesian coordinates.\n";
      throw e;
    }

    MPI_Bcast((void*)(&(model.g.dim)),1,MPI_INT,0,cm);
#endif

    /* set boundary flags */
    IASN(asgpars->lbc, IPNT_0); /* default left bc flag */ 
    IASN(asgpars->rbc, IPNT_0); /* default right bc flag */ 

    for (int i=0;i<model.g.dim;i++) {
      if (crank[i]==0) asgpars->lbc[i]=1;
      if (crank[i]==cdims[i]-1) asgpars->rbc[i]=1;
    }
    /************                          END                      *************/	    

    /* decode half-order - with version 2.0, deprecated syntax "scheme_phys" etc. is dropped */
    asgpars->k = valparse<int>(pars,"order",1);
#ifdef IWAVE_VERBOSE
    fprintf(stream,"NOTE: initializing ASG with half-order = %d\n",asgpars->k);
#endif
    // MAXK set in sgcoeffs.h
    if (asgpars->k<1 || asgpars->k>MAXK) {
      RVLException e;
      e<<"ERROR: asgmodelinit\n";
      e<<"finite difference half-order requested = "<<asgpars->k<<" out of available range [1,"<<MAXK<<"]\n";
      throw e;
    }

    /* set model dimn par */
    asgpars->ndim = model.g.dim;

    /* extract grid steps from grid */
    get_d(dxs, model.g);

    /* initialize scaled Courant arrays */
    for (int idim = 0;idim < asgpars->ndim;idim ++) {
      if (dxs[idim] <= 0.0) {
	RVLException e;
	e<<"ERROR: asg_modelinit\n";
	e<<"  bad input: wrong grid space step, dim="<<idim<<", step="<<dxs[idim]<<"\n";
	throw e;
      }
      lam = model.tsind.dt / dxs[idim];
      for (int i=0;i<asgpars->k;i++) (asgpars->coeffs)[idim][i]=lam*SCHEME_COEFFS[asgpars->k][i];
    }

    /* initialize bound check params */
    asgpars->cmax = valparse<ireal>(pars,"cmax");
    asgpars->cmin = valparse<ireal>(pars,"cmin");

    /* reserve a copy of dt for use in source scaling */
    asgpars->dt = model.tsind.dt;

    /* initialize pml stuff */
    asgpars->amp = valparse<ireal>(pars,"pmlampl",REAL_ZERO);

    for (int i=0;i<RARR_MAX_NDIM; i++) {
      asgpars->ep[i]=NULL;
      asgpars->epp[i]=NULL;
      asgpars->ev[i]=NULL;
      asgpars->evp[i]=NULL;
    }

    /* identify active fields */
    // cerr<<"resize\n";
    model.active.resize(8);
    // cerr<<"assign active\n";
    model.active[0]="bulkmod";
    model.active[1]="buoyancy";
    model.active[2]="p0";
    model.active[5]="v0";
    // cerr<<"dim>1\n";
    if (asgpars->ndim > 1) {
      model.active[3]="p1";
      model.active[6]="v1";
    }
    else {
      model.active[3]="fruitcake";
      model.active[6]="fruitcake";
    }
    // cerr<<"dim>2\n";
    if (asgpars->ndim > 2) {
      model.active[4]="p2";
      model.active[7]="v2";
    }
    else {
      model.active[4]="fruitcake";
      model.active[7]="fruitcake";
    }

    // cerr<<"specs\n";
    /* assign param object pointer */
    model.specs = (void*)asgpars;

    // cerr<<"return\n";
    return 0;
  }
  catch (RVLException & e) {
    e<<"\ncalled from asg_modelinit\n";
    throw e;
  }
}

/*----------------------------------------------------------------------------*/
void asg_modeldest(void ** fdpars) {
  // cerr<<"destroy\n";
    /* destroy asgpars - all data members allocated on stack */
    userfree_(*fdpars);
}

/*----------------------------------------------------------------------------*/
int asg_timegrid(PARARRAY * pars, 
		 FILE * stream, 
		 grid const & g, 
		 ireal & dt) {

  try {

    ireal cmax = valparse<ireal>(*pars,"cmax");
    ireal cfl  = valparse<ireal>(*pars,"cfl");
      
    /* compute max stable step, optionally scaled by cfl from table */
    /* THIS BELONGS IN CHECK
	if ((err=sg_readgustime(dom,stream,&dtgus,g,cflgus,par))) {  
	fprintf(stream,"\ncalled from sg_readtimegrid\n");
	return err;
    }
    */

    /* inline
    if ((err=sg_readcfltime(dom,stream,dt,g,cfl,cmax,par))) {  
	fprintf(stream,"\ncalled from sg_readtimegrid\n");
	return err;
    }
    */
    ireal a = g.axes[0].d;  
    for (int i = 1; i < g.dim; ++i ) a = iwave_min(a, g.axes[i].d);
  
    if ( (a < REAL_EPS) || (cfl < REAL_EPS) ) {
      fprintf(stream,"ERROR: sg_readcfltime - either min dx=%e or cfl=%e "
	      " too small\n", a, cfl);
      return E_BADINPUT;
    }
    dt = a*cfl/(cmax*sqrt((float)(g.dim)));

    /*  MOVE TO CHECK
    if (*dt > dtgus) {
	fprintf(stream,"ERROR: sg_readtimegrid - max_step option\n");
	fprintf(stream,"computed dt based on CFL criterion only = %e\n",*dt);
	fprintf(stream,"exceeds max stable step = %e\n",dtgus);
	return E_BADINPUT;
    }
    */

    /* however if dt is in the parameters, use it and throw the rest away */
    ireal dttmp = dt;
    dt = valparse<ireal>(*pars,"dt",dttmp);

    return 0;
  }
  catch (RVLException & e) {
    e<<"\ncallef from asg_timegrid\n";
    throw e;
  }
}

void asg_pmlaxis(int n0, int nl, int nr, 
		 ireal amp, ireal dt, int gtype, 
		 ireal ** ep, ireal ** epp) {
  
  if (*ep || *epp) return;
  nl=iwave_max(0,nl);
  nr=iwave_max(0,nr);
  *ep  = (ireal*)usermalloc_(sizeof(ireal)*n0);
  *epp = (ireal*)usermalloc_(sizeof(ireal)*n0);
  for (int i=0;i<nl;i++) {
    ireal p = (REAL_ONE - (ireal(i)+0.5*gtype)/ireal(nl));
    p = amp*p*p*p;
    (*ep)[i] = REAL_ONE - 0.5*dt*p;
    (*epp)[i] = REAL_ONE/(REAL_ONE + 0.5*dt*p);
  }
  for (int i=nl; i<n0-nr;i++) {
    (*ep)[i] = REAL_ZERO;
    (*epp)[i]= REAL_ZERO;
  }
  for (int i=n0-nr; i<n0;i++) {
    ireal p =  (ireal(i-n0+nr+1)+0.5*gtype)/ireal(nr);
    p = amp*p*p*p;
    (*ep)[i] = REAL_ONE - 0.5*dt*p;
    (*epp)[i] = REAL_ONE/(REAL_ONE + 0.5*dt*p);
  }
}    

void asg_timestep(std::vector<RDOM *> dom, 
		  bool fwd, 
		  int iv, 
		  void* fdpars) {

  IPNT np;
  IPNT nv;

  /* fd parameter struct */
  ASG_TS_PARS * asgpars = (ASG_TS_PARS *) fdpars;  
    
  /* field indices - from FIELDS struct */
  IPNT i_p;
  i_p[0] = 2;
  i_p[1] = 3;
  i_p[2] = 4;
  IPNT i_v;
  i_v[0] = 5;
  i_v[1] = 6;
  i_v[2] = 7;

  /* fill in pml arrays; these are no-ops after the first call */
  for (int idim=0;idim<asgpars->ndim;idim++) {
    rd_a_size(dom[0], i_p[idim], np);
    rd_a_size(dom[0], i_v[idim], nv);
    asg_pmlaxis(np[idim], asgpars->nls[idim], asgpars->nrs[idim], 
		asgpars->amp, asgpars->dt, 0, 
		&((asgpars->ep)[idim]), &((asgpars->epp)[idim]));
    asg_pmlaxis(nv[idim], asgpars->nls[idim], asgpars->nrs[idim], 
		asgpars->amp, asgpars->dt, 1, 
		&((asgpars->ev)[idim]), &((asgpars->evp)[idim]));
  }
  
  // address list - 2D
  register ireal ** restrict bulk2;
  register ireal ** restrict buoy2;
  register ireal ** restrict p02;
  register ireal ** restrict p12;
  register ireal ** restrict v02;
  register ireal ** restrict v12;

  // address list - 3D
  register ireal *** restrict bulk3;
  register ireal *** restrict buoy3;
  register ireal *** restrict p03;
  register ireal *** restrict p13;
  register ireal *** restrict p23;
  register ireal *** restrict v03;
  register ireal *** restrict v13;
  register ireal *** restrict v23;

  // address list - both
  register ireal ** restrict ep   = asgpars->ep;
  register ireal ** restrict epp  = asgpars->epp;
  register ireal ** restrict ev   = asgpars->ev; 
  register ireal ** restrict evp  = asgpars->evp;
  
  /* sizes of computational domain - P0 same as P1, P2 */
  IPNT gsc_p;
  IPNT gec_p;
  int gsc_v[RARR_MAX_NDIM][RARR_MAX_NDIM];
  int gec_v[RARR_MAX_NDIM][RARR_MAX_NDIM];
  rd_gse(dom[0], i_p[0], gsc_p, gec_p);    
  rd_gse(dom[0], i_v[0], &(gsc_v[0][0]), &(gec_v[0][0]));    
  rd_gse(dom[0], i_v[1], &(gsc_v[1][0]), &(gec_v[1][0]));    
  rd_gse(dom[0], i_v[2], &(gsc_v[2][0]), &(gec_v[2][0]));    

  // divergence workspace
  register ireal * sdiv = (ireal *)usermalloc_((gec_p[0]-gsc_p[0]+1)*sizeof(ireal));
  register ireal ** gradp = (ireal **)usermalloc_(RARR_MAX_NDIM * sizeof(ireal*));
  int ndim;
  ra_ndim(&(dom[0]->_s[i_p[0]]),&ndim);
  for (int i=0;i<ndim; i++) 
    gradp[i] = (ireal *)usermalloc_((gec_v[i][0]-gsc_v[i][0]+1)*sizeof(ireal));

  if (ndim == 2) {
    bulk2 = (dom[0]->_s)[0]._s2;
    buoy2 = (dom[0]->_s)[1]._s2;
    p02 = (dom[0]->_s)[i_p[0]]._s2;
    p12 = (dom[0]->_s)[i_p[1]]._s2;
    v02 = (dom[0]->_s)[i_v[0]]._s2;
    v12 = (dom[0]->_s)[i_v[1]]._s2;
  }
  else if (ndim == 3) {
    bulk3 = (dom[0]->_s)[0]._s3;
    buoy3 = (dom[0]->_s)[1]._s3;
    p03 = (dom[0]->_s)[i_p[0]]._s3;
    p13 = (dom[0]->_s)[i_p[1]]._s3;
    p23 = (dom[0]->_s)[i_p[2]]._s3;
    v03 = (dom[0]->_s)[i_v[0]]._s3;
    v13 = (dom[0]->_s)[i_v[1]]._s3;
    v23 = (dom[0]->_s)[i_v[2]]._s3;
  }

  int n  = dom.size();

  switch (n) {
  case 1: {
    if(fwd == true) {
      if (ndim == 2) {
	if (iv == 0) {
	  asg_pstep2d(bulk2,
		      p02,p12,
		      v02,v12,
		      ep, epp,
		      sdiv,
		      gsc_p,gec_p,
		      asgpars->k,asgpars->coeffs);
	}
	if (iv == 1) {
	  asg_vstep2d(buoy2,
		      p02,p12,
		      v02,v12,
		      ev, evp,
		      gradp,
		      &(gsc_v[0][0]),&(gec_v[0][0]),
		      &(gsc_v[1][0]),&(gec_v[1][0]),
		      asgpars->k,asgpars->coeffs);
	}
      }
      if (ndim == 3) {
	if (iv == 0) {
	  asg_pstep3d(bulk3,
		      p03,p13,p23,
		      v03,v13,v23,
		      ep, epp,
		      sdiv,
		      gsc_p,gec_p,
		      asgpars->k,asgpars->coeffs);
	}
	if (iv == 1) {
	  asg_vstep3d(buoy3,
		      p03,p13,p23,
		      v03,v13,v23,
		      ev, evp,
		      gradp,
		      &(gsc_v[0][0]),&(gec_v[0][0]),
		      &(gsc_v[1][0]),&(gec_v[1][0]),
		      &(gsc_v[2][0]),&(gec_v[2][0]),
		      asgpars->k,asgpars->coeffs);
	}
      }
    }
    if(fwd == false){
      RVLException e;
      e<<"Error: asg_timestep().  iw.size() = "<< n << " has no adjoint!";
      throw e;
    } 
  }
    /*
  case 2: {
    if(fwd == true) {
      return;
    }
    else {
      return;
    }
  }
  case 4: {
    if(fwd == true) {
      return;
    }
    else {
      return;
    }
  }
    */
  default:
    RVLException e;
    e<<"Error: asg_timestep(). Only zero, 1st & 2nd derivatives are implemented.\n";
    throw e;
  }

  userfree_(sdiv);
  for (int i=0; i<ndim; i++) userfree_(gradp[i]);
  userfree_(gradp);
}

/**
 * all stencils are of the same type: for each axis idim, 
 * - source v_idim, target p, offset range -k,...,k-1 on axis idim 
 * - source p, target v_idim, offset range --k+1,...k on axis idim
 * so there are 2*ndim masks in total, each mask storing 2k offset vectors
 */  
int asg_create_sten(void * specs, 
		    FILE * stream, 
		    int ndim, 
		    IPNT gtype[RDOM_MAX_NARR], 
		    STENCIL * sten) {

  // cerr<<"begin stencil\n";

  ASG_TS_PARS * asgpars = (ASG_TS_PARS *)(specs);
  STENCIL_MASK mask;/* workspace */
  int err = 0;

  /* number of masks */
  int nmask = (asgpars->ndim+1)*(asgpars->ndim);

  /* initialize STENCIL to null stencil */
  sten_setnull(sten);
  
  // cerr<<"asg_stencil_create stencil nmask="<<nmask<<endl;
  /* nontrivial STENCIL initialization */
  if ((err = sten_create(sten,nmask))) {
    fprintf(stream,"ERROR: asg_create_sten - failed to create stencil\n");
    return err;
  }

  /* field indices - from FIELDS struct */
  IPNT ip;
  ip[0] = 2;
  ip[1] = 3;
  ip[2] = 4;
  IPNT iv;
  iv[0] = 5;
  iv[1] = 6;
  iv[2] = 7;

  /* every mask stores 2k offsets */
  int len = 2*asgpars->k;

  /* loop over axes */

  int imask=0; // mask counter

  for (int idim=0; idim<asgpars->ndim; idim++) {
    
    for (int jdim=0; jdim<asgpars->ndim; jdim++) {
      // cerr<<"asg_stencil_create imask="<<imask<<" idim="<<idim<<" jdim="<<jdim<<" mask workspace\n";
      /* source = v[idim], target = p[jdim] */
      if ((err=mask_create(&mask, iv[idim], ip[jdim], len))) {
	fprintf(stream,"ERROR: asg_create_sten from mask_create\n");
	sten_destroy(sten);
	return err;
      }
      // cerr<<"asg_stencil_create: idim="<<idim<<endl;
      for (int i=0; i<len; i++) {
	IPNT offs;
	IASN(offs,IPNT_0);
	offs[idim]=i-asgpars->k;
	if ((err=mask_set(&mask,i,offs))) {
	  fprintf(stream,"ERROR: asg_create_sten from mask_set\n");	
	  sten_destroy(sten);
	  return err;
	}	
      }
      /* insert mask */
      if ((err=sten_set(sten,imask, &mask))) {
	fprintf(stream,"ERROR: asg_create_sten from sten_set\n");	
	sten_destroy(sten);
	return err;
      }
      
      imask++;

    }

    /* source = p, target = v[idim] */
    // cerr<<"asg_stencil_create mask workspace\n";
      // cerr<<"asg_stencil_create imask="<<imask<<" idim="<<idim<<" mask workspace\n";
    if ((err=mask_create(&mask, ip[idim], iv[idim], len))) {
      fprintf(stream,"ERROR: asg_create_sten from mask_create\n");
      sten_destroy(sten);
      return err;
    }
    for (int i=0; i<len; i++) {
      IPNT offs;
      IASN(offs,IPNT_0);
      offs[idim]=i-asgpars->k+1;
      if ((err=mask_set(&mask,i,offs))) {
	fprintf(stream,"ERROR: asg_create_sten from mask_set\n");	
	sten_destroy(sten);
	return err;
      }	
    }
    /* insert mask */
    if ((err=sten_set(sten,imask, &mask))) {
      fprintf(stream,"ERROR: asg_create_sten from sten_set\n");	
      sten_destroy(sten);
      return err;
    }
    imask++;
  }
  /* end loop over axes */

  // cerr<<"end stencil\n";
  return 0;
}

void asg_check(RDOM * dom,
	       void * specs,
	       FILE * stream) {}


