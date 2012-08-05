#include "sgn.h"
#include "iwave.h"

/* 
   Functions to convert between array names and indices.
*/

/*------------------- private data declarations ----------------------------*/

static int m_ndim = 0; /* dimension - need only compute once */
static int m_size = 17; /* 16; */
static const char* m_names[] = {"p0", "kappa", "v0","rho0","v1","rho1","v2","rho2",
				"eta_p0", "eta_v0", "p1", "eta_p1", "eta_v1",
				"p2", "eta_p2", "eta_v2","NONE","NONE","NONE","NONE","NONE"};

/*----------------------------------------------------------------------------*/

/*--- private function declarations - assigned to FD_MODEL pointers ----------*/

int asg_isarr(int i);
int asg_numsubsteps();
int asg_update();
int asg_readschemeinfo(PARARRAY *, FILE *, IMODEL *);
int asg_set_grid_type(FILE *, int, IPNT[RDOM_MAX_NARR]);
int asg_build_sten_dep(FILE *, int, int[RDOM_MAX_NARR][RDOM_MAX_NARR]);
int asg_create_sten(FD_MODEL *, FILE *, int, IPNT[RDOM_MAX_NARR], int[RDOM_MAX_NARR][RDOM_MAX_NARR], STENCIL *);
int asg_alter_dom(int, IPNT, IPNT);
const char * asg_ind2str(int);  
int asg_modeldest(IMODEL * model);
void sgn_ts_parcopy(void * tgt, const void * src);

extern int asg_step(RDOM*, int, void *);

/*----------------------------------------------------------------------------*/

const char* asg_ind2str(int ind) {
  if (ind>-1 && ind<m_size) return m_names[ind];
  return m_names[m_size];
}
/*----------------------------------------------------------------------------*/

int asg_modelinit(PARARRAY *pars, 
		  FILE *stream, 
		  IMODEL * im) {

  int err=0;           /* return value */
  SGN_TS_PARS *sgnm;   /* model pars */
  FD_MODEL * fd;       /* function struct */

  IPNT cdims;          /* workspace for cartesian grid dim info */
  IPNT crank;          /* workspace for cartesian grid rank */
#ifdef IWAVE_USE_MPI
  IPNT cpers;          /* workspace for periodic wrap info  - currently not used */
#endif

  int i;               /* counter */

  fd=(FD_MODEL *)usermalloc_(sizeof(FD_MODEL));
  if (fd==NULL) return E_ALLOC;
  /* first the "base class" behaviour */
  //  fprintf(stream,"top of asg_modelinit: im->specs=%x\n",im->specs);
  im->specs=(void *)fd;
  //  fprintf(stream,"top of asg_modelinit: im->specs=%x after assignment\n",im->specs);

  /* allocate sgn model ----------------------------------------------------*/   
  sgnm = (SGN_TS_PARS*)usermalloc_(sizeof(SGN_TS_PARS));
  if ( sgnm == NULL ) { 
    err=E_ALLOC;
    fprintf(stream,"ERROR: asg_modelinit\n");
    fprintf(stream,"failed to allocate SGN_TS_PARS object\n");
    return err;
  }
	    
  // assign heritage pars
  // 11.04.12: remove these
  //  sgnm->psingle = 0; /* TODO: single pressures in Step 3 */
  //  IASN(sgnm->eflags,IPNT_1); /* this is related - with code from old sg added, could try using this device */

  // decode dimensions, parallel rank - read grid dimn on rank 0, broadcast

  IASN(cdims, IPNT_1); /* default grid size */ 
  IASN(crank, IPNT_0); /* default cartisian ranks */ 

#ifdef IWAVE_USE_MPI
  int rk=retrieveRank();
  MPI_Comm cm=retrieveComm();

  if ( MPI_Cart_get(cm, RARR_MAX_NDIM, cdims, cpers, crank) != MPI_SUCCESS )  {
    fprintf(stream, "ERROR. Internal: cannot get Cartesian coordinates.\n");
    return E_INTERNAL;
  }

  if (rk==0) {
#endif
    m_ndim=0;
    
    grid g;
    if (!init_acoustic_geom_par(&g,*pars,stream)) m_ndim=g.dim;
    else {
      err=E_INTERNAL;
      fprintf(stream,"ERROR: in asg_modelinit - failed to read spatial geometry\n");
      fflush(stream);
      abortexit(err,pars,&stream);
    }
    
#ifdef IWAVE_USE_MPI
  }
  MPI_Bcast(&m_ndim,1,MPI_INT,0,cm);
#endif

  if (m_ndim<1 || m_ndim>RARR_MAX_NDIM) {
    err=E_INTERNAL;
    fprintf(stream,"ERROR: in asg_modelinit - failed to read dim=%d\n",m_ndim);
    return err;
  }

  /* set boundary flags */
  IASN(sgnm->lbc, IPNT_0); /* default left bc flag */ 
  IASN(sgnm->rbc, IPNT_0); /* default right bc flag */ 

  for (i=0;i<m_ndim;i++) {
    if (crank[i]==0) sgnm->lbc[i]=1;
    if (crank[i]==cdims[i]-1) sgnm->rbc[i]=1;
  }

  /* decode order - with version 2.0, deprecated syntax "scheme_phys" etc. is dropped */
  sgnm->k=1;
  ps_flint(*pars,"order",&(sgnm->k));

  // initialize scaled Courant arrays
  RASN(sgnm->c11,RPNT_0);
  RASN(sgnm->c12,RPNT_0);
  RASN(sgnm->c22,RPNT_0);
  RASN(sgnm->c14,RPNT_0);
  RASN(sgnm->c24,RPNT_0);
  RASN(sgnm->c34,RPNT_0);
  RASN(sgnm->c44,RPNT_0);

  // initialize aux pml arrays
  sgnm->nep0=0;
  sgnm->nev0=0;
  sgnm->nep1=0;
  sgnm->nev1=0;
  sgnm->nep2=0;
  sgnm->nev2=0;
  sgnm->ep0_p=NULL;
  sgnm->ev0_p=NULL;
  sgnm->ep1_p=NULL;
  sgnm->ev1_p=NULL;
  sgnm->ep2_p=NULL;
  sgnm->ev2_p=NULL;
  sgnm->ep0_pp=NULL;
  sgnm->ev0_pp=NULL;
  sgnm->ep1_pp=NULL;
  sgnm->ev1_pp=NULL;
  sgnm->ep2_pp=NULL;
  sgnm->ev2_pp=NULL;

  // assign param object pointer
  fd->fdpars = (void*)sgnm;

  fd->isarr  = asg_isarr;
  fd->numsubsteps = asg_numsubsteps;
  fd->update = asg_update;
  fd->readgrid = asg_readgrid;
  fd->readtimegrid = asg_readtimegrid; 
  fd->readschemeinfo = asg_readschemeinfo;
  fd->set_grid_type = asg_set_grid_type;
  fd->build_sten_dep = asg_build_sten_dep;
  fd->create_sten = asg_create_sten;
  fd->ind2str = asg_ind2str;
  fd->alter_dom = asg_alter_dom;
  fd->readmedia = asg_readmedia;
  fd->parcopy = sgn_ts_parcopy;
  fd->fd_model_init = asg_modelinit;
  fd->fd_model_dest = asg_modeldest;

  // choose time step
  fd->tsf=asg_step;

  //  fprintf(stream,"\n********** asg indices *************\n");
  /*
    fprintf(stream,"D_P0=%d\n",D_P0);
    fprintf(stream,"D_MP0=%d\n",D_MP0);
    fprintf(stream,"D_V0=%d\n",D_V0);
    fprintf(stream,"D_MV0=%d\n",D_MV0);
    fprintf(stream,"D_V1=%d\n",D_V1);
    fprintf(stream,"D_MV1=%d\n",D_MV1);
    fprintf(stream,"D_V2=%d\n",D_V2);
    fprintf(stream,"D_MV2=%d\n",D_MV2);
    fprintf(stream,"D_EP0=%d\n",D_EP0);
    fprintf(stream,"D_EV0=%d\n",D_EV0);
    fprintf(stream,"D_P1=%d\n",D_P1);
    fprintf(stream,"D_EP1=%d\n",D_EP1);
    fprintf(stream,"D_EV1=%d\n",D_EV1);
    fprintf(stream,"D_P2=%d\n",D_P2);
    fprintf(stream,"D_EP2=%d\n",D_EP2);
    fprintf(stream,"D_EV2=%d\n",D_EV2);
    fprintf(stream,"\n********** asg indices *************\n");  
  */

  return 0;
}

int asg_numsubsteps() { return 2; }

int asg_modeldest(IMODEL * model) {
  FD_MODEL * fdm = (FD_MODEL *)(model->specs);
  SGN_TS_PARS * sgnm;
  if (fdm) {
    if (fdm->fdpars) {
      sgnm=(SGN_TS_PARS *)(fdm->fdpars);
      if (sgnm->ep0_p) userfree_(sgnm->ep0_p);
      if (sgnm->ep0_pp) userfree_(sgnm->ep0_pp);
      if (sgnm->ev0_p) userfree_(sgnm->ev0_p);
      if (sgnm->ev0_pp) userfree_(sgnm->ev0_pp);
      if (sgnm->ep1_p) userfree_(sgnm->ep1_p);
      if (sgnm->ep1_pp) userfree_(sgnm->ep1_pp);
      if (sgnm->ev1_p) userfree_(sgnm->ev1_p);
      if (sgnm->ev1_pp) userfree_(sgnm->ev1_pp);
      if (sgnm->ep2_p) userfree_(sgnm->ep2_p);
      if (sgnm->ep2_pp) userfree_(sgnm->ep2_pp);
      if (sgnm->ev2_p) userfree_(sgnm->ev2_p);
      if (sgnm->ev2_pp) userfree_(sgnm->ev2_pp);
      userfree_(sgnm);
    }
    // since the FD_MODEL is allocated here, destroy here
    userfree_(fdm);
  }
  return 0;
}

void sgn_ts_parcopy(void * tgt, const void * src) {

  int i;

  SGN_TS_PARS * ptgt = (SGN_TS_PARS *)tgt;
  const SGN_TS_PARS * psrc = (const SGN_TS_PARS *)src;

  ptgt->dt=psrc->dt;
  for (i=0;i<RARR_MAX_NDIM;i++) 
    (ptgt->lam)[i]=(psrc->lam)[i];
  ptgt->ndim=psrc->ndim;
  ptgt->k = psrc->k;
  //ptgt->psingle=psrc->psingle;

  //  IASN(ptgt->eflags,psrc->eflags);
  IASN(ptgt->lbc,psrc->lbc);
  IASN(ptgt->rbc,psrc->rbc);   

  if ((ptgt->nep0=psrc->nep0)) {
    if (!(ptgt->ep0_p)) ptgt->ep0_p = (ireal *)usermalloc_(ptgt->nep0*sizeof(ireal));
    if (!(ptgt->ep0_pp)) ptgt->ep0_pp = (ireal *)usermalloc_(ptgt->nep0*sizeof(ireal));
    for (i=0;i<ptgt->nep0;i++) {
      ptgt->ep0_p[i]=psrc->ep0_p[i];
      ptgt->ep0_pp[i]=psrc->ep0_pp[i];
    }
  }
  if ((ptgt->nev0=psrc->nev0)) {
    if (!(ptgt->ev0_p)) ptgt->ev0_p = (ireal *)usermalloc_(ptgt->nev0*sizeof(ireal));
    if (!(ptgt->ev0_pp)) ptgt->ev0_pp = (ireal *)usermalloc_(ptgt->nev0*sizeof(ireal));
    for (i=0;i<ptgt->nev0;i++) {
      ptgt->ev0_p[i]=psrc->ev0_p[i];
      ptgt->ev0_pp[i]=psrc->ev0_pp[i];
    }
  }
  if ((ptgt->nep1=psrc->nep1)) {
    if (!(ptgt->ep1_p)) ptgt->ep1_p = (ireal *)usermalloc_(ptgt->nep1*sizeof(ireal));
    if (!(ptgt->ep1_pp)) ptgt->ep1_pp = (ireal *)usermalloc_(ptgt->nep1*sizeof(ireal));
    for (i=0;i<ptgt->nep1;i++) {
      ptgt->ep1_p[i]=psrc->ep1_p[i];
      ptgt->ep1_pp[i]=psrc->ep1_pp[i];
    }
  }
  if ((ptgt->nev1=psrc->nev1)) {
    if (!(ptgt->ev1_p)) ptgt->ev1_p = (ireal *)usermalloc_(ptgt->nev1*sizeof(ireal));
    if (!(ptgt->ev1_pp)) ptgt->ev1_pp = (ireal *)usermalloc_(ptgt->nev1*sizeof(ireal));
    for (i=0;i<ptgt->nev1;i++) {
      ptgt->ev1_p[i]=psrc->ev1_p[i];
      ptgt->ev1_pp[i]=psrc->ev1_pp[i];
    }
  }
  if ((ptgt->nep2=psrc->nep2)) {
    if (!(ptgt->ep2_p)) ptgt->ep2_p = (ireal *)usermalloc_(ptgt->nep2*sizeof(ireal));
    if (!(ptgt->ep2_pp)) ptgt->ep2_pp = (ireal *)usermalloc_(ptgt->nep2*sizeof(ireal));
    for (i=0;i<ptgt->nep2;i++) {
      ptgt->ep2_p[i]=psrc->ep2_p[i];
      ptgt->ep2_pp[i]=psrc->ep2_pp[i];
    }
  }
  if ((ptgt->nev2=psrc->nev2)) {
    if (!(ptgt->ev2_p)) ptgt->ev2_p = (ireal *)usermalloc_(ptgt->nev2*sizeof(ireal));
    if (!(ptgt->ev2_pp)) ptgt->ev2_pp = (ireal *)usermalloc_(ptgt->nev2*sizeof(ireal));
    for (i=0;i<ptgt->nev2;i++) {
      ptgt->ev2_p[i]=psrc->ev2_p[i];
      ptgt->ev2_pp[i]=psrc->ev2_pp[i];
    }
  }

}

int asg_isarr(int i) {
  //  fprintf(stderr,"asg_isarr: ndim=%d i=%d\n",m_ndim,i);
  if (m_ndim<4) {
    if (m_ndim>0) {
      if (i==D_P0) return 1;
      if (i==D_V0) return 1;
      if (i==D_MP0) return 1;
      if (i==D_MV0) return 1;
      if (i==D_EP0) return 1;
      if (i==D_EV0) return 1;
    }
    if (m_ndim>1) {
      if (i==D_P1) return 1;
      if (i==D_V1) return 1;
      if (i==D_MV1) return 1;
      if (i==D_EP1) return 1;
      if (i==D_EV1) return 1;
    }
    if (m_ndim>2) {
      if (i==D_P2) return 1;
      if (i==D_V2) return 1;
      if (i==D_MV2) return 1;
      if (i==D_EP2) return 1;
      if (i==D_EV2) return 1;
    }
  }
  return 0;
}

int asg_update(int ia, int iv) {
  // update pressure arrays on iv = 0
  if (iv==0) {
    if (m_ndim>2 && ia==D_P2) return 1;
    if (m_ndim>1 && ia==D_P1) return 1;
    if (m_ndim>0 && ia==D_P0) return 1;
  }

  if (iv==1) {
    if (m_ndim>2 && ia==D_V2) return 1;
    if (m_ndim>1 && ia==D_V1) return 1;
    if (m_ndim>0 && ia==D_V0) return 1;
  }

  return 0;
}

int asg_alter_dom(int iv, IPNT gs, IPNT ge) {
  int idim,i;

  /* one-dimensionalize the eta arrays */
  for ( idim = 0; idim < m_ndim ; ++idim ) {
    if (iv==D_EP[idim] || iv==D_EV[idim]) {
      gs[0]=gs[idim];
      ge[0]=ge[idim];
      //      fprintf(stderr,"iv=%d idim=%d\n",iv,idim);
      for ( i = 1; i < m_ndim; ++i ) 
	gs[i] = ge[i] = 0;
    }
  }
  return 0;
}

int asg_set_grid_type(FILE *stream, int ndim, IPNT gtype[RDOM_MAX_NARR] ) {
  
  int iv, i; 
  if ( ndim < 1 || ndim > RARR_MAX_NDIM ) return E_BADINPUT;
  /*< px, py, pz, mp00, mp01 on the primal grid */
  for (iv = 0;iv < RDOM_MAX_NARR;iv ++)  IASN(gtype[iv], IPNT_0);
  
  for (i = 0;i < ndim;i ++) {
    gtype[D_V[i]][i]  = 1;  
    gtype[D_MV[i]][i] = 1;
    gtype[D_EV[i]][i] = 1;
  }

  return 0;
}

/**
 * vary with different PDEs
 * define the stencil dependent matrix 
 * relation type: dependent of value -> DEP_F
 *                dependent of 1st derivative wrt axis-z -> DEP_DFDZ
 *                dependent of 1st derivative wrt axis-x -> DEP_DFDX
 *                dependent of 1st derivative wrt axis-y -> DEP_DFDY
 */
/*< user must provide sten_dep_mat information */
int asg_build_sten_dep(FILE *stream, int ndim, int sten_dep_mat[RDOM_MAX_NARR][RDOM_MAX_NARR]) {

  int i, j;

  if (ndim < 1 || ndim > RARR_MAX_NDIM) return E_BADINPUT;

  for (i = 0;i < RDOM_MAX_NARR;i ++)
    for (j = 0;j < RDOM_MAX_NARR;j ++)
      sten_dep_mat[i][j] = 0;

  /*< participating arrays for normal stresses */
  /*< DEP_DFDZ: 1  DEP_DFDX: 2  DEP_DFDY: 3 */
  for (i = 0;i < ndim;i ++) {
    for (j = 0;j < ndim;j ++) {
      sten_dep_mat[D_P[i]][D_V[j]] = j+1;
    }
  }
  
  for (i = 0;i < ndim;i ++) {
    sten_dep_mat[D_P[i]][D_MP0] = DEP_F;
  }
  
  /*< participating arrays for velocity */
  for (i = 0;i < ndim;i ++) {
    sten_dep_mat[D_V[i]][D_MV[i]] = DEP_F;
    sten_dep_mat[D_V[i]][D_P[i]] = i+1;
  }

  return 0;
}

/*----------------------------------------------------------------------------*/

int sgn_modelinfo(FILE *stream, IMODEL *model)
{
  FD_MODEL *sgnm;
  int ndim, i, iv, ia;
  IPNT n;
  RDOM *dc; /* computational */

  sgnm = (FD_MODEL*)(model->specs);
  ndim = model->g.dim;

  fprintf(stream, "SGN model %dD.\n", ndim);

  dc = &(model->ld_c); /* get computational domain */
  rd_size(dc, D_P0, n);

  fprintf(stream, "  P0 size = [%d", n[0]);
  for ( i = 1; i < ndim; ++i ) fprintf(stream, " %d", n[i]);
  fprintf(stream, "]\n");

  for ( iv = 0; iv < ndim; ++iv )
    {
      ia = D_V[iv];
      rd_size(dc, ia, n);
      fprintf(stream, "  V%d size = [%d", iv, n[0]);
      for ( i = 1; i < ndim; ++i ) fprintf(stream, " %d", n[i]);
      fprintf(stream, "]\n");
    }

  return 0;
}

int asg_create_sten(FD_MODEL * fdm,
		    FILE * stream, 
		    int ndim,
		    IPNT gtype[RDOM_MAX_NARR], 
		    int sten_dep_mat[RDOM_MAX_NARR][RDOM_MAX_NARR], 
		    STENCIL * sten) {

  SGN_TS_PARS * sgnp = (SGN_TS_PARS *)(fdm->fdpars);
  //  return create_sten2_2k(stream, sgnp->k, ndim, m_size, gtype, sten_dep_mat, asg_isdyn, //asg_getindices, 
  return create_sten2_2k(fdm,stream, sgnp->k, ndim, gtype, sten_dep_mat, sten);
}

int asg_readschemeinfo(PARARRAY * par, 
		       FILE * stream, 
		       IMODEL * model) {
  
  RPNT dxs;       // grid steps
  int idim, i;    // counters
  IPNT gsa, gea;  // workspace for aux PML calc
  ireal * tmp;    // "
  ireal dt2;      // "

  RDOM * dom = &(model->ld_a);
  FD_MODEL * fdm = (FD_MODEL *)(model->specs);
  SGN_TS_PARS * sgnp = (SGN_TS_PARS *)(fdm->fdpars);

  get_d(dxs, model->g);
  // set model dimn par
  sgnp->ndim = (model->g).dim;

  for (idim = 0;idim < sgnp->ndim;idim ++) {

    if (dxs[idim] <= 0.0) {
      fprintf(stream, "Error: bad input: wrong grid space step, dim=%d, step=%g\n",
	      idim, dxs[idim]);
      return E_BADINPUT;
    }
    sgnp->lam[idim] = (model->tsind).dt / dxs[idim];
#ifdef VERBOSE
    fprintf(stderr, "lam[%d] = %g\n", idim, sgnp->lam[idim]);
#endif

    // assign scaled Courant numbers for orders 2, 4, and 8 - these are the only
    // choices implemented
    if (sgnp->k==1) {
      sgnp->c11[idim] = sgnp->lam[idim]*COEFF1[0];
    }
    else if (sgnp->k==2) {
      sgnp->c12[idim] = sgnp->lam[idim]*COEFF2[0];
      sgnp->c22[idim] = sgnp->lam[idim]*COEFF2[1];
    }
    else if (sgnp->k==4) {
      sgnp->c14[idim] = sgnp->lam[idim]*COEFF4[0];
      sgnp->c24[idim] = sgnp->lam[idim]*COEFF4[1];
      sgnp->c34[idim] = sgnp->lam[idim]*COEFF4[2];
      sgnp->c44[idim] = sgnp->lam[idim]*COEFF4[3];
    }
    else {
      fprintf(stream,"ERROR: asg_readschemeinfo\n");
      fprintf(stream,"assigned scheme half-order = %d not defined\n",sgnp->k);
      fprintf(stream,"currently defined schemes: half-orders 1, 2, and 4\n");
      fflush(stream);
      return E_BADINPUT;
    }
  }

  // reserve a copy of dt for use in source scaling
  sgnp->dt = (model->tsind).dt;

  fprintf(stream,"NOTE: asg_readschemeinfo: initialize aux PML arrays\n");
  fflush(stream);

  if (sgnp->nep0 || sgnp->ep0_p || sgnp->ep0_pp ||
      sgnp->nev0 || sgnp->ev0_p || sgnp->ev0_pp ||
      sgnp->nep1 || sgnp->ep1_p || sgnp->ep1_pp ||
      sgnp->nev1 || sgnp->ev1_p || sgnp->ev1_pp ||
      sgnp->nep2 || sgnp->ep2_p || sgnp->ep2_pp ||
      sgnp->nev2 || sgnp->ev2_p || sgnp->ev2_pp) {

    fprintf(stream,"ERROR: sgn_setetas \n");
    fprintf(stream,"-- aux pml arrays already initialized, total confusion!\n");
    fflush(stream);
    return E_BADINPUT;
  }

  dt2 = sgnp->dt / ((ireal)2.0);

#if RARR_MAX_NDIM > 0
  // p0
  rd_gse(dom, D_EP[0], gsa, gea);  
  sgnp->nep0 = gea[0]-gsa[0]+1;
  if (sgnp->nep0) {
    sgnp->ep0_p = (ireal *)usermalloc_(sizeof(ireal)*(sgnp->nep0));
    sgnp->ep0_pp = (ireal *)usermalloc_(sizeof(ireal)*(sgnp->nep0));
    tmp = dom->_s[D_EP[0]]._s0;
    for (i=0;i<sgnp->nep0;i++) { 
      sgnp->ep0_p[i]=REAL_ONE - tmp[i]*dt2;
      sgnp->ep0_pp[i]=REAL_ONE/(REAL_ONE + tmp[i]*dt2);
    }
  }

  // v0
  rd_gse(dom, D_EV[0], gsa, gea);  
  sgnp->nev0 = gea[0]-gsa[0]+1;
  if (sgnp->nev0) {
    sgnp->ev0_p = (ireal *)usermalloc_(sizeof(ireal)*(sgnp->nev0));
    sgnp->ev0_pp = (ireal *)usermalloc_(sizeof(ireal)*(sgnp->nev0));
    tmp = dom->_s[D_EV[0]]._s0;
    for (i=0;i<sgnp->nev0;i++) { 
      sgnp->ev0_p[i]=REAL_ONE/(REAL_ONE + tmp[i]*dt2);
      sgnp->ev0_pp[i]=(REAL_ONE-tmp[i]*dt2)*sgnp->ev0_p[i];
    }
  }
#endif

#if RARR_MAX_NDIM > 1
  // p1
  rd_gse(dom, D_EP[1], gsa, gea);  
  sgnp->nep1 = gea[0]-gsa[0]+1;
  if (sgnp->nep1) {
    sgnp->ep1_p = (ireal *)usermalloc_(sizeof(ireal)*(sgnp->nep1));
    sgnp->ep1_pp = (ireal *)usermalloc_(sizeof(ireal)*(sgnp->nep1));
    tmp = dom->_s[D_EP[1]]._s0;
    for (i=0;i<sgnp->nep1;i++) { 
      sgnp->ep1_p[i]=REAL_ONE - tmp[i]*dt2;
      sgnp->ep1_pp[i]=REAL_ONE/(REAL_ONE + tmp[i]*dt2);
    }
  }
  // v1
  rd_gse(dom, D_EV[1], gsa, gea);  
  sgnp->nev1 = gea[0]-gsa[0]+1;
  if (sgnp->nev1) {
    sgnp->ev1_p = (ireal *)usermalloc_(sizeof(ireal)*(sgnp->nev1));
    sgnp->ev1_pp = (ireal *)usermalloc_(sizeof(ireal)*(sgnp->nev1));
    tmp = dom->_s[D_EV[1]]._s0;
    for (i=0;i<sgnp->nev1;i++) { 
      sgnp->ev1_p[i]=REAL_ONE/(REAL_ONE + tmp[i]*dt2);
      sgnp->ev1_pp[i]=(REAL_ONE-tmp[i]*dt2)*(sgnp->ev1_p[i]);
    }
  }
#endif

#if RARR_MAX_NDIM > 2
  // p2
  rd_gse(dom, D_EP[2], gsa, gea);  
  sgnp->nep2 = gea[0]-gsa[0]+1;
  if (sgnp->nep2) {
    sgnp->ep2_p = (ireal *)usermalloc_(sizeof(ireal)*(sgnp->nep2));
    sgnp->ep2_pp = (ireal *)usermalloc_(sizeof(ireal)*(sgnp->nep2));
    tmp = dom->_s[D_EP[2]]._s0;
    for (i=0;i<sgnp->nep2;i++) { 
      sgnp->ep2_p[i]=REAL_ONE - tmp[i]*dt2;
      sgnp->ep2_pp[i]=REAL_ONE/(REAL_ONE + tmp[i]*dt2);
    }
  }

  // v2
  rd_gse(dom, D_EV[2], gsa, gea);  
  sgnp->nev2 = gea[0]-gsa[0]+1;
  if (sgnp->nev2) {
    sgnp->ev2_p = (ireal *)usermalloc_(sizeof(ireal)*(sgnp->nev2));
    sgnp->ev2_pp = (ireal *)usermalloc_(sizeof(ireal)*(sgnp->nev2));
    tmp = dom->_s[D_EV[2]]._s0;
    for (i=0;i<sgnp->nev2;i++) { 
      sgnp->ev2_p[i]=REAL_ONE/(REAL_ONE + tmp[i]*dt2);
      sgnp->ev2_pp[i]=(REAL_ONE-tmp[i]*dt2)*(sgnp->ev2_p[i]);
    }
  }
#endif

  fprintf(stream,"NOTE: asg_schemeinfo\n");
  fprintf(stream," -- initialized aux PML arrays with dimensions\n");
  fprintf(stream," -- nep0=%d nev0=%d\n",sgnp->nep0,sgnp->nev0);
  fprintf(stream," -- nep1=%d nev1=%d\n",sgnp->nep1,sgnp->nev1);
  fprintf(stream," -- nep2=%d nev2=%d\n",sgnp->nep2,sgnp->nev2);
  fflush(stream);

  return 0;

}


