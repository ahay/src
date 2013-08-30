#include "iwave.h"
#include "fd.h"
#include "esgn_read.h"
#include "esgn.h"

#define IWAVE_VERBOSE


/*------------------- private data declarations ----------------------------*/

static int m_ndim = 0; /* dimension - need only compute once */
static int m_size = 24;

static const char *m_names[] = {"pz", "mp00", "mp01", "vz", "rhoz", "szx", "ms0", 
                         "px", "vx", "rhox", "sxy", "ms1",
                         "py", "vy", "rhoy", "szy", "ms2",
                         "eta_p0", "eta_v0", "eta_p1", "eta_v1",
                         "eta_p2", "eta_v2", "NONE"};

//MB: moved to esgn.h
// /*--- private function declarations - assigned to FD_MODEL pointers ----------*/
// 
// int esg_isarr(int i);
// int esg_numsubsteps();
// int esg_update();
// int esg_readschemeinfo(PARARRAY *, FILE *, IMODEL *);
// int esg_set_grid_type(FILE *, int, IPNT[RDOM_MAX_NARR]);
// int esg_build_sten_dep(FILE *, int, int[RDOM_MAX_NARR][RDOM_MAX_NARR]);
// int esg_create_sten(FD_MODEL *, FILE *, int, IPNT[RDOM_MAX_NARR], int[RDOM_MAX_NARR][RDOM_MAX_NARR], STENCIL *);
// int esg_alter_dom(int, IPNT, IPNT);
// const char * asg_ind2str(int);  
// int esg_modeldest(IMODEL * model);
// void esgn_ts_parcopy(void * tgt, const void * src);



/*----------------------------------------------------------------------------*/
const char* esg_ind2str(int ind) {
/*----------------------------------------------------------------------------*/
  if (ind>-1 && ind<m_size) return m_names[ind];
  return m_names[m_size];
}

/*----------------------------------------------------------------------------*/
int esg_modelinit(PARARRAY *pars, FILE *stream, IMODEL * im) {
/*----------------------------------------------------------------------------*/
#ifdef IWAVE_VERBOSE
  fprintf(stream,">>> Inside esg_modelinit\n");
#endif
  
  int err=0;           /* return value */
  ESGN_TS_PARS *esgnp; /* model pars */
  FD_MODEL * fd;       /* function struct */

  IPNT cdims;
  IPNT crank;

#ifdef IWAVE_USE_MPI
    IPNT cpers;          /* workspace for periodic wrap info  - currently not used */
#endif

//   char * jnk;          /* workspace */



	    
  /* decode dimension - only on rk 0 */
#ifdef IWAVE_USE_MPI
  int rk; //=retrieveRank(); RN
  MPI_Comm cm; //=retrieveComm(); RN
  //  if (rk==0) { RN

#endif
    grid g;
//MB	replaced: fd=(FD_MODEL *)malloc(sizeof(FD_MODEL));
  fd=(FD_MODEL *)usermalloc_(sizeof(FD_MODEL));

  if (fd==NULL) return E_ALLOC;

  im->specs=(void *)fd;


  /* allocate sgn model */ 
  //MB replaced: esgnp = (ESGN_TS_PARS*)malloc(sizeof(ESGN_TS_PARS));
  esgnp = (ESGN_TS_PARS*)usermalloc_(sizeof(ESGN_TS_PARS));


  if ( esgnp == NULL ) { 
    err=E_ALLOC;
    fprintf(stream,"ERROR: esg_modelinit\n");
    fprintf(stream,"failed to allocate ESGN_TS_PARS object\n");
    return err;
  }
  m_ndim=0; //RN

  IASN(cdims,IPNT_1);
  IASN(crank,IPNT_0);

#ifdef IWAVE_USE_MPI
    rk=retrieveRank();
    cm=retrieveComm();

    if ( MPI_Cart_get(cm, RARR_MAX_NDIM, cdims, cpers, crank) != MPI_SUCCESS )  {
	fprintf(stream, "ERROR. Internal: cannot get Cartesian coordinates.\n");
	return E_INTERNAL;
    }
    if (rk==0) {
#endif
      // m_ndim=0;

    if (!init_elastic_geom_par(&g,*pars,stream)) m_ndim=g.dim;   
	else {
	    err=E_INTERNAL;
	    fprintf(stream,"ERROR: in esgn_modelinit - failed to read spatial geometry\n");
	    fflush(stream);
	    abortexit(err,pars,&stream);
	}

#ifdef IWAVE_USE_MPI
    }
  

  MPI_Bcast(&m_ndim,1,MPI_INT,0,cm);

#endif
  if (m_ndim<1 || m_ndim>RARR_MAX_NDIM) {
    err=E_INTERNAL;
    fprintf(stream,"ERROR: in esg_modelinit - failed to read dim=%d\n",m_ndim);
    return err;
  }



  /* decode order - take care of deprecated cases */
  esgnp->k=1;
 ps_flint(*pars,"order",&(esgnp->k));

#ifdef IWAVE_VERBOSE
    fprintf(stream,"NOTE: initializing ESG with half-order = %d\n",esgnp->k);
#endif

/*   /\* main branch *\/ */
/*   if (ps_ffint(*pars,"order",&(esgnp->k))) { */
/*     fprintf(stream," RN: read k and is %d \n",esgnp->k); */
/* //MB	comment: possible deprication */
/*     /\* old-style iwave specs - deprecated *\/ */
/*       if ((err=ps_ffcstring(*pars,"scheme_phys",&jnk))) { */
/* 	  fprintf(stream,"ERROR: asg_modelinit from ps_ffstring\n"); */
/* 	  fprintf(stream,"failed to read deprecated order param scheme_phys\n"); */
/* 	  return err; */
/*       } */
/*     else { */
/*       if (!strcmp(jnk,"22")) esgnp->k=1; */
/*       else if (!strcmp(jnk,"24")) esgnp->k=2; */
/*       else if (!strcmp(jnk,"210")) esgnp->k=5; */
/*       else if (!strcmp(jnk,"2k")) { */
/* 	  if ((err=ps_ffint(*pars,"k_phys",&(esgnp->k)))) { */
/* 	      fprintf(stream,"ERROR: esg_modelinit from ps_ffint\n"); */
/* 	      fprintf(stream,"failed to read deprecated order param k_phys\n"); */
/* 	      return err; */
/* 	  } */
/*       } */
/*       else { */
/* 	err=E_BADINPUT; */
/* 	fprintf(stream,"ERROR: asg_modelinit\n");	 */
/* 	fprintf(stream,"invalid value for deprecated param scheme_phys\n"); */
/* 	fprintf(stream,"legal values are 22, 24, 210 (3D only), and 2k\n"); */
/* 	return err; */
/*       } */
/*     } */
/*   } */
/*   /\* sanity-check order *\/ */
/*   if (esgnp->k<1 || esgnp->k>8) { */
/*     err=E_BADINPUT; */
/*     fprintf(stream,"ERROR: asg_modelinit\n");	 */
/*     fprintf(stream,"half-order param = %d not in allowed range [1,7]\n",esgnp->k); */
/*     return err; */
/*   } */

  /* assign param object pointer */
  fd->fdpars 		= (void*)esgnp;
  fd->isarr 		= esg_isarr;
  fd->numsubsteps 	= esg_numsubsteps;
  fd->update 		= esg_update;
  fd->readgrid 		= esg_readgrid;
  fd->readtimegrid 	= esg_readtimegrid; 
  fd->readschemeinfo 	= esg_readschemeinfo;
  fd->set_grid_type 	= esg_set_grid_type;
  fd->build_sten_dep 	= esg_build_sten_dep;
  fd->create_sten 	= esg_create_sten;
  fd->ind2str 		= esg_ind2str;
  fd->alter_dom 	= esg_alter_dom;
  fd->readmedia		= esgn_readmedia;
  fd->parcopy 		= esgn_ts_parcopy;
  fd->fd_model_init 	= esg_modelinit;
  fd->fd_model_dest 	= esg_modeldest;


  /* choose time step */
  fd->tsf = esg_step;

#ifdef IWAVE_VERBOSE
  fprintf(stream,"\n********** esg indices *************\n");
  fprintf(stream,"D_P0=%d\n",D_P0);
  fprintf(stream,"D_MP00=%d\n",D_MP00);
  fprintf(stream,"D_MP01=%d\n",D_MP01);
  fprintf(stream,"D_MP02=%d\n",D_MP02);
  fprintf(stream,"D_V0=%d\n",D_V0);
  fprintf(stream,"D_MV0=%d\n",D_MV0);
  fprintf(stream,"D_S0=%d\n",D_S0);
  fprintf(stream,"D_MS0=%d\n",D_MS0);
  fprintf(stream,"D_E_PRIME0=%d\n",D_E_PRIME0);
  fprintf(stream,"D_E_DUAL0=%d\n",D_E_DUAL0);
  
  if (m_ndim >= 2) {
    fprintf(stream,"D_P1=%d\n",D_P1);
    fprintf(stream,"D_MP10=%d\n",D_MP10);
    fprintf(stream,"D_MP11=%d\n",D_MP11);
    fprintf(stream,"D_MP12=%d\n",D_MP12);
    fprintf(stream,"D_V1=%d\n",D_V1);
    fprintf(stream,"D_MV1=%d\n",D_MV1);
    fprintf(stream,"D_S1=%d\n",D_S1);
    fprintf(stream,"D_MS1=%d\n",D_MS1);
    fprintf(stream,"D_E_PRIME1=%d\n",D_E_PRIME1);
    fprintf(stream,"D_E_DUAL1=%d\n",D_E_DUAL1);
  }
  
  if (m_ndim >=3) {
    fprintf(stream,"D_P2=%d\n",D_P2);
    fprintf(stream,"D_MP20=%d\n",D_MP20);
    fprintf(stream,"D_MP21=%d\n",D_MP21);
    fprintf(stream,"D_MP22=%d\n",D_MP22);
    fprintf(stream,"D_V2=%d\n",D_V2);
    fprintf(stream,"D_MV2=%d\n",D_MV2);
    fprintf(stream,"D_S2=%d\n",D_S2);
    fprintf(stream,"D_MS2=%d\n",D_MS2);
    fprintf(stream,"D_E_PRIME2=%d\n",D_E_PRIME2);
    fprintf(stream,"D_E_DUAL2=%d\n",D_E_DUAL2);
  }
  
  fprintf(stream,"\n********** esg indices *************\n");  
#endif  
  return 0;
}

/*----------------------------------------------------------------------------*/
int esg_numsubsteps() {
/*----------------------------------------------------------------------------*/
  return 2; 
}

/*----------------------------------------------------------------------------*/
int esg_modeldest(IMODEL * model) {
/*----------------------------------------------------------------------------*/
//MB	comment: may have to come back to alter.  
  FD_MODEL * fdm = (FD_MODEL *)(model->specs);
  ESGN_TS_PARS * esgnp;
  RDOM *ld_pml;
  int iv;
  if (fdm) {
    if (fdm->fdpars) {
      esgnp=(ESGN_TS_PARS *)(fdm->fdpars);
      ld_pml = esgnp->ld_pml;
      for (iv = 0;iv < 2*m_ndim*m_ndim;iv ++) {
        /*< prerequest: rd_a_setnull are doing nothing, 
          otherwise need to reverse the loop order */
        rd_a_destroy(ld_pml+iv);
      }
      if (esgnp->ld_pml) userfree_(ld_pml); //MB replaced: free(ld_pml);
      userfree_(esgnp); //MB replaced: free(esgnp);
    }
    userfree_(fdm); //MB added
  }
  /* fdm is freed in im_destroy */
  return 0;
}

/*----------------------------------------------------------------------------*/
void esgn_ts_parcopy(void * tgt, const void * src) {
/*----------------------------------------------------------------------------*/
  int i;

  ESGN_TS_PARS * ptgt = (ESGN_TS_PARS *)tgt;
  const ESGN_TS_PARS * psrc = (const ESGN_TS_PARS *)src;

  ptgt->dt=psrc->dt;
  for (i=0;i<RARR_MAX_NDIM;i++) 
    (ptgt->lam)[i]=(psrc->lam)[i];
  ptgt->ndim=psrc->ndim;
  ptgt->k = psrc->k;

  /* since I put PML auxillary variables in ESGN_TS_PARS */
  fprintf(stderr, "WARNING: copy RDOM * ld_pml in esgn_ts_parcopy\n");
  ptgt->ld_pml = psrc->ld_pml;
}

//MB	possible deprication
/*----------------------------------------------------------------------------*/
int esg_isarr(int i) {
/*----------------------------------------------------------------------------*/
  if (m_ndim<4) {
    if (m_ndim>0) {
      if (i==D_P0) return 1;
      if (i==D_V0) return 1;
      if (i==D_S0) return 1;
      if (i==D_MP00) return 1;
      if (i==D_MP01) return 1;
      if (i==D_MP02) return 1;
      if (i==D_MV0) return 1;
      if (i==D_MS0) return 1;
      if (i==D_E_PRIME0) return 1;
      if (i==D_E_DUAL0) return 1;
    }
    if (m_ndim>1) {
      if (i==D_P1) return 1;
      if (i==D_V1) return 1;
      if (i==D_S1) return 1;
      if (i==D_MP10) return 1;
      if (i==D_MP11) return 1;
      if (i==D_MP12) return 1;
      if (i==D_MV1) return 1;
      if (i==D_MS1) return 1;
      if (i==D_E_PRIME1) return 1;
      if (i==D_E_DUAL1) return 1;
    }
    if (m_ndim>2) {
      if (i==D_P2) return 1;
      if (i==D_V2) return 1;
      if (i==D_S2) return 1;
      if (i==D_MP20) return 1;
      if (i==D_MP21) return 1;
      if (i==D_MP22) return 1;
      if (i==D_MV2) return 1;
      if (i==D_MS2) return 1;
      if (i==D_E_PRIME2) return 1;
      if (i==D_E_DUAL2) return 1;
    }
  }
  return 0;
}

/*----------------------------------------------------------------------------*/
int esg_update(int ia, int iv) {
/*----------------------------------------------------------------------------*/
    /* update pressure and shear stress arrays on iv = 0 */
    if (iv==0) {
	if (m_ndim>2){
		if (ia==D_P2) return 1;
		if (ia==D_S2) return 1;
		if (ia==D_S1) return 1;
	}
	if (m_ndim>1){
		if (ia==D_P1) return 1;
		if (ia==D_S0) return 1;
	}
	if (m_ndim>0 && ia==D_P0 ) return 1;
    }

    if (iv==1) {
	if (m_ndim>2 && ia==D_V2) return 1;
	if (m_ndim>1 && ia==D_V1) return 1;
	if (m_ndim>0 && ia==D_V0) return 1;
    }
    return 0;
}


/*----------------------------------------------------------------------------*/
int esg_alter_dom(int iv, IPNT gs, IPNT ge) {
/*----------------------------------------------------------------------------*/
  int idim,i;

  /* one-dimensionalize the eta arrays */
  for ( idim = 0; idim < m_ndim ; ++idim ) {
    if (iv==D_EP[idim] || iv==D_EV[idim]) {
      gs[0]=gs[idim];
      ge[0]=ge[idim];
      /*      fprintf(stderr,"iv=%d idim=%d\n",iv,idim); */
      for ( i = 1; i < m_ndim; ++i ) 
	gs[i] = ge[i] = 0;
    }
  }
  return 0;
}

//MB	comment: the following is an auxilary function
/*----------------------------------------------------------------------------*/
void esg_fprint_grid_type(FILE * stream, int ndim, IPNT gtype[RDOM_MAX_NARR]) {
/*----------------------------------------------------------------------------*/
  int iv, i;
  fprintf(stream, "============= FD grid type =============\n");
  for (iv = 0;iv < m_size;iv ++) {
    fprintf(stream, "%s:\n", m_names[iv]);
    for (i = 0;i < ndim;i ++) {
      if (gtype[iv][i] == 0) 
        fprintf(stream, "using PRIMAL grid on axis-%d\n", i);
      else 
        fprintf(stream, "using DUAL   grid on axis-%d\n", i);
    }
  }
}

/*----------------------------------------------------------------------------*/
int esg_set_grid_type(FILE *stream, int ndim, IPNT gtype[RDOM_MAX_NARR] ) {
/*----------------------------------------------------------------------------*/
  int iv, i, nss; /*< nss: # shear stress component(s) */
  if ( ndim < 1 || ndim > RARR_MAX_NDIM ) return E_BADINPUT;
  /*< px, py, pz, mp00, mp01 on the primal grid */
  for (iv = 0;iv < RDOM_MAX_NARR;iv ++)  IASN(gtype[iv], IPNT_0);
  
  for (i = 0;i < ndim;i ++) {
    gtype[D_V[i]][i]  = 1;  
    gtype[D_MV[i]][i] = 1;
    gtype[D_EV[i]][i] = 1;
  }
  nss = ndim*(ndim-1)/2;
  for (i = 0;i < nss;i ++) {
    gtype[D_S[i]][i]  = 1;  gtype[D_S[i]][(i+1)%3]  = 1;
    gtype[D_MS[i]][i] = 1;  gtype[D_MS[i]][(i+1)%3] = 1;
  }
#ifdef IWAVE_VERBOSE
  esg_fprint_grid_type(stream, ndim, gtype);
#endif
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

/*----------------------------------------------------------------------------*/
int esg_build_sten_dep(FILE *stream, int ndim, int sten_dep_mat[RDOM_MAX_NARR][RDOM_MAX_NARR]) {
/*----------------------------------------------------------------------------*/
  int nss;
  int i, j;

  if (ndim < 1 || ndim > RARR_MAX_NDIM) return E_BADINPUT;
  nss = ndim*(ndim-1)/2;
  
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
    sten_dep_mat[D_P[i]][D_MP00] = DEP_F;
    sten_dep_mat[D_P[i]][D_MP01] = DEP_F;
  }
  
  
  /*< participating arrays for shear stresses */
 
  for (i = 0;i < nss;i ++) {
    sten_dep_mat[D_S[i]][D_MS[i]] = DEP_F;
    sten_dep_mat[D_S[i]][D_V[i]] = (i+1)%3+1;
    sten_dep_mat[D_S[i]][D_V[(i+1)%3]] = i+1;
  }
  
  /*< participating arrays for velocity */
  for (i = 0;i < ndim;i ++) {
    sten_dep_mat[D_V[i]][D_MV[i]] = DEP_F;
    sten_dep_mat[D_V[i]][D_P[i]] = i+1;
  }
  
  if (ndim >= 2) {
    sten_dep_mat[D_V[0]][D_S[0]] = DEP_DFDX;
    sten_dep_mat[D_V[1]][D_S[0]] = DEP_DFDZ;
  }
  if (ndim >= 3) {
    sten_dep_mat[D_V[0]][D_S[2]] = DEP_DFDY;
    sten_dep_mat[D_V[2]][D_S[2]] = DEP_DFDZ;

    sten_dep_mat[D_V[1]][D_S[1]] = DEP_DFDY;
    sten_dep_mat[D_V[2]][D_S[1]] = DEP_DFDX;
  }
  
/* #ifdef IWAVE_VERBOSE */
/*   /\*< try to print the equations *\/ */
/*   fprintf(stream, "\nWAVE EQUATIONS:\n"); */
/*   for (i = 0;i < m_size;i ++) { */
/*     int ir = i; */
/* //MB	deleted: if (!esg_isdyn(ir)) continue; */
/*     fprintf(stream, "d %4s/ dt    =    ", m_names[i]); */
/*     for (j = 0;j < m_size;j ++) { */
/*       int ip = j; */
/*       int op = fdm->sten_dep_mat[ir][ip]; */
/*       switch (op) { */
/*       case DEP_F: */
/*         fprintf(stream, "%10s    ", m_names[j]); */
/*         break; */
/*       case DEP_DFDZ: */
/*         fprintf(stream, "d %4s/ dz    ", m_names[j]); */
/*         break; */
/*       case DEP_DFDX: */
/*         fprintf(stream, "d %4s/ dx    ", m_names[j]); */
/*         break; */
/*       case DEP_DFDY: */
/*         fprintf(stream, "d %4s/ dy    ", m_names[j]); */
/*         break; */
/*       default: */
/*         break; */
/*       } */
/*     } */
/*     fprintf(stream,"\n"); */
/*   }   */
/* #endif */

  return 0;
}
  

/*----------------------------------------------------------------------------*/
int esgn_modelinfo(FILE *stream, IMODEL *model){
/*----------------------------------------------------------------------------*/
    /* FD_MODEL *fdm; */
  int ndim, i, iv, ia;
  IPNT n;
  int nss;
  RDOM *dc; /* computational */

  /* fdm = (FD_MODEL*)(model->specs); */
  ndim = model->g.dim;
  nss = ndim*(ndim-1)/2;

  fprintf(stream, "ESGN model %dD.\n", ndim);

  dc = &(model->ld_c); /* get computational domain */

  for (iv = 0;iv < ndim; ++iv) {
    ia = D_P[iv];
    rd_size(dc, ia, n);
    fprintf(stream, "  P%d size = [%d", iv, n[0]);
    for ( i = 1; i < ndim; ++i ) fprintf(stream, " %d", n[i]);
    fprintf(stream, "]\n");
  }
    
  for ( iv = 0; iv < ndim; ++iv ) {
    ia = D_V[iv];
    rd_size(dc, ia, n);
    fprintf(stream, "  V%d size = [%d", iv, n[0]);
    for ( i = 1; i < ndim; ++i ) fprintf(stream, " %d", n[i]);
    fprintf(stream, "]\n");
  }

  for (iv = 0;iv < nss; ++iv) {
    ia = D_S[iv];
    rd_size(dc, ia, n);
    fprintf(stream, " S%d size = [%d", iv, n[0]);
    for (i = 1;i < ndim; ++i) fprintf(stream, " %d", n[i]);
    fprintf(stream, "]\n");
  }

  return 0;
}

/*----------------------------------------------------------------------------*/
//MB	modified
int esg_create_sten(FD_MODEL * fdm,
		    FILE * stream, 
		    int ndim, 
		    IPNT gtype[RDOM_MAX_NARR], 
		    int sten_dep_mat[RDOM_MAX_NARR][RDOM_MAX_NARR], 
		    STENCIL * sten) {
/*----------------------------------------------------------------------------*/
#ifdef IWAVE_VERBOSE
	fprintf(stream,">>> Inside esg_create_sten\n");
	fflush(stream);
#endif

  ESGN_TS_PARS * esgnp = (ESGN_TS_PARS *)(fdm->fdpars);
  return create_sten2_2k(fdm,stream, esgnp->k, ndim, gtype, sten_dep_mat,sten);
}  		

/*----------------------------------------------------------------------------*/
int esg_readschemeinfo(PARARRAY * par, FILE * stream, IMODEL * model) {
/*----------------------------------------------------------------------------*/
/*----------------------------------------------------------------------------*/
#ifdef IWAVE_VERBOSE
  fprintf(stream,">>> Inside esg_readschemeinfo\n");
  fflush(stream);
#endif

  RPNT dxs;  /* grid stepsn */
  int idim;  /* counter */
  int ndim;
  RDOM *ld_pml;
  IPNT dgsa[RDOM_MAX_NARR], dgea[RDOM_MAX_NARR], gsc, gec, ns;
  /* IPNT gtype[RDOM_MAX_NARR]; */
  int iv, i;
  FD_MODEL * fdm = (FD_MODEL *)(model->specs); 
  ESGN_TS_PARS * esgnp = (ESGN_TS_PARS *)(fdm->fdpars);
  
  int err;

  get_d(dxs, model->g);
  get_n(ns, model->g);
  /* set model dimn par */
  esgnp->ndim = ndim = (model->g).dim;
  for (idim = 0;idim < esgnp->ndim;idim ++) {
    if (dxs[idim] <= 0.0) {
      fprintf(stream, "Error: bad input: wrong grid space step, dim=%d, step=%g\n",
              idim, dxs[idim]);
      return E_BADINPUT;
    }
    esgnp->lam[idim] = (model->tsind).dt / dxs[idim];
#ifdef IWAVE_VERBOSE
    fprintf(stderr, "lam[%d] = %g\n", idim, esgnp->lam[idim]);
#endif 
  }
  esgnp->dt = (model->tsind).dt;
  
  if (ndim < 2 || ndim > 3) {
    fprintf(stderr, "Error: bad input in esg_build_scheme: ndim = %d\n",
            ndim);
    fprintf(stderr, "       esg only provides 2d/3d elastic wave solvers\n");
    return E_BADINPUT;
  }
  
  /*< initialize ld_pml here */
  /* esg_set_grid_type(stderr, ndim, gtype); */
  //MB repalced: esgnp->ld_pml = (RDOM *)malloc(2*ndim*ndim*sizeof(RDOM));
  esgnp->ld_pml = (RDOM *)usermalloc_(2*ndim*ndim*sizeof(RDOM));
  
  ld_pml = esgnp->ld_pml;
  
  if (ndim == 2) {
    /******* PML region I *******************************************/
    for ( iv = 0;iv < RDOM_MAX_NARR; ++iv ) {
      IASN(dgsa[iv], IPNT_1);
      IASN(dgea[iv], IPNT_0);
//MB	replaced: if (!esg_isdyn(iv))  continue;
      if (!isdyn(fdm,iv)) continue;
      rd_gse(&(model->ld_c), iv, gsc, gec);
      if (gsc[1] >= 0)  continue;
      for (idim = 0;idim < ndim;idim ++) {
        dgsa[iv][idim] = gsc[idim];
        dgea[iv][idim] = gec[idim];
      }
      dgea[iv][1] = iwave_min(gec[1], -1);
    }
    /** tilde U_z in region I */
    err = rd_a_create(ld_pml++, RDOM_MAX_NARR, ndim, dgsa, dgea);
    /** tilde U_x in region I */
    err = rd_a_create(ld_pml++, RDOM_MAX_NARR, ndim, dgsa, dgea);
    /***************************************************************/

    /******* PML region II *****************************************/
    for ( iv = 0;iv < RDOM_MAX_NARR; ++iv ) {
      IASN(dgsa[iv], IPNT_1);
      IASN(dgea[iv], IPNT_0);
//MB	replaced: if (!esg_isdyn(iv))  continue;
      if (!isdyn(fdm,iv)) continue;
      rd_gse(&(model->ld_c), iv, gsc, gec);
      if (gec[1] <= ns[1]-1)  continue;
      for (idim = 0;idim < ndim;idim ++) {
        dgsa[iv][idim] = gsc[idim];
        dgea[iv][idim] = gec[idim];
      }
      dgsa[iv][1] = iwave_max(gsc[1], ns[1]-1);
    }
    /** tilde U_z in region II */
    err = rd_a_create(ld_pml++, RDOM_MAX_NARR, ndim, dgsa, dgea);
    /** tilde U_x in region II */
    err = rd_a_create(ld_pml++, RDOM_MAX_NARR, ndim, dgsa, dgea);
    /***************************************************************/

    /******* PML region III ****************************************/
    for ( iv = 0;iv < RDOM_MAX_NARR; ++iv ) {
      IASN(dgsa[iv], IPNT_1);
      IASN(dgea[iv], IPNT_0);
//MB	replaced: if (!esg_isdyn(iv))  continue;
      if (!isdyn(fdm,iv)) continue;
      rd_gse(&(model->ld_c), iv, gsc, gec);
      if (gsc[0] >= 0)  continue;
      for (idim = 0;idim < ndim;idim ++) {
        dgsa[iv][idim] = gsc[idim];
        dgea[iv][idim] = gec[idim];
      }
      dgea[iv][0] = iwave_min(gec[0], -1);  
    
    }
    /** tilde U_z in region III */
    err = rd_a_create(ld_pml++, RDOM_MAX_NARR, ndim, dgsa, dgea);
    /** tilde U_x in region III */
    err = rd_a_create(ld_pml++, RDOM_MAX_NARR, ndim, dgsa, dgea);   
    /***************************************************************/

    /******* PML region IV *****************************************/
    for ( iv = 0;iv < RDOM_MAX_NARR; ++iv ) {
      IASN(dgsa[iv], IPNT_1);
      IASN(dgea[iv], IPNT_0);
//MB	replaced: if (!esg_isdyn(iv))  continue;
      if (!isdyn(fdm,iv)) continue;
      rd_gse(&(model->ld_c), iv, gsc, gec);
      if (gec[0] <= ns[0]-1)  continue;
      for (idim = 0;idim < ndim;idim ++) {
        dgsa[iv][idim] = gsc[idim];
        dgea[iv][idim] = gec[idim];
      }
      dgsa[iv][0] = iwave_max(gsc[0], ns[0]-1);
    }
    /** tilde U_z in region IV */
    err = rd_a_create(ld_pml++, RDOM_MAX_NARR, ndim, dgsa, dgea);
    /** tilde U_x in region IV */
    err = rd_a_create(ld_pml++, RDOM_MAX_NARR, ndim, dgsa, dgea);   
    /***************************************************************/
  }
  
  if (ndim == 3) { 
    /******* PML region I *******************************************/
    for ( iv = 0;iv < RDOM_MAX_NARR; ++iv ) {
      IASN(dgsa[iv], IPNT_1);
      IASN(dgea[iv], IPNT_0);
//MB	replaced: if (!esg_isdyn(iv))  continue;
      if (!isdyn(fdm,iv)) continue;
      rd_gse(&(model->ld_c), iv, gsc, gec);
      if (gsc[2] >= 0)  continue;
      for (idim = 0;idim < ndim;idim ++) {
        dgsa[iv][idim] = gsc[idim];
        dgea[iv][idim] = gec[idim];
      }
      dgea[iv][2] = iwave_min(gec[2], -1);
    }
    /** tilde U_z in region I */
    err = rd_a_create(ld_pml++, RDOM_MAX_NARR, ndim, dgsa, dgea);
    /** tilde U_x in region I */
    err = rd_a_create(ld_pml++, RDOM_MAX_NARR, ndim, dgsa, dgea);
    /** tilde U_y in region I */
    err = rd_a_create(ld_pml++, RDOM_MAX_NARR, ndim, dgsa, dgea);
    /***************************************************************/

    /******* PML region II *****************************************/
    for ( iv = 0;iv < RDOM_MAX_NARR; ++iv ) {
      IASN(dgsa[iv], IPNT_1);
      IASN(dgea[iv], IPNT_0);
//MB	replaced: if (!esg_isdyn(iv))  continue;
      if (!isdyn(fdm,iv)) continue;
      rd_gse(&(model->ld_c), iv, gsc, gec);
      if (gec[2] <= ns[2]-1)  continue;
      for (idim = 0;idim < ndim;idim ++) {
        dgsa[iv][idim] = gsc[idim];
        dgea[iv][idim] = gec[idim];
      }
      dgsa[iv][2] = iwave_max(gsc[2], ns[2]-1);
    }
    /** tilde U_z in region II */
    err = rd_a_create(ld_pml++, RDOM_MAX_NARR, ndim, dgsa, dgea);
    /** tilde U_x in region II */
    err = rd_a_create(ld_pml++, RDOM_MAX_NARR, ndim, dgsa, dgea);
    /** tilde U_y in region II */
    err = rd_a_create(ld_pml++, RDOM_MAX_NARR, ndim, dgsa, dgea);
    /***************************************************************/

    /******* PML region III ****************************************/
    for ( iv = 0;iv < RDOM_MAX_NARR; ++iv ) {
      IASN(dgsa[iv], IPNT_1);
      IASN(dgea[iv], IPNT_0);
//MB	replaced: if (!esg_isdyn(iv))  continue;
      if (!isdyn(fdm,iv)) continue;
      rd_gse(&(model->ld_c), iv, gsc, gec);
      if (gsc[1] >= 0)  continue;
      for (idim = 0;idim < ndim;idim ++) {
        dgsa[iv][idim] = gsc[idim];
        dgea[iv][idim] = gec[idim];
      }
      dgea[iv][1] = iwave_min(gec[1], -1);  
    
    }
    /** tilde U_z in region III */
    err = rd_a_create(ld_pml++, RDOM_MAX_NARR, ndim, dgsa, dgea);
    /** tilde U_x in region III */
    err = rd_a_create(ld_pml++, RDOM_MAX_NARR, ndim, dgsa, dgea);   
    /** tilde U_y in region III */
    err = rd_a_create(ld_pml++, RDOM_MAX_NARR, ndim, dgsa, dgea);
    /***************************************************************/

    /******* PML region IV *****************************************/
    for ( iv = 0;iv < RDOM_MAX_NARR; ++iv ) {
      IASN(dgsa[iv], IPNT_1);
      IASN(dgea[iv], IPNT_0);
//MB	replaced: if (!esg_isdyn(iv))  continue;
      if (!isdyn(fdm,iv)) continue;
      rd_gse(&(model->ld_c), iv, gsc, gec);
      if (gec[1] <= ns[1]-1)  continue;
      for (idim = 0;idim < ndim;idim ++) {
        dgsa[iv][idim] = gsc[idim];
        dgea[iv][idim] = gec[idim];
      }
      dgsa[iv][1] = iwave_max(gsc[1], ns[1]-1);
    }
    /** tilde U_z in region IV */
    err = rd_a_create(ld_pml++, RDOM_MAX_NARR, ndim, dgsa, dgea);
    /** tilde U_x in region IV */
    err = rd_a_create(ld_pml++, RDOM_MAX_NARR, ndim, dgsa, dgea);
    /** tilde U_y in region IV */
    err = rd_a_create(ld_pml++, RDOM_MAX_NARR, ndim, dgsa, dgea);
    /***************************************************************/
   
    /******* PML region V ******************************************/
    for ( iv = 0;iv < RDOM_MAX_NARR; ++iv ) {
      IASN(dgsa[iv], IPNT_1);
      IASN(dgea[iv], IPNT_0);
//MB	replaced: if (!esg_isdyn(iv))  continue;
      if (!isdyn(fdm,iv)) continue;
      rd_gse(&(model->ld_c), iv, gsc, gec);
      if (gsc[0] >= 0)  continue;
      for (idim = 0;idim < ndim;idim ++) {
        dgsa[iv][idim] = gsc[idim];
        dgea[iv][idim] = gec[idim];
      }
      dgea[iv][0] = iwave_min(gec[0], -1);
    }
    /** tilde U_z in region V */
    err = rd_a_create(ld_pml++, RDOM_MAX_NARR, ndim, dgsa, dgea);
    /** tilde U_x in region V */
    err = rd_a_create(ld_pml++, RDOM_MAX_NARR, ndim, dgsa, dgea);   
    /** tilde U_y in region V */
    err = rd_a_create(ld_pml++, RDOM_MAX_NARR, ndim, dgsa, dgea);
    /***************************************************************/

    /******* PML region VI *****************************************/
    for ( iv = 0;iv < RDOM_MAX_NARR; ++iv ) {
      IASN(dgsa[iv], IPNT_1);
      IASN(dgea[iv], IPNT_0);
//MB	replaced: if (!esg_isdyn(iv))  continue;
      if (!isdyn(fdm,iv)) continue;
      rd_gse(&(model->ld_c), iv, gsc, gec);
      if (gec[0] <= ns[0]-1)  continue;
      for (idim = 0;idim < ndim;idim ++) {
        dgsa[iv][idim] = gsc[idim];
        dgea[iv][idim] = gec[idim];
      }
      dgsa[iv][0] = iwave_max(gsc[0], ns[0]-1);
    }
    /** tilde U_z in region VI */
    err = rd_a_create(ld_pml++, RDOM_MAX_NARR, ndim, dgsa, dgea);
    /** tilde U_x in region VI */
    err = rd_a_create(ld_pml++, RDOM_MAX_NARR, ndim, dgsa, dgea);
    /** tilde U_y in region VI */
    err = rd_a_create(ld_pml++, RDOM_MAX_NARR, ndim, dgsa, dgea);
    /***************************************************************/

  }

//MB	comment: zero out RARRs in ld_pml domain
  for (i = 0;i < 2*ndim*ndim;i ++) {
    for (iv = 0;iv < RDOM_MAX_NARR;iv ++) {
      ra_zero(&(esgnp->ld_pml[i]._s[iv]));
    }
  }
  return 0;
}

//MB	comment: auxilary function
/*----------------------------------------------------------------------------*/
int esg_dynamic_init(IMODEL * model) {
/*----------------------------------------------------------------------------*/
	int iv, i;
	FD_MODEL * fdm = (FD_MODEL *)(model->specs); 
  	ESGN_TS_PARS * esgnp = (ESGN_TS_PARS *)(fdm->fdpars);
  	int ndim = (model->g).dim;  
  	for (i = 0;i < 2*ndim*ndim;i ++) {
    		for (iv = 0;iv < RDOM_MAX_NARR;iv ++) {
      			ra_zero(&(esgnp->ld_pml[i]._s[iv]));
    		}
  	}
  	return 0;
}

//MB	comment: auxilary function
/*----------------------------------------------------------------------------*/
void esg_printpml(PARARRAY *pars, FILE *stream, IMODEL *model){
/*----------------------------------------------------------------------------*/
	FD_MODEL * fdm = (FD_MODEL *)(model->specs);
    	ESGN_TS_PARS * esgnp = (ESGN_TS_PARS *)(fdm->fdpars);

    	fprintf(stream, "PML RDOM dump\n");
	int tot = (model->g.dim) * (model->g.dim) * 2;
	int i;
	for (i=0; i<tot ;i++){
		fprintf(stream, "dumping ld_pml[%d]\n",i);
    		rd_a_dump( &(esgnp->ld_pml[i]), stream );
	}

    	fflush(stream);
}
