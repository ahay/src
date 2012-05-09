#include "sgn.h"

//#define NEW
//#define ASG

/* 
Functions to convert between array names and indices.
*/

/*------------------- private data declarations ----------------------------*/

static int m_ndim = 0; /* dimension - need only compute once */
static int m_size = 17; /* 16; */
static const char* m_names[] = {"p0", "kappa", "v0","rho0","v1","rho1","v2","rho2",
				"eta_p0", "eta_v0", "p1", "eta_p1", "eta_v1",
				"p2", "eta_p2", "eta_v2", "dyn_rho","NONE","NONE","NONE","NONE"};

/*----------------------------------------------------------------------------*/

/*------------------- private function declarations --------------------------*/

int asg_isdyn(int i);
int asg_isarr(int i);

/*----------------------------------------------------------------------------*/

//int asg_getnarr() { return m_size; }
int asg_getnarr() { return RDOM_MAX_NARR; }

/*
const char * asg_getnames(int i) { 
  if (i>-1 && i<m_size) return m_names[i];
  return m_names[m_size];
}
*/

/*int asg_getindices(int i) {
  if (i>-1 && i<m_size) return i;
  return m_size;
}
*/

const char* asg_ind2str(int ind) {
  if (ind>-1 && ind<m_size) return m_names[ind];
  return m_names[m_size];
}
/*----------------------------------------------------------------------------*/

//int asg_getndim() { return m_ndim; }

/*----------------------------------------------------------------------------*/
int asg_modelinit(PARARRAY *pars, 
		  FILE *stream, 
		  IMODEL * im) {

  int err=0;           /* return value */
  char * jnk;          /* workspace */
  SGN_TS_PARS *sgnm;   /* model pars */
  FD_MODEL * fd;       /* function struct */

  IPNT cdims;          /* workspace for cartesian grid dim info */
  IPNT crank;          /* workspace for cartesian grid rank */
#ifdef IWAVE_USE_MPI
  IPNT cpers;          /* workspace for periodic wrap info  - currently not used */
#endif

  int i;               /* counter */

  fd=(FD_MODEL *)malloc(sizeof(FD_MODEL));
  if (fd==NULL) return E_ALLOC;
  /* first the "base class" behaviour */
  //  fprintf(stream,"top of asg_modelinit: im->specs=%x\n",im->specs);
  im->specs=(void *)fd;
  //  fprintf(stream,"top of asg_modelinit: im->specs=%x after assignment\n",im->specs);

  /* allocate sgn model ----------------------------------------------------*/   
  sgnm = (SGN_TS_PARS*)malloc(sizeof(SGN_TS_PARS));
  if ( sgnm == NULL ) { 
    err=E_ALLOC;
    fprintf(stream,"ERROR: asg_modelinit\n");
    fprintf(stream,"failed to allocate SGN_TS_PARS object\n");
    return err;
  }
	    
  // assign heritage pars
  sgnm->psingle = 0; /* TODO: single pressures in Step 3 */
  IASN(sgnm->eflags,IPNT_1); /* this is related - with code from old sg added, could try using this device */

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

#ifdef IWAVE_USE_MPI
  fprintf(stream,"boundary flags: rk=%d\n",rk);
  for (i=0;i<m_ndim;i++) 
    fprintf(stream,"    lbc[%d]=%d rbc[%d]=%d\n",i,sgnm->lbc[i],i,sgnm->rbc[i]);
  
#endif

  /* decode order - take care of deprecated cases */
  sgnm->k=1;

  /* main branch */
  if (ps_ffint(*pars,"order",&(sgnm->k))) {
    /* old-style iwave specs - deprecated */
    if (err=ps_ffcstring(*pars,"scheme_phys",&jnk)) {
      fprintf(stream,"ERROR: asg_modelinit from ps_ffstring\n");
      fprintf(stream,"failed to read deprecated order param scheme_phys\n");
      return err;
    }
    else {
      if (!strcmp(jnk,"22")) sgnm->k=1;
      else if (!strcmp(jnk,"24")) sgnm->k=2;
      else if (!strcmp(jnk,"210")) sgnm->k=5;
      else if (!strcmp(jnk,"2k")) {
        if (err=ps_ffint(*pars,"k_phys",&(sgnm->k))) {
	  fprintf(stream,"ERROR: asg_modelinit from ps_ffint\n");
	  fprintf(stream,"failed to read deprecated order param k_phys\n");
	  return err;
	}
      }
      else {
	err=E_BADINPUT;
	fprintf(stream,"ERROR: asg_modelinit\n");	
	fprintf(stream,"invalid value for deprecated param scheme_phys\n");
	fprintf(stream,"legal values are 22, 24, 210 (3D only), and 2k\n");
	return err;
      }
    }
  }
  /* sanity-check order */
  if (sgnm->k<1 || sgnm->k>8) {
    err=E_BADINPUT;
    fprintf(stream,"ERROR: asg_modelinit\n");	
    fprintf(stream,"half-order param = %d not in allowed range [1,7]\n",sgnm->k);
    return err;
  }
    
  fd->fdpars = (void*)sgnm;

  //  fd->getndim  = asg_getndim;
  //  fd->getnarr = asg_getnarr;
  fd->isarr  = asg_isarr;
  //fd->getnames = asg_getnames;
  //fd->getindices = asg_getindices;
  fd->isdyn  = asg_isdyn;
  fd->readgrid = asg_readgrid;
  fd->readtimegrid = asg_readtimegrid; 
  fd->readschemeinfo = asg_readschemeinfo;
  fd->set_grid_type = asg_set_grid_type;
  fd->build_sten_dep = asg_build_sten_dep;
  fd->create_sten = asg_create_sten;
  fd->ind2str = asg_ind2str;
  fd->assign_action_array = asg_assign_action_array;
  fd->alter_dom = asg_alter_dom;
  fd->readmedia = asg_readmedia;
  // begin migration to new time-step model
#ifdef NEW
  if ((sgnm->k == 2) && (m_ndim == 2)) fd->ptsf = asg_noop;
  else fd->ptsf = asg_modelpostts;
#else
  fd->ptsf = asg_modelpostts;
#endif
  fd->parcopy = sgn_ts_parcopy;
  fd->fd_model_init = asg_modelinit;
  fd->fd_model_dest = asg_modeldest;

  // choose time step
  if (sgnm->k==1) {
    if (m_ndim==1) fd->tsf=sgn_ts1d_22;
    if (m_ndim==2) fd->tsf=sgn_ts2d_22;
    if (m_ndim==3) fd->tsf=sgn_ts3d_22;
  }
  else if (sgnm->k==2) {
    if (m_ndim==1) fd->tsf=sgn_ts1d_24;
#ifdef NEW
#ifdef ASG
    if (m_ndim==2) fd->tsf=asg_fts2d_24;
#else
    if (m_ndim==2) fd->tsf=duhasg24_2d;
#endif
#else
    if (m_ndim==2) fd->tsf=sgn_gts2d_24;
#endif
    if (m_ndim==3) fd->tsf=sgn_ts3d_24;
  }
  else {
    if (sgnm->k<8) {
      if (m_ndim==1) {
	if (sgnm->k==5) fd->tsf=sgn_ts1d_210;
	else {
	  fprintf(stream,"ERROR: asg_modelinit\n");
	  fprintf(stream,"order = %d not implemented for 1D\n",2*sgnm->k);
	  return E_BADINPUT;
	}
      }
      if (m_ndim==2) fd->tsf=sgn_ts2d_2K;
      if (m_ndim==3) {
	if (sgnm->k==5) fd->tsf=sgn_ts3d_210;
	else fd->tsf=sgn_ts3d_2K;
      }
    }
    else {
      fprintf(stream,"ERROR: asg_modelinit\n");
      fprintf(stream,"order = %d not implemented\n",2*sgnm->k);
      return E_BADINPUT;
    }
  }


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
  fprintf(stream,"D_MV_DYN=%d\n",D_MV_DYN);
  fprintf(stream,"\n********** asg indices *************\n");  
  */

  return 0;
}

int asg_modeldest(IMODEL * model) {
  FD_MODEL * fdm = (FD_MODEL *)(model->specs);
  SGN_TS_PARS * sgnm;
  if (fdm) {
    if (fdm->fdpars) {
      sgnm=(SGN_TS_PARS *)(fdm->fdpars);
      free(sgnm);
    }
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
  ptgt->psingle=psrc->psingle;

  IASN(ptgt->eflags,psrc->eflags);
  IASN(ptgt->lbc,psrc->lbc);
  IASN(ptgt->rbc,psrc->rbc);   
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
    if (i==D_MV_DYN) return 1;
  }
  return 0;
}

int asg_isdyn(int i) {

  if (m_ndim<4) {
    if (m_ndim>0) {
      if (i==D_P0) return 1;
      if (i==D_V0) return 1;
    }
    if (m_ndim>1) {
      if (i==D_P1) return 1;
      if (i==D_V1) return 1;
    }
    if (m_ndim>2) {
      if (i==D_P2) return 1;
      if (i==D_V2) return 1;
    }
    if (i==D_MV_DYN) return 1;
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

  /*< participating arrays for dynamic velocity-multiplier */
  for (i = 0;i < ndim;i ++) {
    sten_dep_mat[D_MV_DYN][D_MV_DYN] += (i+1);
    if (i < ndim -1) sten_dep_mat[D_MV_DYN][D_MV_DYN] *= 10;
  }
  return 0;
}

/*----------------------------------------------------------------------------*/

/* this no-op is part of a migration to to new (02.12) time step structure */

int asg_noop(int iarr, IMODEL * model) { return 0; }

/* NOTE: this construction ASSUMES that pars is a pointer to IMODEL
   whose domains include dom (first arg) - in other words, dom is
   passed through BOTH args. This will lead to disaster, money back
   guaranteed!!!

   WWS 13.02.11: eliminated first arg & ambiguity
*/
    
//int asg_modelpostts(RDOM *dom, int iarr, IMODEL * model) {
int asg_modelpostts(int iarr, IMODEL * model) {

  RDOM *dom; /* allocated */
  RDOM *dc; /* computational */
  IPNT gsa, gea, gsc, gec, ip, ns;
  int ndim, shift, i;
  ireal sgn, v;
  
  get_n(ns, model->g);
  
  dom = &(model->ld_a); /* get allocated domain */
  dc = &(model->ld_c); /* get computational domain */
    
  rd_gse(dom, iarr, gsa, gea);    
  rd_gse(dc, iarr, gsc, gec);

  /* WWS 22.02.08 ndim = model->ndim;*/
  ndim = model->g.dim;

  /* we mirror only outside phys+npml */
  for ( i = 0; i < ndim; ++i )
    {
      shift = 1;
      if ( iarr == D_V[i] ) shift = 0;
      if ( gsc[i] != -model->nls[i] + shift ) gsa[i] = gsc[i];
        
      shift = -2;
      if ( gec[i] != ns[i] + model->nrs[i] + shift ) gea[i] = gec[i];
    }

  if ( ndim == 1 )
    {
      if      ( iarr == D_P0 ) 
	{
	  sgn = -1.0; /* odd pressure mirroring */
	  shift = -2;
	}
      else if ( iarr == D_V0 )
	{
	  sgn = 1.0; /* even velocity mirroring */
	  shift = -1;
	}
      else return E_INTERNAL;
      for ( i = gsa[0]; i <= gsc[0] - 1; ++i ) /* left */
	{
	  ip[0] = 2 * gsc[0] - i + shift;
	  if (ip[0]==i && iarr==D_P0) {
	    rd_gset(dom,iarr,ip,REAL_ZERO);
	  }
	  else {
	    v = rd_gget(dc, iarr, ip);
	    ip[0] = i;
	    rd_gset(dom, iarr, ip, sgn * v);
	  }
	}
      for ( i = gea[0]; i >= gec[0] + 1; --i ) /* right */
	{
	  ip[0] = 2 * gec[0] - i - shift;
	  if (ip[0]==i && iarr==D_P0) {
	    rd_gset(dom,iarr,ip,REAL_ZERO);
	  }
	  else {
	    v = rd_gget(dc, iarr, ip);
	    ip[0] = i;
	    rd_gset(dom, iarr, ip, sgn * v);
	  }
	}
        
      return 0;
    }
    
  /* 2D---------------------------------------------------------------------*/
  if ( ndim == 2 ) {

    if  ( (iarr == D_P0) || (iarr == D_P1)  || (iarr == D_V1) ) 
      {
	shift = -2;
	sgn = -1.0;
      }
    else if ( iarr == D_V0 ) 
      {
	shift = -1;
	sgn = 1.0;
      }

    else return E_INTERNAL;

    for ( ip[1] = gsc[1]; ip[1] <= gec[1]; ++ip[1] ) 
      {
	for ( i = gsa[0]; i <= gsc[0] - 1; ++i ) 
	  {
	    ip[0] = 2 * gsc[0] - i + shift;
	    if (ip[0]==i && (iarr==D_P0 || iarr==D_P1)) 
	      {
		rd_gset(dom,iarr,ip,REAL_ZERO);
	      }
	    else 
	      {
		v = rd_gget(dc, iarr, ip);
		ip[0] = i;
		rd_gset(dom, iarr, ip, sgn * v);
	      }
	  }
				
	for ( i = gea[0]; i >= gec[0] + 1; --i ) 
	  {
	    ip[0] = 2 * gec[0] - i - shift;
	    if (ip[0]==i && (iarr==D_P0 || iarr==D_P1)) 
	      {
		rd_gset(dom,iarr,ip,REAL_ZERO);
	      }   
	    else 
	      {
		v = rd_gget(dc, iarr, ip);
		ip[0] = i;
		rd_gset(dom, iarr, ip, sgn * v);
	      }
	  }
      }
        
    if  ( (iarr == D_P0) || (iarr == D_P1)  || (iarr == D_V0) ) 
      {
	shift = -2;
	sgn = -1.0;
      }
    else if ( iarr == D_V1 ) 
      {
	shift = -1;
	sgn = 1.0;
      }
    else return E_INTERNAL;
			

    for ( ip[0] = gsa[0]; ip[0] <= gea[0]; ++ip[0] ) 
      {
	for ( i = gsa[1]; i <= gsc[1] - 1; ++i ) 
	  {
	    ip[1] = 2 * gsc[1] - i + shift;
	    if (ip[1]==i && (iarr==D_P0 || iarr==D_P1)) 
	      {
		rd_gset(dom,iarr,ip,REAL_ZERO); 
	      }   
	    else 
	      {
		v = rd_gget(dc, iarr, ip);
		ip[1] = i;
		rd_gset(dom, iarr, ip, sgn * v);
	      }
	  }
				
	for ( i = gea[1]; i >= gec[1] + 1; --i ) 
	  {
	    ip[1] = 2 * gec[1] - i - shift;
	    if (ip[1]==i && (iarr==D_P0 || iarr==D_P1)) 
	      {
		rd_gset(dom,iarr,ip,REAL_ZERO);
	      }   
	    else 
	      {
		v = rd_gget(dc, iarr, ip);
		ip[1] = i;
		rd_gset(dom, iarr, ip, sgn * v);
	      }
	  }
      }
    return 0;
  }

        
  /* 3D---------------------------------------------------------------------*/
  if ( ndim == 3 )
    {
      if      ( (iarr == D_P0) || (iarr == D_P1) || (iarr == D_P2) ||
		(iarr == D_V1) || (iarr == D_V2) ) 
	{
	  sgn = -1.0;
	  shift = -2;
	}
      else if ( iarr == D_V0 )
	{
	  shift = -1;
	  sgn = 1.0;
	}
      else return E_INTERNAL;
      for ( ip[2] = gsc[2]; ip[2] <= gec[2]; ++ip[2] )
	for ( ip[1] = gsc[1]; ip[1] <= gec[1]; ++ip[1] )
	  {
	    for ( i = gsa[0]; i <= gsc[0] - 1; ++i )
	      {
		ip[0] = 2 * gsc[0] - i + shift;
		if (ip[0]==i && (iarr==D_P0 || iarr==D_P1 || iarr==D_P2)) {
		  rd_gset(dom,iarr,ip,REAL_ZERO);
		}
		else {
		  v = rd_gget(dc, iarr, ip);
		  ip[0] = i;
		  rd_gset(dom, iarr, ip, sgn * v);
		}
	      }
	    for ( i = gea[0]; i >= gec[0] + 1; --i )
	      {
		ip[0] = 2 * gec[0] - i - shift;
		if (ip[0]==i && (iarr==D_P0 || iarr==D_P1 || iarr==D_P2)) {
		  rd_gset(dom,iarr,ip,REAL_ZERO);
		}
		else {
		  v = rd_gget(dc, iarr, ip);
		  ip[0] = i;
		  rd_gset(dom, iarr, ip, sgn * v);
		}
	      }
	  }
        
      if      ( (iarr == D_P0) || (iarr == D_P1) || (iarr == D_P2) ||
		(iarr == D_V0) || (iarr == D_V2) ) 
	{
	  sgn = -1.0;
	  shift = -2;
	}
      else if ( iarr == D_V1 )
	{
	  shift = -1;
	  sgn = 1.0;
	}
      else return E_INTERNAL;
      for ( ip[2] = gsc[2]; ip[2] <= gec[2]; ++ip[2] )
	for ( ip[0] = gsa[0]; ip[0] <= gea[0]; ++ip[0] )
	  {
	    for ( i = gsa[1]; i <= gsc[1] - 1; ++i )
	      {
		ip[1] = 2 * gsc[1] - i + shift;
		if (ip[1]==i && (iarr==D_P0 || iarr==D_P1 || iarr==D_P2)) {
		  rd_gset(dom,iarr,ip,REAL_ZERO);
		}
		else {
		  v = rd_gget(dc, iarr, ip);
		  ip[1] = i;
		  rd_gset(dom, iarr, ip, sgn * v);
		}
	      }
	    for ( i = gea[1]; i >= gec[1] + 1; --i )
	      {
		ip[1] = 2 * gec[1] - i - shift;
		if (ip[1]==i && (iarr==D_P0 || iarr==D_P1 || iarr==D_P2)) {
		  rd_gset(dom,iarr,ip,REAL_ZERO);
		}
		else {
		  v = rd_gget(dc, iarr, ip);
		  ip[1] = i;
		  rd_gset(dom, iarr, ip, sgn * v);
		}
	      }
	  }
        
      if      ( (iarr == D_P0) || (iarr == D_P1) || (iarr == D_P2) ||
		(iarr == D_V0) || (iarr == D_V1) ) 
	{
	  sgn = -1.0;
	  shift = -2;
	}
      else if ( iarr == D_V2 )
	{
	  shift = -1;
	  sgn = 1.0;
	}
      else return E_INTERNAL;
      for ( ip[1] = gsa[1]; ip[1] <= gea[1]; ++ip[1] )
	for ( ip[0] = gsa[0]; ip[0] <= gea[0]; ++ip[0] )
	  {
	    for ( i = gsa[2]; i <= gsc[2] - 1; ++i )
	      {
		ip[2] = 2 * gsc[2] - i + shift;
		if (ip[2]==i && (iarr==D_P0 || iarr==D_P1 || iarr==D_P2)) {
		  rd_gset(dom,iarr,ip,REAL_ZERO);
		}
		else {
		  v = rd_gget(dc, iarr, ip);
		  ip[2] = i;
		  rd_gset(dom, iarr, ip, sgn * v);
		}
	      }
	    for ( i = gea[2]; i >= gec[2] + 1; --i )
	      {
		ip[2] = 2 * gec[2] - i - shift;
		if (ip[2]==i && (iarr==D_P0 || iarr==D_P1 || iarr==D_P2)) {
		  rd_gset(dom,iarr,ip,REAL_ZERO);
		}
		else {
		  v = rd_gget(dc, iarr, ip);
		  ip[2] = i;
		  rd_gset(dom, iarr, ip, sgn * v);
		}
	      }
	  }

      return 0;
    } 
    
  return E_BADINPUT;
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

int asg_create_sten(FILE * stream, 
		    int ndim,
		    IPNT gtype[RDOM_MAX_NARR], 
		    int sten_dep_mat[RDOM_MAX_NARR][RDOM_MAX_NARR], 
		    STENCIL * sten,
		    void * fdpars) {

  SGN_TS_PARS * sgnp = (SGN_TS_PARS *)(fdpars);
  //  return create_sten2_2k(stream, sgnp->k, ndim, m_size, gtype, sten_dep_mat, asg_isdyn, //asg_getindices, 
  return create_sten2_2k(stream, sgnp->k, ndim, gtype, sten_dep_mat, asg_isdyn, sten);
}

int asg_assign_action_array(FILE * stream, IMODEL * model) {

  int idim;
  int i;
  int ndim = (model->g).dim;

  /*-----------------------------------------------------------------*  
   * structure timestep by assigning action array
   *-----------------------------------------------------------------*/
  
  /* set timestep order ----------------------------------------------------*/
  model->tsinfo.narr = 2 * ndim;
  for ( idim = 0; idim < ndim; ++idim ) {
    model->tsinfo.arrs[idim] = D_P[idim];
    model->tsinfo.arrs[idim + ndim] = D_V[idim];
  }
    
  /* ORDER: n + 2*n + n =: 4*n

  px - COMPUTE (includes py, pz)
  py - COMPUTE (for pts)
  pz - COMPUTE (for pts)
        
  px - EXCHANGE
  vx - COMPUTE
        
  py - EXCHANGE
  vy - COMPUTE
        
  pz - EXCHANGE
  vz - COMPUTE
        
  vx - EXCHANGE
  vy - EXCHANGE
  vz - EXCHANGE
  */    
    
  model->tsinfo.npair = ndim + 2 * ndim + ndim;
  i = 0;
    
  for ( idim = 0; idim < ndim; ++idim ) { 
    model->tsinfo.pairs[i].arr = D_P[idim];
    model->tsinfo.pairs[i++].action = ACTION_COMPUTE;
  }    


  for ( idim = 0; idim < ndim; ++idim) { 
    model->tsinfo.pairs[i].arr = D_P[idim];
    model->tsinfo.pairs[i++].action = ACTION_EXCHANGE;
  }

  for ( idim = 0; idim < ndim; ++idim) { 
    model->tsinfo.pairs[i].arr = D_V[idim];
    model->tsinfo.pairs[i++].action = ACTION_COMPUTE;
  } 
        
  for ( idim = 0; idim < ndim; ++idim ) { 
    model->tsinfo.pairs[i].arr = D_V[idim];
    model->tsinfo.pairs[i++].action = ACTION_EXCHANGE;
  }    

  
  if ( i != model->tsinfo.npair ) {
    fprintf(stream,"ERROR: wrong number of timestep actions: %d, "
            "instead of %d.\n", i-1, model->tsinfo.npair);
    return E_INTERNAL;
  }
  return 0;
}

int asg_readschemeinfo(PARARRAY * par, 
		       FILE * stream, 
		       IMODEL * model) {
  
  RPNT dxs;  // grid steps
  int idim;  // counter

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
  }
  sgnp->dt = (model->tsind).dt;
  return 0;
}

void asg_refsubstep(int* itout,int* ivout, const IMODEL * m) {
  
  int ndim;

  ndim = (m->g).dim;
  *itout=(m->tsind).it;
  *ivout= 2*ndim-1;

}
