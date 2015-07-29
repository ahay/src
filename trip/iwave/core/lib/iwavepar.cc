/* 
WWS 29.06.08: combination of iwavepar.c and iwaveser.c, by Igor
Terentyev. Doing it this wave gets rid of annoying warning messages
from the linker.
********************************************************************************
MPI-enabled version.
*/

/*============================================================================*/

/*
#define VERBOSE
*/

/*----------------------------------------------------------------------------*/

#include "utils.h"
#include "iwave.h"
/*----------------------------------------------------------------------------*/

/* added 10.04.11 WWS - cartesian grid info, init in initparallel_local */

/*----------------------------------------------------------------------------*/

#ifdef IWAVE_USE_MPI /* this code for MPI version only */

/* int destroyparallel(PARALLELINFO *pinfo, int fabort) { // name change 10.04.11 */
int destroypinfo(PARALLELINFO *pinfo, int fabort) {
    int iv, i;
    for ( iv = 0; iv < RDOM_MAX_NARR; ++iv )
        for ( i = 0; i < IWAVE_NNEI; ++i ) 
        {
            ei_destroy(&(pinfo->seinfo[iv][i]));
            ei_destroy(&(pinfo->reinfo[iv][i]));
        }

    /* donnot free the communicator 
      if ( ( pinfo->ccomm != MPI_COMM_NULL ) && ( pinfo->ccomm != MPI_COMM_WORLD ) ) MPI_Comm_free(&(pinfo->ccomm));

      pinfo->ccomm = MPI_COMM_NULL;
    */
    
    /*
      deprecated - MPI_Finalize should be called in driver
    MPI_Barrier(MPI_COMM_WORLD); 
    if ( fabort ) MPI_Abort(MPI_COMM_WORLD, 0);
    else MPI_Finalize();
    */
    return 0;
}
/*----------------------------- DEPRECATED -----------------------------------------------*/
int initparallel(int ts)
{
  int rk;
  int sz;
  /*
    Moved to driver - symm with MPI_Finalize
  MPI_Init_thread(argc, argv, MPI_THREAD_FUNNELED, &ts);
  */
  MPI_Comm_size(MPI_COMM_WORLD, &sz);
  MPI_Comm_rank(MPI_COMM_WORLD, &rk);
  /* global params */
  storeGlobalComm(MPI_COMM_WORLD);
  storeGlobalRank(rk);
  storeGlobalSize(sz);
  /* local params initially same as global - can be changed later */
  storeComm(MPI_COMM_WORLD);
  storeRank(rk);
  storeSize(sz);
  /* thread support level */
  storeThreadSupp(ts);
  return 0;
} 

/*----------------------------------------------------------------------------*/
void initparallel_global(int ts) {
  /* global params */
  int rk;
  int sz;

  MPI_Comm_size(MPI_COMM_WORLD, &sz);
  MPI_Comm_rank(MPI_COMM_WORLD, &rk);

  storeGlobalComm(MPI_COMM_WORLD);
  storeGlobalRank(rk);
  storeGlobalSize(sz);
  storeThreadSupp(ts);
}

/*----------------------------------------------------------------------------*/
int initparallel_local(PARARRAY par, FILE * stream) {

  int rkw;      /* global rank */
  int szw;      /* global size */
  int rk;       /* local rank */
  int sz;       /* local size */
  int ng;       /* group number */
  int g;        /* group id (color) */
  int tg;       /* trimmed group id (= 1 or MPI_UNDEFINED) */
  int d;        /* (domain) index within group */
  int pt;       /* task parallel flag */
  div_t res;    /* workspace */
  int idim;     /* workspace */

  int cartdim;  /* dim of cart comm */
  IPNT cdims;   /* lengths of axes in cart comm */
  IPNT pers;    /* periodicity flags for cartesian comm */

  MPI_Comm cmw; /* initial global comm (WORLD) */
  MPI_Comm cmg; /* trimmed global comm (does not included undef processes) */
  MPI_Comm cm;  /* local (intragroup Cartesian) comm */
  MPI_Comm cmr; /* remote (intergroup) comm */
  MPI_Comm cmtmp; /* workspace */

  //  cerr<<"initparallel_local -> enter\n";
  /* global params */
  cmw=retrieveGlobalComm();
  rkw=retrieveGlobalRank();
  szw=retrieveGlobalSize();

  //  cerr<<"initparallel_local -> defaults\n";
  /* defaults (non-task-parallel) for local params */
  cm=cmw;
  cmtmp=cmw;
  cmr=0;
  cmg=cmw;
  rk=rkw;
  sz=1;
  ng=1;
  g =0;
  d =rkw;
  
  //  cerr<<"initparallel_local -> detect task parallelization\n";
  /* detect task parallelization */
  pt=0;
  ps_flint(par,"partask",&pt);

  /* read cartesian grid dimn */
  /*
  ndim = 0;

  if (ps_flint(par,"pardim",&ndim) || ndim<1 || ndim>RARR_MAX_NDIM) {
    fprintf(stream,"ERROR: must define pardim in parameter table\n");
    fprintf(stream,"- should be in range [1,%d] and same as spatial dimension\n",
	    RARR_MAX_NDIM);
    return E_BADINPUT;
  }
  */
  /* expt!!! */
  cartdim=1;
  IASN(cdims, IPNT_1); /* default grid size */
    
  //  cerr<<"initparallel_local -> parse cartesian grid\n";
  /* parse out cartesian grid for local size */
  sz=1;
  for ( idim = 0; idim < RARR_MAX_NDIM; ++idim ) {
    stringstream pname;
    pname<<"mpi_np";
    pname<<idim+1;
    if ( ps_flint(par, pname.str().c_str(), &(cdims[idim]))) {
#ifdef IWAVE_VERBOSE
      fprintf(stream, "NOTE. Cannot read %s. Default = 1.\n", pname.str().c_str());
#endif
      //      cerr<<"NOTE. Cannot read "<<pname.str()<<" pnames for idim="<<idim<<". Default = 1.\n";
    }
    else {
      if ( cdims[idim] < 1 ) {
	fprintf(stream, "ERROR. Bad input: number of processors %s = %d.\n", pname.str().c_str(), cdims[idim]);
	return E_BADINPUT;
      }
      //      cerr<<"cdims["<<idim<<"]="<<cdims[idim]<<"\n";
      sz *= cdims[idim];
      /*      fprintf(stderr,"cdims[%d]=%d\n",idim,cdims[idim]); */
      if (cdims[idim] > 1) cartdim=idim+1;
     
    }
  }
  //  cerr<<"cartdim="<<cartdim<<" sz="<<sz<<"\n";
  if (sz<0 ||sz>szw) {
    fprintf(stream,"Error: initparallel_local\n");
    fprintf(stream,"sz=%d computed from parameter array (total\n",sz);
    fprintf(stream,"number of domains) not in range [0, world size=%d]\n",
	    szw);
#ifdef IWAVE_USE_MPI
    MPI_Abort(cmw,E_BADINPUT);
#else
    return E_BADINPUT;
#endif
  }
  
  /* task-parallel branch - the only difference
     between the task-parallel and task-serial 
     branches is that with task parallelism, the 
     number of groups and group ids for live ranks
     are computed, rather than assigned. */

  //  cerr<<"initparallel_local -> task-parallel\n";
  if (pt) {
    res=div(szw,sz);
    /* allow user control of number of groups - since
    // number of simulations not know at this stage of
    // construction, must assign responsibility to user */
    ng=res.quot;
    ng=iwave_min(ng,pt);
    res=div(rkw,sz);
    g=res.quot;
    d=res.rem;
  }

  /* in all cases, world ranks not belonging to a group (of which
     there is only one in the task-serial case) are assigned the
     "undefined" color, which means that they do not participate in
     any comm 
  */
  tg=1;
  if (rkw>sz*ng-1) {
    g=MPI_UNDEFINED;
    tg=g;
    d=MPI_UNDEFINED;
  }
   
  //  cerr<<"initparallel_local -> create comms\n";
  /* create intragroup comm - all domains, one task */
  MPI_Comm_split(cmw,g,rkw,&cmtmp);
  /* create intergroup comm - all tasks, one domain */
  MPI_Comm_split(cmw,d,rkw,&cmr);
  /* create new global communicator, without dead processes */
  MPI_Comm_split(cmw,tg,rkw,&cmg);

  //  cerr<<"initparallel_local -> create cartesian grid\n";
  /* create cartesian grid, identify this process --------------------------*/

  IASN(pers, IPNT_0);
  
  /* ONLY CREATE CARTESIAN GRID FOR DEFINED COMMUNICATORS */
  if ( g != MPI_UNDEFINED ) {
    if (MPI_Cart_create(cmtmp,cartdim, cdims, pers, 1, &cm) != MPI_SUCCESS ) {
      fprintf(stream, "ERROR. initparallel_local\n");
      fprintf(stream, "Internal: cannot create Cartesian grid.\n");
      return E_INTERNAL;
    }
    if ( cm == MPI_COMM_NULL ) {
#ifdef IWAVE_VERBOSE
      fprintf(stream, "NOTE. Processor %d not in the Cartesian grid.\n", rkw);
#endif
      return E_NOTINGRID;
    }
  }
  else {
    cm=cmtmp;
  }

  //  cerr<<"initparallel_local -> find local, global rank, size\n";
  /* find local, global rank, size */
  if ((g != MPI_UNDEFINED) && (cm != MPI_COMM_NULL)) {
    MPI_Comm_rank(cmg,&rkw);
    MPI_Comm_size(cmg,&szw);
    MPI_Comm_rank(cm,&rk);
    MPI_Comm_size(cm,&sz);
#ifdef IWAVE_VERBOSE
    fprintf(stream,"NOTE: initparallel_local\n");
    fprintf(stream,"NOTE: global rank = %d global size = %d number of groups = %d group = %d group size = %d local rank = %d\n",rkw,szw,ng,g,sz,rk);
    fflush(stream);
#endif
  }
  else {
#ifdef IWAVE_VERBOSE
    fprintf(stream,"NOTE: initparallel_local\n");
    fprintf(stream,"NOTE: world rank = %d not in any group\n",rkw);
    fflush(stream);
#endif
  }  

  //  cerr<<"initparallel_local -> store params\n";
  /* store local params, group info */
  storeComm(cm);
  storeRank(rk); /* should be same as d - index in local comm */
  storeSize(sz);
  storeGroupID(g);
  storeNumGroups(ng);

  /* WWS 17.01.11 */
  storeRemComm(cmr);
    
  /* WWS 24.08.11
     redefine global communicator without undef processes */
  storeGlobalComm(cmg);
  storeGlobalRank(rkw);
  storeGlobalSize(szw);
  
  //  cerr<<"initparallel_local -> exit\n";
  return 0;
} 
  
/*----------------------------------------------------------------------------*/
int initpinfo(PARALLELINFO *pinfo, FILE *stream) {
  int iv, i;
  IPNT pers;
  
  pinfo->threadsupp=retrieveThreadSupp();
  pinfo->wsize=retrieveSize();
  pinfo->wrank=retrieveRank();
  pinfo->ccomm = retrieveComm();
  pinfo->lrank = retrieveRank();

  IASN(pinfo->cdims, IPNT_1); /* default grid size */ 
  IASN(pinfo->crank, IPNT_0); /* default cartisian ranks */ 
  if ( MPI_Cart_get(pinfo->ccomm, RARR_MAX_NDIM, pinfo->cdims, pers, pinfo->crank) != MPI_SUCCESS )  {
    fprintf(stream, "ERROR. Internal: cannot get Cartesian coordinates.\n");
    return E_INTERNAL;
  }

  /* cannot set exch info until problem dimension is known */
  for ( iv = 0; iv < RDOM_MAX_NARR; ++iv ) {
    for ( i = 0; i < IWAVE_NNEI; ++i ) {
      ei_setnull(&(pinfo->seinfo[iv][i]));
      ei_setnull(&(pinfo->reinfo[iv][i]));
    }
  }
  return 0;
}
/*----------------------------------------------------------------------------*/

/* name change 10.04.11
int createparallel(PARALLELINFO *pinfo, int ndim, FILE *stream) {
*/
int initexch(PARALLELINFO *pinfo, int ndim, FILE *stream) {
  int rbproc, sbproc, cartdim, idim, i; /* grid size, boundary flag, cartesian dim, iterators */
  IPNT roffs, soffs; /* neighbor index recv/send offset */

#ifdef VERBOSE
  fprintf(stderr,"In createparallel:\n");
  fprintf(stderr," retrieve infomation related to cartesian comm\n");
#endif

  /* first task: check that cart comm dim is <= physical grid dim  */
  if (MPI_SUCCESS != MPI_Cartdim_get(pinfo->ccomm,&cartdim)) {
    fprintf(stream,"Error: createparallel\n");
    fprintf(stream,"- failed to extract cartesian dims from cartesian comm\n");
    return E_INTERNAL;
  }
  if (cartdim > ndim) {
    fprintf(stream,"ERROR: initexch\n");
    fprintf(stream,"- cartesian communicator dimn = %d\n",cartdim);
    fprintf(stream,"- in this version of IWAVE, must be less than physical grid dimn = %d\n",ndim);
    fprintf(stream,"- (we could make it work, but you are using too many processes!!)\n");
    return E_BADINPUT;
  }

#ifdef VERBOSE
  fprintf(stderr,"call gen_3n1\n");
#endif
  /* set physical dimn, get number of neighbors */
  pinfo->ndim=ndim;
  if ( gen_3n1(ndim, &(pinfo->nnei)) ) {
    fprintf(stream, "ERROR. Internal: cannot compute number of neighbors for ndim = %d.\n", pinfo->ndim);
    return E_INTERNAL;
  }

#ifdef VERBOSE
  fprintf(stderr,"nnei=%d ndim=%d\n",pinfo->nnei,pinfo->ndim);
#endif

  /* get neighbor offsets ---------------------------------------------------------*/
  for ( i = 0; i < pinfo->nnei; ++i )  /* neighbor ranks */ {
    if ( gen_i2pnt(pinfo->ndim, i, roffs) ) {
      fprintf(stream, "ERROR. Internal: cannot get neighbor [%d] coordinates.\n", i);
      return E_INTERNAL;
    }
    IASN(soffs, roffs);
    rbproc = sbproc = 0;
    for ( idim = 0; idim < pinfo->ndim; ++idim ) {
      roffs[idim] = pinfo->crank[idim] + roffs[idim]; /* Cartesian coordinate */
      soffs[idim] = pinfo->crank[idim] - soffs[idim]; /* Cartesian coordinate */
      if ( (roffs[idim] < 0) || (roffs[idim] >= pinfo->cdims[idim]) ) rbproc = 1;
      if ( (soffs[idim] < 0) || (soffs[idim] >= pinfo->cdims[idim]) ) sbproc = 1;
    }
    if ( rbproc ) pinfo->rranks[i] = MPI_PROC_NULL;
    else {
      if ( MPI_Cart_rank(pinfo->ccomm, roffs, pinfo->rranks + i) != MPI_SUCCESS ) {
	fprintf(stream, "ERROR. Internal: cannot get neighbor rank.\n");
	return E_INTERNAL;
      }
    }
    
    if ( sbproc ) pinfo->sranks[i] = MPI_PROC_NULL;
    else {
      if ( MPI_Cart_rank(pinfo->ccomm, soffs, pinfo->sranks + i) != MPI_SUCCESS ) {
	fprintf(stream, "ERROR. Internal: cannot get neighbor rank.\n");
	return E_INTERNAL;
      }
    }
  }

  fflush(stream);
  return 0;
}
/*---------------------- END PARALLEL ----------------------------------------*/
#else
/*-------------------- BEGIN SERIAL ------------------------------------------*/

/* int destroyparallel(PARALLELINFO *pinfo, int fabort) renamed 10.04.11 */
int destroypinfo(PARALLELINFO *pinfo, int fabort)
{
    return 0;
}
/*----------------------------------------------------------------------------*/

int initparallel(int ts) {
  return 0;
}

/*----------------------------------------------------------------------------*/

void initparallel_global(int ts) {}

/*----------------------------------------------------------------------------*/

int initparallel_local(PARARRAY par, FILE * stream) {
  int idim;
  IPNT cdims; /* local declaration added 02.07.11 WWS */
  IASN(cdims, IPNT_1); /* default grid size */
  for ( idim = 0; idim < IWAVE_NDIM; ++idim ) {
    stringstream pname;
    pname<<"mpi_np";
    pname<<idim+1;
    if ( ps_flint(par, pname.str().c_str(), cdims + idim) ) {
#ifdef IWAVE_VERBOSE
      fprintf(stream, "NOTE. Cannot read %s. Default = %d.\n", pname.str().c_str(), cdims[idim]);
#endif
    }
    if (cdims[idim] != 1) {
      fprintf(stream, "ERROR: nontrivial domain decomposition not possible\n");
      fprintf(stream, "       in serial mode of IWAVE; recompile with MPI\n");
      fprintf(stream, "       or use %s=1\n",pname.str().c_str());
      fprintf(stderr, "ERROR: nontrivial domain decomposition not possible\n");
      fprintf(stderr, "       in serial mode of IWAVE; recompile with MPI\n");
      fprintf(stderr, "       or use %s=1\n",pname.str().c_str());
      return E_BADINPUT;
    }
  }
  return 0;
}

/*----------------------------------------------------------------------------*/

int initpinfo(PARALLELINFO *pinfo, FILE *stream) {

    int iv, i;
    
    pinfo->ccomm = 0;

    for ( iv = 0; iv < RDOM_MAX_NARR; ++iv )
      for ( i = 0; i < IWAVE_NNEI; ++i ) {
	ei_setnull(&(pinfo->seinfo[iv][i]));
	ei_setnull(&(pinfo->reinfo[iv][i]));
      }
    
    IASN(pinfo->cdims, IPNT_1); /* default grid size */
    pinfo->ndim = 1;//RARR_MAX_NDIM;
      
    if ( gen_3n1(pinfo->ndim, &(pinfo->nnei)) )
    {
      fprintf(stream, "ERROR. Internal: cannot compute number of neighbors for ndim = %d.\n", pinfo->ndim);
      return E_INTERNAL;
    }

    /* This test already performed in initparallel_local
    
    for ( idim = 0; idim < RARR_MAX_NDIM; ++idim )
    {
      if ( ps_flint(*pars, PNAMES_NP[idim], cdims + idim) )
	fprintf(stream, "NOTE. Cannot read %s. Default = %d.\n", PNAMES_NP[idim], cdims[idim]);
      if (cdims[idim] != 1) {
	fprintf(stream, "ERROR: from createparallel -- serial mode\n");
	fprintf(stream, "       nontrivial domain decomposition not possible\n");
	fprintf(stream, "       in serial mode of IWAVE; recompile with MPI\n");
	fprintf(stream, "       or use %s=1\n",PNAMES_NP[idim]);
	fprintf(stderr, "ERROR: from createparallel -- serial mode\n");
	fprintf(stderr, "       nontrivial domain decomposition not possible\n");
	fprintf(stderr, "       in serial mode of IWAVE; recompile with MPI\n");
	fprintf(stderr, "       or use %s=1\n",PNAMES_NP[idim]);
	return E_BADINPUT;
      }
    }
    */    
    /* set ranks -------------------------------------------------------------*/
    /*    for ( i = 0; i < pinfo->nnei; ++i ) rranks[i] = sranks[i] = MPI_PROC_NULL; */
    for ( i = 0; i < pinfo->nnei; ++i ) pinfo->rranks[i] = pinfo->sranks[i] = 0;
    pinfo->lrank = 0;
    IASN(pinfo->crank, IPNT_0);

    return 0;
}
/*----------------------------------------------------------------------------
No longer needed 
int readparallel(PARALLELINFO *pinfo, PARARRAY *pars, FILE *stream)
{

    }

    return 0;
}
*/
/*----------------------------------------------------------------------------*/
/* name change 10.04.11
int createparallel(PARALLELINFO *pinfo, FILE *stream)
*/
int initexch(PARALLELINFO *pinfo, int ndim, FILE *stream) {
    /* int initexch(PARALLELINFO *pinfo, FILE *stream) */
    return 0;
}
/*----------------------------------------------------------------------------*/

#endif /* IWAVE_USE_MPI*/

int setrecvexchange(IMODEL * model, PARALLELINFO * pinfo, FILE * stream, IWaveInfo const & ic) {

  int err=0;
#ifdef IWAVE_USE_MPI
  int i, iv, j, ia, overlap, ndim, tmp_narr, index;
  RDOM *ld_r, *ld_s, *ld_c;

  MPI_Status status;
  IPNT gs, ge, gs2, ge2;

  /* narr, {array #, ndim, gs[0], ..., gs[ndim-1], ge[0], ..., ge[ndim-1]} */
  int recvinfo[1 + RDOM_MAX_NARR * (2 + 2 * RARR_MAX_NDIM)]; /* info about recv arrays - it is sent */
  int sendinfo[1 + RDOM_MAX_NARR * (2 + 2 * RARR_MAX_NDIM)]; /* info about send arrays - it is received */
  int infocount = 1 + RDOM_MAX_NARR * (2 + 2 * RARR_MAX_NDIM);

  //  FD_MODEL *fdm = (FD_MODEL *)model->specs;
  
  ld_r = model->ld_r;
  ld_s = model->ld_s;
  ld_c = &(model->ld_c);


  /* 12.03.12 WWS */
  /*  narr = model->tsinfo.narr; */
  
  /* Set all receives to empty ---------------------------------------------*/
  for ( ia = 0; ia < RDOM_MAX_NARR; ++ia ) {
    for ( i = 0; i < model->nnei; ++i ) {
      ei_destroy(&(pinfo->reinfo[ia][i]));
      ei_destroy(&(pinfo->seinfo[ia][i]));
    }
  }
        
  /* Process recomputed arrays ---------------------------------------------*/
  for ( ia = 0; ia < RDOM_MAX_NARR; ++ia ) {
    /* only process dynamic arrays */
    //    if (isdyn(fdm,ia)) {
    if (fd_isdyn(ia,ic)) {

      /* parse neighbors */
      for ( i = 0; i < model->nnei; ++i ) {
	/* overlap check */
	/* parse another neighbor */
	for ( j = 0; j < i; ++j ) {
	  err = rd_overlap(ld_r + i, ia, ld_r + j, ia, &overlap);
	  if ( err ) {
	    fprintf(stream, 
		    "ERROR. Internal: setrecvexchange overlap function error #%d. ABORT.\n", 
		    err);
	    return err;
	  }                
	  if ( overlap ) {
	    fprintf(stream, 
		    "ERROR. Internal: setrecvexchange receive ares overlap. ABORT.\n");
	    return E_DOMAINDECOMP;
	  }                
	}
	/* create type */
	err = rd_setexchangeinfo(ld_r + i, ia, &(pinfo->reinfo[ia][i]));
	if ( err ) {
	  fprintf(stream, 
		  "ERROR. Internal: setrecvexchange set exchange information error #%d. ABORT.\n", 
		  err);
	  
	  return err;
	}                
      }
    } 
  }
    
  /* Exchange recv array info ----------------------------------------------*/
  
  /* parse neighbors */
  for ( i = 0; i < model->nnei; ++i ) {
    /* fill recvinfo */
    index = 0;
    // compute number of dynamic fields
    int ndyn = 0;
    for (int kk=0; kk<ic.get_num_fields(); kk++) 
      if (fd_isdyn(kk,ic)) ndyn++;
    //    recvinfo[index++] = narr(fdm);
    recvinfo[index++] = ndyn;
    /* loop over dynamic arrays */
    for ( ia = 0; ia < RDOM_MAX_NARR; ++ia ) {
	/*      ia = model->tsinfo.arrs[iv];  recomputed array index */
      //     if (isdyn(fdm,ia)) {
      if (fd_isdyn(ia,ic)) {
	recvinfo[index++] = ia;
	
	rd_ndim(ld_r + i, ia, &ndim);/* number of dimensions in the array */
	recvinfo[index++] = ndim;
	
	rd_gse(ld_r + i, ia, gs, ge);
	for ( j = 0; j < ndim; ++j ) recvinfo[index++] = gs[j];
	for ( j = 0; j < ndim; ++j ) recvinfo[index++] = ge[j];
      }
    }

    for ( ; index < infocount; ) recvinfo[index++] = 0;
    
    /* exchange */
    err = MPI_Sendrecv(recvinfo, infocount, MPI_INT, pinfo->rranks[i], 0, 
		       /* send this */
		       sendinfo, infocount, MPI_INT, pinfo->sranks[i], 0, 
		       /* recv this */
		       pinfo->ccomm, &status);
    if ( err != MPI_SUCCESS ) {
      fprintf(stream, "ERROR. MPI_Sendrecv for recv info exchange. ABORT.\n");
      return err;
    }
    /* read sends */
    if ( pinfo->sranks[i] == MPI_PROC_NULL ) {
      rd_a_narr(ld_s + i, &j);
      for ( iv = 0; iv < j; ++iv ) rd_setnull(ld_s + i, iv);
      continue;
    }
        
    index = 0;
    tmp_narr = sendinfo[index++];
    for ( iv = 0; iv < tmp_narr; ++iv ) {
      ia = sendinfo[index++];
      ndim = sendinfo[index++];
      for ( j = 0; j < ndim; ++j ) gs[j] = sendinfo[index++];
      for ( j = 0; j < ndim; ++j ) ge[j] = sendinfo[index++];
      
      for ( j = 0; j < ndim; ++j ) if ( gs[j] > ge[j] ) break;
      /* do not greset empty receives */ 
      if ( j < ndim ) {
	rd_setnull(ld_s + i, ia);
	continue;
      }
            
      err = rd_greset(ld_s + i, ia, gs, ge);
      if ( err ) {
	fprintf(stream, 
		"ERROR. Internal: send array (%d) of domain (%d) greset error #%d. ABORT.\n", 
		ia, i, err);
	return err;
      }
      err = rd_setexchangeinfo(ld_s + i, ia, &(pinfo->seinfo[ia][i]));
      if ( err ) {
	fprintf(stream, 
		"ERROR. Internal: set exchange information error #%d. ABORT.\n", 
		err);
	return err;
      }                
    }
  }
    
  /* Process send arrays ---------------------------------------------*/
  for ( ia = 0; ia < RDOM_MAX_NARR; ++ia ) {
    //    if (isdyn(fdm,ia)) {        /* recomputed array index */
    if (fd_isdyn(ia,ic)) {        /* recomputed array index */
      /* parse neighbors */
      for ( i = 0; i < model->nnei; ++i ) {
	/* overlap check */
	/* parse another neighbor */
	for ( j = 0; j < i; ++j ) {
	  rd_empty(ld_s + i, ia, &index);
	  if ( index ) continue;
	  rd_empty(ld_s + j, ia, &index);
	  if ( index ) continue;
	  err = rd_overlap(ld_s + i, ia, ld_s + j, ia, &overlap);
	  if ( err ) {
	    fprintf(stream, 
		    "ERROR. Internal: setrecvexchange overlap function error #%d, arr = %d. ABORT.\n", 
		    err, ia);
	    return err;;
	  }                
	  if ( overlap ) {
	    fprintf(stream, 
		    "ERROR. Internal: setrecvexchange send ares overlap, arr = %d. ABORT.\n", 
		    ia);
	    rd_gse(ld_s + i, ia, gs, ge);
	    rd_gse(ld_s + j, ia, gs2, ge2);
	    for ( index = 0; index < ndim; ++index )
	      fprintf(stream, 
		      "       [%d %d] - [%d %d]\n", gs[index], ge[index], gs2[index], ge2[index]);
	    return err;
	  }                
	}
      }
    }
  }
#endif
  return err;
}
