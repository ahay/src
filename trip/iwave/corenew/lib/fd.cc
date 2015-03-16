/**
 * XW 01-18-2010
 * wangxin.tom@gmail.com
 * fd.c
 */

#include "fd.h"

/* helper functions */

int fd_isarr(int i, IMODEL & model, IWaveInfo const & ic) {
  if (i<0 || i >= ic.get_num_fields()) {
    return 0;
  }
  for (size_t j=0; j<model.active.size();j++) {
    if (ic.iwave_fields[i].field == model.active[j]) {
      return 1;
    }
  }
  return 0;
}

int fd_isdyn(int i, IWaveInfo const & ic) {
  if (i>-1 || i < ic.get_num_fields()) {
    return ic.get_iwave_fields()[i].dynamic;
  }
  return 0;
}

const char* fd_ind2str(int i, IWaveInfo const & ic) {
  if (i>-1 && i<=ic.get_num_fields()) return ic.get_iwave_fields()[i].field.c_str();
  return NULL;
}

int fd_numsubsteps(IWaveInfo const & ic) {
  int ns = true;
  for (int i=0;i<ic.get_num_fields(); i++) 
    ns = iwave_max(ic.get_iwave_fields()[i].substep+1, ns);
  return ns;
}

bool fd_update(int ia, int iv, IWaveInfo const & ic) {
  if (ia<0 || ia>=ic.get_num_fields()) return false;
  if (iv==ic.get_iwave_fields()[ia].substep) return true;
  return false;
}
      
int fd_set_grid_type(FILE * stream, IPNT gtype[RDOM_MAX_NARR], IWaveInfo const & ic) {
  for (int i=0;i<ic.get_num_fields();i++) IASN(gtype[i],ic.get_iwave_fields()[i].gtype);
  for (int i=ic.get_num_fields();i<RDOM_MAX_NARR;i++) IASN(gtype[i],IPNT_0);
  return 0;
}

int fd_readgrid(PARARRAY * pars, 
		FILE *stream, 
		IMODEL * model,
		IWaveInfo const & ic) {
		//		std::string gridkey) {

  /* logic:

     page through fields:

     for each static field rarray index (position in array of fields),
     find i/o key for this rarray index - there should be only one per
     (check!).

     Read grid from file - extract physical spatial grid, should be same for
     each input static field (check!)

     Record grid in array of static field grids (data member of IMODEL).
  */

  // cerr<<"in fd_readgrid\n";
  std::string gridfile="";    /* workspace for filename */
  int err=0;
  
  // simplest method to page through fields is to page through iokeys 
  for (int io=0; io < ic.get_num_iokeys(); io++) {

    // only look at static fields - should all be readable
    // since iokeys only lists the base-level (deriv=0) fields, they 
    // should be unique
    if (!(ic.iwave_fields[ic.iwave_iokeys[io].rarrindex].dynamic)) {

      // extract filename from entry in field table
      gridfile="";
      if (!parse(*pars,ic.iwave_iokeys[io].keyword,gridfile)) {
	//    cerr<<"error return err="<<err<<endl;
	fprintf(stream,"Error: fd_readgrid from ps_flcstring\n");
	fprintf(stream,"  failed to parse key = %s\n",ic.iwave_iokeys[io].keyword.c_str());
	ps_printall(*pars,stream);
	fflush(stream);
	return E_BADINPUT;
      }
      else {
	// mod 30.11.13: read only spatial axes for primal grid - ignore others
	// mod 07.03.15: read all static grids - set primal grid from first one,
	// then check others against it (just spatial axes)
	grid gtmp;
	init_default_grid(&gtmp);
	//	init_default_grid(&(model->g)); should have been done in im_construct
	//	cerr<<"fd_readgrid: reading gridfile="<<gridfile<<endl;
	err=read_grid(&gtmp,gridfile.c_str(),stream);
	if (err==0) {
	  // use model->g.gdim==0 as proxy for "first time" - copy spatial 
	  // axes onto prototype spatial grid 
	  if ((model->g).gdim == 0) {
	    (model->g).gdim = gtmp.dim;
	    (model->g).dim  = gtmp.dim;
	    // this loop copies spatial axes in the order encountered, regardless
	    // of index in source grid - result is appropriate for dynamic field
	    // alloc, but not necessarily for static, which may interpose other axes
	    int j = 0;
	    for (int i=0;i<gtmp.gdim;i++) {
	      if ((gtmp.axes[i].id >= 0) && (gtmp.axes[i].id < gtmp.dim)) {
		copy_axis(&((model->g).axes[j]),&(gtmp.axes[i]));
		j++;
	      }
	    }
	    if (j != gtmp.dim) {
	      fprintf(stream,"Error: fd_readgrid\n");
	      fprintf(stream,"  failed to find correct number of axes with ids between 0 and dim\n");
	      return E_BADINPUT;
	    }
	    
	  }
	  // otherwise compare with prototype axes - should be same
	  // per id's.
	  else {
	    if (((model->g).dim  != gtmp.dim)) {
	      fprintf(stream,"Error: fd_readgrid\n");
	      fprintf(stream,"  grid dim = %d from io keyword %s with differs from prototype grid dim = %d\n",gtmp.dim, ic.iwave_iokeys[io].keyword.c_str(), (model->g).dim);
	      return E_BADINPUT;
	    }
	    // loop over id's of spatial axes
	    for (int j=0; j<(model->g).dim; j++) {
	      // loop over comparison axes
	      for (int i=0;i<RARR_MAX_NDIM;i++) {
		if (gtmp.axes[i].id == j) {
		  if (compare_axis((model->g).axes[j],gtmp.axes[i])) {
		    fprintf(stream,"Error: fd_readgrid\n");
		    fprintf(stream,"  axis %d from io keyword %s with id %d\n",
			     i, ic.iwave_iokeys[io].keyword.c_str(), gtmp.axes[i].id);
		    fprintf(stream,"  differs from axis %d of prototype grid with same id\n",
			   j);
		    return E_BADINPUT;
		  }
		}
	      }
	    }
	  }
	  // copy whole thing into appropriate slot in sfg table
	  //	  cerr<<"fd_readgrid: copying grid onto sfg["<<ic.iwave_iokeys[io].rarrindex<<"]\n";
	  //	  fprint_grid(stderr,gtmp);
	  copy_grid(&((model->sfg)[ic.iwave_iokeys[io].rarrindex]),&gtmp);
	}
	else {
	  fprintf(stream,"Error: fd_readgrid from read_grid, failed to read grid on keyword = %s\n",
		  ic.iwave_iokeys[io].keyword.c_str());
	  fprintf(stream,"  err=%d\n",err);
	  return err;
	}
      }
    }
  }

  /* moved to modelinit, where it belongs - 05.01.15 WWS 
  ireal dx;
  for (int idim=0; idim<model->g.dim; idim++ ) {
    if ((dx=model->g.axes[idim].d)<REAL_EPS) {
      fprintf(stream, "ERROR: fd_readgrid, dx[%d] = %f\n", idim,dx);
      return E_BADINPUT;
    }
    {
      stringstream key;
      key<<"nl"<<idim+1<<'\0';
      //    snprintf(key,kl,"nl%d",idim+1);
      ireal tmp=0.0;
      //      fprintf(stderr,"read key=%s\n"key.str().c_str());
      ps_flreal(*pars, key.str().c_str(), &tmp);
      //      fprintf(stderr,"result=%g\n",tmp);
      model->nls[idim]=iwave_max(0,(int)((ceil)(tmp/dx)));
    }
    //    snprintf(key,kl,"nr%d",idim+1);
    {
      stringstream key;
      key<<"nr"<<idim+1<<'\0';
      ireal tmp=0.0;
      //      fprintf(stderr,"read key=%s\n",key.str().c_str());
      ps_flreal(*pars,key.str().c_str(),&tmp);
      //      fprintf(stderr,"result=%g\n",tmp);
      model->nrs[idim]=iwave_max(0,(int)((ceil)(tmp/dx)));
    }
    //    fprintf(stderr,"fd_readgrid:: nls[%d]=%d, nrs[%d]=%d\n",idim,model->nls[idim],idim,model->nrs[idim]);
    //    ps_printall(*pars,stderr);
  }
  */
  return err;
}

int fd_setcompdom(FILE * stream, IPNT cdims, IPNT crank, 
                  IMODEL * model, IPNT dgs[], IPNT dge[],
		  IPNT gtype[RDOM_MAX_NARR],
		  IWaveInfo const & ic) {
    /**
     * virtual start and end: physical + artificial,  
     * assume physical start at IPNT_0
     */
    IPNT gs_pa, ge_pa;  
    IPNT gn_pa;
    IPNT ls, le; 
    IPNT ns;
    int i, idim;
    int ndim = (model->g).dim;

    for (i = 0;i < RDOM_MAX_NARR;i++) {
      IASN(dgs[i], IPNT_1);
      IASN(dge[i], IPNT_0);
    }
 
    get_n(ns, model->g);
    IASN(gs_pa, IPNT_1);
    IASN(ge_pa, IPNT_0);
    IASN(ls, IPNT_1);
    IASN(le, IPNT_0);

    /*< ls, le: local start, end for primal grid */
    for (idim = 0;idim < ndim;idim ++) {
      gs_pa[idim] = -model->nls[idim];
      ge_pa[idim] = ns[idim] + model->nrs[idim] - 1;
      gn_pa[idim] = ns[idim] + model->nls[idim] + model->nrs[idim];
      ls[idim] = gn_pa[idim] * (long)crank[idim];
      ls[idim] = gs_pa[idim] + ls[idim] / (long)cdims[idim];
      le[idim] = gn_pa[idim] * (long)(crank[idim] + 1);
      le[idim] = gs_pa[idim] + le[idim] / (long)cdims[idim]-1;
      if (le[idim] < ls[idim]) {
	fprintf(stream, "Error: in fd_setcompdom: le[%d] = %d < ls[%d] = %d\n",
		idim, le[idim], idim, ls[idim]);
	return E_INTERNAL;
      }
    } 

#if INCLUDE_BOUNDARY_PNTS
    /*< include left and right boundary pnts in all arrays */   
    for (i = 0;i < RDOM_MAX_NARR;i ++) {
      if (fd_isarr(i,*model,ic)) {
	if (fd_isdyn(i,ic)) {
	  for (idim = 0;idim < ndim;idim ++) {
	    dgs[i][idim] = ls[idim];
	    dge[i][idim] = le[idim];
	    if ( crank[idim] == cdims[idim]-1 && 
		 gtype[i][idim] == DUAL_GRID )
	      dge[i][idim] --;
	  }
	}
	/** for static fields, assign for internal extended axes.
	    Note that, being internal, these axes are not localized by
	    domain decomposition 
	*/
	else {
	  IPNT g_gs;
	  IPNT g_ge;
	  get_gs(g_gs,(model->sfg)[i]);
	  get_ge(g_ge,(model->sfg)[i]);
	  // separate counter for spatial axes - these have
	  // already been initialized, just need placement
	  int j=0;
	  for (idim=0; idim<(model->sfg)[i].gdim; idim++) {
	    if ((model->sfg)[i].axes[idim].id > EXTINT-1) {
	      dgs[i][idim]=g_gs[idim];
	      dge[i][idim]=g_ge[idim];
	    }
	    else {
	      if (((model->sfg)[i].axes[idim].id >-1) && 
		  ((model->sfg)[i].axes[idim].id < (model->sfg)[i].dim)) { 
		dgs[i][idim] = ls[j];
		dge[i][idim] = le[j];
		if ( crank[j] == cdims[j]-1 && 
		     gtype[i][j] == DUAL_GRID )
		  dge[i][idim] --;
		j++;
	      }
	    }
	  }
	}
	
      }
    }
#else
    for (i = 0;i < RDOM_MAX_NARR;i ++) {
      if (fd_isarr(i,*model,ic)) {
	if (fd_isdyn(i,ic)) {
	  /** first, assign base grid limits for spatial axes */
	  for (idim = 0;idim < ndim;idim ++) {
	    dgs[i][idim] = ls[idim]+fd_isdyn(i,ic);
	    dge[i][idim] = le[idim]-fd_isdyn(i,ic);
	    if (crank[idim] == 0 && gtype[i][idim] == DUAL_GRID)
	      dgs[i][idim] --;
	  }
	}
	else {
	  /** for static fields, assign for internal extended axes.
	      Note that, being internal, these axes are not localized by
	      domain decomposition 
	  */
	  IPNT g_gs;
	  IPNT g_ge;
	  get_gs(g_gs,(model->sfg)[i]);
	  get_ge(g_ge,(model->sfg)[i]);
	  //	  cerr<<"fd_setcompdom: static array "<<i<<endl;
	  //	  cerr<<"  dim="<<ndim<<" gdim="<<(model->sfg)[i].gdim<<endl;
	  // separate counter for spatial axes - these have
	  // already been initialized, just need placement
	  // presumes that all axes are either spatial or extended
	  int j=0;
	  for (idim=0; idim<(model->sfg)[i].gdim; idim++) {
	    //	    cerr<<"  axis "<<idim<<" id="<< (model->sfg)[i].axes[idim].id<<endl;
	    // internal extended axes - must be added to RARR
	    if ((model->sfg)[i].axes[idim].id > EXTINT-1) {
	      dgs[i][idim]=g_gs[idim];
	      dge[i][idim]=g_ge[idim];
	    }
	    else {
	      // spatial axes
	      if (((model->sfg)[i].axes[idim].id >-1) && 
		  ((model->sfg)[i].axes[idim].id < (model->sfg)[i].dim)) { 
		dgs[i][idim] = ls[j];
		dge[i][idim] = le[j];
		if (crank[idim] == 0 && gtype[i][idim] == DUAL_GRID)
		  dge[i][idim] --;
		j++;
	      }
	      // all other axes (neither spatial nor interal) are
	      // external extended axes and are not represented in the
	      // RARR - simply result in repeated loads
	    }
	    //	    cerr<<"  dgs["<<i<<"]["<<idim<<"]="<<dgs[i][idim]<<endl;
	    //	    cerr<<"  dge["<<i<<"]["<<idim<<"]="<<dge[i][idim]<<endl;
	  }
	}
      }  
    
      /*
      if (fd_isarr(i,*model,ic)) {
	fprintf(stderr,"ON EXIT FROM fd_setcompdom array index = %d:\n",i);
	fprintf(stderr,"  gs: ");
	for (int iv=0;iv<RARR_MAX_NDIM;iv++) fprintf(stderr,"%6d ",dgs[i][iv]);
	fprintf(stderr,"\n");
	fprintf(stderr,"  ge: ");
	for (int iv=0;iv<RARR_MAX_NDIM;iv++) fprintf(stderr,"%6d ",dge[i][iv]);
	fprintf(stderr,"\n");
      }
      */
    }
#endif

    return 0;

}

int fd_modelcrea(IPNT cdims, IPNT crank, PARARRAY * par, FILE * stream, IMODEL * model, IWaveInfo const & ic) {
    int err;
    int ndim, nnei, inei, idim, iv, i;
    IPNT ns;
    IPNT dgs[RDOM_MAX_NARR], dge[RDOM_MAX_NARR];    /*< computational domain */
    IPNT dgsa[RDOM_MAX_NARR], dgea[RDOM_MAX_NARR];  /*< allocated domain */
    IPNT gs, ge;
    IPNT dgsr[IWAVE_NNEI], dger[IWAVE_NNEI];                        /*< receive domains */
    IPNT dgsrs[IWAVE_NNEI][RDOM_MAX_NARR], dgers[IWAVE_NNEI][RDOM_MAX_NARR];   /*< all receive domains */
    int frcvempty[IWAVE_NNEI], rcvne;                       /*< empty receive flag */
  
    /*< grid type in each dimension: primal grid (=0) or dual grid
      (=1) */
    IPNT gtype[RDOM_MAX_NARR];         
    /** 
     * stencil dependence matrix:
     * if variable ir is dependent of variable ip, sten_dep_mat[ir][ip] is 
     * assigned a positive number which indicates the dependency type, e.g.,
     * DEP_F, DEP_DFDZ,DEP_DFDX,DEP_DFDY, etc.
     */
    //int  sten_dep_mat[RDOM_MAX_NARR][RDOM_MAX_NARR];
    /** Stencil - defines size and shape of all FD timestep updates */
    STENCIL sten;
   
    /*--------------------------------------------------------------------------*/
    /*-assign indices for send and receive domains------------------------------*/
    for (iv = 0;iv < IWAVE_NNEI;iv ++) {
	IASN(dgsr[iv], IPNT_1);
	IASN(dger[iv], IPNT_0);
	frcvempty[iv] = 0;
	for (i = 0;i < RDOM_MAX_NARR;i ++) {
	    IASN(dgsrs[iv][i], IPNT_1);
	    IASN(dgers[iv][i], IPNT_0);
	}
    }
    /*--------------------------------------------------------------------------*/
    /*-read grid info-----------------------------------------------------------*/
    /* moved to iwave_construct - 28.01.11 WWS 
       if ( (err=fdm->readgrid(par, stream, model)) ) {
       return err;
       }
    */
    // set spatial dimension
    ndim = model->g.dim;

    //    fprintf(stderr,"in fd_modelcrea: ndim=%d\n",ndim); 

    /*--------------------------------------------------------------------------*/
    /*-set nnei (num of neighbors in cart grid --*/
    if ( (err=im_setndim(model)) ) {
	fprintf(stream,"ERROR: fd_modelcrea from im_setndim, err=%d\n",err);
	fflush(stream);
	return err;
    }
    nnei = model->nnei;
    //    fprintf(stderr,"in fd_modelcrea: nnei=%d\n",nnei); 
    /*--------------------------------------------------------------------------*/
    /*-create stencil-----------------------------------------------------------*/
    if ( (err=fd_set_grid_type(stream,gtype,ic)) ) {
	fprintf(stream,"ERROR: fd_modelcrea from set_grid_type, err=%d\n",err);
	fflush(stream);
	return err;
    }

    //    cerr<<"->build_sten_dep"<<endl;
    /*
    if ( (err=fdm->build_sten_dep(stream,ndim,sten_dep_mat)) ) {
      fprintf(stream,"ERROR: fd_modelcrea from build_sten_dep, err=%d\n",err);
      fflush(stream);
      return err;
    }
    */
    //    cerr<<"fd_modelcrea->create_sten"<<endl;
    //    if ( (err=fdm->create_sten(fdm,stream,ndim,gtype,sten_dep_mat,&sten)) )  {
    //    if ( (err=ic.get_stencil()(model->specs,stream,ic,ndim,gtype,&sten)) )  {
    if ( (err=ic.get_stencil()(model->specs,stream,ndim,gtype,&sten)) )  {
      fprintf(stream,"ERROR: fd_modelcrea from create_sten, err=%d\n",err);
      fflush(stream);
      return err;
    }

    // cerr<<"finished with fdm\n";
    /* print out stencil if desired */
#ifdef VERBOSE
    //    sten_out(sten, stream, fd_ind2str);
    //    fprint_weqn(fdm,stream,sten_dep_mat);
#endif

    /*--------------------------------------------------------------------------*/
    /*-compute local computational grid size------------------------------------*/
    for ( idim = 0; idim < RDOM_MAX_NARR; ++idim ) {
      IASN(dgs[idim], IPNT_1);
      IASN(dge[idim], IPNT_0);
    }

    // cerr<<"setcompdom\n";
    if ((err=fd_setcompdom(stream, cdims, crank, model, dgs, dge, gtype, ic)) ) {
      fprintf(stream,"ERROR: fd_modelcrea from fd_setcompdom, err=%d\n",err);
      fflush(stream);
      return err;
    }
    /*--------------------------------------------------------------------------*/
    /*-declare computational domain---------------------------------------------*/
    // cerr<<"declare\n";
    err = rd_a_declare(&(model->ld_c), ic.get_num_fields(), dgs, dge);
    if ( err ) {
      fprintf(stream, "ERROR. fd_modelcrea from rd_a_declare err=%d\n", 
	      err);
      fflush(stream);
      return E_BADINPUT;
    }
    /*--------------------------------------------------------------------------*/
    /* ex_compute does two different things, which should probably be separated:*/
    /* 1. sets bounds for allocated arrays                                      */
    /* 2. computes bounds for send/receive virtual subarrays                    */
    /* both start with declared computational array bounds. stencil info used   */
    /* to set send/recv                                                         */
    /* this stuff works right if ex_compute is a no-op for static fields        */
    /*--------------------------------------------------------------------------*/
    
    for ( iv = 0; iv < RDOM_MAX_NARR; ++iv ) {
      IASN(dgsa[iv], IPNT_1);
      IASN(dgea[iv], IPNT_0);
      /*    fprintf(stderr,"iv=%d\n",iv); */
      if (fd_isarr(iv,*model,ic)) {
	if (fd_isdyn(iv,ic)) {
	  err = ex_compute(iv, &sten, &(model->ld_c), 
			   dgsa[iv], dgea[iv], dgsr, dger, frcvempty);
	  if ( err ) {
	    fprintf(stream, "ERROR. fd_modelcrea from ex_compute err=%d for array %d.\n", 
		    err, iv);
	    fflush(stream);
	    return E_INTERNAL;
	  }

	  rcvne = 0;
	  for ( i = 0; i < nnei; ++i ) if ( !(frcvempty[i]) ) ++rcvne;
	  
	  for ( i = 0; i < nnei; ++i ) {
	    /*      fprintf(stream,"iv=%d i=%d dgsr[i][0]=%d dger[i][0]=%d dgsr[i][1]=%d dger[i][1]=%d\n", */
	    /*	      iv, i, dgsr[i][0],dger[i][0],dgsr[i][1],dger[i][1]); */
	    IASN(dgsrs[i][iv], dgsr[i]);
	    IASN(dgers[i][iv], dger[i]);
	  }
	}
	else {
	  IASN(dgsa[iv], dgs[iv]);
	  IASN(dgea[iv], dge[iv]);
	}
      }

      /* TODO: change process receives below */
      
      /* model-specific alterations of any arrays */
      /*    fprintf(stderr,"fd_modelcrea->alter_dom,iv=%d\n",iv); */
      /* eliminated 11.13 - all non-spatial domain arrays in fdspecs
	 err = fdm->alter_dom(iv,dgsa[iv],dgea[iv]);
	 if (err) {
	 fprintf(stream, "ERROR. fd_modelcrea from alter_dom, err=%d\n", err);
	 fflush(stream);
	 return err;
	 }      
      */
      /*    fprintf(stderr,"bottom iv=%d\n",iv); */
    }
    
    /*--------------------------------------------------------------------------*/
    /*-allocate main domain, create computational domain------------------------*/
    // cerr<<"allocate\n";
    err = rd_a_create(&(model->ld_a), ic.get_num_fields(), dgsa, dgea);
    if ( err ) {
      fprintf(stream, "ERROR. fd_modelcrea from rd_a_create allocated domain err=%d.\n", err);
      fflush(stream);
      return E_INTERNAL;
    }

#ifdef IWAVE_VERBOSE
    fprintf(stream,"in modelcrea\n");
    fflush(stream);
#endif

    model->ld_c = model->ld_a;
    for (i = 0;i < RDOM_MAX_NARR;i ++) {
      if ( !(fd_isarr(i,*model,ic) && fd_isdyn(i,ic)) ) continue;
      err = rd_greset(&(model->ld_c), i, dgs[i], dge[i]);
      if ( err ) {
	fprintf(stream,"in modelcrea: reset failed for array %s err=%d\n", 
		fd_ind2str(i,ic), err);
	fflush(stream);
	return E_INTERNAL;
      }
    }
    /*--------------------------------------------------------------------------*/
    /*-create physical domain, here it includes the 2 bnd pnts------------------*/
    /*-not for sure, need advise from Dr. Symes!!!!!!---------------------------*/
    get_n(ns, model->g);
    model->ld_p = model->ld_c;
    /*  for (i = 0;i < fdm->getnarr();i ++) { */
    for (i = 0;i < ic.get_num_fields(); i ++) {
      if ( (err=rd_gse(&(model->ld_p), i, gs, ge)) ) {
	fprintf(stream, 
		"ERROR. fd_modelcrea from rd_gse physical array [%s] err=%d.\n", 
		/*              fdm->ind2str(iv), err); */
		fd_ind2str(i,ic), err);
	fflush(stream);
	return E_INTERNAL;
      }
      for (idim = 0;idim < ndim;idim ++) {
	/*      if (gtype[iv][idim] == PRIMAL_GRID) { */
	if (gtype[i][idim] == PRIMAL_GRID) {
#if INCLUDE_BOUNDARY_PNTS
	  /*< include left and right boundary pnts */
	  if ( gs[idim] < 0 )            gs[idim] = iwave_min(0, ge[idim] + 1);
	  if ( ge[idim] > ns[idim] - 1 ) ge[idim] = ns[idim] - 1;
#else
	  /*< not include left and right boundary pnts */
	  if ( gs[idim] < 1 )            gs[idim] = iwave_min(1, ge[idim] + 1);
	  if ( ge[idim] > ns[idim] - 2 ) ge[idim] = ns[idim] - 2;
#endif     
	  
	  if ( ge[idim] < gs[idim] )     ge[idim] = gs[idim] - 1;
	}
	/*      else if (gtype[iv][idim] == DUAL_GRID) { */
	else if (gtype[i][idim] == DUAL_GRID) {
	  if ( gs[idim] < 0 )            gs[idim] = iwave_min(0, ge[idim] + 1); 
	  if ( ge[idim] > ns[idim] - 2 ) ge[idim] = ns[idim] - 2;
	  if ( ge[idim] < gs[idim] )     ge[idim] = gs[idim] - 1;
	}
	else {
	  fprintf(stream, "ERROR. fd_modelcrea: undefined grid type: %d\n", 
		  /*                gtype[iv][idim]); */
		  gtype[i][idim]);
	  fflush(stream);
	  return E_INTERNAL;
	}
      }
      /*    if ( (err=rd_greset(&(model->ld_p), iv, gs, ge)) ) { */
      if ( (err=rd_greset(&(model->ld_p), i, gs, ge)) ) {
	fprintf(stream, 
		"ERROR. fd_modelcrea from rd_greset physical array [%s] err=%d.\n", 
		/*              fdm->ind2str(iv), err); */
		fd_ind2str(i,ic), err);
	fflush(stream);
	return E_INTERNAL;
      }
    }
    /*--------------------------------------------------------------------------*/
    /*-set virtual receive domains ---------------------------------------------*/
    for (inei = 0;inei < nnei;inei ++) {
      model->ld_r[inei] = model->ld_s[inei] = model->ld_a;
      for (i = 0;i < RDOM_MAX_NARR;i ++) {
	if (!(fd_isarr(i,*model,ic) && fd_isdyn(i,ic)))  continue;
	err = rd_greset(model->ld_r + inei, i, dgsrs[inei][i], dgers[inei][i]);
	if ( err ) {
	  fprintf(stream, 
		  "ERROR. fd_modelcrea from rd_greset\n receive array [%s] if domain (%d) greset error #%d.\n iv=%d dgsrs[0]=%d dgers[0]=%d dgsrs[1]=%d dgers[1]=%d\n",
		  /*		fdm->ind2str(iv), inei, err,iv,dgsrs[inei][iv][0],dgers[inei][iv][0],dgsrs[inei][iv][1],dgers[inei][iv][1]); */
		  fd_ind2str(i,ic), inei, err,i,dgsrs[inei][i][0],dgers[inei][i][0],dgsrs[inei][i][1],dgers[inei][i][1]);
	  fflush(stream);
	  return E_INTERNAL;
	}
      }
    }
    // cerr<<"set local grid\n";
    /*--------------------------------------------------------------------------*/
    /*-set local grid ----------------------------------------------------------*/
    /* must have a variable defined on PRIMAL grid in every axis */
    /*  for (i = 0;i < fdm->getnarr();i ++) { */
    for (i = 0;i < ic.get_num_fields();i ++) {
      /*    iv = fdm->getindices(i); */
      for (idim = 0;idim < ndim;idim ++) {
	/*      if (gtype[iv][idim] == DUAL_GRID)  break; */
	if (gtype[i][idim] == DUAL_GRID)  break;
      }
      if (idim == ndim){
	/*      err = rd_gse(&(model->ld_c), iv, gs, ge); */
	err = rd_gse(&(model->ld_c), i, gs, ge);
	
	if ( err ) {
	  fprintf(stream, "ERROR. fd_modelcrea from rd_gse allocated array [%s] err=%d.\n",
		  /*                fdm->ind2str(iv), err); */
		  fd_ind2str(i,ic), err);
	  fflush(stream);
	  return E_INTERNAL;
	}
	break;
      }
    }
    init_grid(&(model->gl),ndim,ndim);
    for ( idim = 0; idim < ndim; ++idim ) {
      model->gl.axes[idim].d = model->g.axes[idim].d;
      model->gl.axes[idim].o = model->g.axes[idim].o + model->gl.axes[idim].d * gs[idim];
      model->gl.axes[idim].n = ge[idim]-gs[idim]+1;
      model->gl.axes[idim].id = model->g.axes[idim].id;
    }
#ifdef IWAVE_VERBOSE
    fprintf(stream,"NOTE: local grid used to determine trace sampling:\n");
    fprint_grid(stream,model->gl);
#endif
    /*--------------------------------------------------------------------------*/
    
    //    cerr<<"fd_modelcrea->deallocate stencil\n";
    sten_destroy(&sten);
    
    //    cerr<<"leaving fd_modelcrea\n";
    return 0;
}




