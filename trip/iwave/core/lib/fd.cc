/**
 * XW 01-18-2010
 * wangxin.tom@gmail.com
 * fd.c
 */

#include "fd.h"

/* helper functions */

int fd_isdyn(int i, IWaveInfo const & ic) {
  if (i<0 || i>=ic.get_num_fields()) return 0;
  return ic.get_iwave_fields()[i].dynamic;
}

int fd_isarr(int i, IWaveInfo const & ic) {
  if (i >= ic.get_num_fields()) return false;
  return true;
}

const char* fd_ind2str(int i, IWaveInfo const & ic) {
  if (fd_isarr(i,ic)) return ic.get_iwave_fields()[i].field.c_str();
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
      
int fd_set_grid_type(FILE * stream, int ndim, IPNT gtype[RDOM_MAX_NARR], IWaveInfo const & ic) {
  for (int i=0;i<ic.get_num_fields();i++) IASN(gtype[i],ic.get_iwave_fields()[i].gtype);
  for (int i=ic.get_num_fields();i<RDOM_MAX_NARR;i++) IASN(gtype[i],IPNT_0);
  return 0;
}

int fd_readgrid(PARARRAY * pars, 
		FILE *stream, 
		IMODEL * model,
		std::string gridkey) {

  std::string gridfile="";    /* workspace for filename */
  int err=0;
  
  // extract filename from 0 entry in field table
  if (!parse(*pars,gridkey,gridfile)) {
    //    cerr<<"error return err="<<err<<endl;
    fprintf(stream,"Error: fd_readgrid from ps_flcstring\n");
    fprintf(stream,"  failed to parse key = %s\n",gridkey.c_str());
    ps_printall(*pars,stream);
    fflush(stream);
    return E_BADINPUT;
  }
  else {
    // mod 30.11.13: read only spatial axes for primal grid - ignore others
    grid gtmp;
    init_default_grid(&gtmp);
    init_default_grid(&(model->g));
    err=read_grid(&gtmp,gridfile.c_str(),stream);
    if (err==0) {
      (model->g).gdim = gtmp.gdim;
      (model->g).dim  = gtmp.dim;
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
    else {
      fprintf(stream,"Error: fd_readgrid from read_grid\n");
      fprintf(stream,"  err=%d\n",err);
      return err;
    }
  }

  /* until such time as someone implements PML...*/
  IASN(model->nls,IPNT_0);
  IASN(model->nrs,IPNT_0);

  ireal dx;
  for (int idim=0; idim<model->g.dim; idim++ ) {
    if ((dx=model->g.axes[idim].d)<REAL_EPS) {
      fprintf(stream, "dx[%d] = %f\n", idim,dx);
      return E_BADINPUT;
    }
    stringstream key;
    key.str()="nl";
    key<<idim+1;
    //    snprintf(key,kl,"nl%d",idim+1);
    ireal tmp=0.0;
    ps_flreal(*pars, key.str().c_str(), &tmp);
    model->nls[idim]=iwave_max(0,(int)((ceil)(tmp/dx)));
    //    snprintf(key,kl,"nr%d",idim+1);
    key.str()="nr";
    key<<idim+1;
    tmp=0.0;
    ps_flreal(*pars,key.str().c_str(),&tmp);
    model->nrs[idim]=iwave_max(0,(int)((ceil)(tmp/dx)));
  }
  
  return err;
}

int fd_setcompdom(FILE * stream, IPNT cdims, IPNT crank, 
                  IMODEL * model, IPNT dgs[], IPNT dge[],
		  /*		  int m_size, */
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
    /*  int iv, i, idim; */
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
#if INCLUDE_BOUNDARY_PNTS
    /*< include left and right boundary pnts */   
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

    /*  for (i = 0;i < m_size;i ++) { */
    for (i = 0;i < RDOM_MAX_NARR;i ++) {
	/*    iv = getindices(i); */
	/*    if (!fd_isarr(iv)) continue; */
      if (!fd_isarr(i,ic)) continue;
	for (idim = 0;idim < ndim;idim ++) {
	    dgs[i][idim] = ls[idim];
	    dge[i][idim] = le[idim];
	    if ( crank[idim] == cdims[idim]-1 && 
		 gtype[i][idim] == DUAL_GRID )
		dge[i][idim] --;
	}
    }
#else
    /*< not include left and right boundary pnts */
    /*< ls, le: local start, end for primal grid */
    for (idim = 0;idim < ndim;idim ++) {
	gs_pa[idim] = -model->nls[idim]+1;
	ge_pa[idim] = ns[idim] + model->nrs[idim] - 2;
	gn_pa[idim] = ns[idim] + model->nls[idim] + model->nrs[idim] - 2;
	ls[idim] = gn_pa[idim] * (long)crank[idim];
	ls[idim] = gs_pa[idim] + ls[idim] / (long)cdims[idim];
	le[idim] = gn_pa[idim] * (long)(crank[idim] + 1);
	le[idim] = gs_pa[idim] + le[idim] / (long)cdims[idim]-1;
	if (le[idim] < ls[idim]) {
	    fprintf(stream, "Error: in fd_setcompdom: le[%d] = %d < ls[%d] = %d\n",
		    idim, le[idim], idim, ls[idim]);
	    return E_INTERNAL;
	}
	//	cerr<<"fd_setcompdom: ls["<<idim<<"]="<<ls[idim]<<" le["<<idim<<"]="<<le[idim]<<"\n";
    } 

    for (i = 0;i < RDOM_MAX_NARR;i ++) {
	/*    iv = getindices(i); */
	/*    fprintf(stream,"fd_setcompdom - array %d index %d\n",i,iv); */
	/*    if (isarr(iv)) {  */
      if (fd_isarr(i,ic)) {
	    /*      fprintf(stream,"fd_setcompdom - assigning array %d ",iv); */
	    /*      fprintf(stream,"fd_setcompdom - assigning array %d ",i);  */
	    for (idim = 0;idim < ndim;idim ++) {
		/*	dgs[iv][idim] = ls[idim]; */
		dgs[i][idim] = ls[idim];
		/*	dge[iv][idim] = le[idim]; */
		dge[i][idim] = le[idim];
		if (crank[idim] == 0 && 
		    /*	    gtype[iv][idim] == DUAL_GRID) */
		    gtype[i][idim] == DUAL_GRID)
		    /*	  dgs[iv][idim] --; */
		    dgs[i][idim] --;
		/*	fprintf(stream,"dgs[%d]=%d, dge[%d]=%d",idim,dgs[iv][idim],idim,dge[iv][idim]); */
		/*	fprintf(stream,"dgs[%d]=%d, dge[%d]=%d",idim,dgs[i][idim],idim,dge[i][idim]); */
	    }
	    /*      fprintf(stream,"\n"); */
	}   
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
    if ( (err=fd_set_grid_type(stream,ndim,gtype,ic)) ) {
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
    //    cerr<<"->create_sten"<<endl;
    //    if ( (err=fdm->create_sten(fdm,stream,ndim,gtype,sten_dep_mat,&sten)) )  {
    //    if ( (err=ic.get_stencil()(model->specs,stream,ic,ndim,gtype,&sten)) )  {
    if ( (err=ic.get_stencil()(model->specs,stream,ndim,gtype,&sten)) )  {
      fprintf(stream,"ERROR: fd_modelcrea from create_sten, err=%d\n",err);
      fflush(stream);
      return err;
    }

    //    cerr<<"finished with fdm\n";
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

    if ((err=fd_setcompdom(stream, cdims, crank, model, dgs, dge, gtype, ic)) ) {
	fprintf(stream,"ERROR: fd_modelcrea from fd_setcompdom, err=%d\n",err);
	fflush(stream);
	return err;
    }
    /*--------------------------------------------------------------------------*/
    /*-declare computational domain---------------------------------------------*/
    err = rd_a_declare(&(model->ld_c), RDOM_MAX_NARR, ndim, dgs, dge);
    if ( err ) {
	fprintf(stream, "ERROR. fd_modelcrea from rd_a_declare err=%d\n", 
		err);
	fflush(stream);
	return E_BADINPUT;
    }
    /*--------------------------------------------------------------------------*/
    /*-compute send and receive domains------------------------------------*/
    for ( idim = 0; idim < RDOM_MAX_NARR; ++idim ) {
	IASN(dgsa[idim], IPNT_1);
	IASN(dgea[idim], IPNT_0);
    }
    for ( iv = 0; iv < RDOM_MAX_NARR; ++iv ) {
	/*    fprintf(stderr,"iv=%d\n",iv); */
	err = ex_compute(iv, &sten, &(model->ld_c), 
			 dgsa[iv], dgea[iv], dgsr, dger, frcvempty);
	if ( err ) {
	    fprintf(stream, "ERROR. fd_modelcrea from ex_compute err=%d for array [%s].\n", 
		    err, (fd_ind2str)(iv,ic));
	    fflush(stream);
	    return E_INTERNAL;
	}
	/* TODO: change process receives below */
	/* First, check that only P and V are received */
	rcvne = 0;
	for ( i = 0; i < nnei; ++i ) if ( !(frcvempty[i]) ) ++rcvne;
            
	/* Second, store receive domains */
	for ( i = 0; i < nnei; ++i ) {
	    /*      fprintf(stream,"iv=%d i=%d dgsr[i][0]=%d dger[i][0]=%d dgsr[i][1]=%d dger[i][1]=%d\n", */
	    /*	      iv, i, dgsr[i][0],dger[i][0],dgsr[i][1],dger[i][1]); */
	    IASN(dgsrs[i][iv], dgsr[i]);
	    IASN(dgers[i][iv], dger[i]);
	}

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
    err = rd_a_create(&(model->ld_a), RDOM_MAX_NARR, ndim, dgsa, dgea);
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
      //	if ( !(isdyn(fdm,i)) ) continue;
      if ( !(fd_isdyn(i,ic)) ) continue;
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
    for (i = 0;i < RDOM_MAX_NARR;i ++) {
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
	  //	    if (!(isdyn(fdm,i)))  continue;
	  if (!(fd_isdyn(i,ic)))  continue;
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
    /*--------------------------------------------------------------------------*/
    /*-set local grid ----------------------------------------------------------*/
    /* must have a variable defined on PRIMAL grid in every axis */
    /*  for (i = 0;i < fdm->getnarr();i ++) { */
    for (i = 0;i < RDOM_MAX_NARR;i ++) {
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

    /* deallocate stencil */
    sten_destroy(&sten);

    //    cerr<<"leaving modelcrea\n";
    return 0;
}




