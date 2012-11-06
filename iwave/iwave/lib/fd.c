/**
 * XW 01-18-2010
 * wangxin.tom@gmail.com
 * fd.c
 */

#include "fd.h"

/* helper functions */

/* returns 1 if i is index of dynamic array in FD_MODEL, else 0 */
int isdyn(FD_MODEL *fd, int i) {
    int iv;
    int niv=fd->numsubsteps();
    for (iv=0;iv<niv;iv++)
	if (fd->update(i,iv)) return 1;
    return 0;
}

/* returns number of dynamic arrays */
int narr(FD_MODEL * fd) {
    int n=0;
    int nia;
    int ia, iv;
    int niv=fd->numsubsteps();
    for (ia=0;ia<RDOM_MAX_NARR;ia++) {
	nia=0;
	for (iv=0;iv<niv;iv++) nia = nia || fd->update(ia,iv);
	n+=nia;
    }
    return n;
}

int fd_setcompdom(FILE * stream, IPNT cdims, IPNT crank, 
                  IMODEL * model, IPNT dgs[], IPNT dge[],
		  /*		  int m_size, */
		  IPNT gtype[RDOM_MAX_NARR],
		  int (*isarr)(int)) {
    /*		  int (*isarr)(int), */
    /*		  int (*getindices)(int)) { */
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
	/*    if (!isarr(iv)) continue; */
	if (!isarr(i)) continue;
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
    } 

    for (i = 0;i < RDOM_MAX_NARR;i ++) {
	/*    iv = getindices(i); */
	/*    fprintf(stream,"fd_setcompdom - array %d index %d\n",i,iv); */
	/*    if (isarr(iv)) {  */
	if (isarr(i)) {
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

/* int fprint_weqn(FILE * stream, int sten_dep_mat[RDOM_MAX_NARR][RDOM_MAX_NARR], char * (*getnames)(int), int (*getindices)(int), int m_size) { */ 
int fprint_weqn(FD_MODEL * fd,
		FILE * stream, 
		int sten_dep_mat[RDOM_MAX_NARR][RDOM_MAX_NARR]) {

    /*< try to print the equations */

    /*  int i, j, ir, ip; */
    int i, j;
    int op;

    fprintf(stream, "\nWAVE EQUATIONS:\n");
    for (i = 0;i < RDOM_MAX_NARR;i ++) {
	if (!isdyn(fd,i)) continue;
	fprintf(stream, "d %4s/ dt    =    ", fd->ind2str(i));
	for (j = 0;j < RDOM_MAX_NARR;j ++) {
	    op = sten_dep_mat[i][j];
	    switch (op) {
		case DEP_F:
		    fprintf(stream, "%10s    ", fd->ind2str(j));
		    break;
		case DEP_DFDZ:
		    fprintf(stream, "d %4s/ dz    ", fd->ind2str(j));
		    break;
		case DEP_DFDX:
		    fprintf(stream, "d %4s/ dx    ", fd->ind2str(j));
		    break;
		case DEP_DFDY:
		    fprintf(stream, "d %4s/ dy    ", fd->ind2str(j));
		    break;
		default:
		    break;
	    }
	}
	fprintf(stream,"\n");
    }  
    return 0;
}

void fd_model_setnull(FD_MODEL * fd) {
    fd->fdpars=NULL;
    /* leave other (statically allocated) data members uninitialised */
    fd->readgrid=NULL;
    fd->readtimegrid=NULL;
    fd->readschemeinfo=NULL;
    fd->set_grid_type=NULL;
    fd->build_sten_dep=NULL;
    fd->create_sten=NULL;
    fd->readmedia=NULL;
    fd->isarr=NULL;
    fd->ind2str=NULL;
    fd->numsubsteps=NULL;
    fd->update=NULL;
    fd->parcopy=NULL;
    fd->fd_model_init=NULL;
    fd->fd_model_dest=NULL;
}

int fd_modelcrea(IPNT cdims, IPNT crank, PARARRAY * par, FILE * stream, IMODEL * model) {
    FD_MODEL *fdm = (FD_MODEL *)model->specs;
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
    int  sten_dep_mat[RDOM_MAX_NARR][RDOM_MAX_NARR];
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

    /*  fprintf(stderr,"in fd_modelcrea: ndim=%d\n",ndim); */

    /*--------------------------------------------------------------------------*/
    /*-set nnei (num of neighbors in cart grid --*/
    if ( (err=im_setndim(model)) ) {
	fprintf(stream,"ERROR: fd_modelcrea from im_setndim, err=%d\n",err);
	fflush(stream);
	return err;
    }
    nnei = model->nnei;
    /*--------------------------------------------------------------------------*/
    /*-create stencil-----------------------------------------------------------*/
    if ( (err=fdm->set_grid_type(stream,ndim,gtype)) ) {
	fprintf(stream,"ERROR: fd_modelcrea from set_grid_type, err=%d\n",err);
	fflush(stream);
	return err;
    }

    if ( (err=fdm->build_sten_dep(stream,ndim,sten_dep_mat)) ) {
	fprintf(stream,"ERROR: fd_modelcrea from build_sten_dep, err=%d\n",err);
	fflush(stream);
	return err;
    }
  
    if ( (err=fdm->create_sten(fdm,stream,ndim,gtype,sten_dep_mat,&sten)) )  {
	fprintf(stream,"ERROR: fd_modelcrea from create_sten, err=%d\n",err);
	fflush(stream);
	return err;
    }

    /* print out stencil if desired */
#ifdef VERBOSE
    sten_out(sten, stream, fdm->ind2str);
    fprint_weqn(fdm,stream,sten_dep_mat);
#endif

    /*--------------------------------------------------------------------------*/
    /*-compute local computational grid size------------------------------------*/
    for ( idim = 0; idim < RDOM_MAX_NARR; ++idim ) {
	IASN(dgs[idim], IPNT_1);
	IASN(dge[idim], IPNT_0);
    }
    if ( (err=fd_setcompdom(stream, cdims, crank, model, dgs, dge, gtype,fdm->isarr)) ) {
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
		    err, (fdm->ind2str)(iv));
	    fflush(stream);
	    return E_INTERNAL;
	}
	/* TODO: change process receives below */
	/* First, check that only P and V are received */
	rcvne = 0;
	for ( i = 0; i < nnei; ++i ) if ( !(frcvempty[i]) ) ++rcvne;
            
	/* check which array exchanged */
   
	if ( rcvne > 0 ) {
	    for ( idim = 0; idim < ndim; ++idim ) 
		if ( isdyn(fdm,iv) ) break;
	    if ( idim == ndim ) {
		fprintf(stream, "ERROR. fd_modelcrea: wrong array to be exchanged [%s].\n", 
			fdm->ind2str(iv));
		fflush(stream);
		return E_INTERNAL;
	    }
	}

	/* Second, store receive domains */
	for ( i = 0; i < nnei; ++i ) {
	    /*      fprintf(stream,"iv=%d i=%d dgsr[i][0]=%d dger[i][0]=%d dgsr[i][1]=%d dger[i][1]=%d\n", */
	    /*	      iv, i, dgsr[i][0],dger[i][0],dgsr[i][1],dger[i][1]); */
	    IASN(dgsrs[i][iv], dgsr[i]);
	    IASN(dgers[i][iv], dger[i]);
	}

	/* model-specific alterations of any arrays */
	/*    fprintf(stderr,"fd_modelcrea->alter_dom,iv=%d\n",iv); */
	err = fdm->alter_dom(iv,dgsa[iv],dgea[iv]);
	if (err) {
	    fprintf(stream, "ERROR. fd_modelcrea from alter_dom, err=%d\n", err);
	    fflush(stream);
	    return err;
	}      
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

    fprintf(stream,"in modelcrea\n");
    fflush(stream);

    model->ld_c = model->ld_a;
    for (i = 0;i < RDOM_MAX_NARR;i ++) {
	if ( !(isdyn(fdm,i)) ) continue;
	err = rd_greset(&(model->ld_c), i, dgs[i], dge[i]);
	if ( err ) {
	    fprintf(stream, 
		    fdm->ind2str(i), err);
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
		    fdm->ind2str(i), err);
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
		    fdm->ind2str(i), err);
	    fflush(stream);
	    return E_INTERNAL;
	}
    }
    /*--------------------------------------------------------------------------*/
    /*-set virtual receive domains ---------------------------------------------*/
    for (inei = 0;inei < nnei;inei ++) {
	model->ld_r[inei] = model->ld_s[inei] = model->ld_a;
	for (i = 0;i < RDOM_MAX_NARR;i ++) {
	    if (!(isdyn(fdm,i)))  continue;
	    err = rd_greset(model->ld_r + inei, i, dgsrs[inei][i], dgers[inei][i]);
	    if ( err ) {
		fprintf(stream, 
			"ERROR. fd_modelcrea from rd_greset\n receive array [%s] if domain (%d) greset error #%d.\n iv=%d dgsrs[0]=%d dgers[0]=%d dgsrs[1]=%d dgers[1]=%d\n",
			/*		fdm->ind2str(iv), inei, err,iv,dgsrs[inei][iv][0],dgers[inei][iv][0],dgsrs[inei][iv][1],dgers[inei][iv][1]); */
			fdm->ind2str(i), inei, err,i,dgsrs[inei][i][0],dgers[inei][i][0],dgsrs[inei][i][1],dgers[inei][i][1]);
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
			fdm->ind2str(i), err);
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

    return 0;
}

int fd_modeldest(IMODEL * model) {
    FD_MODEL * fdm = (FD_MODEL *)(model->specs);
    return fdm->fd_model_dest(model);
}

int fd_isdyn(IMODEL * model, int i) {
    FD_MODEL * fdm = (FD_MODEL *)(model->specs);
    return isdyn(fdm,i);
}

int fd_modelinfo(FILE * stream, IMODEL * model) {
    return 0;
}

int fd_mread(FILE * stream, IMODEL * model, PARARRAY * par, int panelindex) {

    FD_MODEL *fdm = (FD_MODEL *)model->specs;
    int err;

    err = fdm->readmedia(par, stream, model, panelindex);
    if (err) {
	fprintf(stream, "Error: error in fd_mread after calling readmedia\n");
	return err;
    }

    err = fdm->readtimegrid(par, stream, model);
    if (err) {
	fprintf(stream, "Error: error in fd_mread after calling readtimegrid\n");
	return err;
    }
    if ((model->tsind).dt <= REAL_ZERO) {
	fprintf(stream, "Error: bad input: wrong time step dt=%g\n", 
		(model->tsind).dt);
	return E_BADINPUT;
    }
#ifdef VERBOSE
    fprintf(stderr, "dt = %g\n", (model->tsind).dt);
#endif
  
    /*  fprintf(stderr,"mread->schemeinfo, ndim=%ld\n",(model->g).dim); */
    /*  err = fdm->readschemeinfo(par, stream, fdm, model); */
    /* changed interface 13.02.11 WWS */
    err = fdm->readschemeinfo(par, stream, model);
    if (err) {
	fprintf(stream, "Error: error in fd_mread from readschemeinfo\n");
	return err;
    }
  
    return 0;
}

