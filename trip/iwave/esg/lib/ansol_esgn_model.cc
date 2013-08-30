#include "iwave.h"
#include "fd.h"
#include "esgn_read.h"

#include "esgn.h" /* needed to call some esg functions */ 
#include "ansol_esgn.h"

#include <math.h>

// #define IWAVE_VERBOSE


/*------------------- private data declarations ----------------------------*/

static int m_ndim = 0; /* dimension - need only compute once */
static int m_size = 24;

static const char *m_names[] = {"pz", "mp00", "mp01", "vz", "rhoz", "szx", "ms0", 
                         "px", "vx", "rhox", "sxy", "ms1",
                         "py", "vy", "rhoy", "szy", "ms2",
                         "eta_p0", "eta_v0", "eta_p1", "eta_v1",
                         "eta_p2", "eta_v2", "NONE"};


/*--- private function declarations - assigned to FD_MODEL pointers ----------*/

int  ansol_esg_numsubsteps();
int  ansol_esg_update(int,int);
int  ansol_esg_readschemeinfo(PARARRAY *, FILE *, IMODEL *);
int  ansol_esg_create_sten(FD_MODEL *, FILE *, int, IPNT[RDOM_MAX_NARR], int[RDOM_MAX_NARR][RDOM_MAX_NARR], STENCIL *);
int  ansol_HI_esg_modeldest(IMODEL *);				//HI specific
void ansol_HI_esgn_ts_parcopy(void * tgt, const void * src);	//HI specific


//HI specific!
/*----------------------------------------------------------------------------*/
int ansol_HI_esg_modelinit( PARARRAY *pars, FILE *stream, IMODEL * sol_mdl ) {
/*----------------------------------------------------------------------------*/
#ifdef IWAVE_VERBOSE
  	fprintf(stream,">>> Inside ansol_HI_esg_modelinit\n");
	fflush(stream);
#endif
  
  	int err=0;

  	ANSOL_ESG_PARS * ansolpars;	//model parameters
	RICKER_INFO    * rickerinfo;	//source information for Ricker source
	ESG_HIMED_INFO * himedinfo; 	//medium information for homogeneous isotropic medium
  	FD_MODEL * sol;       		//function struct

  	IPNT cdims;
  	IPNT crank;

#ifdef IWAVE_USE_MPI
    	IPNT cpers;          /* workspace for periodic wrap info  - currently not used */
#endif

  	char * jnk;          /* workspace */
	    
  	/* decode dimension - only on rk 0 */
#ifdef IWAVE_USE_MPI
  	int rk; 
  	MPI_Comm cm;
#endif
    	grid g;

	//allocating for FD_MODEL struct sol
	sol = (FD_MODEL *)usermalloc_(sizeof(FD_MODEL));
  	if (sol==NULL) return E_ALLOC;
	sol_mdl->specs=(void *)sol;

	//allocating for ANSOL_ESG_PARS struct
  	ansolpars = (ANSOL_ESG_PARS*)usermalloc_(sizeof(ANSOL_ESG_PARS));
	if ( ansolpars == NULL ) { 
    		err=E_ALLOC;
    		fprintf(stream,"ERROR: ansol_HI_esg_modelinit\n");
    		fprintf(stream,"failed to allocate ANSOL_ESG_PARS object\n");
    		return err;
  	}
	//fd_mdl is linked in ansol_HI_esg_postinit.

	//allocating for RICKER_INFO
	rickerinfo = (RICKER_INFO *)usermalloc_(sizeof(RICKER_INFO));
	if (rickerinfo==NULL) {
    		err=E_ALLOC;
    		fprintf(stream,"ERROR: ansol_HI_esg_modelinit\n");
    		fprintf(stream,"failed to allocate RICKER_INFO object\n");
    		return err;
  	}
	ansolpars->srcinfo = rickerinfo;
	//allocating for ESG_HIMED_INFO
	himedinfo = (ESG_HIMED_INFO *)usermalloc_(sizeof(ESG_HIMED_INFO));
	if (himedinfo==NULL) {
    		err=E_ALLOC;
    		fprintf(stream,"ERROR: ansol_HI_esg_modelinit\n");
    		fprintf(stream,"failed to allocate ESG_HIMED_INFO object\n");
    		return err;
  	}
	ansolpars->medinfo = himedinfo;

  	m_ndim=0;

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

		if (!init_elastic_geom_par(&g,*pars,stream)) m_ndim=g.dim; //from esgn_read.c  
			else {
				err=E_INTERNAL;
				fprintf(stream,"ERROR: in ansol_HI_esgn_modelinit - failed to read spatial geometry\n");
				fflush(stream);
				abortexit(err,pars,&stream);
		}

#ifdef IWAVE_USE_MPI
    	}
	MPI_Bcast(&m_ndim,1,MPI_INT,0,cm);
#endif

  	if (m_ndim<1 || m_ndim>RARR_MAX_NDIM) {
    		err=E_INTERNAL;
    		fprintf(stream,"ERROR: in ansol_HI_esg_modelinit - failed to read dim=%d\n",m_ndim);
    		return err;
  	}

  	/* decode order - take care of deprecated cases */
  	ansolpars->k=1;
 	ps_flint(*pars,"order",&(ansolpars->k));

#ifdef IWAVE_VERBOSE
    	fprintf(stream,"NOTE: initializing ANSOL_ESG with half-order = %d\n",ansolpars->k);
#endif

  	/* assign param object pointer */
	sol->fdpars 		= (void*)ansolpars;		//from ansol_esgn.h

	sol->isarr 		= esg_isarr;			//from esgn.h
	sol->set_grid_type 	= esg_set_grid_type;		//from esgn.h
	sol->build_sten_dep 	= esg_build_sten_dep;		//from esgn.h
	sol->ind2str 		= esg_ind2str;			//from esgn.h
	sol->alter_dom 		= esg_alter_dom;		//from esgn.h
	sol->readtimegrid 	= esg_readtimegrid; 		//from esgn_read.h
	sol->readmedia		= esgn_readmedia;		//from esgn_read.h
	sol->readgrid 		= esg_readgrid;		 	//from esgn_read.h

	sol->numsubsteps 	= ansol_esg_numsubsteps;
	sol->update 		= ansol_esg_update;
	sol->readschemeinfo 	= ansol_esg_readschemeinfo;
	sol->create_sten 	= ansol_esg_create_sten;
	sol->parcopy 		= ansol_HI_esgn_ts_parcopy; 	//HI specific
	sol->fd_model_init 	= ansol_HI_esg_modelinit;	//HI specific
	sol->fd_model_dest 	= ansol_HI_esg_modeldest;	//HI specific
	
	/* choose time step */
	sol->tsf = ansol_HI_esg_step; //HI spacific

  	return 0;
}

/*----------------------------------------------------------------------------*/
int ansol_esg_numsubsteps() {
/*----------------------------------------------------------------------------*/
  	return 2; //equal to esg_numsubsteps()
}

/*----------------------------------------------------------------------------*/
int ansol_esg_update(int ia, int iv) {
/*----------------------------------------------------------------------------*/
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
		if (m_ndim>0){
			if (ia==D_P0) return 1;
		}
	}
	if (iv==1) {
		if (m_ndim>2 && ia==D_V2) return 1;
		if (m_ndim>1 && ia==D_V1) return 1;
		if (m_ndim>0 && ia==D_V0) return 1;
	}
	return 0;
}

/*----------------------------------------------------------------------------*/
void ansol_HI_esgn_ts_parcopy(void * tgt, const void * src) {
/*----------------------------------------------------------------------------*/
  	ANSOL_ESG_PARS * ptgt = (ANSOL_ESG_PARS *)tgt;
  	const ANSOL_ESG_PARS * psrc = (const ANSOL_ESG_PARS *)src;
	int i,j;

	//copying members of ANSOL_ESG_PARS
	if (!ptgt){
		fprintf(stderr,"ERROR: ansol_HI_esgn_ts_parcopy, copying into empty target. ABORT!\n");
		exit(1);
	}
	if (!psrc){
		ptgt->dt   = psrc->dt;
		ptgt->k    = psrc->k;
		ptgt->ndim = psrc->ndim;
		ptgt->t_off= psrc->t_off;
	
		for (i=0; i<RARR_MAX_NDIM; i++){
			ptgt->lam[i]     = psrc->lam[i];
			ptgt->dx[i]      = psrc->dx[i];
			ptgt->o_coord[i] = psrc->o_coord[i];
			ptgt->o_index[i] = psrc->o_index[i];
			ptgt->src_x[i]   = psrc->src_x[i];
	
			for (j=0; j<RDOM_MAX_NARR; j++){
				ptgt->gtype[j][i] = psrc->gtype[j][i];
			}
		}
		//copying physical domain link
		ptgt->link_mdl = psrc->link_mdl;

		//copying members of RICKER_INFO
		RICKER_INFO * src_ptgt = (RICKER_INFO *)ptgt->srcinfo;
		const RICKER_INFO * src_psrc = (const RICKER_INFO *)psrc->srcinfo;
		
		if (!src_ptgt){
			fprintf(stderr,"ERROR: ansol_HI_esgn_ts_parcopy, copying into empty target. ABORT!\n");
			exit(1);
		}
		if (!src_psrc){
			src_ptgt->fpeak = src_psrc->fpeak;
			src_ptgt->amp   = src_psrc->amp;
			for (i=0; i<RARR_MAX_NDIM; i++)
				src_ptgt->src_d[i] = src_psrc->src_d[i];
		}
		else{
			fprintf(stderr,"WARNING: ansol_HI_esgn_ts_parcopy, no copying due to empty source.\n");
		}
	
		//copying members of ESG_HIMED_INFO
		ESG_HIMED_INFO * med_ptgt = (ESG_HIMED_INFO *)ptgt->medinfo;
		const ESG_HIMED_INFO * med_psrc = (const ESG_HIMED_INFO *)psrc->medinfo;
		if (!med_ptgt){
			fprintf(stderr,"ERROR: ansol_HI_esgn_ts_parcopy, copying into empty target. ABORT!\n");
			exit(1);
		}
		if (!med_psrc){
			med_ptgt->mu     = med_psrc->mu;
			med_ptgt->lambda = med_psrc->lambda;
			med_ptgt->rho    = med_psrc->rho;
			med_ptgt->alpha  = med_psrc->alpha;
			med_ptgt->beta   = med_psrc->beta;
		}
		else{
			fprintf(stderr,"WARNING: ansol_HI_esgn_ts_parcopy, no copying due to empty source.\n");
		}
	}
	else{
		fprintf(stderr,"WARNING: ansol_HI_esgn_ts_parcopy, no copying due to empty source.\n");
	}

	return;
}

/*----------------------------------------------------------------------------*/
int ansol_esg_create_sten( FD_MODEL * sol,
		           FILE     * stream, 
		           int        ndim, 
		           IPNT       gtype[RDOM_MAX_NARR], 
		           int        sten_dep_mat[RDOM_MAX_NARR][RDOM_MAX_NARR], 
		           STENCIL  * sten) {
/*----------------------------------------------------------------------------*/
#ifdef IWAVE_VERBOSE
	fprintf(stream,">>> Inside ansol_esg_create_sten\n");
	fflush(stream);
#endif

  	ANSOL_ESG_PARS * ansolpars = (ANSOL_ESG_PARS *)(sol->fdpars);
	int err=create_sten2_2k(sol,stream, ansolpars->k, ndim, gtype, sten_dep_mat,sten);
	
	if (!err){
		fprintf(stream,"ERROR: in ansol_esg_create_sten from create_sten2_2k\n");
	}
  	return err; 
}  		

/*----------------------------------------------------------------------------*/
int ansol_esg_readschemeinfo(PARARRAY * par, FILE * stream, IMODEL * model) {
/*----------------------------------------------------------------------------*/
	//setting some of the parameters of ANSOL_ESG_PARS struct

	FD_MODEL * sol = (FD_MODEL *)(model->specs); 
	ANSOL_ESG_PARS * ansolpars = (ANSOL_ESG_PARS *)(sol->fdpars);

	int ndim, idim;

	ansolpars->dt = (model->tsind).dt;     //setting dt
	ansolpars->ndim = ndim = (model->g).dim; //setting ndim

	if (ndim < 2 || ndim > 3) {
		fprintf(stderr, "Error: bad input in ansol_esg_build_scheme: ndim = %d\n",ndim);
		fprintf(stderr, "       ansol_esg only provides 2d/3d elastic wave solvers\n");
	return E_BADINPUT;
	}

	get_d (ansolpars->dx,     model->g); //setting dx
	get_o (ansolpars->o_coord,model->g); //setting o_coord
	get_gs(ansolpars->o_index,model->g); //setting o_index
	
	for (idim = 0; idim < ndim; idim ++) {
		if ( ansolpars->dx[idim] <= 0.0) {
			fprintf(stream, "Error: bad input: wrong grid space step, dim=%d, step=%g\n",
				idim, ansolpars->dx[idim]);
			return E_BADINPUT;
		}
		ansolpars->lam[idim] = ansolpars->dt / ansolpars->dx[idim]; //setting lam
#ifdef IWAVE_VERBOSE
		fprintf(stderr, "lam[%d] = %g\n", idim, ansolpars->lam[idim]);
#endif 
	}

	esg_set_grid_type( stream, ndim, ansolpars->gtype ); //setting gtype

	return 0;
}

/*----------------------------------------------------------------------------*/
int ansol_HI_esg_modeldest(IMODEL * model) {
/*----------------------------------------------------------------------------*/
  	FD_MODEL * sol = (FD_MODEL *)(model->specs);

	if (sol){
		if (sol->fdpars){
			ANSOL_ESG_PARS * ansolpars = (ANSOL_ESG_PARS *)(sol->fdpars);
			if (ansolpars->srcinfo){
				RICKER_INFO * rickerinfo = (RICKER_INFO *)(ansolpars->srcinfo);
				userfree_(rickerinfo);
			}
			if (ansolpars->medinfo){
				ESG_HIMED_INFO * himedinfo = (ESG_HIMED_INFO *)(ansolpars->medinfo);
				userfree_(himedinfo);
			}
			userfree_(ansolpars);
		}
		userfree_(sol);
	}

  	/* fdm is freed in im_destroy */
  	return 0;
}

/*----------------------------------------------------------------------------*/
int ansol_HI_esg_medinfo( FILE * stream, IMODEL * model ) {
/*----------------------------------------------------------------------------*/
	int err;

	FD_MODEL * sol = (FD_MODEL *)(model->specs);
  	if (!sol) {
		fprintf(stream,"Error: in ansol_HI_esg_medinfo, IMODEL has null specs!\n");
		return 1;
	}
	if (!sol->fdpars) 
      	{		
		fprintf(stream,"Error: in ansol_HI_esg_medinfo, FD_MODEL null fdpars!\n");
		return 1;
	}
  	ANSOL_ESG_PARS * ansolpars = (ANSOL_ESG_PARS *)(sol->fdpars);

	//Extracting ESG_HIMED_INFO
	ESG_HIMED_INFO * himedinfo = (ESG_HIMED_INFO *)(ansolpars->medinfo); //allocated in ansol_HI_esg_modelinit.
	
	IPNT   n; 	//numb of points per axis
	int    ntot; 	//tot numb of points in field
	ireal  fmax;	//maximum of field
	ireal  fmin;	//minimum of field
	RDOM * dom;	//buffer domain
	int    i, field_i;
	int    ndim;


	dom = &(model->ld_c); 
	ndim = ansolpars->ndim;

	//extracting mu of fd_model, correspoding to field with index D_MS0
	field_i = D_MS0;
	rd_size(dom,field_i,n);
	ntot = 1;
	for (i=0;i<ndim;i++) ntot*=n[i];
	//Computing max and min of mu
	fmax = dom->_s[field_i]._s0[0];
	fmin = fmax;
	for (i=1;i<ntot;i++){
		fmax = iwave_max(fmax,dom->_s[field_i]._s0[i]);
		fmin = iwave_min(fmin,dom->_s[field_i]._s0[i]);
	}
	if (fabs(fmax-fmin)>TOL){
		fprintf(stream,"Error: in ansol_HI_esg_medinfo, field for mu is not constant! ABORT!\n");
		fflush(stream);
		exit(1);
	}
#ifdef IWAVE_VERBOSE
	fprintf(stream,"NOTE: for mu, max-min = %g\n",fmax-fmin);
	fprintf(stream,"mu = %g\n",fmax);
#endif
	//if mu is constant, then copying into himedinfo
	himedinfo->mu = fmax;

	//extracting lambda of fd_model, correspoding to field with index D_MP01
	field_i = D_MP01;
	rd_size(dom,field_i,n);
	ntot = 1;
	for (i=0;i<ndim;i++) ntot*=n[i];
	//Computing max and min of lambda
	fmax = dom->_s[field_i]._s0[0];
	fmin = fmax;
	for (i=1;i<ntot;i++){
		fmax = iwave_max(fmax,dom->_s[field_i]._s0[i]);
		fmin = iwave_min(fmin,dom->_s[field_i]._s0[i]);
	}
	if (fabs(fmax-fmin)>TOL){
		fprintf(stream,"Error: in ansol_HI_esg_readmedia, field for lambda is not constant! ABORT!\n");
		exit(1);
	}
#ifdef IWAVE_VERBOSE
	fprintf(stream,"NOTE: for lambda, max-min = %g\n",fmax-fmin);
	fprintf(stream,"lambda = %g\n",fmax);
#endif
	//if lambda is constant, then copying into himedinfo
	himedinfo->lambda = fmax;

	//extracting lambda of fd_model, correspoding to field with index D_MV0
	field_i = D_MV0;
	rd_size(dom,field_i,n);
	ntot = 1;
	for (i=0;i<ndim;i++) ntot*=n[i];
	//Computing max and min of 1/rho
	fmax = dom->_s[field_i]._s0[0];
	fmin = fmax;
	for (i=1;i<ntot;i++){
		fmax = iwave_max(fmax,dom->_s[field_i]._s0[i]);
		fmin = iwave_min(fmin,dom->_s[field_i]._s0[i]);
	}
	if (fabs(fmax-fmin)>TOL){
		fprintf(stream,"Error: in ansol_HI_esg_readmedia, field for rho is not constant! ABORT!\n");
		exit(1);
	}
#ifdef IWAVE_VERBOSE
	fprintf(stream,"NOTE: for rho, max-min = %g\n",fmax-fmin);
	fprintf(stream,"rho = %g\n",fmax);
#endif
	//if 1/rho is constant, then copying into himedinfo
	himedinfo->rho = 1.0/fmax;
	
	//computing alpha and beta
	himedinfo->alpha = sqrt( (himedinfo->lambda + 2.0*himedinfo->mu)/himedinfo->rho );
	himedinfo->beta = sqrt( himedinfo->mu/himedinfo->rho );

	return 0;
}

/*----------------------------------------------------------------------------*/
int ansol_HI_esg_srcinfo( FILE * stream, PARARRAY * pars, IMODEL * model, tracegeom * tg, POINTSRC * pntsrc ){
/*----------------------------------------------------------------------------*/
	if( !tg){
		fprintf(stream,"ERROR: ansol_HI_esg_srcinfo, empty tracegeom\n");
		return 1;
	}
	if( !pntsrc){
		fprintf(stream,"ERROR: ansol_HI_esg_srcinfo, empty POINTSRC\n");
		return 1;
	}

	FD_MODEL * sol = (FD_MODEL *)(model->specs);
	if (!sol){
		fprintf(stream,"ERROR: ansol_HI_esg_srcinfo, model has empty specs\n");
		return 1;
	}

	ANSOL_ESG_PARS * ansolpars = (ANSOL_ESG_PARS *)(sol->fdpars);
	if (!ansolpars){
		fprintf(stream,"ERROR: ansol_HI_esg_srcinfo, FD_MODEL has empty fdpars\n");
		return 1;
	}

	RICKER_INFO * rickerinfo = (RICKER_INFO *)(ansolpars->srcinfo);
	if (!rickerinfo){
		fprintf(stream,"ERROR: ansol_HI_esg_srcinfo, ANSOL_ESG_PARS has empty srcinfo\n");
		return 1;
	}

	ESG_HIMED_INFO * himedinfo = (ESG_HIMED_INFO *)(ansolpars->medinfo);
	if (!himedinfo){
		fprintf(stream,"ERROR: ansol_HI_esg_srcinfo, ANSOL_ESG_PARS has empty medinfo\n");
		return 1;
	}

	char * srctype;
	char * wp;
	int    iw;
	ireal  fpeak   = FPEAK_DEF;	//default frequency peak (def in defaults.h) of 10 Hz.
	int    idim, ndim;
	ireal  prod_d;

	ndim = ansolpars->ndim;

	ps_ffcstring(*pars,"srctype",&srctype);
	//Asserting source is a point source
	if ( strcmp(srctype,"point")  ){
		fprintf(stream,"ERROR: ansol_HI_esg_srcinfo, analytical solution flag is set for");
		fprintf(stream,"       homogeneous isotropic case with a non-point source!\nABORT!\n");
		abortexit(1,pars,&stream);
	}
	
	//getting peak frequency
	if (ps_ffreal(*pars,"fpeak", &fpeak)) {
    		fprintf(stream,"NOTE: ansol_HI_esg_srcinfo - using default ");
    		fprintf(stream,"peak frequency (fpeak) = %e\n",fpeak);
  	}

	//setting t_off in ANSOL_ESG_PARS from zerophase or causal option for source
	iw = (pntsrc->n - 1.0)/2.0;
	ansolpars->t_off = (iw + pntsrc->istart) * ansolpars->dt;

	//setting src_x in ANSOL_ESG_PARS
	RASN( ansolpars->src_x, tg->src[0] );

	//setting fpeak in RICKER_INFO
	rickerinfo->fpeak = fpeak;

	prod_d = 1.0;
	for (idim=0; idim<ndim; idim++)
		prod_d *= ansolpars->dx[idim];
	//setting amp in RICKER_INFO
	//assuming srcamp = bou * refdist * refamp * dt / dxdydz
	rickerinfo->amp = pntsrc->scramp * prod_d * himedinfo->rho / ansolpars->dt;

	//setting src_d in RICKER_INFO
	for (idim=0; idim<ndim; idim++)
		rickerinfo->src_d[idim] = pntsrc->src_d[idim];

/*	userfree_(wp);*/
/*	userfree_(srctype);*/
	return 0;
}

/*----------------------------------------------------------------------------*/
int ansol_link_mdl( FILE * stream, IMODEL * sol_mdl, IMODEL * fd_mdl ){
/*----------------------------------------------------------------------------*/
	FD_MODEL * sol = (FD_MODEL *)(sol_mdl->specs);
	if (!sol){
		fprintf(stream,"ERROR: ansol_HI_esg_srcinfo, sol_mdl has empty specs\n");
		return 1;
	}

	ANSOL_ESG_PARS * ansolpars = (ANSOL_ESG_PARS *)(sol->fdpars);
	if (!ansolpars){
		fprintf(stream,"ERROR: ansol_HI_esg_srcinfo, FD_MODEL has empty fdpars\n");
		return 1;
	}
	ansolpars->link_mdl = fd_mdl;
	return 0;
}

/*----------------------------------------------------------------------------*/
void ansol_esg_fprintf( FILE * stream, IMODEL * sol_mdl ){
/*----------------------------------------------------------------------------*/
	int idim;

	fprintf(stream,"\n");
	fprintf(stream,"==================================\n");	
	fprintf(stream,"PRINTING ansol model info\n");
	fprintf(stream,"==================================\n");	


	FD_MODEL * sol = (FD_MODEL *)(sol_mdl->specs);
	if (!sol){
		fprintf(stream,"NOTE: model has empty specs\n");
		fflush(stream);
		return;
	}

	ANSOL_ESG_PARS * ansolpars = (ANSOL_ESG_PARS *)(sol->fdpars);
	if (!ansolpars){
		fprintf(stream,"NOTE: FD_MODEL has empty fdpars\n");
		fflush(stream);
		return;
	}

	fprintf(stream,"Parameter info:\n");
	fprintf(stream,"---------------------\n");
	fprintf(stream,"dt          = %g\n",ansolpars->dt);
	fprintf(stream,"k           = %d\n",ansolpars->k);
	fprintf(stream,"ndim        = %d\n",ansolpars->ndim);
for( idim=0; idim<ansolpars->ndim; idim++)
	fprintf(stream,"lam[%d]     = %g\n",idim,ansolpars->lam[idim]);
if (ansolpars->link_mdl)
	fprintf(stream,"link_mdl    = NON-NULL\n");
else
	fprintf(stream,"link_mdl    = NULL\n");
for( idim=0; idim<ansolpars->ndim; idim++)
	fprintf(stream,"dx[%d]      = %g\n",idim,ansolpars->dx[idim]);
for( idim=0; idim<ansolpars->ndim; idim++)
	fprintf(stream,"o_coord[%d] = %g\n",idim,ansolpars->o_coord[idim]);
for( idim=0; idim<ansolpars->ndim; idim++)
	fprintf(stream,"o_index[%d] = %g\n",idim,ansolpars->o_index[idim]);
for( idim=0; idim<ansolpars->ndim; idim++)
	fprintf(stream,"src_x[%d]   = %g\n",idim,ansolpars->src_x[idim]);
	fprintf(stream,"t_off       = %g\n",ansolpars->t_off);
	fprintf(stream,"\n");


	ESG_HIMED_INFO * himedinfo = (ESG_HIMED_INFO *)(ansolpars->medinfo);
	if (!himedinfo){
		fprintf(stream,"NOTE: ANSOL_ESG_PARS has empty medinfo\n");
		fflush(stream);
		return;
	}

	fprintf(stream,"Medium info:\n");
	fprintf(stream,"---------------------\n");
	fprintf(stream,"mu     = %g\n",himedinfo->mu);
	fprintf(stream,"lambda = %g\n",himedinfo->lambda);
	fprintf(stream,"rho    = %g\n",himedinfo->rho);
	fprintf(stream,"alpha  = %g\n",himedinfo->alpha);
	fprintf(stream,"beta   = %g\n\n",himedinfo->beta);


	RICKER_INFO * rickerinfo = (RICKER_INFO *)(ansolpars->srcinfo);
	if (!rickerinfo){
		fprintf(stream,"NOTE: ANSOL_ESG_PARS has empty srcinfo\n");
		fflush(stream);
		return;	
	}

	fprintf(stream,"Source info:\n");
	fprintf(stream,"---------------------\n");
	fprintf(stream,"fpeak     = %g\n",rickerinfo->fpeak);
	fprintf(stream,"amp       = %g\n",rickerinfo->amp);
for( idim=0; idim<ansolpars->ndim; idim++)
	fprintf(stream,"src_d[%d] = %g\n",idim,rickerinfo->src_d[idim]);
	fprintf(stream,"==================================\n\n");	

	fflush(stream);
	return;
}