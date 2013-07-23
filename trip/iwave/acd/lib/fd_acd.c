#include "acd.h"
#include "iwave.h"

/*#define IWAVE_VERBOSE*/

/*------------------- private data declarations ----------------------------*/

static int m_ndim = 0; /* dimension - need only compute once */
static int m_size = 3; 
static const char* m_names[] = {"uc", "up", "csq"};
/* this static variable is available to acd_movie_select, so can
   be used to control the choice of current updated array in movie
   and in other apps */
static int UC = D_UC;

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

/*--- private function declarations - assigned to FD_MODEL pointers ----------*/

int acd_isarr(int i);
int acd_numsubsteps();
int acd_update();
int acd_readschemeinfo(PARARRAY *, FILE *, IMODEL *);
int acd_set_grid_type(FILE *, int, IPNT[RDOM_MAX_NARR]);
int acd_build_sten_dep(FILE *, int, int[RDOM_MAX_NARR][RDOM_MAX_NARR]);
int acd_create_sten(FD_MODEL *, FILE *, int, IPNT[RDOM_MAX_NARR], int[RDOM_MAX_NARR][RDOM_MAX_NARR], STENCIL *);
int acd_alter_dom(int, IPNT, IPNT);
const char * acd_ind2str(int);  
int acd_modeldest(IMODEL * model);
void acd_ts_parcopy(void * tgt, const void * src);
int acd_readgrid(PARARRAY *, FILE *, IMODEL *);
int acd_readtimegrid(PARARRAY *, FILE *, IMODEL *);
int acd_readmedia(PARARRAY *, FILE *, IMODEL *, int);

extern int acd_step(RDOM*, int, void *);

/*----------------------------------------------------------------------------*/
/* no-ops for this implementation                                             */
/*----------------------------------------------------------------------------*/
int acd_build_sten_dep(FILE * stream, 
		       int ndim, 
		       int stendep[RDOM_MAX_NARR][RDOM_MAX_NARR]) {
    return 0;
}

int acd_alter_dom(int ia, IPNT gs, IPNT ge) {
    return 0;
}

/*----------------------------------------------------------------------------*/
/* working functions                                                          */
/*----------------------------------------------------------------------------*/

const char* acd_ind2str(int ind) {
    if (ind>-1 && ind<m_size) return m_names[ind];
    return "NONE";
}

/*----------------------------------------------------------------------------*/

int acd_modelinit(PARARRAY *pars, 
		  FILE *stream, 
		  IMODEL * im) {

    int err=0;           /* return value */
    ACD_TS_PARS *acdpars;   /* model pars */
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
    /*  fprintf(stream,"top of acd_modelinit: im->specs=%x\n",im->specs); */
    im->specs=(void *)fd;
    /*  fprintf(stream,"top of acd_modelinit: im->specs=%x after assignment\n",im->specs); */

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
    int rk=retrieveRank();
    MPI_Comm cm=retrieveComm();

    if ( MPI_Cart_get(cm, RARR_MAX_NDIM, cdims, cpers, crank) != MPI_SUCCESS )  {
	fprintf(stream, "ERROR. Internal: cannot get Cartesian coordinates.\n");
	return E_INTERNAL;
    }

    if (rk==0) {
#endif
	/* must read grid here, even though it's read again later */
	m_ndim=0;
    
	if (!acd_readgrid(pars, stream, im)) {
	    m_ndim = (im->g).dim;
	}
	else {
	    err=E_INTERNAL;
	    fprintf(stream,"ERROR: in acd_modelinit - failed to read spatial geometry\n");
	    fflush(stream);
	    abortexit(err,pars,&stream);
	}
    
#ifdef IWAVE_USE_MPI
    }
    MPI_Bcast(&m_ndim,1,MPI_INT,0,cm);
#endif

    if (m_ndim<1 || m_ndim>RARR_MAX_NDIM) {
	err=E_INTERNAL;
	fprintf(stream,"ERROR: in acd_modelinit - failed to read dim=%d\n",m_ndim);
	return err;
    }

    /* set boundary flags */
    IASN(acdpars->lbc, IPNT_0); /* default left bc flag */ 
    IASN(acdpars->rbc, IPNT_0); /* default right bc flag */ 

    for (i=0;i<m_ndim;i++) {
	if (crank[i]==0) acdpars->lbc[i]=1;
	if (crank[i]==cdims[i]-1) acdpars->rbc[i]=1;
    }

    /* decode order - with version 2.0, deprecated syntax "scheme_phys" etc. is dropped */
    acdpars->k=1;
    ps_flint(*pars,"order",&(acdpars->k));
#ifdef IWAVE_VERBOSE
    fprintf(stream,"NOTE: initializing ACD with half-order = %d\n",acdpars->k);
#endif


    /* initialize scaled Courant arrays */
    acdpars->c0=REAL_ONE;
    RASN(acdpars->c1,RPNT_0);
    RASN(acdpars->c2,RPNT_0);
    RASN(acdpars->c3,RPNT_0);
    RASN(acdpars->c4,RPNT_0);

    /* assign param object pointer */
    fd->fdpars = (void*)acdpars;

    /* assign function pointers for ACD model */
    fd->isarr  = acd_isarr;
    fd->numsubsteps = acd_numsubsteps;
    fd->update = acd_update;
    fd->readgrid = acd_readgrid;
    fd->readtimegrid = acd_readtimegrid; 
    fd->readschemeinfo = acd_readschemeinfo;
    fd->set_grid_type = acd_set_grid_type;
    fd->build_sten_dep = acd_build_sten_dep;
    fd->create_sten = acd_create_sten;
    fd->ind2str = acd_ind2str;
    fd->alter_dom = acd_alter_dom;
    fd->readmedia = acd_readmedia;
    fd->parcopy = acd_ts_parcopy;
    fd->fd_model_init = acd_modelinit;
    fd->fd_model_dest = acd_modeldest;

    /* choose time step */
    fd->tsf=acd_step;

    return 0;
}

/*----------------------------------------------------------------------------*/
int acd_numsubsteps() { return 1; }

int acd_modeldest(IMODEL * model) {
    FD_MODEL * fdm = (FD_MODEL *)(model->specs);
    /* since the FD_MODEL is allocated here, destroy here */
    if (fdm) userfree_(fdm);
    return 0;
}

void acd_ts_parcopy(void * tgt, const void * src) {

    int i;

    ACD_TS_PARS * ptgt = (ACD_TS_PARS *)tgt;
    const ACD_TS_PARS * psrc = (const ACD_TS_PARS *)src;

    ptgt->dt=psrc->dt;
    for (i=0;i<RARR_MAX_NDIM;i++) 
	(ptgt->lam)[i]=(psrc->lam)[i];
    ptgt->ndim=psrc->ndim;
    ptgt->k = psrc->k;

    IASN(ptgt->lbc,psrc->lbc);
    IASN(ptgt->rbc,psrc->rbc);   
}

/*----------------------------------------------------------------------------*/
int acd_isarr(int i) {
    /* exactly three arrays:
    // uc = array 0
    // up = array 1
    // csq = array 2 */
    if (i>-1 && i<3) return 1;
    return 0;
}

/*----------------------------------------------------------------------------*/
/* only return true for arrays actually updated!!!!! */
int acd_update(int ia, int iv) {
    /* update all dynamical arrays on iv = 0 - this would be correct
    // for any one-step update */
    if (((ia==0) || (ia==1)) && iv==0) return 1;
    return 0;
}

/*----------------------------------------------------------------------------*/
int acd_set_grid_type(FILE *stream, int ndim, IPNT gtype[RDOM_MAX_NARR] ) {
    /* all arrays defined on primal grid for ACD modeling */
    int iv; 
    for (iv = 0;iv < RDOM_MAX_NARR;iv ++)  IASN(gtype[iv], IPNT_0);
    /* just an additional sanity check */
    if ( ndim < 1 || ndim > RARR_MAX_NDIM ) return E_BADINPUT;
    return 0;
}


/*----------------------------------------------------------------------------*/
int acd_modelinfo(FILE *stream, IMODEL *model)
{
    int ndim, i;
    IPNT n;
    RDOM *dc; /* computational */

    ndim = model->g.dim;

    fprintf(stream, "ACD model %dD.\n", ndim);

    dc = &(model->ld_c); /* get computational domain */
    rd_size(dc, D_UC, n);

    fprintf(stream, "  UC size = [%d", n[0]);
    for ( i = 1; i < ndim; ++i ) fprintf(stream, " %d", n[i]);
    fprintf(stream, "]\n");

    return 0;
}

/*----------------------------------------------------------------------------*/
/* this variant so simple that might as well write it from
// scratch, so no need for sten_dep_mat - note that all arrays in this
// app are primal, so no
// need for gtype. Clearly this interface should be refactored so that
// these things can be hidden.
// ndim (from fdpars) gtype and sten_dep_mat should be internal details */

int acd_create_sten(FD_MODEL * fdm,
		    FILE * stream, 
		    int ndim,
		    IPNT gtype[RDOM_MAX_NARR], 
		    int sten_dep_mat[RDOM_MAX_NARR][RDOM_MAX_NARR], 
		    STENCIL * sten) {
    ACD_TS_PARS * acdpars = (ACD_TS_PARS *)(fdm->fdpars);
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
int acd_readschemeinfo(PARARRAY * par, 
		       FILE * stream, 
		       IMODEL * model) {
  
    RPNT dxs;    /* grid steps */
    ireal lam;   /* slownesses dt/dx */
    int idim;    /* counter */

    FD_MODEL * fdm = (FD_MODEL *)(model->specs);
    ACD_TS_PARS * acdpars = (ACD_TS_PARS *)(fdm->fdpars);

    /* extract grid steps from grid */
    get_d(dxs, model->g);

    /* set model dimn par */
    acdpars->ndim = (model->g).dim;

    acdpars->c0 = 0.0;
    for (idim = 0;idim < acdpars->ndim;idim ++) {

	if (dxs[idim] <= 0.0) {
	    fprintf(stream, "Error: bad input: wrong grid space step, dim=%d, step=%g\n",
		    idim, dxs[idim]);
	    return E_BADINPUT;
	}
	lam = (model->tsind).dt / dxs[idim];
#ifdef IWAVE_VERBOSE
	fprintf(stderr, "lam[%d] = %g\n", idim, lam);
#endif
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
    }

    /* reserve a copy of dt for use in source scaling */
    acdpars->dt = (model->tsind).dt;

    return 0;

}

/*----------------------------------------------------------------------------*/
/* implements new time grid logic:
   - by default choose stable time step (max_step set), optionally with 
   reduced cfl
   - if max_step unset, use CFL timestep but test against max step
   WWS 02.10
   NEW INTERFACE 22.11.10 WWS
   old interface
   int sg_readtimegrid(RDOM dom, FILE *stream, ireal *dt, grid g, PARARRAY par) {
   where
   dom = model->ld_a;
   dt  = &((model->tsind).dt);
   g=model->g;
   par=*pars;

*/

int acd_readtimegrid(PARARRAY *pars, FILE * stream, IMODEL * model) {

    ireal * dt;                   /* time step, to be computed */

    ireal cmax;                   /* max velo, computed or read from params */
    ireal csq;                    /* max c-squared, from data */
    ireal cfl = CFL_DEF;          /* default cfl fraction */
    int max_step = 1;             /* flag to take max permitted step */
    RDOM dom;                     /* allocated domain */
    IPNT n;                       /* workspace for size */
    int nsz;                      /* size */
    ireal a;                      /* accumulator for max space step */
    int i;                        /* counter */
    int ndim;                     /* model dimension */

    /* pointer into model->tsind struct */
    dt  = &((model->tsind).dt);

    /* branch on presence of parameter dt - if dt set in param table,
       use it */
    if ( !ps_flreal(*pars,"dt", dt ) ){
#ifdef IWAVE_VERBOSE
	fprintf(stream,"NOTE: sg_readtimegrid - dt=%12.4e read from param table\n", *dt);	
	fprintf(stream,"NOTE: NOT CHECKED FOR STABILITY!\n");
#endif
	return 0;
    }
    	
    /* for either method of computing dt, need cfl = proportion
       of max time step to use */
    cfl=CFL_DEF;
    if (ps_flreal(*pars,"cfl",&cfl)) {
#ifdef IWAVE_VERBOSE
	fprintf(stream,"NOTE: acd_readtimegrid\n");
	fprintf(stream,"  using default cfl fraction %g\n",cfl);;
#endif
    }
    /* branch on max_step */
    ps_flint(*pars,"max_step",&max_step);		

    ndim=(model->g).dim;
  
    if (max_step) {
#ifdef IWAVE_VERBOSE
	fprintf(stream,"NOTE: acd_readtimegrid\n");
	fprintf(stream,"  compute max permissible step\n");
	fflush(stream);
#endif
	/* find max velocity */
	dom = model->ld_a;
	csq=0.0;
	ra_a_size(&(dom._s[D_CSQ]),n);
	nsz=1;
	for (i=0;i<ndim;i++) nsz *= n[i];
	for (i=0;i<nsz;i++) {
	    csq = iwave_max(csq,dom._s[D_CSQ]._s0[i]);
	    if (csq<REAL_EPS) {
		fprintf(stream,"Error: acd_readtimegrid\n");
		fprintf(stream,"  at global index = %d\n",i);
		fprintf(stream,"  max of square velocity < %g\n",REAL_EPS);
		return E_BADINPUT;
	    }
	}
    
	cmax = sqrt(csq);
#ifdef IWAVE_VERBOSE
	fprintf(stream,"NOTE: acd_readtimegrid, computed cmax=%e\n",cmax);
	fflush(stream);
#endif
    }

    else {
  
	if (ps_flreal(*pars,"cmax",&cmax)) { 
	    fprintf(stream,"ERROR: acd_readtimegrid - failed to read cmax from parameters \n");
	    return E_BADINPUT;
	}

    }

    a = (model->g).axes[0].d;
  
    for ( i = 1; i < ndim; ++i ) a = iwave_min(a, (model->g).axes[i].d);
#ifdef IWAVE_VERBOSE
    fprintf(stream,"NOTE: acd_readtimegrid: ndim=%d min dx=%e cfl=%e\n",ndim,a,cfl);
    fflush(stream);
#endif

    if ( (a < REAL_EPS) || (cfl < REAL_EPS) ) {
	fprintf(stream,"ERROR: sg_readcfltime - either min dx=%e or cfl=%e "
		" too small\n", a, cfl);
	return E_BADINPUT;
    }

    *dt = a*cfl/(cmax*sqrt((float)(ndim)));
#ifdef IWAVE_VERBOSE
    fprintf(stream,"NOTE: acd_readtimegrid: on return, dt=%e\n",*dt);
    fflush(stream);
#endif

    return 0;
}

/*----------------------------------------------------------------------------*/
/** Input mechanical parameters. Target: velocity squared on
    computational grid. 

    Includes velocity check.
*/

int acd_readmedia(PARARRAY * pars, 
		  FILE * stream,
		  IMODEL * model,
		  int panelindex) {

    /*==================================================
      ================= DECLARATIONS ===================
      ==================================================*/

    /* workspace for filenames */  
    char * velokey=NULL;

    /* domain */
    RDOM dom;

    /* other workspace */
    int i,j;         /* counters */
    int dim;         /* problem dimension */
    IPNT ran;        /* axis lens */
    IPNT rags;       /* axis starts */
    IPNT rage;       /* axis ends */
    unsigned long ntot;        /* total num of words in phys dom */
    ireal cmax = REAL_MAX;  /* max velo */
    ireal cmin = REAL_EPS;  /* min velo */

    ireal vmax=0.0;  /* working vmax */
    ireal vmin=0.0;  /* working vmin */

    /* error flag */
    int err=0;

    /*==================================================
      =============== END DECLARATIONS =================
      ==================================================*/  

    /* initializations */
    dom = model->ld_a;

    /* axis lens, starts for bulk mod - allocated, not virtual */
    rd_ndim(&dom,D_CSQ,&dim);
    rd_size(&dom,D_CSQ,ran);
    rd_gse(&dom, D_CSQ,rags,rage);
    ntot=1;
    for (i=0;i<dim;i++) ntot*=ran[i];

    /* set velocity limits
    // default values - zero and infinity, not much use! */
    if (ps_flreal(*pars,"cmax",&cmax))  {
#ifdef IWAVE_VERBOSE
	fprintf(stream,"NOTE: asg_readmedia - using default max velocity = %e\n",cmax);
#endif
    }
    if (ps_flreal(*pars,"cmin",&cmin)) {
#ifdef IWAVE_VERBOSE
	fprintf(stream,"NOTE: asg_readmedia - using default min velocity = %e\n",cmin);
#endif
    }

    /* FILE I/O SECTION */
    /* read velocity into D_CSQ - then square, sanity check */
    if (!ps_flcstring(*pars,"velocity",&velokey)) {
#ifdef IWAVE_VERBOSE
	fprintf(stream,"velocity option val=%s\n",velokey);
#endif
	err=rsfread(dom._s[D_CSQ]._s0,rags,ran,velokey,1,stream,panelindex);
	userfree_(velokey);
	if (err) {
	    fprintf(stream,"ERROR: acd_readmedia from ra_rsfread\n");
	    return err;
	}
	vmax=dom._s[D_CSQ]._s0[0];
	vmin=dom._s[D_CSQ]._s0[0];
	dom._s[D_CSQ]._s0[0] = dom._s[D_CSQ]._s0[0] * dom._s[D_CSQ]._s0[0];
	for (j=1;j<ntot;j++) {
	    vmax=iwave_max(vmax,dom._s[D_CSQ]._s0[j]);
	    vmin=iwave_min(vmin,dom._s[D_CSQ]._s0[j]);
	    dom._s[D_CSQ]._s0[j] = dom._s[D_CSQ]._s0[j] * dom._s[D_CSQ]._s0[j];
	}
    }
    /* 09.02.13 - add option to read velocity squared directly */
    else if (!ps_flcstring(*pars,"csq",&velokey)) {
#ifdef IWAVE_VERBOSE
	fprintf(stream,"velocity squared option val=%s\n",velokey);
#endif
	err=rsfread(dom._s[D_CSQ]._s0,rags,ran,velokey,1,stream,panelindex);
	userfree_(velokey);
	if (err) {
	    fprintf(stream,"ERROR: acd_readmedia from ra_rsfread\n");
	    return err;
	}
	vmax=sqrt(iwave_max(REAL_ZERO, dom._s[D_CSQ]._s0[0]));
	vmin=vmax;
	for (j=1;j<ntot;j++) {
	    vmax=iwave_max(vmax,sqrt(iwave_max(REAL_ZERO, dom._s[D_CSQ]._s0[j])));
	    vmin=iwave_min(vmin,sqrt(iwave_max(REAL_ZERO, dom._s[D_CSQ]._s0[j])));
	}
    }
    else {
	fprintf(stream,
		"ERROR: acd_readmedia - filename for velocity field\n");
	fprintf(stream,"  not specified in param table\n");
	err=E_BADINPUT;
    }
  
    err = (vmin < cmin) || (vmax > cmax);
    if (err) {
	fprintf(stream,
		"ERROR: acd_readmedia - velocity field read from file violates bounds [%e, %e]\n",
		cmin,cmax);
	fprintf(stream,"min v = %e max v = %e\n",vmin,vmax);
	err = E_BADINPUT;
    }
#ifdef IWAVE_VERBOSE
    else {
	fprintf(stream,
		"NOTE: acd_readmedia - velocity field read from file within bounds [%e, %e]\n",
		cmin,cmax);
	fprintf(stream,"min v = %e max v = %e\n",vmin,vmax);
    }
#endif

    return err;
}

/*----------------------------------------------------------------------------*/
/* read simulation grid from velocity file */
int acd_readgrid(PARARRAY * pars, 
		 FILE *stream, 
		 IMODEL * model) {

    char * _velokey =NULL;    /* workspace for filename */
    int err=0;

    /* check for permissible velo/bulkmod keywords */  
    if (ps_flcstring(*pars,"velocity",&_velokey)) 
	ps_flcstring(*pars,"csq",&_velokey); 

    if (_velokey) {
	err=read_grid(&(model->g),_velokey,stream);
	if (_velokey) userfree_(_velokey);
	if (err) 
	    fprintf(stream,"Error: acd::read_grid from read_grid\n");
    }
    else {
	fprintf(stream,"Error: acd::read_grid\n");
	fprintf(stream,"  failed to extract velocity filename from par file\n");
	ps_printall(*pars,stderr);
	fflush(stream);
	err = E_BADINPUT;
    }
  
    /* until such time as someone implements PML...*/
    IASN(model->nls,IPNT_0);
    IASN(model->nrs,IPNT_0);
  
    return err;
}

/* returns static UC which at after call to acd_step
   is index of updated current field
*/
int acd_movie_select(const char * key) {
    if (!strcmp(key,"p")) return UC;
    return -1;
}

int acd_movie_construct(MOVIE * mt, FILE * stream) {
    movie_setnull(mt);
    mt->iselect=acd_movie_select;
    return 0;
}

int acd_step(RDOM* dom, int iv, void * tspars) {

    /* pointers for 2D case */
    register ireal ** restrict uc2;
    register ireal ** restrict up2;
    register ireal ** restrict csq2;
    /* pointers for 3D case */
    register ireal *** restrict uc3;
    register ireal *** restrict up3;
    register ireal *** restrict csq3;
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
	    acd_2d_2(uc2, up2, csq2, 
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

	for (i[1]=s0[1];i[1]<=e0[1];i[1]++) {
	    for (i[0]=s0[0];i[0]<=e0[0];i[0]++) {
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
