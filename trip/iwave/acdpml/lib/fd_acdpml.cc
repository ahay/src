#include "acdpml.hh"

using RVL::parse;

/*#define IWAVE_VERBOSE*/

/*--- time step functions ---------------------------------------------------*/
extern void acdpml_2d_2(float **, float **, float **,
                        float **, float **,
                        float *,  float *,
                        float *,  float,
                        int *,    int *,
                        float,    float *,
                        int *,    int *);
extern void acdpml_2d_4(float **, float **, float **,
                        float **, float **,
                        float *,  float *,
                        float *,  float,
                        int *,    int *,
                        float,    float *,  float *,
                        int *,    int *);
extern void acdpml_2d_8(float ** uc, float ** up, float ** csq,
                        float ** phi1, float ** phi0,
                        float * dp1,   float * dp0,
                        float *di, float dt,
                        int * s,   int * e,
                        float c0,  float * c1,  float * c2,
                        float *c3, float * c4,
                        int * lbc, int * rbc);

extern int acdpml_step(RDOM*, int, void *);

/*----------------------------------------------------------------------------*/
/* no-ops for this implementation                                             */
/*----------------------------------------------------------------------------*/
int acdpml_build_sten_dep(FILE * stream,
		       int ndim, 
		       int stendep[RDOM_MAX_NARR][RDOM_MAX_NARR]) {
    return 0;
}

/*----------------------------------------------------------------------------*/
/* working functions                                                          */
/*----------------------------------------------------------------------------*/

int acdpml_modelinit(PARARRAY *pars,
		  FILE *stream, 
		  grid const & g,
		  ireal dt,
		  void ** fdpars) {

    int err=0;           /* return value */
    ACDPML_TS_PARS *acdpmlpars;   /* model pars */

    IPNT cdims;          /* workspace for cartesian grid dim info */
    IPNT crank;          /* workspace for cartesian grid rank */
#ifdef IWAVE_USE_MPI
    IPNT cpers;          /* workspace for periodic wrap info  - currently not used */
#endif

    int i;               /* counter */
    RPNT dxs;    /* grid steps */
    ireal lam;   /* slownesses dt/dx */
    int idim, L0, L1;    /* counter */
    int gd1,gd0,i1,i0;
    float pmlampl;
    float freesurface;
    float pi = 4.0*atan(1.0);

    /* allocate sgn model ----------------------------------------------------*/   
    acdpmlpars = (ACDPML_TS_PARS*)usermalloc_(sizeof(ACDPML_TS_PARS));
    if ( acdpmlpars == NULL ) {
	err=E_ALLOC;
	fprintf(stream,"ERROR: acd_modelinit\n");
	fprintf(stream,"failed to allocate SGN_TS_PARS object\n");
	return err;
    }

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

    MPI_Bcast((void*)(&(g.dim)),1,MPI_INT,0,cm);
#endif

    /* set boundary flags */
    IASN(acdpmlpars->lbc, IPNT_0); /* default left bc flag */
    IASN(acdpmlpars->rbc, IPNT_0); /* default right bc flag */

    for (i=0;i<g.dim;i++) {
	if (crank[i]==0) acdpmlpars->lbc[i]=1;
	if (crank[i]==cdims[i]-1) acdpmlpars->rbc[i]=1;
    }

    /* decode order - with version 2.0, deprecated syntax "scheme_phys" etc. is dropped */
    acdpmlpars->k=1;
    parse(*pars,"order",acdpmlpars->k);
#ifdef IWAVE_VERBOSE
    fprintf(stream,"NOTE: initializing ACDPML with half-order = %d\n",acdpmlpars->k);
#endif

    /* initialize scaled Courant arrays */
    acdpmlpars->c0=REAL_ONE;
    RASN(acdpmlpars->c1,RPNT_0);
    RASN(acdpmlpars->c2,RPNT_0);
    RASN(acdpmlpars->c3,RPNT_0);
    RASN(acdpmlpars->c4,RPNT_0);
    /* initialize dz, dx */
    RASN(acdpmlpars->di,RPNT_1);
    // initialize damping parameters
    acdpmlpars->dp0 = NULL;
    acdpmlpars->dp1 = NULL;
    /* initialize bound check params */
    if (!parse(*pars,"cmax",acdpmlpars->cmax)) {
      RVLException e;
      e<<"Error: acd_modelinit\n";
      e<<"  failed to find cmax in param table\n";
      throw e;
    }
    if (!parse(*pars,"cmin",acdpmlpars->cmin)) {
      RVLException e;
      e<<"Error: acd_modelinit\n";
      e<<"  failed to find cmin in param table\n";
      throw e;
    }

    /* extract grid steps from grid */
    get_d(dxs, g);
    gd0 = g.axes[0].n+1;
    gd1 = g.axes[1].n+1;
    if (!(acdpmlpars->dp0)) {
        //acdpmlpars->dp0 = (ireal *)usermalloc_(g.axes[0].n * sizeof(ireal));
        acdpmlpars->dp0 = (ireal *)malloc((gd0)*sizeof(ireal));
        //fprintf(stderr, "dp0 address = %d\n", acdpmlpars->dp0);
    }
    if (!(acdpmlpars->dp1)) {
        //fprintf(stderr, "allocate dp1 field\n");
        //acdpmlpars->dp1 = (ireal *)usermalloc_(g.axes[1].n * sizeof(ireal));
        acdpmlpars->dp1 = (ireal *)malloc((gd1)*sizeof(ireal));

    }
    // initialize dp0 and dp1 to zero
    for (int i=0; i<gd0; i++)
        acdpmlpars->dp0[i] = 0.0f;
    for (int i=0; i<gd1; i++)
        acdpmlpars->dp1[i] = 0.0f;
    // read in pml info
    if(!parse(*pars,"freesurface", freesurface)){// flag of freesurface
        freesurface = 1;
        fprintf(stream, "NOTE: DEFAULT using free surface\n");
        fprintf(stream, "  set freesurface=0 in par file to use pml boundary on the surface\n");
    }
    
    L0=0;
    L1=0;
    pmlampl=0.0;
    parse(*pars,"LZ",L0); // thickness of PML in z-axis
    parse(*pars,"LX",L1); // thickness of PML in x-axis
    parse(*pars,"pmlampl",pmlampl);
#ifdef IWAVE_VERBOSE
    fprintf(stream,"NOTE: initializing ACDPML with thickness in x-axis = %d\n",L1);
    fprintf(stream,"NOTE: initializing ACDPML with thickness in z-axis = %d\n",L0);
    fprintf(stream,"NOTE: initializing ACDPML with damping constant = %f\n",pmlampl);
#endif
    //L0 = (g.axes[0].n-2)/10;
    //L1 = (g.axes[1].n-2)/10;
    // calculate damping profile
    int PMLtype=1;
    parse(*pars,"PMLtype", PMLtype);
    if (PMLtype==0){
    for (i1=0; i1<=L1; i1++) {
        acdpmlpars->dp1[i1]=pmlampl*( fabs(i1-L1-1.0)/L1-sin(2*pi*fabs(i1-L1-1.0)/L1)/2.0/pi);
    }
    for (i1=gd1-L1-1; i1<gd1; i1++) {
        acdpmlpars->dp1[i1]=pmlampl*(fabs(i1-(gd1-1.0-L1))/L1-sin(2*pi*fabs(i1-(gd1-1.0-L1))/L1)/2.0/pi);
    }
    if (freesurface == 0) {
        for (i0=0; i0<=L0; i0++) {
            acdpmlpars->dp0[i0]=pmlampl*(fabs(i0-L0-1.0)/L0-sin(2*pi*fabs(i0-L0-1.0)/L0)/2.0/pi);
        }
    }
    for (i0=gd0-L0-1; i0<gd0; i0++) {
        acdpmlpars->dp0[i0]=pmlampl*(fabs(i0-(gd0-1.0-L0))/L0-sin(2*pi*fabs(i0-(gd0-1.0-L0))/L0)/2.0/pi);
    }
    }
    else{
        for (i1=0; i1<=L1; i1++) {
            acdpmlpars->dp1[i1]=pmlampl*(fabs((i1-L1-1.0)*(i1-L1-1)*(i1-L1-1))/L1/L1/L1);
        }
        for (i1=gd1-L1; i1<gd1; i1++) {
            acdpmlpars->dp1[i1]=pmlampl*(fabs((i1-(gd1-1.0-L1))*(i1-(gd1-1-L1))*(i1-(gd1-1-L1)))/L1/L1/L1);
        }
        if (freesurface == 0) {
            for (i0=0; i0<=L0; i0++) {
                acdpmlpars->dp0[i0]=pmlampl*(fabs((i0-L0-1.0)*(i0-L0-1)*(i0-L0-1))/L0/L0/L0);
            }
        }
        for (i0=gd0-L0; i0<gd0; i0++) {
            acdpmlpars->dp0[i0]=pmlampl*(fabs((i0-(gd0-1.0-L0))*(i0-(gd0-1-L0))*(i0-(gd0-1-L0)))/L0/L0/L0);
        }
    }
/*    fprintf(stderr, "gd0 = %d\n", gd0);
    fprintf(stderr, "gd1 = %d\n", gd1);
    fprintf(stderr, "L0 = %f\n", L0);
    fprintf(stderr, "L1 = %f\n", L1);
    fprintf(stderr, "pmlampl = %f\n", pmlampl);
    FILE *f0 = fopen("dp0.txt", "w");
    FILE *f1 = fopen("dp1.txt", "w");
    for (i0=0; i0<gd0; i0++)
        fprintf(f0, "%f\n", acdpmlpars->dp0[i0]);
    for (i1=0; i1<gd1; i1++)
        fprintf(f1, "%f\n", acdpmlpars->dp1[i1]);
    fclose(f1);
    fclose(f0);
*/
    acdpmlpars->dp0 = acdpmlpars->dp0+1;
    acdpmlpars->dp1 = acdpmlpars->dp1+1;
   
    
    /* set model dimn par */
    acdpmlpars->ndim = g.dim;

    acdpmlpars->c0 = 0.0;
    for (idim = 0;idim < acdpmlpars->ndim;idim ++) {
        acdpmlpars->di[idim] = dxs[idim];
        if (dxs[idim] <= 0.0) {
            fprintf(stream, "Error: bad input: wrong grid space step, dim=%d, step=%g\n",
		    idim, dxs[idim]);
            return E_BADINPUT;
        }
	lam = dt / dxs[idim];
#ifdef IWAVE_VERBOSE
	fprintf(stderr, "lam[%d] = %g\n", idim, lam);
#endif
	/* assign scaled Courant numbers for orders 2, 4, and 8 - these are the only */
	/* choices implemented */
	if (acdpmlpars->k==1) {
        acdpmlpars->c1[idim]   = lam*lam;
	    acdpmlpars->c0        += lam*lam*(-2.0);
    }
    else if (acdpmlpars->k==2){
	    acdpmlpars->c1[idim]   = lam*lam*(4.0/3.0);
        acdpmlpars->c2[idim]   = lam*lam*(-1.0/12.0);
	    acdpmlpars->c0        += lam*lam*(-5.0/2.0);
	}
    else if (acdpmlpars->k==4) {
	    acdpmlpars->c1[idim]   = lam*lam*(8.0/5.0);
	    acdpmlpars->c2[idim]   = lam*lam*(-1.0/5.0);
	    acdpmlpars->c3[idim]   = lam*lam*(8.0/315.0);
	    acdpmlpars->c4[idim]   = lam*lam*(-1.0/560.0);
	    acdpmlpars->c0        += lam*lam*(-205.0/72.0);
	}
	else {
	    fprintf(stream,"ERROR: acd_readschemeinfo\n");
	    fprintf(stream,"assigned scheme half-order = %d not defined\n",acdpmlpars->k);
	    fprintf(stream,"currently defined schemes: half-orders 1\n");
	    fflush(stream);
	    return E_BADINPUT;
	}
    }

    /* reserve a copy of dt for use in source scaling */
    acdpmlpars->dt = dt;

    /* assign param object pointer */
    *fdpars = (void*)acdpmlpars;
    return 0;
}

/*----------------------------------------------------------------------------*/
void acdpml_modeldest(void ** fdpars) {

    ACDPML_TS_PARS * acdpmlpars = (ACDPML_TS_PARS *)(*fdpars);
    //fprintf(stderr, "acdpmlpars->dt = %f\n", acdpmlpars->dt);
    //fprintf(stderr, "acdpmlpars->dz = %f\n", acdpmlpars->di[0]);
    //fprintf(stderr, "dp0 address = %d\n", acdpmlpars->dp0);

    acdpmlpars->dp1 = acdpmlpars->dp1-1;
    acdpmlpars->dp0 = acdpmlpars->dp0-1;
    //fprintf(stderr, "dp0 address = %d\n", acdpmlpars->dp0);

    if(acdpmlpars){
        if ((acdpmlpars->dp0)!=NULL) {
            free(acdpmlpars->dp0);
        }
        if (acdpmlpars->dp1) {
            free(acdpmlpars->dp1);
        }
        
    }
    /* destroy acdpars - all data members allocated on stack */

    userfree_(*fdpars);
    //fprintf(stderr, " after acdpml_modeldest!\n");

}

/*----------------------------------------------------------------------------*/
/* this variant so simple that might as well write it from
// scratch, so no need for sten_dep_mat - note that all arrays in this
// app are primal, so no
// need for gtype. Clearly this interface should be refactored so that
// these things can be hidden.
// ndim (from fdpars) gtype and sten_dep_mat should be internal details */

int acdpml_create_sten(void * fdm,
		    FILE * stream,
            IWaveInfo const & ic,
		    int ndim,
		    IPNT gtype[RDOM_MAX_NARR], 
		    //		    int sten_dep_mat[RDOM_MAX_NARR][RDOM_MAX_NARR], 
		    STENCIL * sten) {
    ACDPML_TS_PARS * acdpmlpars = (ACDPML_TS_PARS *)(fdm);
    STENCIL_MASK mask;/* workspace */
    int nmask;        /* number of masks - dependent pairs of dynamic arrays */
    int ipair[10][2];  /* index pairs for use in defining mask */
    int imask;        /* mask counter */
    int idim;         /* dim counter */
    int iv;           /* mask entry counter */
    int len;          /* length of mask - number of entries */
    int j;            /* counter */
    int k;            /* scheme order */
    IPNT ind;         /* workspace for mask entry */
    int err = 0;
    
    /* set order variable */
    k = acdpmlpars->k;
    
    /* initialize index pairs */
    /* first - uc->up */
    ipair[0][0] = D_UC;
    ipair[0][1] = D_UP;
    /* first - up->uc */
    ipair[1][0] = D_UP;
    ipair[1][1] = D_UC;
    /* second - uc->phi1 */
    ipair[2][0] = D_UC;
    ipair[2][1] = D_PHI1;
    /* third - up->phi1 */
    ipair[3][0] = D_UP;
    ipair[3][1] = D_PHI1;
    /* forth - uc->phi0 */
    ipair[4][0] = D_UC;
    ipair[4][1] = D_PHI0;
    /* fifth - up->phi0 */
    ipair[5][0] = D_UP;
    ipair[5][1] = D_PHI0;
    /* sixth - phi1->up */
    ipair[6][0] = D_PHI1;
    ipair[6][1] = D_UP;
    /* seventh - phi0->up */
    ipair[7][0] = D_PHI0;
    ipair[7][1] = D_UP;
//    /* eighth - du1->phi1 */
//    ipair[8][0] = D_DU1;
//    ipair[8][1] = D_PHI1;
//    /* nineth - du0->phi0 */
//    ipair[9][0] = D_DU0;
//    ipair[9][1] = D_PHI0;
    
    /* initialize STENCIL to null stencil */
    sten_setnull(sten);
    
    /* sanity check */
    if (k < 1) {
        fprintf(stream,"ERROR: acd_create_sten - illegal value of k = %d\n",k);
        return E_BADINPUT;
    }
    
    /* declare number of masks: one for shape of uc->up stencil, 
     one for shape of up->uc stencil */
    nmask=8;
 
    /* nontrivial STENCIL initialization */
    if ((err = sten_create(sten,nmask))) {
        fprintf(stream,"ERROR: acd_create_sten - failed to create stencil\n");
        return err;
    }
    
    /* length of stencil is 2k+1 in each direction, but origin is common to all
     // directions, so */
    imask=0;
    //k = 2;
    len = 2*k*ndim+1;
    for (imask=0; imask<2; imask++) {
    mask_create(&mask, ipair[imask][0], ipair[imask][1], len);
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
	}
    IASN(ind,IPNT_0);
    iv=2*k*ndim;
    if ((err = mask_set(&mask,iv,ind))) {
        fprintf(stream,"ERROR: acd_create_sten from mask_set\n");
        sten_destroy(sten);
        return err;
    }

    if ((err=sten_set(sten,imask, &mask))) {
        fprintf(stream,"ERROR: acd_create_sten from sten_set\n");
	    sten_destroy(sten);
	    return err;
	}
    }
    k = acdpmlpars->k;
    len = 4;
    for (imask=2; imask<nmask; imask++) {
        IASN(ind, IPNT_0);
        if ((err = mask_create(&mask, ipair[imask][0], ipair[imask][1], len))) {
            fprintf(stream,"ERROR: acd_create_sten from mask_create\n");
            sten_destroy(sten);
            return err;
        }
        if(imask<6){
        for (int i=0; i<len; i++) {
            ind[0] = i % 2;
            ind[1] = i >> 1;
            if ((err = mask_set(&mask,i,ind))) {
                fprintf(stream,"ERROR: acd_create_sten from mask_set\n");
                sten_destroy(sten);
                return err;
            }
        }
        }
        else{
            for (int i=0; i<len; i++) {
                ind[0] = -1*(i % 2);
                ind[1] = -1*(i >> 1);
                if ((err = mask_set(&mask,i,ind))) {
                    fprintf(stream,"ERROR: acd_create_sten from mask_set\n");
                    sten_destroy(sten);
                    return err;
                }
            }
        }
        if ((err=sten_set(sten,imask, &mask))) {
            fprintf(stream,"ERROR: acd_create_sten from sten_set\n");
            sten_destroy(sten);
            return err;
        }
    }
    return 0;
    
}

/*----------------------------------------------------------------------------*/
/* implements new time grid logic: choose stable time step (max_step set), 
   optionally with reduced cfl - to be called in iwave_construct BEFORE
   any data are read, so must depend on max velo, cannot check here.
*/

//int acd_readtimegrid(PARARRAY *pars, FILE * stream, IMODEL * model) {
//int acd_readtimegrid(PARARRAY *pars, FILE * stream, grid const & g, ireal & dt) {
int acdpml_timegrid(PARARRAY *pars, FILE * stream, grid const & g, ireal & dt) {

    ireal cmax;                   /* max velo, computed or read from params */
    ireal cfl = CFL_DEF;          /* default cfl fraction */
    ireal a;                      /* accumulator for max space step */
    int i;                        /* counter */
    
    /* branch on presence of parameter dt - if dt set in param table,
     use it */
    if (parse(*pars,"dt", dt) ){
#ifdef IWAVE_VERBOSE
        fprintf(stream,"NOTE: sg_readtimegrid - dt=%12.4e read from param table\n", dt);
        fprintf(stream,"NOTE: NOT CHECKED FOR STABILITY!\n");
#endif
        return 0;
    }
    
    /* for either method of computing dt, need cfl = proportion
     of max time step to use */
    cfl=CFL_DEF;
    if (!parse(*pars,"cfl",cfl)) {
#ifdef IWAVE_VERBOSE
        fprintf(stream,"NOTE: acd_readtimegrid\n");
        fprintf(stream,"  using default cfl fraction %g\n",cfl);;
#endif
    }
    
    if (!parse(*pars,"cmax",cmax)) {
        fprintf(stream,"ERROR: acd_readtimegrid - failed to read cmax from parameters \n");
        return E_BADINPUT;
    }
    
    a = g.axes[0].d;
    
    for ( i = 1; i < g.dim; ++i ) a = iwave_min(a, g.axes[i].d);
#ifdef IWAVE_VERBOSE
    fprintf(stream,"NOTE: acd_readtimegrid: ndim=%d min dx=%e cfl=%e\n",g.dim,a,cfl);
    fflush(stream);
#endif
    
    if ( (a < REAL_EPS) || (cfl < REAL_EPS) ) {
        fprintf(stream,"ERROR: acd_readtimegrid - either min dx=%e or cfl=%e "
                " too small\n", a, cfl);
        return E_BADINPUT;
    }
    
    dt = a*cfl/(cmax*sqrt((float)(g.dim)));
#ifdef IWAVE_VERBOSE
    fprintf(stream,"NOTE: acd_readtimegrid: on return, dt=%e\n",dt);
    fflush(stream);
#endif
    fprintf(stderr, "dt = %f\n", dt);
    return 0;
}

int acdpml_step(RDOM* dom, int iv, void * tspars) {

    /* pointers for 2D case */
    register ireal ** restrict uc2;
    register ireal ** restrict up2;
    register ireal ** restrict csq2;
    register ireal ** restrict phi1;
    register ireal ** restrict phi0;
    /* pointers for 3D case */
    //register ireal *** restrict uc3;
    //register ireal *** restrict up3;
    //register ireal *** restrict csq3;
    int ndim;                       /* problem dmn */
    IPNT s, s0;                     /* loop starts  */
    IPNT e, e0;                     /* loop ends */

    ireal tmp;
    IPNT i;
    /* acdpml struct */
    ACDPML_TS_PARS * acdpmlpars = (ACDPML_TS_PARS *)tspars;
    
    /* extract dimn info */
    ra_ndim(&(dom->_s[D_UC]),&ndim);
    ra_gse(&(dom->_s[D_UC]),s,e);
    ra_a_gse(&(dom->_s[D_UC]),s0,e0);

    if (ndim == 2) {

	/* 2D computational arrays */
	uc2   = (dom->_s)[D_UC ]._s2;
	up2   = (dom->_s)[D_UP ]._s2;
	csq2  = (dom->_s)[D_CSQ]._s2;
    phi1  = (dom->_s)[D_PHI1]._s2;
    phi0  = (dom->_s)[D_PHI0]._s2;

	/* 2nd order case */
	if (acdpmlpars->k == 1) {
        acdpml_2d_2(uc2,    // current field
                    up2,        // previous field
                    csq2,      // csq
                    phi1,      // phi1
                    phi0,      // phi0
                    acdpmlpars->dp1,        // damping profile zeta_x
                    acdpmlpars->dp0,        // damping profile zeta_x
                    acdpmlpars->di,
                    acdpmlpars->dt,
                    s,            // start index
                    e,            // end index
                    acdpmlpars->c0,
                    acdpmlpars->c1,
                    acdpmlpars->lbc,
                    acdpmlpars->rbc);
	}
    else if(acdpmlpars->k == 2){
        acdpml_2d_4(uc2, up2, csq2,
                    phi1, phi0,
                    acdpmlpars->dp1, acdpmlpars->dp0,
                    acdpmlpars->di,  acdpmlpars->dt,
                    s,   e,
                    acdpmlpars->c0,  acdpmlpars->c1,  acdpmlpars->c2,
                    acdpmlpars->lbc, acdpmlpars->rbc);
    }
    else if(acdpmlpars->k == 4){
        acdpml_2d_8(uc2, up2, csq2,
                    phi1, phi0,
                    acdpmlpars->dp1, acdpmlpars->dp0,
                    acdpmlpars->di,  acdpmlpars->dt,
                    s,   e,
                    acdpmlpars->c0,  acdpmlpars->c1,  acdpmlpars->c2,
                    acdpmlpars->c3,  acdpmlpars->c4,
                    acdpmlpars->lbc, acdpmlpars->rbc);
    }
	else {
	    fprintf(stderr,"ERROR: acd_step\n");
	    fprintf(stderr,"called with half-order != 1\n");
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
    /*else if (ndim == 3) {
    
	uc3   = (dom->_s)[D_UC ]._s3;
	up3   = (dom->_s)[D_UP ]._s3;
	csq3  = (dom->_s)[D_CSQ]._s3;

	cerr<<"timestep order = "<< acdpars->k<<endl;

	// 2nd order case
	if (acdpars->k == 1) {
	    acd_3d_2(uc3, up3, csq3, 
		     s, e, 
		     acdpars->c0, 
		     acdpars->c1);
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
    }*/
    else {
	fprintf(stderr,"ERROR: acdpml_step\n");
	fprintf(stderr,"called with space dim != 2\n");
	return E_BADINPUT;
    }
  
    return 0;
}

/*----------------------------------------------------------------------------*/
void acdpml_check(RDOM * dom, void * specs, FILE * stream) {
  // arrays to test: only csq = rarr[0]
  size_t n = 0; // length of data array
  ireal vmax = -FLT_MAX;
  ireal vmin = FLT_MAX;

  // acd struct
  ACDPML_TS_PARS * acdpmlpars = (ACDPML_TS_PARS *)specs;
  
  // extract dimn info
  ra_a_datasize(&(dom->_s[0]),&n);

  // max & min
  for (int i=0;i<(int)n;i++) {
    vmax=iwave_max((dom->_s[0]._s0)[i],vmax);
    vmin=iwave_min((dom->_s[0]._s0)[i],vmin);
    if (vmax<0.0f || sqrt(vmax)>acdpmlpars->cmax ||
	vmin<0.0f || sqrt(vmin)<acdpmlpars->cmin) {
      RVLException e;
      e<<"Error: input csq at index "<<i<<" = "<<(dom->_s[0]._s0)[i]<<"\n";
      e<<"  out of bounds ["<<acdpmlpars->cmin<<", "<<acdpmlpars->cmax<<"]\n";
      throw e;
    }
  }
 
}
