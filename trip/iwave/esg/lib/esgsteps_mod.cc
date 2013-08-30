#include "esgsteps_mod.h"

#ifdef __ESG_STEPS_MOD__

#ifdef _OPENMP
#include <omp.h>
#endif /*_OPENMP*/

// #define MARIO_VERBOSE
// #define SUPER_MARIO_VERBOSE

/*----------------------------------------------------------------------------*/
int esg_step_s(RDOM * dom, void * pars){
/*----------------------------------------------------------------------------*/
#ifdef MARIO_VERBOSE
	fprintf(stderr,">>>> Inside esg_step_s\n");
#endif

	int ndim = dom->_s[0].ndim;
	int err = 0;
	ESGN_TS_PARS * esgnp = (ESGN_TS_PARS *)pars;
	int iarr, idim;

// #ifdef SUPER_MARIO_VERBOSE
// 	fprintf(stderr,"Printing max of field before updating: ");
// 	IPNT ran; 	    //numb of points per axis
// 	unsigned long ntot; //tot numb of points in field
// 	ireal fieldmax;	    //maximum of P0
// 	ireal coeffmax;
// 	int i;
// 
// 	rd_size(dom,D_S0,ran);
// 	ntot = 1;
// 	for (i=0;i<ndim;i++) ntot*=ran[i];
// 	//Computing max of field
// /*	fieldmax = dom->_s[D_S0]._s0[0];*/
// 	coeffmax = dom->_s[D_MS0]._s0[0];
// 	for (i=1;i<ntot;i++){
// 		fieldmax = iwave_max(fieldmax,dom->_s[D_S0]._s0[i]);
// /*		coeffmax = iwave_max(coeffmax,dom->_s[D_MS0]._s0[i]);*/
// 	}
// 	fprintf(stderr,"fieldmax = %g\n",fieldmax);
// 	fprintf(stderr,"coeffmax = %g\n",coeffmax);
// #endif

	/* Updating pressures, need only to send D_P0 since pressures are all updated at once.*/
	iarr = D_P[0];
	//nothing implemented in 1D or 3D
	if (ndim != 2) {
		fprintf(stderr,"ERROR: time step module for pressure\n");
		fprintf(stderr,"       not implemented for %dD\n",ndim);
		return E_BADINPUT;
	}
	//2D case
	if (ndim==2) err = esg_2d(dom,iarr,pars);
// 	fprintf(stderr,"after ns, err = %d\n",err);


	/* Updating shear stresses */
	for (idim=0; idim<ndim; idim++){
		iarr = D_S[idim];
		//nothing implemented in 1D or 3D
		if (ndim != 2) {
			fprintf(stderr,"ERROR: time step module for shear stress\n");
			fprintf(stderr,"       not implemented for %dD\n",ndim);
			return E_BADINPUT;
		}
		//2D case
		if (ndim==2) err = esg_2d(dom,iarr,pars);
	}
// 	fprintf(stderr,"after ss, err = %d\n",err);


//MB added to check that step does not yield 0.
// #ifdef SUPER_MARIO_VERBOSE
// 	fprintf(stderr,"Printing max field after updating: ");
// 	//Computing max of P0 field
// 	fieldmax = dom->_s[D_S0]._s0[0];
// 	for (i=1;i<ntot;i++) fieldmax = iwave_max(fieldmax,dom->_s[D_S0]._s0[i]);
// 	fprintf(stderr,"fieldmax = %g\n",fieldmax);
// #endif


	return err;
}

/*----------------------------------------------------------------------------*/
int esg_step_v(RDOM * dom, void * pars){
/*----------------------------------------------------------------------------*/
#ifdef MARIO_VERBOSE
	fprintf(stderr,">>>> Inside esg_step_v\n");
#endif
	int ndim = dom->_s[0].ndim;
	int err = 0;
	ESGN_TS_PARS * esgnp = (ESGN_TS_PARS *)pars;
	int iarr, idim;

// #ifdef SUPER_MARIO_VERBOSE
// 	fprintf(stderr,"Printing max of V0 field before updating: ");
// 	IPNT ran; 	    //numb of points per axis
// 	unsigned long ntot; //tot numb of points in field
// 	ireal vmax;	    //maximum of V0
// 	int i;
// 
// 	rd_size(dom,D_V0,ran);
// 	ntot = 1;
// 	for (i=0;i<ndim;i++) ntot*=ran[i];
// 	//Computing max of V0 field
// 	vmax = dom->_s[D_V0]._s0[0];
// 	for (i=1;i<ntot;i++) vmax = iwave_max(vmax,dom->_s[D_V0]._s0[i]);
// 	fprintf(stderr,"vmax = %g\n",vmax);
// #endif

	// Updating velocities
	for (idim=0; idim<ndim; idim++){
		iarr = D_V[idim];
		//nothing implemented in 1D or 3D
		if (ndim != 2) {
			fprintf(stderr,"ERROR: time step module for velocities\n");
			fprintf(stderr,"       not implemented for %dD\n",ndim);
			return E_BADINPUT;
		}
		//2D case
		if (ndim==2) err = esg_2d(dom,iarr,pars); //esgn_gts2d_24(dom,iarr,pars);
	}

// /*MB added to check that step does not yield 0.*/
// #ifdef SUPER_MARIO_VERBOSE
// 	fprintf(stderr,"Printing max V0 after updating: ");
// 	//Computing max of V0 field
// 	vmax = dom->_s[D_V0]._s0[0];
// 	for (i=1;i<ntot;i++) vmax = iwave_max(vmax,dom->_s[D_V0]._s0[i]);
// 	fprintf(stderr,"vmax = %g\n",vmax);
// #endif

	return err;

}

/*----------------------------------------------------------------------------*/
int esg_step(RDOM * dom, int iv, void * pars){
/*----------------------------------------------------------------------------*/
#ifdef MARIO_VERBOSE
	fprintf(stderr,">>>> Inside esg_step, for iv = %d\n",iv);
#endif
    	if ( iv == 0 ) return esg_step_s(dom, pars);
    	if ( iv == 1 ) return esg_step_v(dom, pars);
    	return 0;
}

#endif /*__ESG_STEPS_MOD__*/