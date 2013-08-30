#include "esgsteps.h"

#ifdef __ESG_STEPS__
/*============================================================================*/

#ifdef _OPENMP
#include <omp.h>
#endif /*_OPENMP*/

// #define MARIO_VERBOSE

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

#ifdef MARIO_VERBOSE
	fprintf(stderr,"Printing max of P0 field before updating: ");
	IPNT ran; 	    //numb of points per axis
	unsigned long ntot; //tot numb of points in field
	ireal pmax;	    //maximum of P0
	int i;

	rd_size(dom,D_P0,ran);
	ntot = 1;
	for (i=0;i<ndim;i++) ntot*=ran[i];
	//Computing max of P0 field
	pmax = dom->_s[D_P0]._s0[0];
	for (i=1;i<ntot;i++) pmax = iwave_max(pmax,dom->_s[D_P0]._s0[i]);
	fprintf(stderr,"pmax = %g\n",pmax);
#endif
	/* Updating pressures, need only to send D_P0 since pressures are all updated at once. */
	iarr = D_P[0];
	//nothing implemented in 1D
	if (ndim == 1) {
		fprintf(stderr,"ERROR: time step module\n");
		fprintf(stderr,"not implemented for 1D\n");
		return E_BADINPUT;
	}
	//nothing implemented for k=1, i.e., 2-2 scheme
	if (esgnp->k==1) {
		fprintf(stderr,"ERROR: 2-%d scheme not implemented for %dD\n",esgnp->k*2,ndim);
		return E_BADINPUT;
	}
	//k=2, i.e., 2-4 scheme case
	if (esgnp->k==2) {
	  //fprintf(stderr," In k=2 RN \n");
		if (ndim==2) err = esgn_gts2d_24(dom,dom,dom,iarr,pars);
		if (ndim==3) fprintf(stderr,"ERROR: 2-%d scheme not implemented for %dD\n",esgnp->k*2,ndim);
		// temp: err = esgn_gts3d_24(dom,dom,dom,iarr,pars);
	}
	//k=5, i.e., 2-10 scheme case
	if (esgnp->k==5) {
	  fprintf(stderr," In k=5 RN \n");
		if (ndim==2) err = esgn_gts2d_210(dom,dom,dom,iarr,pars);
		if (ndim==3) fprintf(stderr,"ERROR: 2-%d scheme not implemented for %dD\n",esgnp->k*2,ndim);
		// temp: err = esgn_gts3d_210(dom,dom,dom,iarr,pars);
	}
	//case where k!= 2 or 5
	if(esgnp->k !=2 && esgnp->k !=5){

		fprintf(stderr,"ERROR: time step module\n");
		fprintf(stderr,"order = %d not implemented for 2D or 3D\n",esgnp->k*2);
	       	return E_BADINPUT;
	}


	// Updating shear stresses
	for (idim=0; idim<ndim; idim++){
		iarr = D_S[idim];
		//nothing implemented in 1D
		if (ndim == 1) {
			fprintf(stderr,"ERROR: time step module\n");
			fprintf(stderr,"not implemented for 1D\n");
			return E_BADINPUT;
		}
		//nothing implemented for k=1, i.e., 2-2 scheme
		if (esgnp->k==1) {
			fprintf(stderr,"ERROR: 2-%d scheme not implemented for %dD\n",esgnp->k*2,ndim);
			return E_BADINPUT;
		}
		//k=2, i.e., 2-4 scheme case
		if (esgnp->k==2) {
			if (ndim==2) err = esgn_gts2d_24(dom,dom,dom,iarr,pars);
			if (ndim==3) err = esgn_gts3d_24(dom,dom,dom,iarr,pars); 
				//fprintf(stderr,"ERROR: 2-%d scheme not implemented for %dD\n",esgnp->k*2,ndim);
		}
		//k=5, i.e., 2-10 scheme case
		if (esgnp->k==5) {
			if (ndim==2) err = esgn_gts2d_210(dom,dom,dom,iarr,pars);
			if (ndim==3) err = esgn_gts3d_210(dom,dom,dom,iarr,pars);
				//fprintf(stderr,"ERROR: 2-%d scheme not implemented for %dD\n",esgnp->k*2,ndim);
		}
		//case where k!= 2 or 5
		if(esgnp->k !=2 && esgnp->k !=5){
			fprintf(stderr,"ERROR: time step module\n");
			fprintf(stderr,"order = %d not implemented for 2D or 3D\n",esgnp->k*2);
			return E_BADINPUT;
		}
	}
//MB added to check that step does not yield 0.
#ifdef MARIO_VERBOSE
	fprintf(stderr,"Printing max P0 after updating: ");
	//Computing max of P0 field
	pmax = dom->_s[D_P0]._s0[0];
	for (i=1;i<ntot;i++) pmax = iwave_max(pmax,dom->_s[D_P0]._s0[i]);
	fprintf(stderr,"pmax = %g\n",pmax);
#endif
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

	// Updating velocities
	for (idim=0; idim<ndim; idim++){
		iarr = D_V[idim];
		//nothing implemented in 1D
		if (ndim == 1) {
			fprintf(stderr,"ERROR: time step module\n");
			fprintf(stderr,"not implemented for 1D\n");
			return E_BADINPUT;
		}
		//nothing implemented for k=1, i.e., 2-2 scheme
		if (esgnp->k==1) {
			fprintf(stderr,"ERROR: 2-%d scheme not implemented for %dD\n",esgnp->k*2,ndim);
			return E_BADINPUT;
		}
		//k=2, i.e., 2-4 scheme case
		if (esgnp->k==2) {
			if (ndim==2) err = esgn_gts2d_24(dom,dom,dom,iarr,pars);
			if (ndim==3) fprintf(stderr,"ERROR: 2-%d scheme not implemented for %dD\n",esgnp->k*2,ndim);
				// temp: err = esgn_gts3d_24(dom,dom,dom,iarr,pars);
		}
		//k=5, i.e., 2-10 scheme case
		if (esgnp->k==5) {
			if (ndim==2) err = esgn_gts2d_210(dom,dom,dom,iarr,pars);
			if (ndim==3) fprintf(stderr,"ERROR: 2-%d scheme not implemented for %dD\n",esgnp->k*2,ndim);
				// temp: err = esgn_gts3d_210(dom,dom,dom,iarr,pars);  
		}
		//case where k!= 2 or 5
		if(esgnp->k !=2 && esgnp->k !=5){
			fprintf(stderr,"ERROR: time step module\n");
			fprintf(stderr,"order = %d not implemented for 2D or 3D\n",esgnp->k*2);
			return E_BADINPUT;
		}
	}
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
#endif /*__ESG_STEPS__*/