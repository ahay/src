#include "ansol_esgsteps.h"

/*============================================================================*/

#ifdef _OPENMP
#include <omp.h>
#endif /*_OPENMP*/

/*#define MARIO_VERBOSE*/
// #define SUPER_MARIO_VERBOSE

/*----------------------------------------------------------------------------*/
int ansol_HI_esg_step(RDOM * dom, int iv, void * pars){
/*----------------------------------------------------------------------------*/
#ifdef MARIO_VERBOSE
	fprintf(stderr,">>>> Inside ansol_HI_esg_step, iv=%d\n",iv);
#endif

	int ndim = dom->_s[0].ndim;
	int err = 0;
	int iarr, idim;

	if (iv==0) {

// #ifdef SUPER_MARIO_VERBOSE
// 		fprintf(stderr,"Printing max of P0 field before updating: ");
// 		IPNT ran; 	    //numb of points per axis
// 		unsigned long ntot; //tot numb of points in field
// 		ireal pmax;	    //maximum of P0
// 		int i;
// 	
// 		//Computing max of P0 field
// 		rd_size(dom,D_P0,ran);
// 		ntot = 1;
// 		for (i=0;i<ndim;i++) 
// 			ntot*=ran[i];
// 		pmax = dom->_s[D_P0]._s0[0];
// 		for (i=1;i<ntot;i++) 
// 			pmax = iwave_max(pmax,dom->_s[D_P0]._s0[i]);
// 		fprintf(stderr,"pmax = %g\n",pmax);
// #endif
		//----------------------------------------------------------------------------//
		//--Updating pressures -------------------------------------------------------//
		for (idim=0; idim<ndim; idim++){
			iarr = D_P[idim];
			//nothing implemented in 1D
			if (ndim == 1) {
				fprintf(stderr,"Error: in ansol_esg_step, not implemented for 1D\n");
				return E_BADINPUT;
			}
			if (ndim ==2) {
				ansol_HI_esg_ker2d(dom,iarr,pars);
			}
			if (ndim == 3){
				fprintf(stderr,"Error: in ansol_HI_esg_step, 3D case not yet implemented! ABORT!\n");
				exit(1);
// 				ansol_HI_esg_ker2d(dom,iarr,pars);
			}
		}
	
		//----------------------------------------------------------------------------//
		//--Updating shear stresses --------------------------------------------------//
		for (idim=0; idim<ndim; idim++){
			iarr = D_S[idim];
			//nothing implemented in 1D
			if (ndim == 1) {
				fprintf(stderr,"Error: in ansol_esg_step, not implemented for 1D\n");
				return E_BADINPUT;
			}
			if (ndim ==2) {
				ansol_HI_esg_ker2d(dom,iarr,pars);
			}
			if (ndim == 3){
				fprintf(stderr,"Error: in ansol_HI_esg_step, 3D case not yet implemented! ABORT!\n");
				exit(1);
// 				ansol_HI_esg_ker2d(dom,iarr,pars);
			}
		}

// #ifdef SUPER_MARIO_VERBOSE
// 		//Computing max of P0 field
// 		fprintf(stderr,"Printing max P0 after updating: ");
// 		pmax = dom->_s[D_P0]._s0[0];
// 		for (i=1;i<ntot;i++) pmax = iwave_max(pmax,dom->_s[D_P0]._s0[i]);
// 		fprintf(stderr,"pmax = %g\n",pmax);
// #endif

	}
	if (iv==1) {

#ifdef SUPER_MARIO_VERBOSE
		fprintf(stderr,"Printing max of V0 field before updating: ");
		IPNT ran; 	    //numb of points per axis
		unsigned long ntot; //tot numb of points in field
		ireal vmax;	    //maximum of V0
		int i;

		rd_size(dom,D_V0,ran);
		ntot = 1;
		for (i=0;i<ndim;i++) ntot*=ran[i];
		//Computing max of V0 field
		vmax = dom->_s[D_V0]._s0[0];
		for (i=1;i<ntot;i++) vmax = iwave_max(vmax,dom->_s[D_V0]._s0[i]);
		fprintf(stderr,"vmax = %g\n",vmax);
#endif

		//----------------------------------------------------------------------------//
		//--Updating velocities ------------------------------------------------------//
		for (idim=0; idim<ndim; idim++){
			iarr = D_V[idim];
			//nothing implemented in 1D
			if (ndim == 1) {
				fprintf(stderr,"Error: in ansol_esg_step, not implemented for 1D\n");
				return E_BADINPUT;
			}
			if (ndim ==2) {
				ansol_HI_esg_ker2d(dom,iarr,pars);
			}
			if (ndim == 3){
				fprintf(stderr,"Error: in ansol_HI_esg_step, 3D case not yet implemented! ABORT!\n");
				exit(1);
// 				ansol_HI_esg_ker2d(dom,iarr,pars);
			}
		}

/*MB added to check that step does not yield 0.*/
#ifdef SUPER_MARIO_VERBOSE
		fprintf(stderr,"Printing max V0 after updating: ");
		//Computing max of V0 field
		vmax = dom->_s[D_V0]._s0[0];
		for (i=1;i<ntot;i++) vmax = iwave_max(vmax,dom->_s[D_V0]._s0[i]);
		fprintf(stderr,"vmax = %g\n",vmax);
#endif

	}
	return err;
}


