#include "esgsteps.h"
#include "esgn.h"

/*============================================================================*/

#ifdef _OPENMP
#include <omp.h>
#endif /*_OPENMP*/
/*----------------------------------------------------------------------------*/

int old_tsf_test(RDOM * unext, RDOM * ucurr, RDOM * coef) {
  if ((unext != ucurr)||(unext != coef)) {
    fprintf(stderr,"ERROR: time step module\n");
    fprintf(stderr,"you are trying to call a generalized time step function\n");
    fprintf(stderr,"(probably from an inversion application) which has not been defined\n");
    fprintf(stderr,"Check input parameters and timestep codes\n");
    return E_BADINPUT;
  }
  return 0;
}

int esg_ts(RDOM * unext, RDOM * ucurr, RDOM * coef, int iarr, int _fwd, void *pars) {
  
  int ndim = unext->_s[iarr].ndim;
  int err=0;
  ESGN_TS_PARS * esgnp = (ESGN_TS_PARS *)pars;

  if (ndim == 1) {
    fprintf(stderr,"ERROR: time step module\n");
    fprintf(stderr,"not implemented for 1D\n");
    return E_BADINPUT;
  }
  if (esgnp->k==1) {
    err = old_tsf_test(unext,ucurr,coef);
    if (err) { 
      fprintf(stderr,"RETURNED from 22 scheme\n");
      return err;
    }
    fprintf(stderr,"ERROR: time step module\n");
    fprintf(stderr,"order = %d not implemented for 1D\n",2*esgnp->k);
    return E_BADINPUT;
  }
  else if (esgnp->k==2) {
    err = old_tsf_test(unext,ucurr,coef);
    if (err) { 
      fprintf(stderr,"RETURNED from 24 scheme\n");
      return err;
    }
    if (ndim==2) return esgn_gts2d_24(unext,ucurr,coef,iarr,pars,_fwd);
    if (ndim==3) return esgn_gts3d_24(unext,ucurr,coef,iarr,pars,_fwd);
  }
  else if (esgnp->k<8) {
    err = old_tsf_test(unext,ucurr,coef);
    if (err) { 
      fprintf(stderr,"RETURNED from 2k scheme\n");
      return err;
    }    
    if (ndim==2) {
	/* if (esgnp->k==5) return esgn_gts2d_210(unext,ucurr,coef,iarr,pars,_fwd); */
	/* else { */
        fprintf(stderr,"ERROR: time step module\n");
        fprintf(stderr,"order = %d not implemented for 2D\n",2*esgnp->k);
        return E_BADINPUT;
        /* } */
    }
    if (ndim==3) {
      if (esgnp->k==5) return esgn_gts3d_210(unext,ucurr,coef,iarr,pars,_fwd);
      else {
        fprintf(stderr,"ERROR: time step module\n");
        fprintf(stderr,"order = %d not implemented for 3D\n",2*esgnp->k);
        return E_BADINPUT;
      }
    }
  }
  
  return E_NOTIMESTEP;
}

int esgts_fwd(RDOM * u, int iarr, void *pars) {
  return esg_ts(u,u,u,iarr,1,pars);
}

int esgts_adj(RDOM * u, int iarr, void *pars) {
  return esg_ts(u,u,u,iarr,0,pars);
}

int esgtsm_fwd(RDOM * unext, RDOM * ucurr, RDOM * coef, int iarr, void *pars) {
  return esg_ts(unext,ucurr,coef,iarr,1,pars);
}

int esgtsm_adj(RDOM * unext, RDOM * ucurr, RDOM * coef, int iarr, void *pars) {
  return esg_ts(unext,ucurr,coef,iarr,0,pars);
}

