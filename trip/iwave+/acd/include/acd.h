#ifndef __IWAVE_ACD__
#define __IWAVE_ACD__

#include "fd.h"
#include "movie.h"

/**----------------------------------------------------------------------------
 * \f$ x=y \nabla u \f$ 
 *	
 * Parameters for time step function.
 *----------------------------------------------------------------------------
 */
typedef struct {
    ireal dt;      /* time step - copied from IMODEL.tsinfo */
    RPNT lam;      /* courant params */
    int k;         /* scheme order */
    int ndim;      /* dimension, copied from IMODEL.grid */
    IPNT lbc;      /* flag left boundary conditions */
    IPNT rbc;      /* flag right boundary conditions */
    /* coefficient arrays for FD schemes - set in readschemeinfo as they
    // are data-dependent
    // encoding (as RPNT): c[diff index][half-order]
    // thus for 3D, c23[0] is the coefficient of u[.][.][.+2]-u[.][.][.-1]  
    // or u[.][.][.+1]-u[.][.][.-2] in a sixth-order
    // staggered difference formula */
    ireal c0;
    RPNT c1;
    RPNT c2;
    RPNT c3; 
    RPNT c4;
} ACD_TS_PARS;  

/*
  Indices of arrays - const density acoustic. 
  D_UC  :  current acoustic potential
  D_UP  :  past acoustic potential
  D_CSQ :  square velocity (multiplier of Laplacian)
*/
#define D_UC  0
#define D_UP  1
#define D_CSQ 2

/* default fraction of max time step */
#define CFL_DEF 0.95

/* this declaration must be publicly available to external drivers:
   modelinit is passed as function pointer argument to iwave_construct. 
   On return, all other FD_MODEL function pointers are initialized and
   point to functions defined in this package - this is the only 
   package function which must be used before being assigned to an 
   FD_MODEL data member.
*/
int acd_modelinit(PARARRAY *pars, FILE *stream, IMODEL *model);

/** Movie constructor - must be public for similar reasons */
int acd_movie_construct(MOVIE * mt, FILE * stream);

#endif /*__ACD_H_*/
