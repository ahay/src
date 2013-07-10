#ifndef __GIWAVE_H_
#define __GIWAVE_H_

#include "iwave.h"
/* included to enable access to MAXPATHLEN */
#include <sys/param.h>

/*---------------*/
#ifdef __cplusplus
extern "C" {
#endif
  
  /**
     Generalized time step function, adequate to express time step and
     all first derivatives thereof. 

     @param[out] pd - RDOM storing perturbation dynamic and static fields to be upated, 
     @param[in] rd - RDOM storing reference dynamic and static fields
     @param[in] ia - index of array to be updated
     @param[in] fdpars - opaque FD_MODEL parameter object 
     @return 0 on success, else error code
     
     called in giwave_run. 
     
  */
  typedef int (*GEN_TIMESTEP_FUN)(RDOM * pd, RDOM * rd, int ia, void* fdpars);
  
  /**
     Generalized time step function to advance 2nd derivative, either fwd or adjoint mode

     @param[out] dd - RDOM storing 2nd order perturbation of dynamic fields to be upated, 
     @param[in] d1 - RDOM storing 1st perturbation dynamic and static fields
     @param[in] d2 - RDOM storing 2nd perturbation dynamic and static fields
     @param[in] rd - RDOM storing reference dynamic and static fields
     @param[in] ia - index of array to be updated
     @param[in] fdpars - opaque FD_MODEL parameter object 
     @return 0 on success, else error code
     
     called in giwave_run. 
     
  */
  typedef int (*GEN_TIMESTEP_FUN2)(RDOM * dd, RDOM * d1, RDOM * d2, RDOM * rd, int ia, void* fdpars);

  /* FD_MODEL initialization function 
     
  @param[in] pars - parameter array, assumed initialized.
  
  @param[in] stream - verbose output stream
  
  @param[out] mdl - \ref IMODEL whose static arrays are being initialized
  
  @return - 0 on success, else error code.
  */

  typedef int (*FDM_INIT_FUN)(PARARRAY * pars, 
			      FILE * stream, 
			      IMODEL * mdl);
  
  /** generalized update function - returns true 
      if field ia is updated at substep iv
  */
  typedef int (*GEN_UPDATE_FUN)(int ia, int iv, const IMODEL * m);

  /* --------------------------------------------------------------- */

  /** Abstract time-stepping simulator with derivative and adjoint derivative.
      Must be defined to implement inversion module. 

      ABSTRACT: must be defined in any instantiable subtype

      Definition via concrete definition of the gfd_model_init
      function, which assigns all data members (function pointers).
  */  
  typedef struct GFD_MODEL {

    /** main constructor - should do at least these three things: - 

    <ul> 
    <li> assign all function pointers in GFD_MODEL, including this
    one </li>
    
    <li> call fd_model_init function in initialize mdl->FD_MODEL
    
    </ul>
    
    */
    int (*gfd_model_init)(struct GFD_MODEL * gfdm);

    /** FD_MODEL constructor - passed to iwave_construct */
    FDM_INIT_FUN fd_model_init;

    /** timestep function - forward step, only one needed for pure
	modeling. called in iwave_run: initialize by (1) calling
	fd_model_init to initialize FD_MODEL, (2) copying
	FD_MODEL.tsf */
    TIMESTEP_FUN tsf;

    /** generalized timestep function - defines model derivative of
	step */
    GEN_TIMESTEP_FUN tsfm;

    /** generalized adjoint timestep function - defines model adjoint
	derivative of step */
    GEN_TIMESTEP_FUN tsam;

    /** generalized timestep function, 2nd order - defines model 2nd derivative of
	step */
    GEN_TIMESTEP_FUN2 tsfm2;

    /** generalized adjoint timestep function, 2nd order - defines model adjoint
	2nd derivative of step */
    GEN_TIMESTEP_FUN2 tsam2;

    /** update flag for lin */
    GEN_UPDATE_FUN udfm;

    /** update flag for adj */
    GEN_UPDATE_FUN udam;

    /** for forward linearized time loop, specifies which step (it)
	and substep (iv) should be used as the target for a given
	substep of the perturbational time step
	@param[out] it - target time step for reference field
	@param[out] iv - target time stubstep for reference field
	@param[in] m - perturbational model, tsind data member 
	contains necessary time info - it, iv outputs are 
	functions of m->tsind.it, m->tsind.iv
    */
    void (*linrefsubstep)(int* it,int* iv, const IMODEL* m);

    /** for adjoint linearized time loop, specifies which step (it)
	and substep (iv) should be used as the target for a given
	substep of the perturbational time step
	@param[out] it - target time step for reference field
	@param[out] iv - target time stubstep for reference field
	@param[in] m - perturbational model, tsind data member
	contains necessary time info - it, iv outputs are 
	functions of m->tsind.it, m->tsind.iv
    */
    void (*adjrefsubstep)(int* it,int* iv, const IMODEL* m);

  } GFD_MODEL;
  
  /** GFD_MODEL initialization function 

      CONCRETE: defined by reference to GFD_MODEL data members 
      
      @param[in] gfdm - generalized fd model struct, includes initialization
      function for fd model struct and generalized time step functions
      
      @return - 0 on success, else error code.
  */  
  typedef int (*GFDM_INIT_FUN)(GFD_MODEL * gfdm);

  /** Linearized and adjoint linearized modeling
      _fwd = true:   linearized;
      _fwd = false:  adjoint adjoint linearized

      CONCRETE: defined by reference to GFD_MODEL data members 

      @param[out] pstate - perturbation IWAVE object 
      @param[in] rstate - reference IWAVE object 
      @param[in] ts - generalized time step, either lin or adj
      @param[in] ud - flags fields updated in current substep
      @param[out] stream - verbose output
      
      @return - 0 on success, else error code.
  */
  int giwave_dmod(IWAVE * pstate, 
		  IWAVE * rstate, 
		  GEN_TIMESTEP_FUN ts,
		  //		  GEN_UPDATE_FUN ud,
		  FILE * stream);
  
/** 2nd derivative and and adjoint 2nd deriv  modeling
      _fwd = true:   linearized;
      _fwd = false:  adjoint;

      CONCRETE: defined by reference to GFD_MODEL data members 

      @param[out] ddstate - 2nd deriv IWAVE object
      @param[in] d1state - 1st linearized IWAVE object 
      @param[in] d2state - 2nd linearized IWAVE object  
      @param[in] rstate - reference IWAVE object 
      @param[in] ts - generalized time step, either lin or adj
      @param[in] ud - flags fields updated in current substep
      @param[out] stream - verbose output
      
      @return - 0 on success, else error code.
  */

  int giwave_dmod2(IWAVE * ddstate, 
		   IWAVE * d1state,
		   IWAVE * d2state,
		   IWAVE * rstate, 
		   GEN_TIMESTEP_FUN2 ts,
		   //		   GEN_UPDATE_FUN ud,
		   FILE * stream);

  /** zero receive buffers in preparation for adjoint
      accumulation. Should be called before giwave_dmod, and update
      called before giwave_synch

      CONCRETE: defined by reference to GFD_MODEL data members 

      @param pstate - perturbation IWAVE object 
      @param ud - update rule (for adjoint)
      @param stream - verbose output
      @return - 0 on success
  */
  int giwave_zero_recv(IWAVE * pstate,
		       GEN_UPDATE_FUN ud,
		       FILE * stream);

  /** synchronize perturbation fields after update 

      CONCRETE: defined by reference to GFD_MODEL data members 

      @param[out] pstate - perturbation IWAVE object 
      @param[in] ud - flags fields updated in current substep
      @param[in] fwd - 1 for fwd lin, 0 for adj lin
      @param[out] stream - verbose output
      
      @return - 0 on success, else error code.
  */
  int giwave_synch(IWAVE * pstate, 
		   GEN_UPDATE_FUN ud,
		   int fwd,
		   FILE * stream);
    
  /** initialize dynamic fields at time(it) from data files 
      or zero out dynamic fields at time(itstart) 

      CONCRETE: defined by reference to GFD_MODEL data members 
  */ 
  int giwave_dynamic_init(IWAVE * state,
			  int it,         /* time to initialize dynamic fields*/
			  int itoff);     /* it - itstart */
  
  /** store dynamic fields at time (itcheck = itstart + itoff) 
      to data files, such as 
      $DATAPATH/statecheckpoint_itcheck_proc_?_array_?.bin

      CONCRETE: defined by reference to GFD_MODEL data members 
    */
  int iwave_dynamic_takeshot(IWAVE * state, 
			     int itoff ); /* itoff = itcheck - itstart */
  
  /** remove unused data files hoding dynamic fields at time 
      (itcheck = itstart + itoff),
      such as $DATAPATH/statecheckpoint_itcheck_proc_?_array_?.bin

      CONCRETE: defined by reference to GFD_MODEL data members 
    */
  int iwave_remove_checkfile(IWAVE * state, int itoff);

  /** generate checkfile names at time (itcheck = itstart + itoff) 
      for array (iarr)

      CONCRETE: defined by reference to GFD_MODEL data members 
 */
  char * iwave_get_a_checkfilename(int iarr,int itoff); 

#ifdef __cplusplus
}
#endif

#endif

