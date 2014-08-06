#ifndef __GIWAVE_H_
#define __GIWAVE_H_

#include "iwave.h"
/* included to enable access to MAXPATHLEN */
#include <sys/param.h>

/*---------------*/
#ifdef __cplusplus
extern "C" {
#endif
  
  /** generalized update function - returns true 
      if field ia is updated at substep iv
  */
  typedef int (*GEN_UPDATE_FUN)(int ia, int iv, const IMODEL * m);

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

