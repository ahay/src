/* 
exchangeinfo.h
Igor Terentyev.
*/
/*============================================================================*/
/**
 * \file exchangeinfo.h
 * Data exchange information
 */
#ifndef __EXCHANGEINFO_H_
#define __EXCHANGEINFO_H_
/*----------------------------------------------------------------------------*/

#include "utils.h"
#include "usempi.h"

/*----------------------------------------------------------------------------*/
/**
 * Data exchange information.
 *
 * Set by \ref ra_setexchangeinfo. buf points to _s in an \ref RARR.
 */
typedef struct s_EXCHANGEINFO
{
  /** pointer to a buffer for data exchange.*/
  void *buf;
  /** MPI_Datatype of data exchange */
  MPI_Datatype type;
  /** 
   * intermedia MPI_Datatype of data exchange for 3D problems. 
   * TODO: make it a local variable in \ref ra_setexchangeinfo ?!
   */
  MPI_Datatype type2;
} EXCHANGEINFO;
/*----------------------------------------------------------------------------*/
/**
 * Set data exchange information defaults.
 *
 * @param[in,out] einfo - (EXCHANGEINFO *) date exchange information
 * 
 * @return 0
 */
int ei_setnull(EXCHANGEINFO *einfo);
/*----------------------------------------------------------------------------*/
/**
 * Destroy date exchange info.
 * Free MPI_Datatype, then set defaults
 *
 * @param[in] einfo - (EXCHANGEINFO *) date exchange information
 *
 * @return 0
 */
int ei_destroy(EXCHANGEINFO *einfo);
/*----------------------------------------------------------------------------*/

#endif /*__EXCHANGEINFO_H_*/

