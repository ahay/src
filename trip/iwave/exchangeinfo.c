/* Data exchange information */
/*************************************************************************

Copyright Rice University, 2008.
All rights reserved.

Permission is hereby granted, free of charge, to any person obtaining a
copy of this software and associated documentation files (the "Software"),
to deal in the Software without restriction, including without limitation
the rights to use, copy, modify, merge, publish, distribute, and/or sell
copies of the Software, and to permit persons to whom the Software is
furnished to do so, provided that the above copyright notice(s) and this
permission notice appear in all copies of the Software and that both the
above copyright notice(s) and this permission notice appear in supporting
documentation.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT OF THIRD PARTY
RIGHTS. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR HOLDERS INCLUDED IN THIS
NOTICE BE LIABLE FOR ANY CLAIM, OR ANY SPECIAL INDIRECT OR CONSEQUENTIAL
DAMAGES, OR ANY DAMAGES WHATSOEVER RESULTING FROM LOSS OF USE, DATA OR
PROFITS, WHETHER IN AN ACTION OF CONTRACT, NEGLIGENCE OR OTHER TORTIOUS
ACTION, ARISING OUT OF OR IN CONNECTION WITH THE USE OR PERFORMANCE OF
THIS SOFTWARE.

Except as contained in this notice, the name of a copyright holder shall
not be used in advertising or otherwise to promote the sale, use or other
dealings in this Software without prior written authorization of the
copyright holder.

**************************************************************************/
#include <stdlib.h>

#include "exchangeinfo.h"

#ifndef _sf_exchangeinfo_h

typedef struct
{
  /** pointer to a buffer for data exchange.*/
  void *buf;
#ifdef IWAVE_USE_MPI
  /** MPI_Datatype of data exchange */
  MPI_Datatype type;
  /** 
   * intermedia MPI_Datatype of data exchange for 3D problems. 
   * TODO: make it a local variable in \ref ra_setexchangeinfo ?!
   */
  MPI_Datatype type2;
#endif
} EXCHANGEINFO;
/*^*/

#endif

/* 
exchangeinfo.c
Igor Terentyev.
*/
/*============================================================================*/

int ei_setnull(EXCHANGEINFO *einfo)
/*< Set data exchange information defaults. >*/
{
	einfo->buf = NULL;
#ifdef IWAVE_USE_MPI
	einfo->type = einfo->type2 = MPI_DATATYPE_NULL;
#endif	
	return 0;
}
/*----------------------------------------------------------------------------*/

int ei_destroy(EXCHANGEINFO *einfo)
/*< Destroy date exchange info. >*/
{
#ifdef IWAVE_USE_MPI
    if ( einfo->type != MPI_DATATYPE_NULL ) MPI_Type_free(&(einfo->type));
    if ( einfo->type2 != MPI_DATATYPE_NULL ) MPI_Type_free(&(einfo->type2));
#endif    
    return ei_setnull(einfo);
}
/*----------------------------------------------------------------------------*/

