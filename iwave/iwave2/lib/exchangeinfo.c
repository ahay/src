/* 
exchangeinfo.c
Igor Terentyev.
*/
/*============================================================================*/

#include "utils.h"
#include "exchangeinfo.h"
/*----------------------------------------------------------------------------*/

int ei_setnull(EXCHANGEINFO *einfo)
{
	einfo->buf = NULL;
	einfo->type = einfo->type2 = MPI_DATATYPE_NULL;
	
	return 0;
}
/*----------------------------------------------------------------------------*/

int ei_destroy(EXCHANGEINFO *einfo)
{
	if ( einfo->type != MPI_DATATYPE_NULL ) MPI_Type_free(&(einfo->type));
	if ( einfo->type2 != MPI_DATATYPE_NULL ) MPI_Type_free(&(einfo->type2));
	
	return ei_setnull(einfo);
}
/*----------------------------------------------------------------------------*/

