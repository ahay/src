/* 
usempi.h
Igor Terentyev.
********************************************************************************
Flag for iwave to use MPI.
*/
/*============================================================================*/

#ifndef __USEMPI_H_
#define __USEMPI_H_
/*----------------------------------------------------------------------------*/

/*#define IWAVE_USE_MPI*/
/*#undef IWAVE_USE_MPI*/
/*----------------------------------------------------------------------------*/
#define MPICH_SKIP_MPICXX
#ifdef IWAVE_USE_MPI
	#include "mpi.h"
	
#else
/*---------------*/
	#define MPI_Comm int
	#define MPI_PROC_NULL -1
	#define MPI_Datatype int
	#define MPI_DATATYPE_NULL 0
	
/*	static inline int MPI_Type_free(MPI_Datatype *t) { *t = MPI_DATATYPE_NULL; return 0; } */
#endif
/*----------------------------------------------------------------------------*/

#endif /*__USEMPI_H_*/
