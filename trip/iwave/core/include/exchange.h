/* 
exchange.h
Igor Terentyev.
********************************************************************************
Function to compute allocated and exchange domains' coordinates.
Domains are declared, not allocated.
*/
/*============================================================================*/

#ifndef __EXCHANGE_H_
#define __EXCHANGE_H_
/*----------------------------------------------------------------------------*/

#include "utils.h"
#include "stencil.h"
#include "rdomain.h"

/*----------------------------------------------------------------------------*/
/*
Computes allocated and receive areas' global coordinates.
Number of receive areas is given by ex_rcvn(...) function.
Allocated area is a minimum cube containing the union of the computational
and receive areas.
In case the computational area is empty, all receive areas will be empty and
allocated area will be as usual (I did it this way because it is useful 
to compute arrays that are not exchanged, e.g. density, bulk modulus).

int iarr            :  array number.
STENCIL *sten       :  stencil pointer.
RDOM *dom           :  input domain (declared/allocated memory).
IPNT a_gs, a_ge     :  allocated domain global indices.
IPNT r_gs[], r_ge[] :  receive domains' global indices.

int return          :  error code.
*/
int ex_compute(int iarr, STENCIL *sten, RDOM *dom,
               IPNT a_gs, IPNT a_ge, IPNT r_gs[], IPNT r_ge[], int frcvempty[]);
/*----------------------------------------------------------------------------*/

#endif /*__EXCHANGE_H_*/
