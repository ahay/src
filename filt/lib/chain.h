#ifndef _sf_chain_h
#define _sf_chain_h

#include "bigsolver.h"

/* chain
   -----
   Chains two operators, computing oper1{oper2{mod}} or its adjoint.
   The tmp[nt] array is used for temporary storage. */
void sf_chain( sf_operator oper1, sf_operator oper2, bool adj, bool add, 
	       int nm, int nd, int nt, float* mod, float* dat, float* tmp);

void sf_chain3 (sf_operator oper1, sf_operator oper2, sf_operator oper3, 
		bool adj, bool add, int nm, int nt1, int nt2, int nd, 
		float* mod, float* dat, float* tmp1, float* tmp2);

/* array
   -----
   Constructs an array of two operators, computing {oper1{mod},oper2{mod}} 
   or its adjoint. */
void sf_array( sf_operator oper1, sf_operator oper2, bool adj, bool add, 
	       int nm, int nd1, int nd2, float* mod, float* dat1, float* dat2);

void sf_normal (sf_operator oper, bool add, 
		int nm, int nd, float *mod, float *dat, float *tmp);
#endif

/* 	$Id: chain.h,v 1.1 2003/10/21 15:12:39 fomels Exp $	 */
