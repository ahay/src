#ifndef _chain_h
#define _chain_h

#include "bigsolver.h"

/* chain
   -----
   Chains two operators, computing oper1{oper2{mod}} or its adjoint.
   The tmp[nt] array is used for temporary storage. */
void chain( operator oper1, operator oper2, bool adj, bool add, 
	    int nm, int nd, int nt, float* mod, float* dat, float* tmp);

void chain3 (operator oper1, operator oper2, operator oper3, 
	     bool adj, bool add, int nm, int nt1, int nt2, int nd, 
	     float* mod, float* dat, float* tmp1, float* tmp2);

/* array
   -----
   Constructs an array of two operators, computing {oper1{mod},oper2{mod}} 
   or its adjoint. */
void array( operator oper1, operator oper2, bool adj, bool add, 
	    int nm, int nd1, int nd2, float* mod, float* dat1, float* dat2);

void normal (operator oper, bool add, 
	     int nm, int nd, float *mod, float *dat, float *tmp);
#endif
