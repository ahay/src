#include "bigsolver.h"
#include "chain.h"

/* chain
   -----
   Chains two operators, computing oper1{oper2{mod}} or its adjoint.
   The tmp[nt] array is used for temporary storage. */
void sf_chain( sf_operator oper1, sf_operator oper2, bool adj, bool add, 
	    int nm, int nd, int nt, float* mod, float* dat, float* tmp) {
    if (adj) {
	oper1 (true, false, nt, nd, tmp, dat);
	oper2 (true, add, nm, nt, mod, tmp);
    } else {
	oper2 (false, false, nm, nt, mod, tmp);
	oper1 (false, add, nt, nd, tmp, dat);
    }
}

/* array
   -----
   Constructs an array of two operators, computing {oper1{mod},oper2{mod}} 
   or its adjoint. */
void sf_array( sf_operator oper1, sf_operator oper2, bool adj, bool add, 
	       int nm, int nd1, int nd2, float* mod, float* dat1, float* dat2) {
    if (adj) {
	oper1 (true, add,  nm, nd1, mod, dat1);
	oper2 (true, true, nm, nd2, mod, dat2);
    } else {
	oper1 (false, add, nm, nd1, mod, dat1);
	oper2 (false, add, nm, nd2, mod, dat2);
    }
}

void sf_normal (sf_operator oper, bool add, 
		int nm, int nd, float *mod, float *dat, float *tmp)
{
    oper (true, false, nm, nd, tmp, mod);
    oper (false, add,  nm, nd, tmp, dat);
}

void sf_chain3 (sf_operator oper1, sf_operator oper2, sf_operator oper3, 
		bool adj, bool add, int nm, int nt1, int nt2, int nd, 
		float* mod, float* dat, float* tmp1, float* tmp2)
{
    if (adj) {
	oper1 (true, false, nt2, nd, tmp2, dat);
	oper2 (true, false, nt1, nt2, tmp1, tmp2);
	oper3 (true, add,   nm, nt1, mod, tmp1);
    } else {
	oper3 (false, false, nm, nt1, mod, tmp1);
	oper2 (false, false, nt1, nt2, tmp1, tmp2);
	oper1 (false, add, nt2, nd, tmp2, dat);
    }
}

/* 	$Id: chain.c,v 1.1 2003/10/21 15:12:39 fomels Exp $	 */

