/* Chaining linear operators. */
/*
  Copyright (C) 2004 University of Texas at Austin
  
  This program is free software; you can redistribute it and/or modify
  it under the terms of the GNU General Public License as published by
  the Free Software Foundation; either version 2 of the License, or
  (at your option) any later version.
  
  This program is distributed in the hope that it will be useful,
  but WITHOUT ANY WARRANTY; without even the implied warranty of
  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
  GNU General Public License for more details.
  
  You should have received a copy of the GNU General Public License
  along with this program; if not, write to the Free Software
  Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA
*/

#include "_bool.h"
#include "_solver.h"
/*^*/

#include "chain.h"

void sf_chain( sf_operator oper1     /* outer operator */, 
	       sf_operator oper2     /* inner operator */, 
	       bool adj              /* adjoint flag */, 
	       bool add              /* addition flag */, 
	       int nm                /* model size */, 
	       int nd                /* data size */, 
	       int nt                /* intermediate size */, 
	       /*@out@*/ float* mod  /* [nm] model */, 
	       /*@out@*/ float* dat  /* [nd] data */, 
	       float* tmp            /* [nt] intermediate */) 
/*< Chains two operators, computing oper1{oper2{mod}} 
  or its adjoint. The tmp array is used for temporary storage. >*/
{
    if (adj) {
	oper1 (true, false, nt, nd, tmp, dat);
	oper2 (true, add, nm, nt, mod, tmp);
    } else {
	oper2 (false, false, nm, nt, mod, tmp);
	oper1 (false, add, nt, nd, tmp, dat);
    }
}

void sf_cchain( sf_coperator oper1         /* outer operator */, 
		sf_coperator oper2         /* inner operator */, 
		bool adj                   /* adjoint flag */, 
		bool add                   /* addition flag */, 
		int nm                     /* model size */, 
		int nd                     /* data size */, 
		int nt                     /* intermediate size */, 
		/*@out@*/ sf_complex* mod  /* [nm] model */, 
		/*@out@*/ sf_complex* dat  /* [nd] data */, 
		sf_complex* tmp            /* [nt] intermediate */) 
/*< Chains two complex operators, computing oper1{oper2{mod}} 
  or its adjoint. The tmp array is used for temporary storage. >*/
{
    if (adj) {
	oper1 (true, false, nt, nd, tmp, dat);
	oper2 (true, add, nm, nt, mod, tmp);
    } else {
	oper2 (false, false, nm, nt, mod, tmp);
	oper1 (false, add, nt, nd, tmp, dat);
    }
}

void sf_array( sf_operator oper1     /* top operator */, 
	       sf_operator oper2     /* bottom operator */, 
	       bool adj              /* adjoint flag */, 
	       bool add              /* addition flag */, 
	       int nm                /* model size */, 
	       int nd1               /* top data size */, 
	       int nd2               /* bottom data size */, 
	       /*@out@*/ float* mod  /* [nm] model */, 
	       /*@out@*/ float* dat1 /* [nd1] top data */, 
	       /*@out@*/ float* dat2 /* [nd2] bottom data */) 
/*< Constructs an array of two operators, 
  computing {oper1{mod},oper2{mod}} or its adjoint. >*/
{
    if (adj) {
	oper1 (true, add,  nm, nd1, mod, dat1);
	oper2 (true, true, nm, nd2, mod, dat2);
    } else {
	oper1 (false, add, nm, nd1, mod, dat1);
	oper2 (false, add, nm, nd2, mod, dat2);
    }
}

void sf_carray( sf_coperator oper1     /* top operator */, 
	       sf_coperator oper2     /* bottom operator */, 
	       bool adj              /* adjoint flag */, 
	       bool add              /* addition flag */, 
	       int nm                /* model size */, 
	       int nd1               /* top data size */, 
	       int nd2               /* bottom data size */, 
	       /*@out@*/ sf_complex* mod  /* [nm] model */, 
	       /*@out@*/ sf_complex* dat1 /* [nd1] top data */, 
	       /*@out@*/ sf_complex* dat2 /* [nd2] bottom data */) 
/*< Constructs an array of two complex operators, 
  computing {oper1{mod},oper2{mod}} or its adjoint. >*/
{
    if (adj) {
	oper1 (true, add,  nm, nd1, mod, dat1);
	oper2 (true, true, nm, nd2, mod, dat2);
    } else {
	oper1 (false, add, nm, nd1, mod, dat1);
	oper2 (false, add, nm, nd2, mod, dat2);
    }
}


void sf_normal (sf_operator oper /* operator */, 
		bool add         /* addition flag */, 
		int nm           /* model size */, 
		int nd           /* data size */, 
		float *mod       /* [nd] model */, 
		float *dat       /* [nd] data */, 
		float *tmp       /* [nm] intermediate */)
/*< Applies a normal operator (self-adjoint) >*/
{
    oper (true, false, nm, nd, tmp, mod);
    oper (false, add,  nm, nd, tmp, dat);
}

void sf_chain3 (sf_operator oper1 /* outer operator */, 
		sf_operator oper2 /* middle operator */, 
		sf_operator oper3 /* inner operator */, 
		bool adj          /* adjoint flag */, 
		bool add          /* addition flag */, 
		int nm            /* model size */, 
		int nt1           /* inner intermediate size */, 
		int nt2           /* outer intermediate size */, 
		int nd            /* data size */, 
		float* mod        /* [nm] model */, 
		float* dat        /* [nd] data */, 
		float* tmp1       /* [nt1] inner intermediate */, 
		float* tmp2       /* [nt2] outer intermediate */)
/*< Chains three operators, computing oper1{oper2{poer3{{mod}}} or its adjoint.
  The tmp1 and tmp2 arrays are used for temporary storage. >*/
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

/* 	$Id: chain.c 7624 2011-08-24 20:27:36Z sfomel $	 */
