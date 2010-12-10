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
#include <rsf.h>

/*^*/

#include "carray.h"


void carray( sf_coperator oper1     /* top operator */, 
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


