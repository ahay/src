/* Operator types for linear solvers. */
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

#ifndef _sf__solver_h
#define _sf__solver_h

#include "_bool.h"
#include "c99.h"

typedef void (*sf_operator)(bool,bool,int,int,float*,float*);
typedef void (*sf_solverstep)(bool,int,int,float*,
			   const float*,float*,const float*);
typedef void (*sf_weight)(int,const float*,float*);
/*^*/

typedef void (*sf_coperator)(bool,bool,int,int,sf_complex*,sf_complex*);
typedef void (*sf_csolverstep)(bool,int,int,sf_complex*,
			       const sf_complex*,sf_complex*,
			       const sf_complex*);
typedef void (*sf_cweight)(int,const sf_complex*,float*);
/*^*/

#endif
