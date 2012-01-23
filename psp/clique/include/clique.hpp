/*
   Clique: a scalable implementation of the multifrontal algorithm

   Copyright (C) 2010-2011 Jack Poulson <jack.poulson@gmail.com>
   Copyright (C) 2011 Jack Poulson, Lexing Ying, and 
   The University of Texas at Austin
 
   This program is free software: you can redistribute it and/or modify
   it under the terms of the GNU General Public License as published by
   the Free Software Foundation, either version 3 of the License, or
   (at your option) any later version.
 
   This program is distributed in the hope that it will be useful,
   but WITHOUT ANY WARRANTY; without even the implied warranty of
   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
   GNU General Public License for more details.
 
   You should have received a copy of the GNU General Public License
   along with this program.  If not, see <http://www.gnu.org/licenses/>.
*/
#ifndef CLIQUE_HPP
#define CLIQUE_HPP 1

#include "clique/core/environment.hpp"
#include "clique/symbolic/symmetric_factorization.hpp"

#include "clique/numeric/front_diagonal_solve.hpp"
#include "clique/numeric/front_lower_solve.hpp"
#include "clique/numeric/front_lower_multiply.hpp"

#include "clique/numeric/symm_front_tree.hpp"

#include "clique/numeric/diagonal_solve.hpp"
#include "clique/numeric/lower_multiply.hpp"
#include "clique/numeric/lower_solve.hpp"

#include "clique/numeric/front_ldl.hpp"
#include "clique/numeric/ldl.hpp"
#include "clique/numeric/ldl_solve.hpp"

#endif /* CLIQUE_HPP */
