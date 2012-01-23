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
#ifndef CLIQUE_NUMERIC_FRONT_LDL_HPP
#define CLIQUE_NUMERIC_FRONT_LDL_HPP 1

namespace clique {
namespace numeric {

template<typename F>
void LocalFrontLDL
( Orientation orientation, Matrix<F>& AL, Matrix<F>& ABR );

template<typename F>
void DistFrontLDL
( Orientation orientation, DistMatrix<F,MC,MR>& AL, DistMatrix<F,MC,MR>& ABR );

namespace internal {
template<typename F>
void DistFrontLDLGeneral
( Orientation orientation, DistMatrix<F,MC,MR>& AL, DistMatrix<F,MC,MR>& ABR );

template<typename F>
void DistFrontLDLSquare
( Orientation orientation, DistMatrix<F,MC,MR>& AL, DistMatrix<F,MC,MR>& ABR );
} // namespace internal

} // namespace numeric
} // namespace clique

#endif /* CLIQUE_NUMERIC_FRONT_LDL_HPP */

