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
#ifndef CLIQUE_NUMERIC_FRONT_LOWER_MULTIPLY_HPP
#define CLIQUE_NUMERIC_FRONT_LOWER_MULTIPLY_HPP 1

namespace clique {
namespace numeric {

template<typename F>
void LocalFrontLowerMultiply
( Orientation orientation, Diagonal diag, int diagOffset,
  const Matrix<F>& L, Matrix<F>& X );
template<typename F>
void DistFrontLowerMultiply
( Orientation orientation, Diagonal diag, int diagOffset,
  const DistMatrix<F,VC,STAR>& L, DistMatrix<F,VC,STAR>& X );

template<typename F>
void LocalFrontLowerMultiplyNormal
( Diagonal diag, int diagOffset, const Matrix<F>& L, Matrix<F>& X );
template<typename F>
void DistFrontLowerMultiplyNormal
( Diagonal diag, int diagOffset, 
  const DistMatrix<F,VC,STAR>& L, DistMatrix<F,VC,STAR>& X );

// Handles the TRANSPOSE and ADJOINT cases
template<typename F>
void LocalFrontLowerMultiplyTranspose
( Orientation orientation, Diagonal diag, int diagOffset, 
  const Matrix<F>& L, Matrix<F>& X );
template<typename F>
void DistFrontLowerMultiplyTranspose
( Orientation orientation, Diagonal diag, int diagOffset,
  const DistMatrix<F,VC,STAR>& L, DistMatrix<F,VC,STAR>& X );

//----------------------------------------------------------------------------//
// Implementation begins here                                                 //
//----------------------------------------------------------------------------//

template<typename F>
void LocalFrontLowerMultiply
( Orientation orientation, Diagonal diag, int diagOffset,
  const Matrix<F>& L, Matrix<F>& X )
{
#ifndef RELEASE
    PushCallStack("numeric::LocalFrontLowerMultiply");
#endif
    if( orientation == NORMAL )
        LocalFrontLowerMultiplyNormal( diag, diagOffset, L, X );
    else
        LocalFrontLowerMultiplyTranspose( orientation, diag, diagOffset, L, X );
#ifndef RELEASE
    PopCallStack();
#endif
}

template<typename F>
void DistFrontLowerMultiply
( Orientation orientation, Diagonal diag, int diagOffset, 
  const DistMatrix<F,VC,STAR>& L, DistMatrix<F,VC,STAR>& X )
{
#ifndef RELEASE
    PushCallStack("numeric::DistFrontLowerMultiply");
#endif
    if( orientation == NORMAL )
        DistFrontLowerMultiplyNormal( diag, diagOffset, L, X );
    else
        DistFrontLowerMultiplyTranspose( orientation, diag, diagOffset, L, X );
#ifndef RELEASE
    PopCallStack();
#endif
}

} // namespace numeric
} // namespace clique

#endif /* CLIQUE_NUMERIC_FRONT_LOWER_MULTIPLY_HPP */

