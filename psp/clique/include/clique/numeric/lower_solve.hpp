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
#ifndef CLIQUE_NUMERIC_LOWER_SOLVE_HPP
#define CLIQUE_NUMERIC_LOWER_SOLVE_HPP 1

namespace clique {
namespace numeric {

template<typename F>
void LowerSolve
( Orientation orientation, Diagonal diag,
  const symbolic::SymmFact& S,
  const numeric::SymmFrontTree<F>& L,
        Matrix<F>& localX );

// Helpers

template<typename F>
void LocalLowerForwardSolve
( Diagonal diag,
  const symbolic::SymmFact& S, 
  const numeric::SymmFrontTree<F>& L, 
        Matrix<F>& localX );

template<typename F>
void DistLowerForwardSolve
( Diagonal diag,
  const symbolic::SymmFact& S,
  const numeric::SymmFrontTree<F>& L,
        Matrix<F>& localX );

template<typename F>
void LocalLowerBackwardSolve
( Orientation orientation, Diagonal diag,
  const symbolic::SymmFact& S, 
  const numeric::SymmFrontTree<F>& L, 
        Matrix<F>& localX );

template<typename F>
void DistLowerBackwardSolve
( Orientation orientation, Diagonal diag,
  const symbolic::SymmFact& S,
  const numeric::SymmFrontTree<F>& L,
        Matrix<F>& localX );

//----------------------------------------------------------------------------//
// Implementation begins here                                                 //
//----------------------------------------------------------------------------//

template<typename F>
void LowerSolve
( Orientation orientation, Diagonal diag, 
  const symbolic::SymmFact& S,
  const numeric::SymmFrontTree<F>& L,
        Matrix<F>& localX )
{
#ifndef RELEASE
    PushCallStack("numeric::LowerSolve");
#endif
    if( orientation == NORMAL )
    {
        LocalLowerForwardSolve( diag, S, L, localX );
        DistLowerForwardSolve( diag, S, L, localX );
    }
    else
    {
        DistLowerBackwardSolve( orientation, diag, S, L, localX );
        LocalLowerBackwardSolve( orientation, diag, S, L, localX );
    }
#ifndef RELEASE
    PopCallStack();
#endif
}

} // namespace numeric
} // namespace clique

#endif /* CLIQUE_NUMERIC_LOWER_SOLVE_HPP */

