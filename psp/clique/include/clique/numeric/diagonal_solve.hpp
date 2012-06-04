/*
   Clique: a scalable implementation of the multifrontal algorithm

   Copyright (C) 2011-2012 Jack Poulson, Lexing Ying, and 
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
#ifndef CLIQUE_NUMERIC_DIAGONAL_SOLVE_HPP
#define CLIQUE_NUMERIC_DIAGONAL_SOLVE_HPP 1

namespace cliq {
namespace numeric {

template<typename F>
void DiagonalSolve
( const symbolic::SymmFact& S,
  const numeric::SymmFrontTree<F>& L,
        Matrix<F>& localX );

// Helpers

template<typename F>
void LocalDiagonalSolve
( const symbolic::SymmFact& S,
  const numeric::SymmFrontTree<F>& L, 
        Matrix<F>& localX );

template<typename F>
void DistDiagonalSolve
( const symbolic::SymmFact& S,
  const numeric::SymmFrontTree<F>& L,
        Matrix<F>& localX );

//----------------------------------------------------------------------------//
// Implementation begins here                                                 //
//----------------------------------------------------------------------------//

template<typename F>
inline void LocalDiagonalSolve
( const symbolic::SymmFact& S,
  const numeric::SymmFrontTree<F>& L,
        Matrix<F>& X )
{
    using namespace symbolic;
#ifndef RELEASE
    PushCallStack("numeric::LocalDiagonalSolve");
#endif
    const int numLocalSupernodes = S.local.supernodes.size();
    const int width = X.Width();
    Matrix<F> XSub;
    for( int s=0; s<numLocalSupernodes; ++s )
    {
        const LocalSymmFactSupernode& sn = S.local.supernodes[s];
        const Matrix<F>& frontL = L.local.fronts[s].frontL;
        XSub.View( X, sn.myOffset, 0, sn.size, width );

        Matrix<F> frontTL;
        frontTL.LockedView( frontL, 0, 0, sn.size, sn.size );
        Matrix<F> d;
        frontTL.GetDiagonal( d );
        elem::DiagonalSolve( LEFT, NORMAL, d, XSub, true );
    }
#ifndef RELEASE
    PopCallStack();
#endif
}

template<typename F> 
void DistDiagonalSolve
( const symbolic::SymmFact& S,
  const numeric::SymmFrontTree<F>& L,
        Matrix<F>& localX )
{
    using namespace symbolic;
#ifndef RELEASE
    PushCallStack("numeric::DistDiagonalSolve");
#endif
    const int numDistSupernodes = S.dist.supernodes.size();
    const int width = localX.Width();

    for( int s=1; s<numDistSupernodes; ++s )
    {
        const DistSymmFactSupernode& sn = S.dist.supernodes[s];
        const DistSymmFront<F>& front = L.dist.fronts[s];

        Matrix<F> localXT;
        localXT.View( localX, sn.localOffset1d, 0, sn.localSize1d, width );

        elem::DiagonalSolve
        ( LEFT, NORMAL, front.diag.LockedLocalMatrix(), localXT, true );
    }
#ifndef RELEASE
    PopCallStack();
#endif
}

template<typename F>
inline void DiagonalSolve
( const symbolic::SymmFact& S,
  const numeric::SymmFrontTree<F>& L,
        Matrix<F>& localX )
{
#ifndef RELEASE
    PushCallStack("numeric::DiagonalSolve");
#endif
    cliq::numeric::LocalDiagonalSolve( S, L, localX );
    cliq::numeric::DistDiagonalSolve( S, L, localX );
#ifndef RELEASE
    PopCallStack();
#endif
}

} // namespace numeric
} // namespace cliq

#endif // CLIQUE_NUMERIC_DIAGONAL_SOLVE_HPP 
