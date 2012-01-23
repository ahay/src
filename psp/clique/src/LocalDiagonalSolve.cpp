/*
   Clique: a scalable implementation of the multifrontal algorithm

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
#include "clique.hpp"

template<typename F> // F represents a real or complex field
void clique::numeric::LocalDiagonalSolve
( const symbolic::SymmFact& S,
  const numeric::SymmFrontTree<F>& L,
        Matrix<F>& X )
{
    using namespace clique::symbolic;
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
        LocalFrontDiagonalSolve( d, XSub );
    }
#ifndef RELEASE
    PopCallStack();
#endif
}

template void clique::numeric::LocalDiagonalSolve
( const symbolic::SymmFact& S,
  const numeric::SymmFrontTree<float>& L,
        Matrix<float>& X );

template void clique::numeric::LocalDiagonalSolve
( const symbolic::SymmFact& S,
  const numeric::SymmFrontTree<double>& L,
        Matrix<double>& X );

template void clique::numeric::LocalDiagonalSolve
( const symbolic::SymmFact& S,
  const numeric::SymmFrontTree<std::complex<float> >& L,
        Matrix<std::complex<float> >& X );

template void clique::numeric::LocalDiagonalSolve
( const symbolic::SymmFact& S,
  const numeric::SymmFrontTree<std::complex<double> >& L,
        Matrix<std::complex<double> >& X );
