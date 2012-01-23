/*
   Modification of include/elemental/basic/level3/Trsm/TrsmLLN.hpp 
   from Elemental.
   Copyright (c) 2009-2011, Jack Poulson
   All rights reserved.

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

template<typename F>
void clique::numeric::LocalFrontDiagonalSolve
( const Matrix<F>& d, Matrix<F>& X )
{
#ifndef RELEASE
    clique::PushCallStack("numeric::LocalFrontDiagonalSolve");
    if( d.Width() != 1 )
        throw std::logic_error("d must be a column vector");
    if( d.Height() != X.Height() )
        throw std::logic_error("Invalid height of X");
#endif
    elemental::DiagonalSolve( LEFT, NORMAL, d, X, true );
#ifndef RELEASE
    clique::PopCallStack();
#endif
}

template void clique::numeric::LocalFrontDiagonalSolve
( const Matrix<float>& d, Matrix<float>& X );

template void clique::numeric::LocalFrontDiagonalSolve
( const Matrix<double>& d, Matrix<double>& X );

template void clique::numeric::LocalFrontDiagonalSolve
( const Matrix<std::complex<float> >& d, Matrix<std::complex<float> >& X );

template void clique::numeric::LocalFrontDiagonalSolve
( const Matrix<std::complex<double> >& d, Matrix<std::complex<double> >& X );
