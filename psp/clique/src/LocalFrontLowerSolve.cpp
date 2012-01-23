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
void clique::numeric::LocalFrontLowerForwardSolve
( Diagonal diag, const Matrix<F>& L, Matrix<F>& X )
{
#ifndef RELEASE
    clique::PushCallStack("numeric::LocalFrontLowerForwardSolve");
    if( L.Height() < L.Width() || L.Height() != X.Height() )
    {
        std::ostringstream msg;
        msg << "Nonconformal solve:\n"
            << "  L ~ " << L.Height() << " x " << L.Width() << "\n"
            << "  X ~ " << X.Height() << " x " << X.Width() << "\n";
        throw std::logic_error( msg.str().c_str() );
    }
#endif
    Matrix<F> LT,
              LB;
    elemental::LockedPartitionDown
    ( L, LT,
         LB, L.Width() );

    Matrix<F> XT, 
              XB;
    elemental::PartitionDown
    ( X, XT,
         XB, L.Width() );

    elemental::Trsm( LEFT, LOWER, NORMAL, diag, (F)1, LT, XT, true );
    elemental::Gemm( NORMAL, NORMAL, (F)-1, LB, XT, (F)1, XB );
#ifndef RELEASE
    clique::PopCallStack();
#endif
}

template<typename F>
void clique::numeric::LocalFrontLowerBackwardSolve
( Orientation orientation, Diagonal diag, 
  const Matrix<F>& L, Matrix<F>& X )
{
#ifndef RELEASE
    clique::PushCallStack("numeric::LocalFrontLowerBackwardSolve");
    if( L.Height() < L.Width() || L.Height() != X.Height() )
    {
        std::ostringstream msg;
        msg << "Nonconformal solve:\n"
            << "  L ~ " << L.Height() << " x " << L.Width() << "\n"
            << "  X ~ " << X.Height() << " x " << X.Width() << "\n";
        throw std::logic_error( msg.str().c_str() );
    }
    if( orientation == NORMAL )
        throw std::logic_error("This solve must be (conjugate-)transposed");
#endif
    Matrix<F> LT,
              LB;
    elemental::LockedPartitionDown
    ( L, LT,
         LB, L.Width() );

    Matrix<F> XT,
              XB;
    elemental::PartitionDown
    ( X, XT,
         XB, L.Width() );

    elemental::Gemm( orientation, NORMAL, (F)-1, LB, XB, (F)1, XT );
    elemental::Trsm( LEFT, LOWER, orientation, diag, (F)1, LT, XT, true );
#ifndef RELEASE
    clique::PopCallStack();
#endif
}

template void clique::numeric::LocalFrontLowerForwardSolve
( Diagonal diag, 
  const Matrix<float>& L, Matrix<float>& X );
template void clique::numeric::LocalFrontLowerBackwardSolve
( Orientation orientation, Diagonal diag,
  const Matrix<float>& L, Matrix<float>& X );

template void clique::numeric::LocalFrontLowerForwardSolve
( Diagonal diag, 
  const Matrix<double>& L, Matrix<double>& X );
template void clique::numeric::LocalFrontLowerBackwardSolve
( Orientation orientation, Diagonal diag,
  const Matrix<double>& L, Matrix<double>& X );

template void clique::numeric::LocalFrontLowerForwardSolve
( Diagonal diag,
  const Matrix<std::complex<float> >& L, Matrix<std::complex<float> >& X );
template void clique::numeric::LocalFrontLowerBackwardSolve
( Orientation orientation, Diagonal diag,
  const Matrix<std::complex<float> >& L, Matrix<std::complex<float> >& X );

template void clique::numeric::LocalFrontLowerForwardSolve
( Diagonal diag,
  const Matrix<std::complex<double> >& L, Matrix<std::complex<double> >& X );
template void clique::numeric::LocalFrontLowerBackwardSolve
( Orientation orientation, Diagonal diag,
  const Matrix<std::complex<double> >& L, Matrix<std::complex<double> >& X );
