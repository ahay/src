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

namespace {
using namespace elemental;
template<typename F> // represents a real or complex ring
void ModifyForTrmm
( Matrix<F>& D, Diagonal diag, int diagOffset, 
  std::vector<Matrix<F> >& diagonals )
{
#ifndef RELEASE
    PushCallStack("ModifyForTrmm");
#endif
    if( diag == UNIT )
    {
        diagonals.resize( 1-diagOffset );
        for( int i=0; i<-diagOffset; ++i )
            diagonals[i].ResizeTo( D.DiagonalLength(-i), 1 );
        diagonals[-diagOffset].ResizeTo( D.DiagonalLength(-diagOffset), 1 );
    }
    else
    {
        diagonals.resize( -diagOffset );
        for( int i=0; i<-diagOffset; ++i )
            diagonals[i].ResizeTo( D.DiagonalLength(-i), 1 );
    }

    const int height = D.Height();
    for( int j=0; j<height; ++j )
    {
        const int length = std::min(-diagOffset,height-j);
        for( int i=0; i<length; ++i )    
        {
            diagonals[i].Set( j, 0, D.Get(j+i,j) );
            D.Set( j+i, j, (F)0 );
        }
        if( diag == UNIT && j-diagOffset < height )
        {
            diagonals[-diagOffset].Set( j, 0, D.Get(j-diagOffset,j) );
            D.Set( j-diagOffset, j, (F)1 );
        }
    }
#ifndef RELEASE
    PopCallStack();
#endif
}

template<typename F> // represents a real or complex ring
void ReplaceAfterTrmm
( Matrix<F>& D, Diagonal diag, int diagOffset, 
  const std::vector<Matrix<F> >& diagonals )
{
#ifndef RELEASE
    PushCallStack("ReplaceAfterTrmm");
#endif
    const int height = D.Height();
    for( int j=0; j<height; ++j )
    {
        const int length = std::min(-diagOffset,height-j);
        for( int i=0; i<length; ++i )    
            D.Set( j+i, j, diagonals[i].Get(j,0) );
        if( diag == UNIT && j-diagOffset < height )
            D.Set( j-diagOffset, j, diagonals[-diagOffset].Get(j,0) );
    }
#ifndef RELEASE
    PopCallStack();
#endif
}
} // anonymous namespace

template<typename F>
void clique::numeric::LocalFrontLowerMultiplyNormal
( Diagonal diag, int diagOffset, const Matrix<F>& L, Matrix<F>& X )
{
#ifndef RELEASE
    clique::PushCallStack("numeric::LocalFrontLowerMultiplyNormal");
    if( L.Height() < L.Width() || L.Height() != X.Height() )
    {
        std::ostringstream msg;
        msg << "Nonconformal multiply:\n"
            << "  L ~ " << L.Height() << " x " << L.Width() << "\n"
            << "  X ~ " << X.Height() << " x " << X.Width() << "\n";
        throw std::logic_error( msg.str().c_str() );
    }
    if( diagOffset > 0 )
        throw std::logic_error("Diagonal offsets cannot be positive");
#endif
    Matrix<F>* LMod = const_cast<Matrix<F>*>(&L);
    Matrix<F> LT,
              LB;
    elemental::PartitionDown
    ( *LMod, LT,
             LB, L.Width() );

    Matrix<F> XT, 
              XB;
    elemental::PartitionDown
    ( X, XT,
         XB, L.Width() );

    elemental::Gemm( NORMAL, NORMAL, (F)1, LB, XT, (F)1, XB );

    if( diagOffset == 0 )
    {
        elemental::Trmm( LEFT, LOWER, NORMAL, diag, (F)1, LT, XT );
    }
    else
    {
        std::vector<Matrix<F> > diagonals;
        ModifyForTrmm( LT, diag, diagOffset, diagonals );
        elemental::Trmm( LEFT, LOWER, NORMAL, NON_UNIT, (F)1, LT, XT );
        ReplaceAfterTrmm( LT, diag, diagOffset, diagonals );
    }
#ifndef RELEASE
    clique::PopCallStack();
#endif
}

template<typename F>
void clique::numeric::LocalFrontLowerMultiplyTranspose
( Orientation orientation, Diagonal diag, int diagOffset,
  const Matrix<F>& L, Matrix<F>& X )
{
#ifndef RELEASE
    clique::PushCallStack("numeric::LocalFrontLowerMultiplyTranspose");
    if( L.Height() < L.Width() || L.Height() != X.Height() )
    {
        std::ostringstream msg;
        msg << "Nonconformal solve:\n"
            << "  L ~ " << L.Height() << " x " << L.Width() << "\n"
            << "  X ~ " << X.Height() << " x " << X.Width() << "\n";
        throw std::logic_error( msg.str().c_str() );
    }
    if( orientation == NORMAL )
        throw std::logic_error("Orientation must be (conjugate-)transposed");
    if( diagOffset > 0 )
        throw std::logic_error("Diagonal offsets cannot be positive");
#endif
    Matrix<F>* LMod = const_cast<Matrix<F>*>(&L);
    Matrix<F> LT,
              LB;
    elemental::PartitionDown
    ( *LMod, LT,
             LB, L.Width() );

    Matrix<F> XT, 
              XB;
    elemental::PartitionDown
    ( X, XT,
         XB, L.Width() );

    if( diagOffset == 0 )
    {
        elemental::Trmm( LEFT, LOWER, orientation, diag, (F)1, LT, XT );
    }
    else
    {
        std::vector<Matrix<F> > diagonals;
        ModifyForTrmm( LT, diag, diagOffset, diagonals );
        elemental::Trmm( LEFT, LOWER, orientation, NON_UNIT, (F)1, LT, XT );
        ReplaceAfterTrmm( LT, diag, diagOffset, diagonals );
    }

    elemental::Gemm( orientation, NORMAL, (F)1, LB, XB, (F)1, XT );
#ifndef RELEASE
    clique::PopCallStack();
#endif
}

template void clique::numeric::LocalFrontLowerMultiplyNormal
( Diagonal diag, int diagOffset, 
  const Matrix<float>& L, Matrix<float>& X );
template void clique::numeric::LocalFrontLowerMultiplyTranspose
( Orientation orientation, Diagonal diag, int diagOffset, 
  const Matrix<float>& L, Matrix<float>& X );

template void clique::numeric::LocalFrontLowerMultiplyNormal
( Diagonal diag, int diagOffset,
  const Matrix<double>& L, Matrix<double>& X );
template void clique::numeric::LocalFrontLowerMultiplyTranspose
( Orientation orientation, Diagonal diag, int diagOffset, 
  const Matrix<double>& L, Matrix<double>& X );

template void clique::numeric::LocalFrontLowerMultiplyNormal
( Diagonal diag, int diagOffset, 
  const Matrix<std::complex<float> >& L, Matrix<std::complex<float> >& X );
template void clique::numeric::LocalFrontLowerMultiplyTranspose
( Orientation orientation, Diagonal diag, int diagOffset, 
  const Matrix<std::complex<float> >& L, Matrix<std::complex<float> >& X );

template void clique::numeric::LocalFrontLowerMultiplyNormal
( Diagonal diag, int diagOffset, 
  const Matrix<std::complex<double> >& L, Matrix<std::complex<double> >& X );
template void clique::numeric::LocalFrontLowerMultiplyTranspose
( Orientation orientation, Diagonal diag, int diagOffset, 
  const Matrix<std::complex<double> >& L, Matrix<std::complex<double> >& X );
