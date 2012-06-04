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
#ifndef CLIQUE_NUMERIC_LOCAL_FRONT_LOWER_MULTIPLY_HPP
#define CLIQUE_NUMERIC_LOCAL_FRONT_LOWER_MULTIPLY_HPP 1

namespace cliq {
namespace numeric {

template<typename F>
void LocalFrontLowerMultiply
( Orientation orientation, UnitOrNonUnit diag, int diagOffset,
  const Matrix<F>& L, Matrix<F>& X );

template<typename F>
void LocalFrontLowerMultiplyNormal
( UnitOrNonUnit diag, int diagOffset, const Matrix<F>& L, Matrix<F>& X );

template<typename F>
void LocalFrontLowerMultiplyTranspose
( Orientation orientation, UnitOrNonUnit diag, int diagOffset,
  const Matrix<F>& L, Matrix<F>& X );

//----------------------------------------------------------------------------//
// Implementation begins here                                                 //
//----------------------------------------------------------------------------//

namespace internal {
using namespace elem;

template<typename F> 
void ModifyForTrmm
( Matrix<F>& D, UnitOrNonUnit diag, int diagOffset, 
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
( Matrix<F>& D, UnitOrNonUnit diag, int diagOffset, 
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

} // namespace internal

template<typename F>
inline void LocalFrontLowerMultiply
( Orientation orientation, UnitOrNonUnit diag, int diagOffset,
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
inline void LocalFrontLowerMultiplyNormal
( UnitOrNonUnit diag, int diagOffset, const Matrix<F>& L, Matrix<F>& X )
{
#ifndef RELEASE
    PushCallStack("numeric::LocalFrontLowerMultiplyNormal");
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
    elem::PartitionDown
    ( *LMod, LT,
             LB, L.Width() );

    Matrix<F> XT, 
              XB;
    elem::PartitionDown
    ( X, XT,
         XB, L.Width() );

    elem::Gemm( NORMAL, NORMAL, (F)1, LB, XT, (F)1, XB );

    if( diagOffset == 0 )
    {
        elem::Trmm( LEFT, LOWER, NORMAL, diag, (F)1, LT, XT );
    }
    else
    {
        std::vector<Matrix<F> > diagonals;
        internal::ModifyForTrmm( LT, diag, diagOffset, diagonals );
        elem::Trmm( LEFT, LOWER, NORMAL, NON_UNIT, (F)1, LT, XT );
        internal::ReplaceAfterTrmm( LT, diag, diagOffset, diagonals );
    }
#ifndef RELEASE
    PopCallStack();
#endif
}

template<typename F>
inline void LocalFrontLowerMultiplyTranspose
( Orientation orientation, UnitOrNonUnit diag, int diagOffset,
  const Matrix<F>& L, Matrix<F>& X )
{
#ifndef RELEASE
    PushCallStack("numeric::LocalFrontLowerMultiplyTranspose");
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
    elem::PartitionDown
    ( *LMod, LT,
             LB, L.Width() );

    Matrix<F> XT, 
              XB;
    elem::PartitionDown
    ( X, XT,
         XB, L.Width() );

    if( diagOffset == 0 )
    {
        elem::Trmm( LEFT, LOWER, orientation, diag, (F)1, LT, XT );
    }
    else
    {
        std::vector<Matrix<F> > diagonals;
        internal::ModifyForTrmm( LT, diag, diagOffset, diagonals );
        elem::Trmm( LEFT, LOWER, orientation, NON_UNIT, (F)1, LT, XT );
        internal::ReplaceAfterTrmm( LT, diag, diagOffset, diagonals );
    }

    elem::Gemm( orientation, NORMAL, (F)1, LB, XB, (F)1, XT );
#ifndef RELEASE
    PopCallStack();
#endif
}

} // namespace numeric
} // namespace cliq

#endif // CLIQUE_NUMERIC_LOCAL_FRONT_LOWER_MULTIPLY_HPP
