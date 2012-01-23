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
void ModifyForTrmm( DistMatrix<F,STAR,STAR>& D, Diagonal diag, int diagOffset )
{
#ifndef RELEASE
    PushCallStack("ModifyForTrmm");
#endif
    const int height = D.Height();
    for( int j=0; j<height; ++j )
    {
        const int length = std::min(-diagOffset,height-j);
        std::memset( D.LocalBuffer(j,j), 0, length*sizeof(F) );
        if( diag == UNIT && j-diagOffset < height )
            D.SetLocalEntry( j-diagOffset, j, (F)1 );
    }
#ifndef RELEASE
    PopCallStack();
#endif
}
} // anonymous namespace

template<typename F>
void clique::numeric::DistFrontLowerMultiplyNormal
( Diagonal diag, int diagOffset,
  const DistMatrix<F,VC,STAR>& L, DistMatrix<F,VC,STAR>& X )
{
#ifndef RELEASE
    clique::PushCallStack("numeric::DistFrontLowerMultiplyNormal");
    if( L.Grid() != X.Grid() )
        throw std::logic_error
        ("L and X must be distributed over the same grid");
    if( L.Height() < L.Width() || L.Height() != X.Height() )
    {
        std::ostringstream msg;
        msg << "Nonconformal multiply:\n"
            << "  L ~ " << L.Height() << " x " << L.Width() << "\n"
            << "  X ~ " << X.Height() << " x " << X.Width() << "\n";
        throw std::logic_error( msg.str().c_str() );
    }
    if( L.ColAlignment() != X.ColAlignment() )
        throw std::logic_error("L and X are assumed to be aligned");
    if( diagOffset > 0 )
        throw std::logic_error("Diagonal offsets cannot be positive");
#endif
    const Grid& g = L.Grid();

    // Matrix views
    DistMatrix<F,VC,STAR>
        LTL(g), LTR(g),  L00(g), L01(g), L02(g),
        LBL(g), LBR(g),  L10(g), L11(g), L12(g),
                         L20(g), L21(g), L22(g);

    DistMatrix<F,VC,STAR> XT(g),  X0(g),
                          XB(g),  X1(g),
                                  X2(g);

    // Temporary distributions
    DistMatrix<F,STAR,STAR> L11_STAR_STAR(g);
    DistMatrix<F,STAR,STAR> X1_STAR_STAR(g);

    // Start the algorithm
    elemental::LockedPartitionDownDiagonal
    ( L, LTL, LTR,
         LBL, LBR, L.Width() );
    elemental::PartitionDown
    ( X, XT,
         XB, L.Width() );
    while( XT.Height() > 0 )
    {
        elemental::LockedRepartitionUpDiagonal
        ( LTL, /**/ LTR,  L00, L01, /**/ L02,
               /**/       L10, L11, /**/ L12,
         /*************/ /******************/
          LBL, /**/ LBR,  L20, L21, /**/ L22 );

        elemental::RepartitionUp
        ( XT,  X0,
               X1,
         /**/ /**/
          XB,  X2 );

        //--------------------------------------------------------------------//
        X1_STAR_STAR = X1;
        elemental::internal::LocalGemm
        ( NORMAL, NORMAL, (F)1, L21, X1_STAR_STAR, (F)1, X2 );

        if( diagOffset == 0 )
        {
            L11_STAR_STAR = L11;
            elemental::internal::LocalTrmm
            ( LEFT, LOWER, NORMAL, diag, (F)1, L11_STAR_STAR, X1_STAR_STAR );
        }
        else
        {
            L11_STAR_STAR = L11;
            ModifyForTrmm( L11_STAR_STAR, diag, diagOffset );
            elemental::internal::LocalTrmm
            ( LEFT, LOWER, NORMAL, NON_UNIT, 
              (F)1, L11_STAR_STAR, X1_STAR_STAR );
        }
        X1 = X1_STAR_STAR;
        //--------------------------------------------------------------------//

        elemental::SlideLockedPartitionUpDiagonal
        ( LTL, /**/ LTR,  L00, /**/ L01, L02,
         /*************/ /******************/
               /**/       L10, /**/ L11, L12,
          LBL, /**/ LBR,  L20, /**/ L21, L22 );

        elemental::SlidePartitionUp
        ( XT,  X0,
         /**/ /**/
               X1,
          XB,  X2 );
    }
#ifndef RELEASE
    clique::PopCallStack();
#endif
}

template<typename F>
void clique::numeric::DistFrontLowerMultiplyTranspose
( Orientation orientation, Diagonal diag, int diagOffset,
  const DistMatrix<F,VC,STAR>& L, DistMatrix<F,VC,STAR>& X )
{
#ifndef RELEASE
    clique::PushCallStack("numeric::DistFrontLowerMultiplyTranspose");
    if( L.Grid() != X.Grid() )
        throw std::logic_error
        ("L and X must be distributed over the same grid");
    if( L.Height() < L.Width() || L.Height() != X.Height() )
    {
        std::ostringstream msg;
        msg << "Nonconformal multiply:\n"
            << "  L ~ " << L.Height() << " x " << L.Width() << "\n"
            << "  X ~ " << X.Height() << " x " << X.Width() << "\n";
        throw std::logic_error( msg.str().c_str() );
    }
    if( L.ColAlignment() != X.ColAlignment() )
        throw std::logic_error("L and X are assumed to be aligned");
    if( orientation == NORMAL )
        throw std::logic_error("Orientation must be (conjugate-)transposed");
    if( diagOffset > 0 )
        throw std::logic_error("Diagonal offsets cannot be positive");
#endif
    const Grid& g = L.Grid();

    // Matrix views
    DistMatrix<F,VC,STAR>
        LTL(g), LTR(g),  L00(g), L01(g), L02(g),
        LBL(g), LBR(g),  L10(g), L11(g), L12(g),
                         L20(g), L21(g), L22(g);

    DistMatrix<F,VC,STAR> XT(g),  X0(g),
                          XB(g),  X1(g),
                                  X2(g);

    // Temporary distributions
    DistMatrix<F,STAR,STAR> X1_STAR_STAR(g);
    DistMatrix<F,STAR,STAR> L11_STAR_STAR(g);
    DistMatrix<F,STAR,STAR> Z1_STAR_STAR(g);

    LockedPartitionDownDiagonal
    ( L, LTL, LTR,
         LBL, LBR, 0 );
    PartitionDown
    ( X, XT,
         XB, 0 );
    while( LTL.Width() < L.Width() )
    {
        LockedRepartitionDownDiagonal
        ( LTL, /**/ LTR,   L00, /**/ L01, L02,
         /*************/  /******************/
               /**/        L10, /**/ L11, L12,
          LBL, /**/ LBR,   L20, /**/ L21, L22 );

        RepartitionDown
        ( XT,  X0,
         /**/ /**/
               X1,
          XB,  X2, L11.Height() );

        //--------------------------------------------------------------------//
        X1_STAR_STAR = X1; // Can this be avoided?
        L11_STAR_STAR = L11;
        if( diagOffset == 0 )
        {
            elemental::internal::LocalTrmm
            ( LEFT, LOWER, orientation, diag, 
              (F)1, L11_STAR_STAR, X1_STAR_STAR );
        }
        else
        {
            ModifyForTrmm( L11_STAR_STAR, diag, diagOffset );
            elemental::internal::LocalTrmm
            ( LEFT, LOWER, orientation, NON_UNIT, 
              (F)1, L11_STAR_STAR, X1_STAR_STAR );
        }
        X1 = X1_STAR_STAR;

        Z1_STAR_STAR.ResizeTo( X1.Height(), X1.Width() );
        elemental::internal::LocalGemm
        ( orientation, NORMAL, (F)1, L21, X2, (F)0, Z1_STAR_STAR );
        X1.SumScatterUpdate( (F)1, Z1_STAR_STAR );
        //--------------------------------------------------------------------//

        SlideLockedPartitionDownDiagonal
        ( LTL, /**/ LTR,  L00, L01, /**/ L02,
               /**/       L10, L11, /**/ L12,
         /*************/ /******************/
          LBL, /**/ LBR,  L20, L21, /**/ L22 );

        SlidePartitionDown
        ( XT,  X0,
               X1,
         /**/ /**/
          XB,  X2 );
    }
#ifndef RELEASE
    clique::PopCallStack();
#endif
}

template void clique::numeric::DistFrontLowerMultiplyNormal
( Diagonal diag, int diagOffset, 
  const DistMatrix<float,VC,STAR>& L,
        DistMatrix<float,VC,STAR>& X );
template void clique::numeric::DistFrontLowerMultiplyTranspose
( Orientation orientation, Diagonal diag, int diagOffset,
  const DistMatrix<float,VC,STAR>& L,
        DistMatrix<float,VC,STAR>& X );

template void clique::numeric::DistFrontLowerMultiplyNormal
( Diagonal diag, int diagOffset, 
  const DistMatrix<double,VC,STAR>& L, 
        DistMatrix<double,VC,STAR>& X );
template void clique::numeric::DistFrontLowerMultiplyTranspose
( Orientation orientation, Diagonal diag, int diagOffset,
  const DistMatrix<double,VC,STAR>& L,
        DistMatrix<double,VC,STAR>& X );

template void clique::numeric::DistFrontLowerMultiplyNormal
( Diagonal diag, int diagOffset, 
  const DistMatrix<std::complex<float>,VC,STAR>& L, 
        DistMatrix<std::complex<float>,VC,STAR>& X );
template void clique::numeric::DistFrontLowerMultiplyTranspose
( Orientation orientation, Diagonal diag, int diagOffset, 
  const DistMatrix<std::complex<float>,VC,STAR>& L, 
        DistMatrix<std::complex<float>,VC,STAR>& X );

template void clique::numeric::DistFrontLowerMultiplyNormal
( Diagonal diag, int diagOffset, 
  const DistMatrix<std::complex<double>,VC,STAR>& L, 
        DistMatrix<std::complex<double>,VC,STAR>& X );
template void clique::numeric::DistFrontLowerMultiplyTranspose
( Orientation orientation, Diagonal diag, int diagOffset, 
  const DistMatrix<std::complex<double>,VC,STAR>& L,
        DistMatrix<std::complex<double>,VC,STAR>& X );
