/*
   Modification of include/elemental/advanced/LDL.hpp from Elemental.
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

template<typename F> // F represents a real or complex field
void clique::numeric::LocalFrontLDL
( Orientation orientation, Matrix<F>& AL, Matrix<F>& ABR )
{
#ifndef RELEASE
    clique::PushCallStack("numeric::LocalFrontLDL");
    if( ABR.Height() != ABR.Width() )
        throw std::logic_error("ABR must be square");
    if( AL.Height() != AL.Width() + ABR.Width() )
        throw std::logic_error("AL and ABR don't have conformal dimensions");
    if( orientation == NORMAL )
        throw std::logic_error("LocalFrontLDL must be (conjugate-)transposed.");
#endif
    Matrix<F>
        ALTL, ALTR,  AL00, AL01, AL02,
        ALBL, ALBR,  AL10, AL11, AL12,
                     AL20, AL21, AL22;
    Matrix<F> d1;
    Matrix<F> S21;

    Matrix<F> S21T,
              S21B;
    Matrix<F> AL21T,
              AL21B;

    // Start the algorithm
    elemental::PartitionDownDiagonal
    ( AL, ALTL, ALTR,
          ALBL, ALBR, 0 );
    while( ALTL.Width() < AL.Width() )
    {
        elemental::RepartitionDownDiagonal
        ( ALTL, /**/ ALTR,  AL00, /**/ AL01, AL02,
         /***************/ /*********************/
                /**/        AL10, /**/ AL11, AL12,
          ALBL, /**/ ALBR,  AL20, /**/ AL21, AL22 );

        //--------------------------------------------------------------------//
        // This routine is unblocked, hence the need for us to generalize to 
        // an (ideally) faster blocked algorithm.
        elemental::internal::LDLVar3( orientation, AL11, d1 );

        elemental::Trsm( RIGHT, LOWER, orientation, UNIT, (F)1, AL11, AL21 );

        S21 = AL21;
        elemental::DiagonalSolve( RIGHT, NORMAL, d1, AL21 );

        elemental::PartitionDown
        ( S21, S21T,
               S21B, AL22.Width() );
        elemental::PartitionDown
        ( AL21, AL21T,
                AL21B, AL22.Width() );
        elemental::Gemm( NORMAL, orientation, (F)-1, S21, AL21T, (F)1, AL22 );
        elemental::internal::TrrkNT
        ( LOWER, orientation, (F)-1, S21B, AL21B, (F)1, ABR );
        //--------------------------------------------------------------------//

        elemental::SlidePartitionDownDiagonal
        ( ALTL, /**/ ALTR,  AL00, AL01, /**/ AL02,
                /**/        AL10, AL11, /**/ AL12,
         /***************/ /*********************/
          ALBL, /**/ ALBR,  AL20, AL21, /**/ AL22 );
    }
#ifndef RELEASE
    clique::PopCallStack();
#endif
}

template void clique::numeric::LocalFrontLDL
( Orientation orientation, 
  Matrix<float>& AL, Matrix<float>& ABR);

template void clique::numeric::LocalFrontLDL
( Orientation orientation, 
  Matrix<double>& AL, Matrix<double>& ABR );

template void clique::numeric::LocalFrontLDL
( Orientation orientation, 
  Matrix<std::complex<float> >& AL, Matrix<std::complex<float> >& ABR );

template void clique::numeric::LocalFrontLDL
( Orientation orientation, 
  Matrix<std::complex<double> >& AL, Matrix<std::complex<double> >& ABR );
