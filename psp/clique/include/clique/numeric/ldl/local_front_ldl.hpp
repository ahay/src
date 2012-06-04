/*
   Modification of include/elemental/advanced/LDL.hpp from Elemental.
   Copyright (c) 2009-2012, Jack Poulson
   All rights reserved.

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
#ifndef CLIQUE_NUMERIC_LOCAL_FRONT_LDL_HPP
#define CLIQUE_NUMERIC_LOCAL_FRONT_LDL_HPP 1

namespace cliq {
namespace numeric {

template<typename F> 
void LocalFrontLDL( Orientation orientation, Matrix<F>& AL, Matrix<F>& ABR );

//----------------------------------------------------------------------------//
// Implementation begins here                                                 //
//----------------------------------------------------------------------------//

template<typename F> 
inline void LocalFrontLDL
( Orientation orientation, Matrix<F>& AL, Matrix<F>& ABR )
{
#ifndef RELEASE
    PushCallStack("numeric::LocalFrontLDL");
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
    elem::PartitionDownDiagonal
    ( AL, ALTL, ALTR,
          ALBL, ALBR, 0 );
    while( ALTL.Width() < AL.Width() )
    {
        elem::RepartitionDownDiagonal
        ( ALTL, /**/ ALTR,  AL00, /**/ AL01, AL02,
         /***************/ /*********************/
                /**/        AL10, /**/ AL11, AL12,
          ALBL, /**/ ALBR,  AL20, /**/ AL21, AL22 );

        //--------------------------------------------------------------------//
        // This routine is unblocked, hence the need for us to generalize to 
        // an (ideally) faster blocked algorithm.
        elem::internal::LDLVar3( orientation, AL11, d1 );

        elem::Trsm( RIGHT, LOWER, orientation, UNIT, (F)1, AL11, AL21 );

        S21 = AL21;
        elem::DiagonalSolve( RIGHT, NORMAL, d1, AL21 );

        elem::PartitionDown
        ( S21, S21T,
               S21B, AL22.Width() );
        elem::PartitionDown
        ( AL21, AL21T,
                AL21B, AL22.Width() );
        elem::Gemm( NORMAL, orientation, (F)-1, S21, AL21T, (F)1, AL22 );
        elem::internal::TrrkNT
        ( LOWER, orientation, (F)-1, S21B, AL21B, (F)1, ABR );
        //--------------------------------------------------------------------//

        elem::SlidePartitionDownDiagonal
        ( ALTL, /**/ ALTR,  AL00, AL01, /**/ AL02,
                /**/        AL10, AL11, /**/ AL12,
         /***************/ /*********************/
          ALBL, /**/ ALBR,  AL20, AL21, /**/ AL22 );
    }
#ifndef RELEASE
    PopCallStack();
#endif
}

} // namespace numeric
} // namespace cliq

#endif // CLIQUE_NUMERIC_LOCAL_FRONT_LDL_HPP
