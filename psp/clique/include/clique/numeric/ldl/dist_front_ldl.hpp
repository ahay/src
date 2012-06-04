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
#ifndef CLIQUE_NUMERIC_DIST_FRONT_LDL_HPP
#define CLIQUE_NUMERIC_DIST_FRONT_LDL_HPP 1

namespace cliq {
namespace numeric {

template<typename F> 
void DistFrontLDL
( Orientation orientation, DistMatrix<F>& AL, DistMatrix<F>& ABR );

//----------------------------------------------------------------------------//
// Implementation begins here                                                 //
//----------------------------------------------------------------------------//

namespace internal {

template<typename F> 
inline void DistFrontLDLGeneral
( Orientation orientation, DistMatrix<F>& AL, DistMatrix<F>& ABR )
{
#ifndef RELEASE
    PushCallStack("numeric::internal::DistFrontLDLGeneral");
    if( ABR.Height() != ABR.Width() )
        throw std::logic_error("ABR must be square");
    if( AL.Height() != AL.Width()+ABR.Height() )
        throw std::logic_error("AL and ABR must have compatible dimensions");
    if( AL.Grid() != ABR.Grid() )
        throw std::logic_error("AL and ABR must use the same grid");
    if( ABR.ColAlignment() !=
        (AL.ColAlignment()+AL.Width()) % AL.Grid().Height() )
        throw std::logic_error
        ("AL and ABR must have compatible col alignments");
    if( ABR.RowAlignment() != 
        (AL.RowAlignment()+AL.Width()) % AL.Grid().Width() )
        throw std::logic_error
        ("AL and ABR must have compatible row alignments");
    if( orientation == NORMAL )
        throw std::logic_error("DistFrontLDL must be (conjugate-)transposed.");
#endif
    const Grid& g = AL.Grid();

    // Matrix views
    DistMatrix<F>
        ALTL(g), ALTR(g),  AL00(g), AL01(g), AL02(g),
        ALBL(g), ALBR(g),  AL10(g), AL11(g), AL12(g),
                           AL20(g), AL21(g), AL22(g);

    // Temporary matrices
    DistMatrix<F,STAR,STAR> AL11_STAR_STAR(g);
    DistMatrix<F,STAR,STAR> d1_STAR_STAR(g);
    DistMatrix<F,VC,  STAR> AL21_VC_STAR(g);
    DistMatrix<F,VR,  STAR> AL21_VR_STAR(g);
    DistMatrix<F,STAR,MC  > S21Trans_STAR_MC(g);
    DistMatrix<F,STAR,MR  > AL21AdjOrTrans_STAR_MR(g);

    DistMatrix<F,STAR,MC> leftL(g), leftR(g);
    DistMatrix<F,STAR,MR> rightL(g), rightR(g);
    DistMatrix<F> AL22T(g), 
                  AL22B(g);

    // Start the algorithm
    elem::PartitionDownDiagonal
    ( AL, ALTL, ALTR,
          ALBL, ALBR, 0 );
    while( ALTL.Width() < AL.Width() )
    {
        elem::RepartitionDownDiagonal
        ( ALTL, /**/ ALTR,  AL00, /**/ AL01, AL02,
         /***************/ /*********************/
               /**/         AL10, /**/ AL11, AL12,
          ALBL, /**/ ALBR,  AL20, /**/ AL21, AL22 );

        AL21_VC_STAR.AlignWith( AL22 );
        AL21_VR_STAR.AlignWith( AL22 );
        S21Trans_STAR_MC.AlignWith( AL22 );
        AL21AdjOrTrans_STAR_MR.AlignWith( AL22 );
        //--------------------------------------------------------------------//
        AL11_STAR_STAR = AL11; 
        elem::internal::LocalLDL
        ( orientation, AL11_STAR_STAR, d1_STAR_STAR );
        AL11 = AL11_STAR_STAR;

        AL21_VC_STAR = AL21;
        elem::internal::LocalTrsm
        ( RIGHT, LOWER, orientation, UNIT, 
          (F)1, AL11_STAR_STAR, AL21_VC_STAR );

        S21Trans_STAR_MC.TransposeFrom( AL21_VC_STAR );
        elem::DiagonalSolve
        ( RIGHT, NORMAL, d1_STAR_STAR, AL21_VC_STAR );
        AL21_VR_STAR = AL21_VC_STAR;
        if( orientation == ADJOINT )
            AL21AdjOrTrans_STAR_MR.AdjointFrom( AL21_VR_STAR );
        else
            AL21AdjOrTrans_STAR_MR.TransposeFrom( AL21_VR_STAR );

        // Partition the update of the bottom-right corner into three pieces
        PartitionRight
        ( S21Trans_STAR_MC, 
          leftL, leftR, AL22.Width() );
        PartitionRight
        ( AL21AdjOrTrans_STAR_MR,
          rightL, rightR, AL22.Width() );
        PartitionDown
        ( AL22, AL22T,
                AL22B, AL22.Width() );
        elem::internal::LocalTrrk
        ( LOWER, orientation, (F)-1, leftL, rightL, (F)1, AL22T );
        elem::internal::LocalGemm
        ( orientation, NORMAL, (F)-1, leftR, rightL, (F)1, AL22B );
        elem::internal::LocalTrrk
        ( LOWER, orientation, (F)-1, leftR, rightR, (F)1, ABR );

        elem::DiagonalSolve
        ( LEFT, NORMAL, d1_STAR_STAR, S21Trans_STAR_MC );
        AL21.TransposeFrom( S21Trans_STAR_MC );
        //--------------------------------------------------------------------//
        AL21_VC_STAR.FreeAlignments();
        AL21_VR_STAR.FreeAlignments();
        S21Trans_STAR_MC.FreeAlignments();
        AL21AdjOrTrans_STAR_MR.FreeAlignments();

        elem::SlidePartitionDownDiagonal
        ( ALTL, /**/ ALTR,  AL00, AL01, /**/ AL02,
               /**/         AL10, AL11, /**/ AL12,
         /***************/ /*********************/
          ALBL, /**/ ALBR,  AL20, AL21, /**/ AL22 );
    }
#ifndef RELEASE
    PopCallStack();
#endif
}

template<typename F>
inline void DistFrontLDLSquare
( Orientation orientation, DistMatrix<F>& AL, DistMatrix<F>& ABR )
{
#ifndef RELEASE
    PushCallStack("numeric::internal::DistFrontLDLSquare");
    if( ABR.Height() != ABR.Width() )
        throw std::logic_error("ABR must be square");
    if( AL.Height() != AL.Width()+ABR.Height() )
        throw std::logic_error("AL & ABR must have compatible dimensions");
    if( AL.Grid() != ABR.Grid() )
        throw std::logic_error("AL & ABR must use the same grid");
    if( ABR.ColAlignment() !=
        (AL.ColAlignment()+AL.Width()) % AL.Grid().Height() )
        throw std::logic_error("AL & ABR must have compatible col alignments");
    if( ABR.RowAlignment() != 
        (AL.RowAlignment()+AL.Width()) % AL.Grid().Width() )
        throw std::logic_error("AL & ABR must have compatible row alignments");
    if( orientation == NORMAL )
        throw std::logic_error("DistFrontLDL must be (conjugate-)transposed.");
#endif
    const Grid& g = AL.Grid();
#ifndef RELEASE
    if( g.Height() != g.Width() )
        throw std::logic_error("This routine assumes a square process grid");
#endif

    // Find the process holding our transposed data
    const int r = g.Height();
    int transposeRank;
    {
        const int colAlignment = AL.ColAlignment();
        const int rowAlignment = AL.RowAlignment();
        const int colShift = AL.ColShift();
        const int rowShift = AL.RowShift();

        const int transposeRow = (colAlignment+rowShift) % r;
        const int transposeCol = (rowAlignment+colShift) % r;
        transposeRank = transposeRow + r*transposeCol;
    }
    const bool onDiagonal = ( transposeRank == g.VCRank() );

    // Matrix views
    DistMatrix<F>
        ALTL(g), ALTR(g),  AL00(g), AL01(g), AL02(g),
        ALBL(g), ALBR(g),  AL10(g), AL11(g), AL12(g),
                           AL20(g), AL21(g), AL22(g);

    // Temporary matrices
    DistMatrix<F,STAR,STAR> AL11_STAR_STAR(g);
    DistMatrix<F,STAR,STAR> d1_STAR_STAR(g);
    DistMatrix<F,VC,  STAR> AL21_VC_STAR(g);
    DistMatrix<F,STAR,MC  > S21Trans_STAR_MC(g);
    DistMatrix<F,STAR,MR  > AL21AdjOrTrans_STAR_MR(g);

    DistMatrix<F,STAR,MC> leftL(g), leftR(g);
    DistMatrix<F,STAR,MR> rightL(g), rightR(g);
    DistMatrix<F> AL22T(g), 
                  AL22B(g);
    DistMatrix<F> AR2T(g), 
                  AR2B(g);

    // Start the algorithm
    elem::PartitionDownDiagonal
    ( AL, ALTL, ALTR,
          ALBL, ALBR, 0 );
    while( ALTL.Width() < AL.Width() )
    {
        elem::RepartitionDownDiagonal
        ( ALTL, /**/ ALTR,  AL00, /**/ AL01, AL02,
         /***************/ /*********************/
               /**/         AL10, /**/ AL11, AL12,
          ALBL, /**/ ALBR,  AL20, /**/ AL21, AL22 );

        AL21_VC_STAR.AlignWith( AL22 );
        S21Trans_STAR_MC.AlignWith( AL22 );
        AL21AdjOrTrans_STAR_MR.AlignWith( AL22 );
        //--------------------------------------------------------------------//
        AL11_STAR_STAR = AL11; 
        elem::internal::LocalLDL
        ( orientation, AL11_STAR_STAR, d1_STAR_STAR );
        AL11 = AL11_STAR_STAR;

        AL21_VC_STAR = AL21;
        elem::internal::LocalTrsm
        ( RIGHT, LOWER, orientation, UNIT, 
          (F)1, AL11_STAR_STAR, AL21_VC_STAR );

        S21Trans_STAR_MC.TransposeFrom( AL21_VC_STAR );
        // SendRecv to form AL21^T[* ,MR] from S21^T[* ,MC], then conjugate
        // if necessary.
        AL21AdjOrTrans_STAR_MR.ResizeTo( AL21.Width(), AL21.Height() );
        {
            if( onDiagonal )
            {
                const int size = AL21.LocalHeight()*AL11.Width();    
                elem::MemCopy
                ( AL21AdjOrTrans_STAR_MR.LocalBuffer(), 
                  S21Trans_STAR_MC.LocalBuffer(), size );
            }
            else
            {
                const int sendSize = AL21.LocalHeight()*AL11.Width();
                const int recvSize = 
                    (AL22.LocalWidth()+ABR.LocalWidth())*AL11.Height();
                // We know that the ldim is the height since we have manually
                // created both temporary matrices.
                mpi::SendRecv
                ( S21Trans_STAR_MC.LocalBuffer(), sendSize, transposeRank, 0,
                  AL21AdjOrTrans_STAR_MR.LocalBuffer(), 
                  recvSize, transposeRank, 0, g.VCComm() );
            }
            elem::DiagonalSolve
            ( LEFT, NORMAL, d1_STAR_STAR, AL21AdjOrTrans_STAR_MR );
            if( orientation == ADJOINT )
                elem::Conjugate( AL21AdjOrTrans_STAR_MR );
        }

        // Partition the update of the bottom-right corner into three pieces
        PartitionRight
        ( S21Trans_STAR_MC, 
          leftL, leftR, AL22.Width() );
        PartitionRight
        ( AL21AdjOrTrans_STAR_MR,
          rightL, rightR, AL22.Width() );
        PartitionDown
        ( AL22, AL22T,
                AL22B, AL22.Width() );
        elem::internal::LocalTrrk
        ( LOWER, orientation, (F)-1, leftL, rightL, (F)1, AL22T );
        elem::internal::LocalGemm
        ( orientation, NORMAL, (F)-1, leftR, rightL, (F)1, AL22B );
        elem::internal::LocalTrrk
        ( LOWER, orientation, (F)-1, leftR, rightR, (F)1, ABR );

        elem::DiagonalSolve
        ( LEFT, NORMAL, d1_STAR_STAR, S21Trans_STAR_MC );
        AL21.TransposeFrom( S21Trans_STAR_MC );
        //--------------------------------------------------------------------//
        AL21_VC_STAR.FreeAlignments();
        S21Trans_STAR_MC.FreeAlignments();
        AL21AdjOrTrans_STAR_MR.FreeAlignments();

        elem::SlidePartitionDownDiagonal
        ( ALTL, /**/ ALTR,  AL00, AL01, /**/ AL02,
               /**/         AL10, AL11, /**/ AL12,
         /***************/ /*********************/
          ALBL, /**/ ALBR,  AL20, AL21, /**/ AL22 );
    }
#ifndef RELEASE
    PopCallStack();
#endif
}

} // namespace internal

template<typename F> 
inline void DistFrontLDL
( Orientation orientation, DistMatrix<F>& AL, DistMatrix<F>& ABR )
{
#ifndef RELEASE
    PushCallStack("numeric::DistFrontLDL");
#endif
    const Grid& grid = AL.Grid();
    if( grid.Height() == grid.Width() )
        internal::DistFrontLDLSquare( orientation, AL, ABR );
    else
        internal::DistFrontLDLGeneral( orientation, AL, ABR );
#ifndef RELEASE
    PopCallStack();
#endif
}

} // namespace numeric
} // namespace cliq

#endif // CLIQUE_NUMERIC_DIST_FRONT_LDL_HPP
