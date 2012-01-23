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
using namespace elemental;

template<typename F>
void clique::numeric::DistFrontFastLowerForwardSolve
( Diagonal diag, DistMatrix<F,VC,STAR>& L, DistMatrix<F,VC,STAR>& X )
{
#ifndef RELEASE
    clique::PushCallStack("numeric::DistFrontFastLowerForwardSolve");
    if( L.Grid() != X.Grid() )
        throw std::logic_error
        ("L and X must be distributed over the same grid");
    if( L.Height() < L.Width() || L.Height() != X.Height() )
    {
        std::ostringstream msg;
        msg << "Nonconformal solve:\n"
            << "  L ~ " << L.Height() << " x " << L.Width() << "\n"
            << "  X ~ " << X.Height() << " x " << X.Width() << "\n";
        throw std::logic_error( msg.str().c_str() );
    }
    if( L.ColAlignment() != X.ColAlignment() )
        throw std::logic_error("L and X are assumed to be aligned");
#endif
    const Grid& g = L.Grid();
    const int commRank = g.VCRank();
    const int commSize = g.Size();
    if( commSize == 1 )
    {
        clique::numeric::LocalFrontLowerForwardSolve
        ( diag, L.LockedLocalMatrix(), X.LocalMatrix() );
#ifndef RELEASE
        PopCallStack();
#endif
        return;
    }

    // Separate the top and bottom portions of X and L
    const int snSize = L.Width();
    DistMatrix<F,VC,STAR> LT(g),
                          LB(g);
    PartitionDown
    ( L, LT, 
         LB, snSize );
    DistMatrix<F,VC,STAR> XT(g),
                          XB(g);
    PartitionDown
    ( X, XT,
         XB, snSize );

    const int localTopHeight = LT.LocalHeight();
    std::vector<F> localDiag;
    if( diag == UNIT )
    {
        // Extract the diagonal of the top triangle and replace it with ones
        localDiag.resize( localTopHeight );
        F* LTBuffer = LT.LocalBuffer();
        const int LTLDim = LT.LocalLDim();
        for( int iLocal=0; iLocal<localTopHeight; ++iLocal )
        {
            const int i = commRank + iLocal*commSize;
            localDiag[iLocal] = LTBuffer[iLocal+i*LTLDim];
            LTBuffer[iLocal+i*LTLDim] = 1;
        }
    }

    // Get a copy of all of XT
    DistMatrix<F,STAR,STAR> XT_STAR_STAR( XT );

    // XT := LT XT
    elemental::internal::LocalGemm
    ( NORMAL, NORMAL, (F)1, LT, XT_STAR_STAR, (F)0, XT );

    if( diag == UNIT )
    {
        // Put the diagonal back
        F* LTBuffer = LT.LocalBuffer();
        const int LTLDim = LT.LocalLDim();
        for( int iLocal=0; iLocal<localTopHeight; ++iLocal )
        {
            const int i = commRank + iLocal*commSize;
            LTBuffer[iLocal+i*LTLDim] = localDiag[iLocal];
        }
    }

    if( LB.Height() != 0 )
    {
        // Gather all of XT again
        XT_STAR_STAR = XT;

        // XB := XB - LB XT
        elemental::internal::LocalGemm
        ( NORMAL, NORMAL, (F)-1, LB, XT_STAR_STAR, (F)1, XB );
    }
#ifndef RELEASE
    PopCallStack();
#endif
}

template<typename F>
void clique::numeric::DistFrontFastLowerBackwardSolve
( Orientation orientation, Diagonal diag, 
  DistMatrix<F,VC,STAR>& L, DistMatrix<F,VC,STAR>& X )
{
#ifndef RELEASE
    clique::PushCallStack("numeric::DistFrontFastLowerBackwardSolve");
    if( L.Grid() != X.Grid() )
        throw std::logic_error
        ("L and X must be distributed over the same grid");
    if( L.Height() < L.Width() || L.Height() != X.Height() )
    {
        std::ostringstream msg;
        msg << "Nonconformal solve:\n"
            << "  L ~ " << L.Height() << " x " << L.Width() << "\n"
            << "  X ~ " << X.Height() << " x " << X.Width() << "\n";
        throw std::logic_error( msg.str().c_str() );
    }
    if( L.ColAlignment() != X.ColAlignment() )
        throw std::logic_error("L and X are assumed to be aligned");
    if( orientation == NORMAL )
        throw std::logic_error("This solve must be (conjugate-)transposed");
#endif
    const Grid& g = L.Grid();
    const int commSize = g.Size();
    const int commRank = g.VCRank();
    if( commSize == 1 )
    {
        LocalFrontLowerBackwardSolve
        ( orientation, diag, L.LockedLocalMatrix(), X.LocalMatrix() );
#ifndef RELEASE
        PopCallStack();
#endif
        return;
    }

    const int snSize = L.Width();
    DistMatrix<F,VC,STAR> LT(g),
                          LB(g);
    elemental::PartitionDown
    ( L, LT,
         LB, snSize );
    DistMatrix<F,VC,STAR> XT(g),
                          XB(g);
    elemental::PartitionDown
    ( X, XT,
         XB, snSize );

    // XT := XT - LB^{T/H} XB
    DistMatrix<F,STAR,STAR> Z( snSize, XT.Width(), g );
    if( XB.Height() != 0 )
    {
        elemental::internal::LocalGemm
        ( orientation, NORMAL, (F)-1, LB, XB, (F)0, Z );
        XT.SumScatterUpdate( (F)1, Z );
    }

    const int localTopHeight = LT.LocalHeight();
    std::vector<F> localDiag;
    if( diag == UNIT )
    {
        // Extract the diagonal of the top triangle and replace it with ones
        localDiag.resize( localTopHeight );
        F* LTBuffer = LT.LocalBuffer();
        const int LTLDim = LT.LocalLDim();
        for( int iLocal=0; iLocal<localTopHeight; ++iLocal )
        {
            const int i = commRank + iLocal*commSize;
            localDiag[iLocal] = LTBuffer[iLocal+i*LTLDim];
            LTBuffer[iLocal+i*LTLDim] = 1;
        }
    }

    // XT := LT^{T/H} XT
    elemental::internal::LocalGemm
    ( orientation, NORMAL, (F)1, LT, XT, (F)0, Z );
    XT.SumScatterFrom( Z );

    if( diag == UNIT )
    {
        // Put the diagonal back
        F* LTBuffer = LT.LocalBuffer();
        const int LTLDim = LT.LocalLDim();
        for( int iLocal=0; iLocal<localTopHeight; ++iLocal )
        {
            const int i = commRank + iLocal*commSize;
            LTBuffer[iLocal+i*LTLDim] = localDiag[iLocal];
        }
    }
#ifndef RELEASE
    PopCallStack();
#endif
}

template void clique::numeric::DistFrontFastLowerForwardSolve
( Diagonal diag, 
  DistMatrix<float,VC,STAR>& L,
  DistMatrix<float,VC,STAR>& X );
template void clique::numeric::DistFrontFastLowerBackwardSolve
( Orientation orientation, Diagonal diag,
  DistMatrix<float,VC,STAR>& L,
  DistMatrix<float,VC,STAR>& X );

template void clique::numeric::DistFrontFastLowerForwardSolve
( Diagonal diag,
  DistMatrix<double,VC,STAR>& L, 
  DistMatrix<double,VC,STAR>& X );
template void clique::numeric::DistFrontFastLowerBackwardSolve
( Orientation orientation, Diagonal diag,
  DistMatrix<double,VC,STAR>& L,
  DistMatrix<double,VC,STAR>& X );

template void clique::numeric::DistFrontFastLowerForwardSolve
( Diagonal diag,
  DistMatrix<std::complex<float>,VC,STAR>& L, 
  DistMatrix<std::complex<float>,VC,STAR>& X );
template void clique::numeric::DistFrontFastLowerBackwardSolve
( Orientation orientation, Diagonal diag,
  DistMatrix<std::complex<float>,VC,STAR>& L, 
  DistMatrix<std::complex<float>,VC,STAR>& X );

template void clique::numeric::DistFrontFastLowerForwardSolve
( Diagonal diag,
  DistMatrix<std::complex<double>,VC,STAR>& L, 
  DistMatrix<std::complex<double>,VC,STAR>& X );
template void clique::numeric::DistFrontFastLowerBackwardSolve
( Orientation orientation, Diagonal diag,
  DistMatrix<std::complex<double>,VC,STAR>& L,
  DistMatrix<std::complex<double>,VC,STAR>& X );
