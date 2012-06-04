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
#ifndef CLIQUE_NUMERIC_DIST_FRONT_BLOCK_LDL_HPP
#define CLIQUE_NUMERIC_DIST_FRONT_BLOCK_LDL_HPP 1

namespace cliq {
namespace numeric {

template<typename F> 
void DistFrontBlockLDL
( Orientation orientation, DistMatrix<F>& AL, DistMatrix<F>& ABR );

//----------------------------------------------------------------------------//
// Implementation begins here                                                 //
//----------------------------------------------------------------------------//

template<typename F> 
inline void DistFrontBlockLDL
( Orientation orientation, DistMatrix<F>& AL, DistMatrix<F>& ABR )
{
#ifndef RELEASE
    PushCallStack("numeric::internal::DistFrontBlockLDL");
#endif
    const Grid& g = AL.Grid();
    DistMatrix<F> ATL(g),
                        ABL(g);
    elem::PartitionDown
    ( AL, ATL,
          ABL, AL.Width() );

    // Make a copy of the original contents of ABL
    DistMatrix<F> BBL( ABL );

    // Call the standard routine
    DistFrontLDL( orientation, AL, ABR );

    // Copy the original contents of ABL back
    ABL = BBL;

    // Overwrite ATL with inv(L D L^{T/H}) = L^{-T/H} D^{-1} L^{-1}
    // Overwrite ATL with inv(L D L^[T/H]) = L^[-T/H] D^{-1} L^{-1}
    elem::TriangularInverse( LOWER, UNIT, ATL );
    elem::Trdtrmm( orientation, LOWER, ATL );
    elem::MakeTrapezoidal( elem::LEFT, elem::LOWER, 0, ATL );
    if( orientation == TRANSPOSE )
    {
        elem::DistMatrix<F> ATLTrans( g );
        elem::Transpose( ATL, ATLTrans );
        elem::MakeTrapezoidal( elem::LEFT, elem::UPPER, 1, ATLTrans );
        elem::Axpy( (F)1, ATLTrans, ATL );
    }
    else
    {
        elem::DistMatrix<F> ATLAdj( g );
        elem::Adjoint( ATL, ATLAdj );
        elem::MakeTrapezoidal( elem::LEFT, elem::UPPER, 1, ATLAdj );
        elem::Axpy( (F)1, ATLAdj, ATL );
    }
#ifndef RELEASE
    PopCallStack();
#endif
}

} // namespace numeric
} // namespace cliq

#endif // CLIQUE_NUMERIC_DIST_FRONT_BLOCK_LDL_HPP
