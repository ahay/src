/*
   Clique: a scalable implementation of the multifrontal algorithm

   Copyright (C) 2010-2011 Jack Poulson <jack.poulson@gmail.com>
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
#ifndef CLIQUE_NUMERIC_LDL_HPP
#define CLIQUE_NUMERIC_LDL_HPP 1

namespace clique {
namespace numeric {

// All fronts of L are required to be initialized to the expansions of the 
// original sparse matrix before calling the following factorizations.

template<typename F>
void InitializeDistLeaf
( const symbolic::SymmFact& S, numeric::SymmFrontTree<F>& L );

template<typename F>
void LDL
( Orientation orientation, 
  symbolic::SymmFact& S, numeric::SymmFrontTree<F>& L );

template<typename F>
void LocalLDL
( Orientation orientation, 
  symbolic::SymmFact& S, numeric::SymmFrontTree<F>& L );

template<typename F>
void DistLDL
( Orientation orientation, 
  symbolic::SymmFact& S, numeric::SymmFrontTree<F>& L );

//----------------------------------------------------------------------------//
// Implementation begins here                                                 //
//----------------------------------------------------------------------------//

template<typename F>
void InitializeDistLeaf
( const symbolic::SymmFact& S, numeric::SymmFrontTree<F>& L )
{
#ifndef RELEASE
    PushCallStack("numeric::InitializeDistLeaf");
#endif
    const symbolic::DistSymmFactSupernode& sn = S.dist.supernodes[0];
    Matrix<F>& topLocalFrontL = L.local.fronts.back().frontL;
    DistMatrix<F,MC,MR>& front2dL = L.dist.fronts[0].front2dL;

    front2dL.LockedView
    ( topLocalFrontL.Height(), topLocalFrontL.Width(), 0, 0,
      topLocalFrontL.LockedBuffer(), topLocalFrontL.LDim(), 
      *sn.grid );
#ifndef RELEASE
    PopCallStack();
#endif
}

template<typename F>
void LDL
( Orientation orientation, 
  symbolic::SymmFact& S, numeric::SymmFrontTree<F>& L )
{
#ifndef RELEASE
    PushCallStack("numeric::LDL");
#endif
    LocalLDL( orientation, S, L );
    DistLDL( orientation, S, L );
#ifndef RELEASE
    PopCallStack();
#endif
}

} // namespace numeric
} // namespace clique

#endif /* CLIQUE_NUMERIC_LDL_HPP */

