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
#ifndef CLIQUE_NUMERIC_SET_SOLVE_MODE_HPP
#define CLIQUE_NUMERIC_SET_SOLVE_MODE_HPP 1

namespace cliq {
namespace numeric {

template<typename F>
void SetSolveMode( SymmFrontTree<F>& L, SolveMode mode );

//----------------------------------------------------------------------------//
// Implementation begins here                                                 //
//----------------------------------------------------------------------------//

// This routine could be modified later so that it uses much less memory
// by replacing the '=' redistributions with piece-by-piece redistributions.
template<typename F>
inline void SetSolveMode( SymmFrontTree<F>& L, SolveMode mode )
{
#ifndef RELEASE
    PushCallStack("numeric::SetSolveMode");
#endif
    // Check if this call can be a no-op
    if( mode == L.dist.mode ) 
    {
#ifndef RELEASE
        PopCallStack();
#endif
        return;
    }
    const int numDistSupernodes = L.dist.fronts.size();    
    DistSymmFront<F>& leafFront = L.dist.fronts[0];
    const SolveMode oldMode = L.dist.mode;

    if( mode == NORMAL_1D && oldMode == NORMAL_2D )
    {
        leafFront.front1dL.LockedView
        ( leafFront.front2dL.Height(), 
          leafFront.front2dL.Width(), 0,
          leafFront.front2dL.LockedLocalBuffer(), 
          leafFront.front2dL.LocalLDim(),
          leafFront.front2dL.Grid() );
        for( int s=1; s<numDistSupernodes; ++s )
        {
            DistSymmFront<F>& front = L.dist.fronts[s];
            front.front1dL.Empty();
            front.front1dL.SetGrid( front.front2dL.Grid() );
            front.front1dL = front.front2dL;
            front.front2dL.Empty();
        }
    }
    else if( mode == FAST_2D_LDL && oldMode == NORMAL_2D )
    {
        for( int s=1; s<numDistSupernodes; ++s )
        {
            DistSymmFront<F>& front = L.dist.fronts[s];
            const elem::Grid& grid = front.front2dL.Grid();
            const int snSize = front.front2dL.Width();

            // Invert the strictly lower portion of the diagonal block, and
            // then make the strictly upper triangle zero
            DistMatrix<F> LT( grid );
            LT.View( front.front2dL, 0, 0, snSize, snSize );
            elem::TriangularInverse( elem::LOWER, elem::UNIT, LT );
            elem::MakeTrapezoidal( elem::LEFT, elem::LOWER, 0, LT );
        }
    }
    else if( mode == FAST_1D_LDL && oldMode == NORMAL_2D )
    {
        leafFront.front1dL.LockedView
        ( leafFront.front2dL.Height(), 
          leafFront.front2dL.Width(), 0,
          leafFront.front2dL.LockedLocalBuffer(), 
          leafFront.front2dL.LocalLDim(),
          leafFront.front2dL.Grid() );
        for( int s=1; s<numDistSupernodes; ++s )
        {
            DistSymmFront<F>& front = L.dist.fronts[s];
            const elem::Grid& grid = front.front2dL.Grid();
            const int snSize = front.front2dL.Width();

            // Invert the strictly lower portion of the diagonal block 
            DistMatrix<F> LT( grid );
            LT.View( front.front2dL, 0, 0, snSize, snSize );
            elem::TriangularInverse( elem::LOWER, elem::UNIT, LT );

            // Copy the data and make the strictly upper triangle zero
            front.front1dL.Empty();
            front.front1dL.SetGrid( grid );
            front.front1dL = front.front2dL;
            front.front2dL.Empty();
            elem::MakeTrapezoidal( elem::LEFT, elem::LOWER, 0, front.front1dL );
        }
    }
    else if( mode == NORMAL_2D && oldMode == NORMAL_1D )
    {
        leafFront.front2dL.LockedView
        ( leafFront.front1dL.Height(), 
          leafFront.front1dL.Width(), 0, 0,
          leafFront.front1dL.LockedLocalBuffer(), 
          leafFront.front1dL.LocalLDim(),
          leafFront.front1dL.Grid() );
        for( int s=1; s<numDistSupernodes; ++s )
        {
            DistSymmFront<F>& front = L.dist.fronts[s];
            front.front2dL.Empty();
            front.front2dL.SetGrid( front.front1dL.Grid() );
            front.front2dL = front.front1dL;
            front.front1dL.Empty();
        }
    }
    else
        throw std::logic_error("Unavailable solve mode change");
    L.dist.mode = mode;
#ifndef RELEASE
    PopCallStack();
#endif
}

} // namespace numeric
} // namespace cliq

#endif // CLIQUE_NUMERIC_SET_SOLVE_MODE_HPP
