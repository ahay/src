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
#ifndef CLIQUE_NUMERIC_LOCAL_BLOCK_LOWER_SOLVE_HPP
#define CLIQUE_NUMERIC_LOCAL_BLOCK_LOWER_SOLVE_HPP 1

namespace cliq {
namespace numeric {

template<typename F> 
void LocalBlockLowerForwardSolve
( const symbolic::SymmFact& S,
  const numeric::SymmFrontTree<F>& L,
        Matrix<F>& X );

template<typename F> 
void LocalBlockLowerBackwardSolve
( Orientation orientation, 
  const symbolic::SymmFact& S, 
  const numeric::SymmFrontTree<F>& L,
        Matrix<F>& X );

//----------------------------------------------------------------------------//
// Implementation begins here                                                 //
//----------------------------------------------------------------------------//

template<typename F> 
inline void LocalBlockLowerForwardSolve
( const symbolic::SymmFact& S,
  const numeric::SymmFrontTree<F>& L,
        Matrix<F>& X )
{
    using namespace symbolic;
#ifndef RELEASE
    PushCallStack("numeric::LocalBlockLowerForwardSolve");
#endif
    const int numLocalSupernodes = S.local.supernodes.size();
    const int width = X.Width();
    for( int s=0; s<numLocalSupernodes; ++s )
    {
        const LocalSymmFactSupernode& sn = S.local.supernodes[s];
        const Matrix<F>& frontL = L.local.fronts[s].frontL;
        Matrix<F>& W = L.local.fronts[s].work;

        // Set up a workspace
        W.ResizeTo( frontL.Height(), width );
        Matrix<F> WT, WB;
        elem::PartitionDown
        ( W, WT,
             WB, sn.size );

        // Pull in the relevant information from the RHS
        Matrix<F> XT;
        XT.View( X, sn.myOffset, 0, sn.size, width );
        WT = XT;
        elem::MakeZeros( WB );

        // Update using the children (if they exist)
        const int numChildren = sn.children.size();
        if( numChildren == 2 )
        {
            const int leftIndex = sn.children[0];
            const int rightIndex = sn.children[1];
            Matrix<F>& leftWork = L.local.fronts[leftIndex].work;
            Matrix<F>& rightWork = L.local.fronts[rightIndex].work;
            const int leftSupernodeSize = S.local.supernodes[leftIndex].size;
            const int rightSupernodeSize = S.local.supernodes[rightIndex].size;
            const int leftUpdateSize = leftWork.Height()-leftSupernodeSize;
            const int rightUpdateSize = rightWork.Height()-rightSupernodeSize;

            // Add the left child's update onto ours
            Matrix<F> leftUpdate;
            leftUpdate.LockedView
            ( leftWork, leftSupernodeSize, 0, leftUpdateSize, width );
            for( int iChild=0; iChild<leftUpdateSize; ++iChild )
            {
                const int iFront = sn.leftChildRelIndices[iChild]; 
                for( int j=0; j<width; ++j )
                    W.Update( iFront, j, leftUpdate.Get(iChild,j) );
            }
            leftWork.Empty();

            // Add the right child's update onto ours
            Matrix<F> rightUpdate;
            rightUpdate.LockedView
            ( rightWork, rightSupernodeSize, 0, rightUpdateSize, width );
            for( int iChild=0; iChild<rightUpdateSize; ++iChild )
            {
                const int iFront = sn.rightChildRelIndices[iChild];
                for( int j=0; j<width; ++j )
                    W.Update( iFront, j, rightUpdate.Get(iChild,j) );
            }
            rightWork.Empty();
        }
        // else numChildren == 0

        // Solve against this front
        LocalFrontBlockLowerForwardSolve( frontL, W );

        // Store the supernode portion of the result
        XT = WT;
    }
#ifndef RELEASE
    PopCallStack();
#endif
}

template<typename F> 
inline void LocalBlockLowerBackwardSolve
( Orientation orientation, 
  const symbolic::SymmFact& S, 
  const numeric::SymmFrontTree<F>& L,
        Matrix<F>& X )
{
    using namespace symbolic;
#ifndef RELEASE
    PushCallStack("numeric::LocalBlockLowerBackwardSolve");
#endif
    const int numLocalSupernodes = S.local.supernodes.size();
    const int width = X.Width();
    if( width == 0 )
    {
#ifndef RELEASE
        PopCallStack();
#endif
        return;
    }

    // Pull in the top local information from the bottom distributed information
    L.local.fronts.back().work.LockedView
    ( L.dist.fronts[0].work1d.LocalMatrix() );

    for( int s=numLocalSupernodes-2; s>=0; --s )
    {
        const LocalSymmFactSupernode& sn = S.local.supernodes[s];
        const Matrix<F>& frontL = L.local.fronts[s].frontL;
        Matrix<F>& W = L.local.fronts[s].work;

        // Set up a workspace
        W.ResizeTo( frontL.Height(), width );
        Matrix<F> WT, WB;
        elem::PartitionDown
        ( W, WT,
             WB, sn.size );

        // Pull in the relevant information from the RHS
        Matrix<F> XT;
        XT.View( X, sn.myOffset, 0, sn.size, width );
        WT = XT;

        // Update using the parent
        const int parent = sn.parent;
        Matrix<F>& parentWork = L.local.fronts[parent].work;
        const LocalSymmFactSupernode& parentSN = S.local.supernodes[parent];
        const int currentUpdateSize = WB.Height();
        const std::vector<int>& parentRelIndices = 
          ( sn.isLeftChild ? 
            parentSN.leftChildRelIndices :
            parentSN.rightChildRelIndices );
        for( int iCurrent=0; iCurrent<currentUpdateSize; ++iCurrent )
        {
            const int iParent = parentRelIndices[iCurrent];
            for( int j=0; j<width; ++j )
                WB.Set( iCurrent, j, parentWork.Get(iParent,j) );
        }

        // The left child is numbered lower than the right child, so 
        // we can safely free the parent's work if we are the left child
        if( sn.isLeftChild )
        {
            parentWork.Empty();
            if( parent == numLocalSupernodes-1 )
                L.dist.fronts[0].work1d.Empty();
        }

        // Solve against this front
        LocalFrontBlockLowerBackwardSolve( orientation, frontL, W );

        // Store the supernode portion of the result
        XT = WT;
    }

    // Ensure that all of the temporary buffers are freed (this is overkill)
    L.dist.fronts[0].work1d.Empty();
    for( int s=0; s<numLocalSupernodes; ++s )
        L.local.fronts[s].work.Empty();
#ifndef RELEASE
    PopCallStack();
#endif
}

} // namespace numeric
} // namespace cliq

#endif // CLIQUE_NUMERIC_LOCAL_BLOCK_LOWER_SOLVE_HPP
