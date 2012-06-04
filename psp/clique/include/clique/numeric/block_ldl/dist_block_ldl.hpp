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
#ifndef CLIQUE_NUMERIC_DIST_BLOCK_LDL_HPP 
#define CLIQUE_NUMERIC_DIST_BLOCK_LDL_HPP 1

// NOTE: This routine is almost identical to DistLDL, so perhaps it could
//       be partially merged.

namespace cliq {
namespace numeric {

template<typename F> 
void DistBlockLDL
( Orientation orientation, 
  symbolic::SymmFact& S, numeric::SymmFrontTree<F>& L );

//----------------------------------------------------------------------------//
// Implementation begins here                                                 //
//----------------------------------------------------------------------------//

template<typename F> 
inline void DistBlockLDL
( Orientation orientation, symbolic::SymmFact& S, numeric::SymmFrontTree<F>& L )
{
    using namespace symbolic;
#ifndef RELEASE
    PushCallStack("numeric::DistBlockLDL");
    if( orientation == NORMAL )
        throw std::logic_error("LDL must be (conjugate-)transposed");
#endif
    L.dist.mode = NORMAL_2D;

    // The bottom front is already computed, so just view it
    LocalSymmFront<F>& topLocalFront = L.local.fronts.back();
    DistSymmFront<F>& bottomDistFront = L.dist.fronts[0];
    const Grid& bottomGrid = *S.dist.supernodes[0].grid;
    bottomDistFront.front2dL.Empty();
    bottomDistFront.front2dL.LockedView
    ( topLocalFront.frontL.Height(), topLocalFront.frontL.Width(), 0, 0, 
      topLocalFront.frontL.LockedBuffer(), topLocalFront.frontL.LDim(), 
      bottomGrid );
    bottomDistFront.work2d.Empty();
    bottomDistFront.work2d.LockedView
    ( topLocalFront.work.Height(), topLocalFront.work.Width(), 0, 0,
      topLocalFront.work.LockedBuffer(), topLocalFront.work.LDim(),
      bottomGrid );

    // Perform the distributed portion of the factorization
    const unsigned numDistSupernodes = S.dist.supernodes.size();
    for( unsigned s=1; s<numDistSupernodes; ++s )
    {
        const DistSymmFactSupernode& childSN = S.dist.supernodes[s-1];
        const DistSymmFactSupernode& sn = S.dist.supernodes[s];
        const int updateSize = sn.lowerStruct.size();
        DistSymmFront<F>& childFront = L.dist.fronts[s-1];
        DistSymmFront<F>& front = L.dist.fronts[s];
        front.work2d.Empty();
#ifndef RELEASE
        if( front.front2dL.Height() != sn.size+updateSize ||
            front.front2dL.Width() != sn.size )
            throw std::logic_error("Front was not the proper size");
#endif

        const bool computeFactRecvIndices = 
            ( sn.childFactRecvIndices.size() == 0 );

        // Grab this front's grid information
        const Grid& grid = front.front2dL.Grid();
        mpi::Comm comm = grid.VCComm();
        const unsigned commRank = mpi::CommRank( comm );
        const unsigned commSize = mpi::CommSize( comm );
        const unsigned gridHeight = grid.Height();
        const unsigned gridWidth = grid.Width();

        // Grab the child's grid information
        const Grid& childGrid = childFront.front2dL.Grid();
        const unsigned childGridHeight = childGrid.Height();
        const unsigned childGridWidth = childGrid.Width();

        // Pack our child's update
        const DistMatrix<F>& childUpdate = childFront.work2d;
        const bool isLeftChild = ( commRank < commSize/2 );
        std::vector<int> sendCounts(commSize), sendDispls(commSize);
        int sendBufferSize = 0;
        for( unsigned proc=0; proc<commSize; ++proc )
        {
            const int sendSize = sn.numChildFactSendIndices[proc];
            sendCounts[proc] = sendSize;
            sendDispls[proc] = sendBufferSize;
            sendBufferSize += sendSize;
        }
        std::vector<F> sendBuffer( sendBufferSize );

        const std::vector<int>& myChildRelIndices = 
            ( isLeftChild ? sn.leftChildRelIndices
                          : sn.rightChildRelIndices );
        const int updateColShift = childUpdate.ColShift();
        const int updateRowShift = childUpdate.RowShift();
        const int updateLocalHeight = childUpdate.LocalHeight();
        const int updateLocalWidth = childUpdate.LocalWidth();
        std::vector<int> packOffsets = sendDispls;
        for( int jChildLocal=0; jChildLocal<updateLocalWidth; ++jChildLocal )
        {
            const int jChild = updateRowShift + jChildLocal*childGridWidth;
            const int destGridCol = myChildRelIndices[jChild] % gridWidth;

            int localColShift;
            if( updateColShift > jChild )
                localColShift = 0;
            else if( (jChild-updateColShift) % childGridHeight == 0 )
                localColShift = (jChild-updateColShift)/childGridHeight;
            else
                localColShift = (jChild-updateColShift)/childGridHeight + 1;
            for( int iChildLocal=localColShift; 
                     iChildLocal<updateLocalHeight; ++iChildLocal )
            {
                const int iChild = updateColShift + iChildLocal*childGridHeight;
                const int destGridRow = myChildRelIndices[iChild] % gridHeight;

                const int destRank = destGridRow + destGridCol*gridHeight;
                sendBuffer[packOffsets[destRank]++] = 
                    childUpdate.GetLocalEntry(iChildLocal,jChildLocal);
            }
        }
#ifndef RELEASE
        for( unsigned proc=0; proc<commSize; ++proc )
        {
            if( packOffsets[proc]-sendDispls[proc] != 
                sn.numChildFactSendIndices[proc] )
                throw std::logic_error("Error in packing stage");
        }
#endif
        packOffsets.clear();
        childFront.work2d.Empty();
        if( s == 1 )
            topLocalFront.work.Empty();

        // Set up the recv buffer for the AllToAll
        if( computeFactRecvIndices )
            ComputeFactRecvIndices( sn, childSN );
        std::vector<int> recvCounts(commSize), recvDispls(commSize);
        int recvBufferSize=0;
        for( unsigned proc=0; proc<commSize; ++proc )
        {
            const int recvSize = sn.childFactRecvIndices[proc].size()/2;
            recvCounts[proc] = recvSize;
            recvDispls[proc] = recvBufferSize;
            recvBufferSize += recvSize;
        }
        std::vector<F> recvBuffer( recvBufferSize );
#ifndef RELEASE
        // Verify the send and recv counts match
        std::vector<int> actualRecvCounts(commSize);
        mpi::AllToAll
        ( &sendCounts[0],       1, 
          &actualRecvCounts[0], 1, comm );
        for( unsigned proc=0; proc<commSize; ++proc )
        {
            if( actualRecvCounts[proc] != recvCounts[proc] )
            {
                std::ostringstream msg;
                msg << "Expected recv count of " << recvCounts[proc]
                    << " but recv'd " << actualRecvCounts[proc] 
                    << " from process " << proc << " for supernode "
                    << s << "\n";
                throw std::logic_error( msg.str().c_str() );
            }
        }
        actualRecvCounts.clear();
#endif

        // AllToAll to send and receive the child updates
#ifdef USE_CUSTOM_ALLTOALLV_FOR_FACT
        int numSends=0,numRecvs=0;
        for( unsigned proc=0; proc<commSize; ++proc )
        {
            if( sendCounts[proc] != 0 )
                ++numSends;
            if( recvCounts[proc] != 0 )
                ++numRecvs;
        }
        std::vector<mpi::Status> statuses(numSends+numRecvs); 
        std::vector<mpi::Request> requests(numSends+numRecvs);
        int rCount=0;
        for( unsigned proc=0; proc<commSize; ++proc )
        {
            int count = recvCounts[proc];
            int displ = recvDispls[proc];
            if( count != 0 )
                mpi::IRecv
                ( &recvBuffer[displ], count, proc, 0, comm, requests[rCount++] );
        }
        for( unsigned proc=0; proc<commSize; ++proc )
        {
            int count = sendCounts[proc];
            int displ = sendDispls[proc];
            if( count != 0 )
                mpi::ISend
                ( &sendBuffer[displ], count, proc, 0, comm, requests[rCount++] );
        }
        mpi::WaitAll( numSends+numRecvs, &requests[0], &statuses[0] );
        statuses.clear();
        requests.clear();
#else
        mpi::AllToAll
        ( &sendBuffer[0], &sendCounts[0], &sendDispls[0],
          &recvBuffer[0], &recvCounts[0], &recvDispls[0], comm );
#endif
        sendBuffer.clear();
        sendCounts.clear();
        sendDispls.clear();

        // Unpack the child udpates (with an Axpy)
        front.work2d.SetGrid( front.front2dL.Grid() );
        front.work2d.Align( sn.size%grid.Height(), sn.size%grid.Width() );
        elem::Zeros( updateSize, updateSize, front.work2d );
        const int leftLocalWidth = front.front2dL.LocalWidth();
        const int topLocalHeight = 
            elem::LocalLength<int>( sn.size, grid.MCRank(), gridHeight );
        for( unsigned proc=0; proc<commSize; ++proc )
        {
            const F* recvValues = &recvBuffer[recvDispls[proc]];
            const std::deque<int>& recvIndices = sn.childFactRecvIndices[proc];
            const int numRecvIndexPairs = recvIndices.size()/2;
            for( int k=0; k<numRecvIndexPairs; ++k )
            {
                const int iFrontLocal = recvIndices[2*k+0];
                const int jFrontLocal = recvIndices[2*k+1];
                const F value = recvValues[k];
                if( jFrontLocal < leftLocalWidth )
                    front.front2dL.UpdateLocalEntry
                    ( iFrontLocal, jFrontLocal, value );
                else
                    front.work2d.UpdateLocalEntry
                    ( iFrontLocal-topLocalHeight, 
                      jFrontLocal-leftLocalWidth, value );
            }
        }
        recvBuffer.clear();
        recvCounts.clear();
        recvDispls.clear();
        if( computeFactRecvIndices )
            sn.childFactRecvIndices.clear();

        // Now that the frontal matrix is set up, perform the factorization
        DistFrontBlockLDL( orientation, front.front2dL, front.work2d );

        // Store the diagonal in a [VC,* ] distribution
        DistMatrix<F,MD,STAR> diag( grid );
        front.front2dL.GetDiagonal( diag );
        front.diag.SetGrid( grid );
        front.diag = diag;
    }
    L.local.fronts.back().work.Empty();
    L.dist.fronts.back().work2d.Empty();
#ifndef RELEASE
    PopCallStack();
#endif
}

} // namespace numeric
} // namespace cliq

#endif // CLIQUE_NUMERIC_DIST_BLOCK_LDL_HPP
