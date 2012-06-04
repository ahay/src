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
#ifndef CLIQUE_NUMERIC_DIST_LOWER_MULTIPLY
#define CLIQUE_NUMERIC_DIST_LOWER_MULTIPLY 1

namespace cliq {
namespace numeric {

template<typename F> 
void DistLowerMultiplyNormal
( UnitOrNonUnit diag, int diagOffset,
  const symbolic::SymmFact& S,
  const numeric::SymmFrontTree<F>& L,
        Matrix<F>& localX );

template<typename F> 
void DistLowerMultiplyTranspose
( Orientation orientation, UnitOrNonUnit diag, int diagOffset,
  const symbolic::SymmFact& S,
  const numeric::SymmFrontTree<F>& L,
        Matrix<F>& localX );

//----------------------------------------------------------------------------//
// Implementation begins here                                                 //
//----------------------------------------------------------------------------//

template<typename F> 
inline void DistLowerMultiplyNormal
( UnitOrNonUnit diag, int diagOffset,
  const symbolic::SymmFact& S,
  const numeric::SymmFrontTree<F>& L,
        Matrix<F>& localX )
{
    using namespace symbolic;
#ifndef RELEASE
    PushCallStack("numeric::DistLowerMultiplyNormal");
#endif
    const int numDistSupernodes = S.dist.supernodes.size();
    const int width = localX.Width();
    if( L.dist.mode != NORMAL_1D )
        throw std::logic_error("This multiply mode is not yet implemented");
    if( width == 0 )
    {
#ifndef RELEASE
        PopCallStack();
#endif
        return;
    }

    // Copy the information from the local portion into the distributed leaf
    const LocalSymmFront<F>& localRootFront = L.local.fronts.back();
    const DistSymmFront<F>& distLeafFront = L.dist.fronts[0];
    distLeafFront.work1d.LockedView
    ( localRootFront.work.Height(), localRootFront.work.Width(), 0,
      localRootFront.work.LockedBuffer(), localRootFront.work.LDim(),
      distLeafFront.front1dL.Grid() );
    
    // Perform the distributed portion of the forward multiply
    for( int s=1; s<numDistSupernodes; ++s )
    {
        const DistSymmFactSupernode& childSN = S.dist.supernodes[s-1];
        const DistSymmFactSupernode& sn = S.dist.supernodes[s];
        const DistSymmFront<F>& childFront = L.dist.fronts[s-1];
        const DistSymmFront<F>& front = L.dist.fronts[s];
        const Grid& childGrid = childFront.front1dL.Grid();
        const Grid& grid = front.front1dL.Grid();
        mpi::Comm comm = grid.VCComm();
        mpi::Comm childComm = childGrid.VCComm();
        const int commSize = mpi::CommSize( comm );
        const int commRank = mpi::CommRank( comm );
        const int childCommSize = mpi::CommSize( childComm );
        const int childCommRank = mpi::CommRank( childComm );

        // Set up a workspace
        DistMatrix<F,VC,STAR>& W = front.work1d;
        W.SetGrid( grid );
        W.ResizeTo( front.front1dL.Height(), width );
        DistMatrix<F,VC,STAR> WT(grid), WB(grid);
        elem::PartitionDown
        ( W, WT,
             WB, sn.size );

        // Pull in the relevant information from the RHS
        Matrix<F> localXT;
        localXT.View( localX, sn.localOffset1d, 0, sn.localSize1d, width );
        WT.LocalMatrix() = localXT;
        elem::MakeZeros( WB );

        // Now that the right-hand side is set up, perform the multiply
        DistFrontLowerMultiply( NORMAL, diag, diagOffset, front.front1dL, W );

        // Pack our child's update
        DistMatrix<F,VC,STAR>& childW = childFront.work1d;
        const int updateSize = childW.Height()-childSN.size;
        DistMatrix<F,VC,STAR> childUpdate;
        childUpdate.LockedView( childW, childSN.size, 0, updateSize, width );
        int sendBufferSize = 0;
        std::vector<int> sendCounts(commSize), sendDispls(commSize);
        for( int proc=0; proc<commSize; ++proc )
        {
            const int sendSize = sn.numChildSolveSendIndices[proc]*width;
            sendCounts[proc] = sendSize;
            sendDispls[proc] = sendBufferSize;
            sendBufferSize += sendSize;
        }
        std::vector<F> sendBuffer( sendBufferSize );

        const bool isLeftChild = ( commRank < commSize/2 );
        const std::vector<int>& myChildRelIndices = 
            ( isLeftChild ? sn.leftChildRelIndices
                          : sn.rightChildRelIndices );
        const int updateColShift = childUpdate.ColShift();
        const int updateLocalHeight = childUpdate.LocalHeight();
        std::vector<int> packOffsets = sendDispls;
        for( int iChildLocal=0; iChildLocal<updateLocalHeight; ++iChildLocal )
        {
            const int iChild = updateColShift + iChildLocal*childCommSize;
            const int destRank = myChildRelIndices[iChild] % commSize;
            F* packBuf = &sendBuffer[packOffsets[destRank]];
            for( int jChild=0; jChild<width; ++jChild )
                packBuf[jChild] = childUpdate.GetLocalEntry(iChildLocal,jChild);
            packOffsets[destRank] += width;
        }
        packOffsets.clear();
        childW.Empty();
        if( s == 1 )
            L.local.fronts.back().work.Empty();

        // Set up the receive buffer
        int recvBufferSize = 0;
        std::vector<int> recvCounts(commSize), recvDispls(commSize);
        for( int proc=0; proc<commSize; ++proc )
        {
            const int recvSize = sn.childSolveRecvIndices[proc].size()*width;
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
        for( int proc=0; proc<commSize; ++proc )
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
#ifdef USE_CUSTOM_ALLTOALLV_FOR_MULT
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

        // Unpack the child updates (with an Axpy)
        for( int proc=0; proc<commSize; ++proc )
        {
            const F* recvValues = &recvBuffer[recvDispls[proc]];
            const std::deque<int>& recvIndices = sn.childSolveRecvIndices[proc];
            for( int k=0; k<recvIndices.size(); ++k )
            {
                const int iFrontLocal = recvIndices[k];
                const F* recvRow = &recvValues[k*width];
                F* WRow = W.LocalBuffer( iFrontLocal, 0 );
                const int WLDim = W.LocalLDim();
                for( int jFront=0; jFront<width; ++jFront )
                    WRow[jFront*WLDim] += recvRow[jFront];
            }
        }
        recvBuffer.clear();
        recvCounts.clear();
        recvDispls.clear();

        // Store the supernode portion of the result
        localXT = WT.LocalMatrix();
    }
    L.local.fronts.back().work.Empty();
    L.dist.fronts.back().work1d.Empty();
#ifndef RELEASE
    PopCallStack();
#endif
}

template<typename F> // F represents a real or complex field
inline void DistLowerMultiplyTranspose
( Orientation orientation, UnitOrNonUnit diag, int diagOffset,
  const symbolic::SymmFact& S,
  const numeric::SymmFrontTree<F>& L,
        Matrix<F>& localX )
{
    using namespace symbolic;
#ifndef RELEASE
    PushCallStack("numeric::DistLowerMultiplyTranspose");
#endif
    const int numDistSupernodes = S.dist.supernodes.size();
    const int width = localX.Width();
    if( L.dist.mode != NORMAL_1D )
        throw std::logic_error("This multiply mode is not yet implemented");
    if( width == 0 )
    {
#ifndef RELEASE
        PopCallStack();
#endif
        return;
    }

    // Directly operate on the root separator's portion of the right-hand sides
    const DistSymmFactSupernode& rootSN = S.dist.supernodes.back();
    const DistSymmFront<F>& rootFront = L.dist.fronts.back();
    const Grid& rootGrid = rootFront.front1dL.Grid();
    DistMatrix<F,VC,STAR> XRoot(rootGrid);
    XRoot.View
    ( rootSN.size, width, 0,
      localX.Buffer(rootSN.localOffset1d,0), localX.LDim(), rootGrid );
    rootFront.work1d = XRoot; // store the RHS for use by the children
    DistFrontLowerMultiply
    ( orientation, diag, diagOffset, rootFront.front1dL, XRoot );

    std::vector<int>::const_iterator it;
    for( int s=numDistSupernodes-2; s>=0; --s )
    {
        const DistSymmFactSupernode& parentSN = S.dist.supernodes[s+1];
        const DistSymmFactSupernode& sn = S.dist.supernodes[s];
        const DistSymmFront<F>& parentFront = L.dist.fronts[s+1];
        const DistSymmFront<F>& front = L.dist.fronts[s];
        const Grid& grid = front.front1dL.Grid();
        const Grid& parentGrid = parentFront.front1dL.Grid();
        mpi::Comm comm = grid.VCComm(); 
        mpi::Comm parentComm = parentGrid.VCComm();
        const int commSize = mpi::CommSize( comm );
        const int commRank = mpi::CommRank( comm );
        const int parentCommSize = mpi::CommSize( parentComm );
        const int parentCommRank = mpi::CommRank( parentComm );

        // Set up a copy of the RHS in our workspace.
        DistMatrix<F,VC,STAR>& W = front.work1d;
        W.SetGrid( grid );
        W.ResizeTo( front.front1dL.Height(), width );
        DistMatrix<F,VC,STAR> WT(grid), WB(grid);
        elem::PartitionDown
        ( W, WT,
             WB, sn.size );

        // Pull in the relevant information from the RHS
        Matrix<F> localXT;
        localXT.View( localX, sn.localOffset1d, 0, sn.localSize1d, width );
        WT.LocalMatrix() = localXT;

        //
        // Set the bottom from the parent's workspace
        //

        // Pack the relevant portions of the parent's RHS's
        // (which are stored in 'work1d')
        int sendBufferSize = 0;
        std::vector<int> sendCounts(parentCommSize), sendDispls(parentCommSize);
        for( int proc=0; proc<parentCommSize; ++proc )
        {
            const int sendSize = 
                parentSN.childSolveRecvIndices[proc].size()*width;
            sendCounts[proc] = sendSize;
            sendDispls[proc] = sendBufferSize;
            sendBufferSize += sendSize;
        }
        std::vector<F> sendBuffer( sendBufferSize );

        DistMatrix<F,VC,STAR>& parentWork = parentFront.work1d;
        for( int proc=0; proc<parentCommSize; ++proc )
        {
            F* sendValues = &sendBuffer[sendDispls[proc]];
            const std::deque<int>& recvIndices = 
                parentSN.childSolveRecvIndices[proc];
            for( int k=0; k<recvIndices.size(); ++k )
            {
                const int iFrontLocal = recvIndices[k];
                F* packedRow = &sendValues[k*width];
                const F* workRow = 
                    parentWork.LockedLocalBuffer( iFrontLocal, 0 );
                const int workLDim = parentWork.LocalLDim();
                for( int jFront=0; jFront<width; ++jFront )
                    packedRow[jFront] = workRow[jFront*workLDim];
            }
        }
        parentWork.Empty();

        // Set up the receive buffer
        int recvBufferSize = 0;
        std::vector<int> recvCounts(parentCommSize), recvDispls(parentCommSize);
        for( int proc=0; proc<parentCommSize; ++proc )
        {
            const int recvSize = 
                parentSN.numChildSolveSendIndices[proc]*width;
            recvCounts[proc] = recvSize;
            recvDispls[proc] = recvBufferSize;
            recvBufferSize += recvSize;
        }
        std::vector<F> recvBuffer( recvBufferSize );
#ifndef RELEASE
        // Verify the send and recv counts match
        std::vector<int> actualRecvCounts(parentCommSize);
        mpi::AllToAll
        ( &sendCounts[0],       1,
          &actualRecvCounts[0], 1, parentComm );
        for( int proc=0; proc<parentCommSize; ++proc )
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

        // AllToAll to send and recv parent updates
#ifdef USE_CUSTOM_ALLTOALLV_FOR_MULT
        int numSends=0,numRecvs=0;
        for( unsigned proc=0; proc<parentCommSize; ++proc )
        { 
            if( sendCounts[proc] != 0 )
                ++numSends;
            if( recvCounts[proc] != 0 )
                ++numRecvs;
        }
        std::vector<mpi::Status> statuses(numSends+numRecvs);
        std::vector<mpi::Request> requests(numSends+numRecvs);
        int rCount=0;
        for( unsigned proc=0; proc<parentCommSize; ++proc )
        {   
            int count = recvCounts[proc];
            int displ = recvDispls[proc];
            if( count != 0 )
                mpi::IRecv
                ( &recvBuffer[displ], count, proc, 0, parentComm, 
                  requests[rCount++] );
        }
        for( unsigned proc=0; proc<parentCommSize; ++proc )
        {
            int count = sendCounts[proc];
            int displ = sendDispls[proc];
            if( count != 0 )
                mpi::ISend
                ( &sendBuffer[displ], count, proc, 0, parentComm, 
                  requests[rCount++] );
        }
        mpi::WaitAll( numSends+numRecvs, &requests[0], &statuses[0] );
        statuses.clear();
        requests.clear();
#else
        mpi::AllToAll
        ( &sendBuffer[0], &sendCounts[0], &sendDispls[0],
          &recvBuffer[0], &recvCounts[0], &recvDispls[0], parentComm );
#endif
        sendBuffer.clear();
        sendCounts.clear();
        sendDispls.clear();

        // Unpack the updates using the send approach from the forward solve
        const bool isLeftChild = ( parentCommRank < parentCommSize/2 );
        const std::vector<int>& myRelIndices = 
            ( isLeftChild ? parentSN.leftChildRelIndices
                          : parentSN.rightChildRelIndices );
        const int updateColShift = WB.ColShift();
        const int updateLocalHeight = WB.LocalHeight();
        for( int iUpdateLocal=0; 
                 iUpdateLocal<updateLocalHeight; ++iUpdateLocal )
        {
            const int iUpdate = updateColShift + iUpdateLocal*commSize;
            const int startRank = myRelIndices[iUpdate] % parentCommSize;
            const F* recvBuf = &recvBuffer[recvDispls[startRank]];
            for( int jUpdate=0; jUpdate<width; ++jUpdate )
                WB.SetLocalEntry(iUpdateLocal,jUpdate,recvBuf[jUpdate]);
            recvDispls[startRank] += width;
        }
        recvBuffer.clear();
        recvCounts.clear();
        recvDispls.clear();

        // Make a copy of the unmodified RHS
        DistMatrix<F,VC,STAR> XNode( front.work1d );

        // Perform the multiply for this front
        DistFrontLowerMultiply
        ( orientation, diag, diagOffset, front.front1dL, XNode );

        // Store the supernode portion of the result
        DistMatrix<F,VC,STAR> XNodeT(grid), XNodeB(grid);
        elem::PartitionDown
        ( XNode, XNodeT,
                 XNodeB, sn.size );
        localXT = XNodeT.LocalMatrix();
        XNode.Empty();
    }
#ifndef RELEASE
    PopCallStack();
#endif
}

} // namespace numeric
} // namespace cliq

#endif // CLIQUE_NUMERIC_DIST_LOWER_MULTIPLY
