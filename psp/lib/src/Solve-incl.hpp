/*
   Parallel Sweeping Preconditioner (PSP): a distributed-memory implementation
   of a sweeping preconditioner for 3d Helmholtz equations.

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

template<typename R>
void
psp::DistHelmholtz<R>::SolveWithGMRES
( GridData<C>& gridB, int m, R relTol ) const
{
    if( !elemental::mpi::CongruentComms( comm_, gridB.Comm() ) )
        throw std::logic_error("B does not have a congruent comm");
    const int commRank = elemental::mpi::CommRank( comm_ );

    // Convert B into custom nested-dissection based ordering
    elemental::Matrix<C> B;
    {
        if( commRank == 0 )
        {
            std::cout << "  pulling right-hand sides...";
            std::cout.flush();
        }
        const double startTime = elemental::mpi::Time();

        PullRightHandSides( gridB, B );

        const double stopTime = elemental::mpi::Time();
        if( commRank == 0 )
            std::cout << stopTime-startTime << " secs" << std::endl;
    }

    // Solve the systems of equations
    InternalSolveWithGMRES( B, m, relTol );

    // Restore the solutions back into the GridData form
    {
        if( commRank == 0 )
        {
            std::cout << "  pushing right-hand sides...";
            std::cout.flush();
        }
        const double startTime = elemental::mpi::Time();

        PushRightHandSides( gridB, B );

        const double stopTime = elemental::mpi::Time();
        if( commRank == 0 )
            std::cout << stopTime-startTime << " secs" << std::endl;
    }
}

template<typename R>
void
psp::DistHelmholtz<R>::SolveWithSQMR( GridData<C>& gridB, R bcgRelTol ) const
{
    if( !elemental::mpi::CongruentComms( comm_, gridB.Comm() ) )
        throw std::logic_error("B does not have a congruent comm");
    const int commRank = elemental::mpi::CommRank( comm_ );

    // Convert B into custom nested-dissection based ordering
    elemental::Matrix<C> B;
    {
        if( commRank == 0 )
        {
            std::cout << "  pulling right-hand sides...";
            std::cout.flush();
        }
        const double startTime = elemental::mpi::Time();

        PullRightHandSides( gridB, B );

        const double stopTime = elemental::mpi::Time();
        if( commRank == 0 )
            std::cout << stopTime-startTime << " secs" << std::endl;
    }

    // Solve the systems of equations
    InternalSolveWithSQMR( B, bcgRelTol );

    // Restore the solutions back into the GridData form
    {
        if( commRank == 0 )
        {
            std::cout << "  pushing right-hand sides...";
            std::cout.flush();
        }
        const double startTime = elemental::mpi::Time();

        PushRightHandSides( gridB, B );

        const double stopTime = elemental::mpi::Time();
        if( commRank == 0 )
            std::cout << stopTime-startTime << " secs" << std::endl;
    }
}

template<typename R>
void
psp::DistHelmholtz<R>::PullRightHandSides
( const GridData<C>& gridB, elemental::Matrix<C>& B ) const
{
    const int commSize = elemental::mpi::CommSize( comm_ );

    // Pack and send the amount of data that we will need to recv
    std::vector<int> recvCounts( commSize, 0 );
    for( int iLocal=0; iLocal<localHeight_; ++iLocal )
    {
        const int naturalIndex = localToNaturalMap_[iLocal];
        const int proc = gridB.OwningProcess( naturalIndex );
        ++recvCounts[proc];
    }
    std::vector<int> sendCounts( commSize );
    elemental::mpi::AllToAll
    ( &recvCounts[0], 1,
      &sendCounts[0], 1, comm_ );

    // Compute the send and recv offsets and total sizes
    int totalSendCount=0, totalRecvCount=0;
    std::vector<int> sendDispls( commSize ), recvDispls( commSize );
    for( int proc=0; proc<commSize; ++proc )
    {
        sendDispls[proc] = totalSendCount;
        recvDispls[proc] = totalRecvCount;
        totalSendCount += sendCounts[proc];
        totalRecvCount += recvCounts[proc];
    }

    // Pack and send the indices that we will need to recv from each process.
    std::vector<int> offsets = recvDispls;
    std::vector<int> recvIndices( totalRecvCount );
    for( int iLocal=0; iLocal<localHeight_; ++iLocal )
    {
        const int naturalIndex = localToNaturalMap_[iLocal];
        const int proc = gridB.OwningProcess( naturalIndex );
        recvIndices[offsets[proc]++] = naturalIndex;
    }
    std::vector<int> sendIndices( totalSendCount );
    elemental::mpi::AllToAll
    ( &recvIndices[0], &recvCounts[0], &recvDispls[0],
      &sendIndices[0], &sendCounts[0], &sendDispls[0], comm_ );
    recvIndices.clear();

    // Scale the send and recv counts to make up for their being several RHS
    const int numRhs = gridB.NumScalars();
    for( int proc=0; proc<commSize; ++proc )
    {
        sendCounts[proc] *= numRhs;
        recvCounts[proc] *= numRhs;
        sendDispls[proc] *= numRhs;
        recvDispls[proc] *= numRhs;
    }
    totalSendCount *= numRhs;
    totalRecvCount *= numRhs;

    // Pack and send our right-hand side data.
    std::vector<C> sendB( totalSendCount );
    const C* gridBBuffer = gridB.LockedLocalBuffer();
    for( int proc=0; proc<commSize; ++proc )
    {
        C* procB = &sendB[sendDispls[proc]];
        const int* procIndices = &sendIndices[sendDispls[proc]/numRhs];
        const int numLocalIndices = sendCounts[proc]/numRhs;
        for( int s=0; s<numLocalIndices; ++s )
        {
            const int naturalIndex = procIndices[s];
            const int localIndex = gridB.LocalIndex( naturalIndex );
            for( int k=0; k<numRhs; ++k )
                procB[s*numRhs+k] = gridBBuffer[localIndex+k];
        }
    }
    sendIndices.clear();
    std::vector<C> recvB( totalRecvCount );
    elemental::mpi::AllToAll
    ( &sendB[0], &sendCounts[0], &sendDispls[0], 
      &recvB[0], &recvCounts[0], &recvDispls[0], comm_ );
    sendB.clear();

    // Unpack the received right-hand side data
    offsets = recvDispls;
    B.ResizeTo( localHeight_, numRhs, localHeight_ );
    for( int iLocal=0; iLocal<localHeight_; ++iLocal )
    {
        const int naturalIndex = localToNaturalMap_[iLocal];
        const int proc = gridB.OwningProcess( naturalIndex );
        for( int k=0; k<numRhs; ++k )
            B.Set( iLocal, k, recvB[offsets[proc]+k] );
        offsets[proc] += numRhs;
    }
}

template<typename R>
void
psp::DistHelmholtz<R>::PushRightHandSides
( GridData<C>& gridB, const elemental::Matrix<C>& B ) const
{
    const int numRhs = gridB.NumScalars();
    const int commSize = elemental::mpi::CommSize( comm_ );

    // Pack and send the amount of data that we will need to send.
    std::vector<int> sendCounts( commSize, 0 );
    for( int iLocal=0; iLocal<localHeight_; ++iLocal )
    {
        const int naturalIndex = localToNaturalMap_[iLocal];
        const int proc = gridB.OwningProcess( naturalIndex );
        ++sendCounts[proc];
    }
    std::vector<int> recvCounts( commSize );
    elemental::mpi::AllToAll
    ( &sendCounts[0], 1, 
      &recvCounts[0], 1, comm_ );

    // Compute the send and recv offsets and total sizes
    int totalSendCount=0, totalRecvCount=0;
    std::vector<int> sendDispls( commSize ), recvDispls( commSize );
    for( int proc=0; proc<commSize; ++proc )
    {
        sendDispls[proc] = totalSendCount;
        recvDispls[proc] = totalRecvCount;
        totalSendCount += sendCounts[proc];
        totalRecvCount += recvCounts[proc];
    }

    // Pack and send the particular indices that we will need to send to 
    // each process.
    std::vector<int> offsets = sendDispls;
    std::vector<int> sendIndices( totalSendCount );
    for( int iLocal=0; iLocal<localHeight_; ++iLocal )
    {
        const int naturalIndex = localToNaturalMap_[iLocal];
        const int proc = gridB.OwningProcess( naturalIndex );
        sendIndices[offsets[proc]++] = naturalIndex;
    }
    std::vector<int> recvIndices( totalRecvCount );
    elemental::mpi::AllToAll
    ( &sendIndices[0], &sendCounts[0], &sendDispls[0], 
      &recvIndices[0], &recvCounts[0], &recvDispls[0], comm_ );
    sendIndices.clear();

    // Scale the counts and offsets by the number of right-hand sides
    totalSendCount *= numRhs;
    totalRecvCount *= numRhs;
    for( int proc=0; proc<commSize; ++proc )
    {
        sendCounts[proc] *= numRhs;
        recvCounts[proc] *= numRhs;
        sendDispls[proc] *= numRhs;
        recvDispls[proc] *= numRhs;
    }

    // Pack and send the right-hand side data
    offsets = sendDispls;
    std::vector<C> sendB( totalSendCount );
    for( int iLocal=0; iLocal<localHeight_; ++iLocal )
    {
        const int naturalIndex = localToNaturalMap_[iLocal]; 
        const int proc = gridB.OwningProcess( naturalIndex );
        for( int k=0; k<numRhs; ++k )
            sendB[offsets[proc]+k] = B.Get( iLocal, k );
        offsets[proc] += numRhs;
    }
    std::vector<C> recvB( totalRecvCount );
    elemental::mpi::AllToAll
    ( &sendB[0], &sendCounts[0], &sendDispls[0],
      &recvB[0], &recvCounts[0], &recvDispls[0], comm_ );
    sendB.clear();

    // Unpack the right-hand side data
    C* gridBBuffer = gridB.LocalBuffer();
    for( int proc=0; proc<commSize; ++proc )
    {
        const C* procB = &recvB[recvDispls[proc]];
        const int* procIndices = &recvIndices[recvDispls[proc]/numRhs];
        const int numLocalIndices = recvCounts[proc]/numRhs;
        for( int s=0; s<numLocalIndices; ++s )
        {
            const int naturalIndex = procIndices[s];
            const int localIndex = gridB.LocalIndex( naturalIndex );
            for( int k=0; k<numRhs; ++k )
                gridBBuffer[localIndex+k] = procB[s*numRhs+k];
        }
    }
}

template<typename R>
void
psp::DistHelmholtz<R>::InternalSolveWithGMRES
( elemental::Matrix<C>& bList, int m, R relTol ) const
{
    const int numRhs = bList.Width();
    const int localHeight = bList.Height();
    const int commRank = elemental::mpi::CommRank( comm_ );

    elemental::Matrix<C> VInter( localHeight, numRhs*m ), // interwoven
                         x0List( localHeight, numRhs   ), // contiguous
                         xList(  localHeight, numRhs   ), // contiguous
                         wList(  localHeight, numRhs   ), // contiguous
                         zList(  m+1,         numRhs   ), // contiguous
                         HList(  m,           m*numRhs ); // contiguous

    // For storing Givens rotations
    elemental::Matrix<R> csList( m, numRhs );
    elemental::Matrix<C> snList( m, numRhs );

    // Various scalars
    std::vector<C> alphaList( numRhs );
    std::vector<R> betaList( numRhs ), deltaList( numRhs );
    std::vector<R> origResidNormList( numRhs ), residNormList( numRhs ), 
                   relResidNormList( numRhs );

    // x := 0
    xList.SetToZero();

    // w := b (= b - A x_0)
    // origResidNorm := ||w||_2
    wList = bList;
    Norms( wList, origResidNormList );
    const bool origResidHasNaN = CheckForNaN( origResidNormList );
    if( origResidHasNaN )
        throw std::runtime_error("Original residual norms had a NaN");

    int it=0;
    bool converged=false;
    while( !converged )
    {
        if( commRank == 0 )
            std::cout << "  starting iteration " << it << "..." << std::endl;

        // x0 := x
        x0List = xList;

        // w := inv(M) w
        // beta := ||w||_2
        {
#ifndef RELEASE
            elemental::mpi::Barrier( comm_ );
#endif
            if( commRank == 0 )
            {
                std::cout << "  startup preconditioner application...";
                std::cout.flush();
            }
            const double startTime = elemental::mpi::Time();
            Precondition( wList );
#ifndef RELEASE
            elemental::mpi::Barrier( comm_ );
#endif
            const double stopTime = elemental::mpi::Time();
            if( commRank == 0 )
                std::cout << stopTime-startTime << " secs" << std::endl;
        }
        Norms( wList, betaList );
        const bool betaListHasNaN = CheckForNaN( betaList );
        if( betaListHasNaN )
            throw std::runtime_error("beta list had a NaN");

        // v0 := w / beta
        elemental::Matrix<C> v0List;
        v0List.View( VInter, 0, 0, localHeight, numRhs );
        v0List = wList;
        DivideColumns( v0List, betaList );

        // z := beta e1
        zList.SetToZero();
        for( int k=0; k<numRhs; ++k )
            zList.Set(0,k,betaList[k]);

        for( int j=0; j<m; ++j )
        {
            // w := A v_j
            elemental::Matrix<C> vjList;
            vjList.LockedView( VInter, 0, j*numRhs, localHeight, numRhs );
            wList = vjList;
            {
#ifndef RELEASE
                elemental::mpi::Barrier( comm_ );
#endif
                if( commRank == 0 )
                {
                    std::cout << "    multiplying...";
                    std::cout.flush();
                }
                const double startTime = elemental::mpi::Time();
                Multiply( wList );
#ifndef RELEASE
                elemental::mpi::Barrier( comm_ );
#endif
                const double stopTime = elemental::mpi::Time();
                if( commRank == 0 )
                    std::cout << stopTime-startTime << " secs" << std::endl;
            }
#ifndef RELEASE
            Norms( wList, deltaList );
            const bool multiplyHasNaN = CheckForNaN( deltaList );
            if( multiplyHasNaN )
                throw std::runtime_error("multiply had a NaN");
#endif

            // w := inv(M) w
            {
                if( commRank == 0 )
                {
                    std::cout << "    preconditioning...";
                    std::cout.flush();
                }
                const double startTime = elemental::mpi::Time();
                Precondition( wList );
#ifndef RELEASE
                elemental::mpi::Barrier( comm_ );
#endif
                const double stopTime = elemental::mpi::Time();
                if( commRank == 0 )
                    std::cout << stopTime-startTime << " secs" << std::endl;
            }
#ifndef RELEASE
            Norms( wList, deltaList );
            const bool preconditionHasNaN = CheckForNaN( deltaList );
            if( preconditionHasNaN )
                throw std::runtime_error("precondition had a NaN");
#endif

            // Run the j'th step of Arnoldi
            {
                if( commRank == 0 )
                {
                    std::cout << "    Arnoldi step...";
                    std::cout.flush();
                }
                const double startTime = elemental::mpi::Time();
                for( int i=0; i<=j; ++i )
                {
                    // H(i,j) := v_i' w
                    elemental::Matrix<C> viList;
                    viList.LockedView
                    ( VInter, 0, i*numRhs, localHeight, numRhs );
                    InnerProducts( viList, wList, alphaList );
                    for( int k=0; k<numRhs; ++k )
                        HList.Set(i,j+k*m,alphaList[k]);

                    // w := w - H(i,j) v_i
                    SubtractScaledColumns( alphaList, viList, wList );
                }
                Norms( wList, deltaList );
                const bool deltaListHasNaN = CheckForNaN( deltaList );
                if( deltaListHasNaN )
                    throw std::runtime_error("delta list had a NaN");
                // TODO: Handle "lucky breakdown" much more carefully
                const bool zeroDelta = CheckForZero( deltaList );
                if( zeroDelta )
                {
                    if( commRank == 0 ) 
                        std::cout 
                            << "GMRES halted due to a (usually) lucky "
                               "breakdown, but this is trickier for multiple "
                               "right-hand sides." << std::endl;
                    return;
                }
                if( j+1 != m )
                {
                    elemental::Matrix<C> vjp1List;
                    vjp1List.View
                    ( VInter, 0, (j+1)*numRhs, localHeight, numRhs );
                    vjp1List = wList;
                    DivideColumns( vjp1List, deltaList );
                }
#ifndef RELEASE
                elemental::mpi::Barrier( comm_ );
#endif
                const double stopTime = elemental::mpi::Time();
                if( commRank == 0 )
                    std::cout << stopTime-startTime << " secs" << std::endl;
            }

            // Apply the previous rotations to the new column of each H
            {
                if( commRank == 0 )
                {
                    std::cout << "    applying previous rotations...";
                    std::cout.flush();
                }
                const double startTime = elemental::mpi::Time();
                for( int k=0; k<numRhs; ++k )
                {
                    elemental::Matrix<C> H;
                    H.View( HList, 0, k*m, j+1, j+1 );
                    for( int i=0; i<j; ++i )
                    {
                        const R c = csList.Get(i,k);
                        const C s = snList.Get(i,k);
                        const C sConj = elemental::Conj(s);
                        const C eta_i_j = H.Get(i,j);
                        const C eta_ip1_j = H.Get(i+1,j);
                        H.Set( i,   j,  c    *eta_i_j + s*eta_ip1_j );
                        H.Set( i+1, j, -sConj*eta_i_j + c*eta_ip1_j );
                    }
                }
#ifndef RELEASE
                elemental::mpi::Barrier( comm_ );
#endif
                const double stopTime = elemental::mpi::Time();
                if( commRank == 0 )
                    std::cout << stopTime-startTime << " secs" << std::endl;
            }

            // Generate the new rotation and apply it to our current column
            // and to z, the rotated beta*e1, then solve for the residual 
            // minimizer
            {
                if( commRank == 0 )
                {
                    std::cout << "    rotating and minimizing residual...";
                    std::cout.flush();
                }
                const double startTime = elemental::mpi::Time();
                for( int k=0; k<numRhs; ++k )
                {
                    // Apply the rotation to the new column of H
                    elemental::Matrix<C> H;
                    H.View( HList, 0, k*m, j+1, j+1 );
                    const C eta_j_j = H.Get(j,j);
                    const C eta_jp1_j = deltaList[k];
                    if( CheckForNaN(eta_j_j) )
                        throw std::runtime_error("H(j,j) was NaN");
                    if( CheckForNaN(eta_jp1_j) )
                        throw std::runtime_error("H(j+1,j) was NaN");
                    R c;
                    C s, rho;
                    elemental::lapack::ComputeGivens
                    ( eta_j_j, eta_jp1_j, &c, &s, &rho );
                    if( CheckForNaN(c) )
                        throw std::runtime_error("c in Givens was NaN");
                    if( CheckForNaN(s) )
                        throw std::runtime_error("s in Givens was NaN");
                    if( CheckForNaN(rho) )
                        throw std::runtime_error("rho in Givens was NaN");
                    H.Set(j,j,rho);
                    csList.Set(j,k,c);
                    snList.Set(j,k,s);

                    // Apply the rotation to z
                    const C sConj = elemental::Conj(s);
                    const C zeta_j = zList.Get(j,k);
                    const C zeta_jp1 = zList.Get(j+1,k);
                    zList.Set( j,   k,  c    *zeta_j + s*zeta_jp1 );
                    zList.Set( j+1, k, -sConj*zeta_j + c*zeta_jp1 );

                    // Minimize the residual
                    elemental::Matrix<C> y, z;
                    z.LockedView( zList, 0, k, j+1, 1 );
                    y = z;
                    elemental::Trsv
                    ( elemental::UPPER, elemental::NORMAL, elemental::NON_UNIT,
                      H, y );

                    // x := x0 + Vj y
                    elemental::Matrix<C> x, x0, vi;
                    x.View(         xList, 0, k, localHeight, 1 );
                    x0.LockedView( x0List, 0, k, localHeight, 1 );
                    x = x0;
                    for( int i=0; i<=j; ++i )
                    {
                        const C eta_i = y.Get(i,0);
                        vi.LockedView( VInter, 0, i*numRhs+k, localHeight, 1 );
                        elemental::Axpy( eta_i, vi, x );
                    }
                }
#ifndef RELEASE
                elemental::mpi::Barrier( comm_ );
#endif
                const double stopTime = elemental::mpi::Time();
                if( commRank == 0 )
                    std::cout << stopTime-startTime << " secs" << std::endl;
            }

            // w := b - A x
            wList = xList; 
            elemental::Scal( (C)-1, wList );
            {
#ifndef RELEASE
                elemental::mpi::Barrier( comm_ );
#endif
                if( commRank == 0 )
                {
                    std::cout << "    residual multiply...";
                    std::cout.flush();
                }
                const double startTime = elemental::mpi::Time();
                Multiply( wList );
#ifndef RELEASE
                elemental::mpi::Barrier( comm_ );
#endif
                const double stopTime = elemental::mpi::Time();
                if( commRank == 0 )
                    std::cout << stopTime-startTime << " secs" << std::endl;
            }
            elemental::Axpy( (C)1, bList, wList );

            // Residual checks
            Norms( wList, residNormList );
            const bool residNormListHasNaN = CheckForNaN( residNormList );
            if( residNormListHasNaN )
                throw std::runtime_error("resid norm list has NaN");
            for( int k=0; k<numRhs; ++k )
                relResidNormList[k] = residNormList[k]/origResidNormList[k];
            R maxRelResidNorm = 0;
            for( int k=0; k<numRhs; ++k )
                maxRelResidNorm = std::max(maxRelResidNorm,relResidNormList[k]);
            if( maxRelResidNorm < relTol )
            {
                if( commRank == 0 )
                    std::cout << "  converged with relative tolerance: " 
                              << maxRelResidNorm << std::endl;
                converged = true;
                break;
            }
            else
            {
                if( commRank == 0 )
                    std::cout << "  finished iteration " << it << " with "
                              << "maxRelResidNorm=" << maxRelResidNorm 
                              << std::endl;
            }
            ++it;
        }
    }
    bList = xList;
}

// Based on several different papers from R. Freund et al. on preconditioned QMR
// for complex symmetric matrices.
template<typename R>
void
psp::DistHelmholtz<R>::InternalSolveWithSQMR
( elemental::Matrix<C>& bList, R bcgRelTol ) const
{
    const R one = 1;
    const int numRhs = bList.Width();
    const int localHeight = bList.Height();
    const int commRank = elemental::mpi::CommRank( comm_ );

    elemental::Matrix<C> vList( localHeight, numRhs ),
                         tList( localHeight, numRhs ),
                         qList( localHeight, numRhs ),
                         pList( localHeight, numRhs );
    std::vector<C> cList( numRhs ),
                   tauList( numRhs ),
                   rhoList( numRhs ),
                   betaList( numRhs ),
                   thetaList( numRhs ),
                   alphaList( numRhs ),
                   sigmaList( numRhs ),
                   rhoLastList( numRhs ),
                   thetaLastList( numRhs ),
                   tempList( numRhs );
    std::vector<R> origResidNormList( numRhs ), 
                   bcgResidNormList( numRhs ),
                   relBcgResidNormList( numRhs );

    // v := b
    // origResidNorm := ||v||_2
    vList = bList;
    Norms( vList, origResidNormList );
    const bool origResidHasNaN = CheckForNaN( origResidNormList );
    if( origResidHasNaN )
        throw std::runtime_error("Original resid has NaN");

    // t := inv(M) v 
    // tau := sqrt(t^T t)
    tList = vList;
    {
#ifndef RELEASE
        elemental::mpi::Barrier( comm_ );
#endif
        if( commRank == 0 )
        {
            std::cout << "  initial preconditioner application...";
            std::cout.flush();
        }
        const double startTime = elemental::mpi::Time();
        Precondition( tList );
#ifndef RELEASE
        elemental::mpi::Barrier( comm_ );
#endif
        const double stopTime = elemental::mpi::Time();
        if( commRank == 0 )
            std::cout << stopTime-startTime << " secs" << std::endl;
    }
    PseudoNorms( tList, tauList );
    const bool tauListHasNaN = CheckForNaN( tauList );
    if( tauListHasNaN )
        throw std::runtime_error("tau list has NaN");

    // q := t
    // x := 0 (use the B matrix for storing the x vectors)
    // p := 0
    // theta := 0
    // rho := v^T q
    qList = tList;
    bList.SetToZero();
    pList.SetToZero();
    for( int k=0; k<numRhs; ++k )
        thetaList[k] = 0;
    PseudoInnerProducts( vList, qList, rhoList );
    const bool rhoListHasNaN = CheckForNaN( rhoList );
    if( rhoListHasNaN )
        throw std::runtime_error("rho list has NaN");

    int it=0;
    while( true )
    {
        if( commRank == 0 )
            std::cout << "  starting iteration " << it << "..." << std::endl;
        // t := A q
        tList = qList;
        {
#ifndef RELEASE
            elemental::mpi::Barrier( comm_ );
#endif
            if( commRank == 0 )
            {
                std::cout << "  multiplying...";
                std::cout.flush();
            }
            const double startTime = elemental::mpi::Time();
            Multiply( tList );
#ifndef RELEASE
            elemental::mpi::Barrier( comm_ );
#endif
            const double stopTime = elemental::mpi::Time();
            if( commRank == 0 )
                std::cout << stopTime-startTime << " secs" << std::endl;
        }

        // sigma := q^T t
        // alpha := rho / sigma
        // v := v - alpha t
        PseudoInnerProducts( qList, tList, sigmaList );
        const bool sigmaListHasNaN = CheckForNaN( sigmaList );
        if( sigmaListHasNaN )
            throw std::runtime_error("sigma list has NaN");
        const bool zeroSigma = CheckForZero( sigmaList );
        if( zeroSigma )
        {
            if( commRank == 0 ) 
                std::cout << "SQMR stopped due to a zero sigma" << std::endl;
            break;
        }
        for( int k=0; k<numRhs; ++k )
            alphaList[k] = rhoList[k] / sigmaList[k];
        const bool alphaListHasNaN = CheckForNaN( alphaList );
        if( alphaListHasNaN )
            throw std::runtime_error("alpha list has NaN");
        SubtractScaledColumns( alphaList, tList, vList );

        // t         := inv(M) v
        // thetaLast := theta
        // theta     := sqrt(t^T t) / tau
        // c         := 1 / sqrt(1+theta^2)
        // tau       := tau theta c
        // p         := c^2 thetaLast^2 p + c^2 alpha q
        // x         := x + p
        tList = vList;
        {
#ifndef RELEASE
            elemental::mpi::Barrier( comm_ );
#endif
            if( commRank == 0 )
            {
                std::cout << "  preconditioning...";
                std::cout.flush();
            }
            const double startTime = elemental::mpi::Time();
            Precondition( tList );
#ifndef RELEASE
            elemental::mpi::Barrier( comm_ );
#endif
            const double stopTime = elemental::mpi::Time();
            if( commRank == 0 )
                std::cout << stopTime-startTime << " secs" << std::endl;
        }
        thetaLastList = thetaList;
        PseudoNorms( tList, thetaList );
        for( int k=0; k<numRhs; ++k )
        {
            thetaList[k] /= tauList[k];
            cList[k]      = one/sqrt(one+thetaList[k]*thetaList[k]);
            tauList[k]    = tauList[k]*thetaList[k]*cList[k];
        }
        for( int k=0; k<numRhs; ++k ) 
            tempList[k] = cList[k]*cList[k]*thetaLastList[k]*thetaLastList[k];
        MultiplyColumns( pList, tempList );
        for( int k=0; k<numRhs; ++k )
            tempList[k] = cList[k]*cList[k]*alphaList[k];
        AddScaledColumns( tempList, qList, pList );
        elemental::Axpy( (C)1, pList, bList );
        const bool thetaListHasNaN = CheckForNaN( thetaList );
        if( thetaListHasNaN )
            throw std::runtime_error("theta list has NaN");
        const bool zeroTheta = CheckForZero( thetaList );
        if( zeroTheta )
        {
            if( commRank == 0 )
                std::cout << "SQMR stopped due to a zero theta" << std::endl;
            break;
        }

        // Residual checks
        Norms( vList, bcgResidNormList );
        const bool bcgResidHasNaN = CheckForNaN( bcgResidNormList );
        if( bcgResidHasNaN )
            throw std::runtime_error("BCG residuals have NaN");
        for( int k=0; k<numRhs; ++k )
            relBcgResidNormList[k] = 
                bcgResidNormList[k]/origResidNormList[k];
        R maxRelBcgResidNorm = 0;
        for( int k=0; k<numRhs; ++k )
            maxRelBcgResidNorm = 
                std::max(maxRelBcgResidNorm,relBcgResidNormList[k]);
        if( maxRelBcgResidNorm < bcgRelTol )
        {
            if( commRank == 0 )
                std::cout << "  converged with BCG relative tolerance: " 
                          << maxRelBcgResidNorm << std::endl;
            break;
        }
        else
        {
            if( commRank == 0 )
                std::cout << "  finished iteration " << it << " with "
                          << "maxRelBcgResidNorm=" << maxRelBcgResidNorm 
                          << std::endl;
        }

        // rhoLast := rho
        // rho     := v^T t
        // beta    := rho / rhoLast
        // q       := t + beta q
        const bool zeroRho = CheckForZero( rhoList );
        if( zeroRho )
        {
            if( commRank == 0 ) 
                std::cout << "SQMR stopped due to a zero rho" << std::endl;
            break;
        }
        rhoLastList = rhoList;
        PseudoInnerProducts( vList, tList, rhoList );
        for( int k=0; k<numRhs; ++k )
            betaList[k] = rhoList[k] / rhoLastList[k];
        MultiplyColumns( qList, betaList );
        elemental::Axpy( (C)1, tList, qList );

        ++it;
    }
}

template<typename R>
bool
psp::DistHelmholtz<R>::CheckForNaN( R alpha ) const
{
    return alpha != alpha; // hopefully this is not optimized away
}

template<typename R>
bool
psp::DistHelmholtz<R>::CheckForNaN( C alpha ) const
{
    return alpha != alpha; // hopefully this is not optimized away
}

template<typename R>
bool
psp::DistHelmholtz<R>::CheckForNaN( const std::vector<R>& alphaList ) const
{
    bool foundNaN = false;
    for( unsigned k=0; k<alphaList.size(); ++k )
        if( CheckForNaN(alphaList[k]) )
            foundNaN = true;
    return foundNaN;
}

template<typename R>
bool
psp::DistHelmholtz<R>::CheckForNaN( const std::vector<C>& alphaList ) const
{
    bool foundNaN = false;
    for( unsigned k=0; k<alphaList.size(); ++k )
        if( CheckForNaN(alphaList[k]) )
            foundNaN = true;
    return foundNaN;
}

template<typename R>
bool
psp::DistHelmholtz<R>::CheckForZero( const std::vector<R>& alphaList ) const
{
    bool foundZero = false;
    for( unsigned k=0; k<alphaList.size(); ++k )
        if( alphaList[k] == (R)0 ) // think about using a tolerance instead
            foundZero = true;
    return foundZero;
}

template<typename R>
bool
psp::DistHelmholtz<R>::CheckForZero( const std::vector<C>& alphaList ) const
{
    bool foundZero = false;
    for( unsigned k=0; k<alphaList.size(); ++k )
        if( alphaList[k] == (C)0 ) // think about using a tolerance instead
            foundZero = true;
    return foundZero;
}

template<typename R>
void
psp::DistHelmholtz<R>::Norms
( const elemental::Matrix<C>& xList, std::vector<R>& normList ) const
{
    const int numCols = xList.Width();
    const int localHeight = xList.Height();
    const int commSize = elemental::mpi::CommSize( comm_ );
    std::vector<R> localNorms( numCols );
    for( int j=0; j<numCols; ++j )
        localNorms[j] = 
            elemental::blas::Nrm2( localHeight, xList.LockedBuffer(0,j), 1 );
    std::vector<R> allLocalNorms( numCols*commSize );
    elemental::mpi::AllGather
    ( &localNorms[0], numCols, &allLocalNorms[0], numCols, comm_ );
    normList.resize( numCols );
    for( int j=0; j<numCols; ++j )
        normList[j] = 
            elemental::blas::Nrm2( commSize, &allLocalNorms[j], numCols );
}

template<typename R>
void
psp::DistHelmholtz<R>::PseudoNorms
( const elemental::Matrix<C>& xList, std::vector<C>& alphaList ) const
{
    const int numCols = xList.Width();
    const int localHeight = xList.Height();
    std::vector<C> localAlphaList( numCols );
    for( int j=0; j<numCols; ++j )
        localAlphaList[j] = 
            elemental::blas::Dotu
            ( localHeight, xList.LockedBuffer(0,j), 1,
                           xList.LockedBuffer(0,j), 1 );
    alphaList.resize( numCols );
    // TODO: Think about avoiding overflow?
    elemental::mpi::AllReduce
    ( &localAlphaList[0], &alphaList[0], numCols, MPI_SUM, comm_ );
    for( int j=0; j<numCols; ++j )
        alphaList[j] = sqrt(alphaList[j]);
}

template<typename R>
void
psp::DistHelmholtz<R>::InnerProducts
( const elemental::Matrix<C>& xList, const elemental::Matrix<C>& yList,
  std::vector<C>& alphaList ) const
{
    const int numCols = xList.Width();
    const int localHeight = xList.Height();
    std::vector<C> localAlphaList( numCols );
    for( int j=0; j<numCols; ++j )
        localAlphaList[j] = 
            elemental::blas::Dot
            ( localHeight, xList.LockedBuffer(0,j), 1,
                           yList.LockedBuffer(0,j), 1 );
    alphaList.resize( numCols );
    elemental::mpi::AllReduce
    ( &localAlphaList[0], &alphaList[0], numCols, MPI_SUM, comm_ );
}

template<typename R>
void
psp::DistHelmholtz<R>::PseudoInnerProducts
( const elemental::Matrix<C>& xList, const elemental::Matrix<C>& yList,
  std::vector<C>& alphaList ) const
{
    const int numCols = xList.Width();
    const int localHeight = xList.Height();
    std::vector<C> localAlphaList( numCols );
    for( int j=0; j<numCols; ++j )
        localAlphaList[j] = 
            elemental::blas::Dotu
            ( localHeight, xList.LockedBuffer(0,j), 1,
                           yList.LockedBuffer(0,j), 1 );
    alphaList.resize( numCols );
    elemental::mpi::AllReduce
    ( &localAlphaList[0], &alphaList[0], numCols, MPI_SUM, comm_ );
}

template<typename R>
void
psp::DistHelmholtz<R>::DivideColumns
( elemental::Matrix<C>& xList, const std::vector<R>& deltaList ) const
{
    const R one = 1;
    const int numCols = xList.Width();
    const int localHeight = xList.Height();
    for( int j=0; j<numCols; ++j )
    {
        const R invDelta = one/deltaList[j];
        C* x = xList.Buffer(0,j);
        for( int iLocal=0; iLocal<localHeight; ++iLocal )
            x[iLocal] *= invDelta;
    }
}

template<typename R>
void
psp::DistHelmholtz<R>::MultiplyColumns
( elemental::Matrix<C>& xList, const std::vector<C>& deltaList ) const
{
    const int numCols = xList.Width();
    const int localHeight = xList.Height();
    for( int j=0;j<numCols; ++j )
    {
        const C delta = deltaList[j];
        C* x = xList.Buffer(0,j);
        for( int iLocal=0; iLocal<localHeight; ++iLocal )
            x[iLocal] *= delta;
    }
}

template<typename R>
void
psp::DistHelmholtz<R>::AddScaledColumns
( const std::vector<C>& deltaList, 
  const elemental::Matrix<C>& xList, elemental::Matrix<C>& yList ) const
{
    const int numCols = xList.Width();
    const int localHeight = xList.Height();
    for( int j=0; j<numCols; ++j )
    {
        const C delta = deltaList[j];
        const C* x = xList.LockedBuffer(0,j);
        C* y = yList.Buffer(0,j);
        for( int iLocal=0; iLocal<localHeight; ++iLocal )
            y[iLocal] += delta*x[iLocal];
    }
}

template<typename R>
void
psp::DistHelmholtz<R>::SubtractScaledColumns
( const std::vector<C>& deltaList, 
  const elemental::Matrix<C>& xList, elemental::Matrix<C>& yList ) const
{
    const int numCols = xList.Width();
    const int localHeight = xList.Height();
    for( int j=0; j<numCols; ++j )
    {
        const C delta = deltaList[j];
        const C* x = xList.LockedBuffer(0,j);
        C* y = yList.Buffer(0,j);
        for( int iLocal=0; iLocal<localHeight; ++iLocal )
            y[iLocal] -= delta*x[iLocal];
    }
}

// B := A B
template<typename R>
void
psp::DistHelmholtz<R>::Multiply( elemental::Matrix<C>& B ) const
{
    const int numRhs = B.Width();
    const int commSize = elemental::mpi::CommSize( comm_ );

    // Modify the basic send/recv information for the number of right-hand sides
    std::vector<int> sendCounts = globalSendCounts_;
    std::vector<int> sendDispls = globalSendDispls_;
    std::vector<int> recvCounts = globalRecvCounts_;
    std::vector<int> recvDispls = globalRecvDispls_;
    for( int proc=0; proc<commSize; ++proc )
    {
        sendCounts[proc] *= numRhs;
        sendDispls[proc] *= numRhs;
        recvCounts[proc] *= numRhs;
        recvDispls[proc] *= numRhs;
    }
    const int totalSendCount = sendDispls.back() + sendCounts.back();
    const int totalRecvCount = recvDispls.back() + recvCounts.back();

    // Pack and scatter our portion of the right-hand sides
    std::vector<C> sendRhs( totalSendCount );
    for( int proc=0; proc<commSize; ++proc )
    {
        const int sendSize = globalSendCounts_[proc];
        C* procRhs = &sendRhs[sendDispls[proc]];
        const int* procIndices = &globalSendIndices_[globalSendDispls_[proc]];
        for( int s=0; s<sendSize; ++s )
        {
            const int iLocal = procIndices[s];
            for( int k=0; k<numRhs; ++k )
                procRhs[s*numRhs+k] = B.Get( iLocal, k );
        }
    }
    std::vector<C> recvRhs( totalRecvCount );
    elemental::mpi::AllToAll
    ( &sendRhs[0], &sendCounts[0], &sendDispls[0], 
      &recvRhs[0], &recvCounts[0], &recvDispls[0], comm_ );
    sendRhs.clear();

    // Run the local multiplies to form the result
    std::vector<int> offsets = recvDispls;
    C* BBuffer = B.Buffer();
    for( int iLocal=0; iLocal<localHeight_; ++iLocal )
    {
        // Multiply by the diagonal value
        const int rowOffset = localRowOffsets_[iLocal];
        const C diagVal = localEntries_[rowOffset];
        for( int k=0; k<numRhs; ++k )
            BBuffer[iLocal+k*localHeight_] *= diagVal;

        // Multiply by the off-diagonal values
        const int rowSize = localRowOffsets_[iLocal+1]-rowOffset;
        for( int jLocal=1; jLocal<rowSize; ++jLocal )
        {
            const int proc = owningProcesses_[rowOffset+jLocal];
            const C offDiagVal = localEntries_[rowOffset+jLocal];
            for( int k=0; k<numRhs; ++k )
                BBuffer[iLocal+k*localHeight_] +=
                    offDiagVal*recvRhs[offsets[proc]+k];
            offsets[proc] += numRhs;
        }
    }
}

template<typename R>
void
psp::DistHelmholtz<R>::Precondition( elemental::Matrix<C>& B ) const
{
    // Apply the sweeping preconditioner
    //
    // Simple algorithm:
    //   // Solve against L
    //   for i=0,...,m-2
    //     B_{i+1} := B_{i+1} - A_{i+1,i} T_i B_i
    //   end
    //   // Solve against D 
    //   for i=0,...,m-1
    //     B_i := T_i B_i
    //   end
    //   // Solve against L^T
    //   for i=m-2,...,0
    //     B_i := B_i - T_i A_{i,i+1} B_{i+1}
    //   end
    //
    // Practical algorithm:
    //   // Solve against L D
    //   for i=0,...,m-2
    //     B_i := T_i B_i
    //     B_{i+1} := B_{i+1} - A_{i+1,i} B_i
    //   end
    //   B_{m-1} := T_{m-1} B_{m-1}
    //   // Solve against L^T
    //   for i=m-2,...,0
    //     Z := B_i
    //     B_i := -A_{i,i+1} B_{i+1}
    //     B_i := T_i B_i
    //     B_i := B_i + Z
    //   end
    //

    // Solve against L D
    for( int i=0; i<numPanels_-1; ++i )
    {
        SolvePanel( B, i );
        SubdiagonalUpdate( B, i );
    }
    SolvePanel( B, numPanels_-1 );

    // Solve against L^T
    elemental::Matrix<C> Z;
    for( int i=numPanels_-2; i>=0; --i )
    {
        ExtractPanel( B, i, Z );
        MultiplySuperdiagonal( B, i );    
        SolvePanel( B, i );
        UpdatePanel( B, i, Z );
    }
}

// B_i := T_i B_i
template<typename R>
void
psp::DistHelmholtz<R>::SolvePanel( elemental::Matrix<C>& B, int i ) const
{
    const clique::symbolic::SymmFact& symbFact = 
        PanelSymbolicFactorization( i );
    const int numRhs = B.Width();
    const int panelPadding = PanelPadding( i );
    const int panelDepth = PanelDepth( i );
    const int localHeight1d = 
        symbFact.dist.supernodes.back().localOffset1d + 
        symbFact.dist.supernodes.back().localSize1d;

    elemental::Matrix<C> localPanelB( localHeight1d, numRhs );
    localPanelB.SetToZero();

    // For each supernode, pull in each right-hand side with a memcpy
    int BOffset = LocalPanelOffset( i );
    const int numLocalSupernodes = symbFact.local.supernodes.size();
    for( int t=0; t<numLocalSupernodes; ++t )
    {
        const clique::symbolic::LocalSymmFactSupernode& sn = 
            symbFact.local.supernodes[t];
        const int size = sn.size;
        const int myOffset = sn.myOffset;

#ifndef RELEASE
        if( size % (panelPadding+panelDepth) != 0 )
            throw std::logic_error("Local supernode size problem");
#endif
        const int xySize = size/(panelPadding+panelDepth);
        const int paddingSize = xySize*panelPadding;
        const int remainingSize = xySize*panelDepth;

        for( int k=0; k<numRhs; ++k )
            std::memcpy
            ( localPanelB.Buffer(myOffset+paddingSize,k), 
              B.LockedBuffer(BOffset,k),
              remainingSize*sizeof(C) );
        BOffset += remainingSize;
    }
    const int numDistSupernodes = symbFact.dist.supernodes.size();
    for( int t=1; t<numDistSupernodes; ++t )
    {
        const clique::symbolic::DistSymmFactSupernode& sn = 
            symbFact.dist.supernodes[t];
        const int size = sn.size;
        const int localOffset1d = sn.localOffset1d;
        const int localSize1d = sn.localSize1d;

        const elemental::Grid& grid = *sn.grid;
        const int gridSize = grid.Size();
        const int gridRank = grid.VCRank();

#ifndef RELEASE
        if( size % (panelPadding+panelDepth) != 0 )
            throw std::logic_error("Dist supernode size problem");
#endif
        const int xySize = size/(panelPadding+panelDepth);
        const int paddingSize = xySize*panelPadding;
        const int localPaddingSize = 
            elemental::LocalLength( paddingSize, gridRank, gridSize );
        const int localRemainingSize = localSize1d - localPaddingSize;

        for( int k=0; k<numRhs; ++k )
            std::memcpy
            ( localPanelB.Buffer(localOffset1d+localPaddingSize,k),
              B.LockedBuffer(BOffset,k),
              localRemainingSize*sizeof(C) );
        BOffset += localRemainingSize;
    }
#ifndef RELEASE
    if( BOffset != LocalPanelOffset(i)+LocalPanelHeight(i) )
        throw std::logic_error("Invalid BOffset usage in pull");
#endif

    // Solve against the panel
    const clique::numeric::SymmFrontTree<C>& fact = 
        PanelNumericFactorization( i );
    clique::numeric::LDLSolve
    ( elemental::TRANSPOSE, symbFact, fact, localPanelB );

    // For each supernode, extract each right-hand side with memcpy
    BOffset = LocalPanelOffset( i );
    for( int t=0; t<numLocalSupernodes; ++t )
    {
        const clique::symbolic::LocalSymmFactSupernode& sn = 
            symbFact.local.supernodes[t];
        const int size = sn.size;
        const int myOffset = sn.myOffset;

#ifndef RELEASE
        if( size % (panelPadding+panelDepth) != 0 )
            throw std::logic_error("Local supernode size problem");
#endif
        const int xySize = size/(panelPadding+panelDepth);
        const int paddingSize = xySize*panelPadding;
        const int remainingSize = size - paddingSize;

        for( int k=0; k<numRhs; ++k )
            std::memcpy
            ( B.Buffer(BOffset,k),
              localPanelB.LockedBuffer(myOffset+paddingSize,k), 
              remainingSize*sizeof(C) );
        BOffset += remainingSize;
    }
    for( int t=1; t<numDistSupernodes; ++t )
    {
        const clique::symbolic::DistSymmFactSupernode& sn = 
            symbFact.dist.supernodes[t];
        const int size = sn.size;
        const int localOffset1d = sn.localOffset1d;
        const int localSize1d = sn.localSize1d;

        const elemental::Grid& grid = *sn.grid;
        const int gridSize = grid.Size();
        const int gridRank = grid.VCRank();

#ifndef RELEASE
        if( size % (panelPadding+panelDepth) != 0 )
            throw std::logic_error("Dist supernode size problem");
#endif
        const int xySize = size/(panelPadding+panelDepth);
        const int paddingSize = xySize*panelPadding;
        const int localPaddingSize = 
            elemental::LocalLength( paddingSize, gridRank, gridSize );
        const int localRemainingSize = localSize1d - localPaddingSize;

        for( int k=0; k<numRhs; ++k )
            std::memcpy
            ( B.Buffer(BOffset,k),
              localPanelB.LockedBuffer(localOffset1d+localPaddingSize,k),
              localRemainingSize*sizeof(C) );
        BOffset += localRemainingSize;
    }
#ifndef RELEASE
    if( BOffset != LocalPanelOffset(i)+LocalPanelHeight(i) )
        throw std::logic_error("Invalid BOffset usage in push");
#endif
}

// B_{i+1} := B_{i+1} - A_{i+1,i} B_i
template<typename R>
void
psp::DistHelmholtz<R>::SubdiagonalUpdate( elemental::Matrix<C>& B, int i ) const
{
    const int commSize = elemental::mpi::CommSize( comm_ );
    const int numRhs = B.Width();
    const int panelSendCount = subdiagPanelSendCounts_[i];
    const int panelRecvCount = subdiagPanelRecvCounts_[i];
    const int panelSendOffset = subdiagPanelSendDispls_[i];
    const int panelRecvOffset = subdiagPanelRecvDispls_[i];

    // Pack and alltoall the local d.o.f. at the back of B_i
    std::vector<C> sendBuffer( panelSendCount*numRhs );
    for( int s=0; s<panelSendCount; ++s )
    {
        const int iLocal = subdiagSendIndices_[panelSendOffset+s];
        for( int k=0; k<numRhs; ++k ) 
            sendBuffer[s*numRhs+k] = B.Get( iLocal, k );
    }

    // Prepare the send and recv information
    std::vector<int> sendCounts( commSize ), sendDispls( commSize ),
                     recvCounts( commSize ), recvDispls( commSize );
    for( int proc=0; proc<commSize; ++proc )
    {
        const int index = i*commSize + proc;
        sendCounts[proc] = subdiagSendCounts_[index]*numRhs;
        sendDispls[proc] = subdiagSendDispls_[index]*numRhs;
        recvCounts[proc] = subdiagRecvCounts_[index]*numRhs;
        recvDispls[proc] = subdiagRecvDispls_[index]*numRhs;
    }

    std::vector<C> recvBuffer( panelRecvCount*numRhs );
    elemental::mpi::AllToAll
    ( &sendBuffer[0], &sendCounts[0], &sendDispls[0],
      &recvBuffer[0], &recvCounts[0], &recvDispls[0], comm_ );
    sendBuffer.clear();
    sendCounts.clear();
    sendDispls.clear();

    // Perform the local update
    const int* recvLocalRows = &subdiagRecvLocalRows_[panelRecvOffset];
    const int* recvLocalIndices = &subdiagRecvLocalIndices_[panelRecvOffset];
    for( int proc=0; proc<commSize; ++proc )
    {
        const int procSize = recvCounts[proc]/numRhs;
        const int procOffset = recvDispls[proc]/numRhs;

        const C* procValues = &recvBuffer[recvDispls[proc]];
        const int* procLocalRows = &recvLocalRows[procOffset];
        const int* procLocalIndices = &recvLocalIndices[procOffset];
        for( int s=0; s<procSize; ++s )
        {
            const int iLocal = procLocalRows[s];
            const int localIndex = procLocalIndices[s];
            const C alpha = localEntries_[localIndex];
            for( int k=0; k<numRhs; ++k )
            {
                const C beta = procValues[s*numRhs+k];
                B.Update( iLocal, k, -alpha*beta );
            }
        }
    }
}

// Z := B_i
// B_i := 0
template<typename R>
void
psp::DistHelmholtz<R>::ExtractPanel
( elemental::Matrix<C>& B, int i, elemental::Matrix<C>& Z ) const
{
    const int localPanelOffset = LocalPanelOffset( i );
    const int localPanelHeight = LocalPanelHeight( i );
    const int numRhs = B.Width();
    Z.ResizeTo( localPanelHeight, numRhs );

    for( int k=0; k<numRhs; ++k )
    {
        std::memcpy
        ( Z.Buffer(0,k), B.LockedBuffer(localPanelOffset,k),
          localPanelHeight*sizeof(C) );
        std::memset
        ( B.Buffer(localPanelOffset,k), 0, localPanelHeight*sizeof(C) );
    }
}

// B_i := -A_{i,i+1} B_{i+1}
template<typename R>
void
psp::DistHelmholtz<R>::MultiplySuperdiagonal
( elemental::Matrix<C>& B, int i ) const
{
    const int commSize = elemental::mpi::CommSize( comm_ );
    const int numRhs = B.Width();
    const int panelSendCount = supdiagPanelSendCounts_[i];
    const int panelRecvCount = supdiagPanelRecvCounts_[i];
    const int panelSendOffset = supdiagPanelSendDispls_[i];
    const int panelRecvOffset = supdiagPanelRecvDispls_[i];

    // Pack and alltoall the local d.o.f. at the front of B_{i+1}
    std::vector<C> sendBuffer( panelSendCount*numRhs );
    for( int s=0; s<panelSendCount; ++s )
    {
        const int iLocal = supdiagSendIndices_[panelSendOffset+s];
#ifndef RELEASE
        if( iLocal < LocalPanelOffset( i+1 ) )
        {
            std::cout << "s=" << s << "\n"
                      << "offset i+1=" << LocalPanelOffset(i+1) << ", \n"
                      << "iLocal=" << iLocal << std::endl;
            throw std::logic_error("Send index was too small");
        }
        if( iLocal >= LocalPanelOffset(i+1)+LocalPanelHeight(i+1) )
        {
            std::cout << "s=" << s << "\n"
                      << "offset i+1=" << LocalPanelOffset(i+1) << ", \n"
                      << "height i+1=" << LocalPanelHeight(i+1) << ", \n"
                      << "iLocal    =" << iLocal << std::endl;
            throw std::logic_error("Send index was too big");
        }
#endif
        for( int k=0; k<numRhs; ++k ) 
            sendBuffer[s*numRhs+k] = B.Get( iLocal, k );
    }
    std::vector<int> sendCounts( commSize ), sendDispls( commSize ),
                     recvCounts( commSize ), recvDispls( commSize );
    for( int proc=0; proc<commSize; ++proc )
    {
        const int index = i*commSize + proc;
        sendCounts[proc] = supdiagSendCounts_[index]*numRhs;
        sendDispls[proc] = supdiagSendDispls_[index]*numRhs;
        recvCounts[proc] = supdiagRecvCounts_[index]*numRhs;
        recvDispls[proc] = supdiagRecvDispls_[index]*numRhs;
    }
    std::vector<C> recvBuffer( panelRecvCount*numRhs );
    elemental::mpi::AllToAll
    ( &sendBuffer[0], &sendCounts[0], &sendDispls[0],
      &recvBuffer[0], &recvCounts[0], &recvDispls[0], comm_ );
    sendBuffer.clear();
    sendCounts.clear();
    sendDispls.clear();

    // Perform the local multiply
    const int* recvLocalRows = &supdiagRecvLocalRows_[panelRecvOffset];
    const int* recvLocalIndices = &supdiagRecvLocalIndices_[panelRecvOffset];
    for( int proc=0; proc<commSize; ++proc )
    {
        const int procSize = recvCounts[proc]/numRhs;
        const int procOffset = recvDispls[proc]/numRhs;

        const C* procValues = &recvBuffer[recvDispls[proc]];
        const int* procLocalRows = &recvLocalRows[procOffset];
        const int* procLocalIndices = &recvLocalIndices[procOffset];
        for( int s=0; s<procSize; ++s )
        {
            const int iLocal = procLocalRows[s];
            const int localIndex = procLocalIndices[s];
            const C alpha = localEntries_[localIndex];
            for( int k=0; k<numRhs; ++k )
            {
                const C beta = procValues[s*numRhs+k];
                B.Set( iLocal, k, -alpha*beta );
            }
        }
    }
}

// B_i := B_i + Z
template<typename R>
void
psp::DistHelmholtz<R>::UpdatePanel
( elemental::Matrix<C>& B, int i, const elemental::Matrix<C>& Z ) const
{
    const int localPanelOffset = LocalPanelOffset( i );
    const int localPanelHeight = LocalPanelHeight( i );
    const int numRhs = Z.Width();
    for( int k=0; k<numRhs; ++k )
        for( int s=0; s<localPanelHeight; ++s )
            B.Update( localPanelOffset+s, k, Z.Get(s,k) );
}

