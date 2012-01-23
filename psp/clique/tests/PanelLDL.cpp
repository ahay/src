/*
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

void Usage()
{
    std::cout << "PanelLDL <numPanels> <nx> <ny> <nz> <cutoff> <fact blocksize>"
                 " <solve blocksize> <write info?>\n"
              << "<numPanels>: number of panels to factor in memory\n"
              << "<nx>: size of panel in x direction\n"
              << "<ny>: size of panel in y direction\n"
              << "<nz>: size of panel in z direction\n"
              << "<cutoff>: minimum required leaf size\n" 
              << "<fact blocksize>: algorithmic blocksize for factorization\n"
              << "<solve blocksize>: algorithm blocksize for solve\n"
              << "<write info?>: write basic local info to file? [0/1]\n"
              << "<prefix>: prefix for each process's info file\n"
              << std::endl;
}

int ReorderedIndexRecursion
( int x, int y, int z, int nx, int ny, int nz,
  int stepsLeft, int cutoff, int offset );

int ReorderedIndex
( int x, int y, int z, int nx, int ny, int nz,
  int minBalancedDepth, int cutoff );

void CountLocalTreeSize
( int nxSub, int nySub, int nz, int cutoff, int& numSupernodes );

void FillOrigStruct
( int nx, int ny, int nz, int cutoff, clique::mpi::Comm comm, int log2CommSize,
  clique::symbolic::SymmOrig& S );

void FillDistOrigStruct
( int nx, int ny, int nz, int& nxSub, int& nySub, int& xOffset, int& yOffset, 
  int cutoff, clique::mpi::Comm comm, int log2CommSize,
  clique::symbolic::SymmOrig& S );

void FillLocalOrigStruct
( int nx, int ny, int nz, int nxSub, int nySub, int xOffset, int yOffset, 
  int cutoff, int log2CommSize,
  clique::symbolic::SymmOrig& S );

int
main( int argc, char* argv[] )
{
    clique::Initialize( argc, argv );
    clique::mpi::Comm comm = clique::mpi::COMM_WORLD;
    const int commRank = clique::mpi::CommRank( comm );
    const unsigned commSize = clique::mpi::CommSize( comm );
    typedef std::complex<double> F;

    // Ensure that we have a power of two number of processes
    unsigned temp = commSize;
    unsigned log2CommSize = 0;
    while( temp >>= 1 )
        ++log2CommSize;
    if( 1u<<log2CommSize != commSize )
    {
        if( commRank == 0 )
        {
            std::cerr << "Must use a power of two number of processes" 
                      << std::endl;
        }
        clique::Finalize();
        return 0;
    }

    if( argc < 9 )
    {
        if( commRank == 0 )        
            Usage();
        clique::Finalize();
        return 0;
    }

    int argNum = 1;
    const int numPanels = atoi( argv[argNum++] );
    const int nx = atoi( argv[argNum++] );
    const int ny = atoi( argv[argNum++] );
    const int nz = atoi( argv[argNum++] );
    const int cutoff = atoi( argv[argNum++] );
    const int factBlocksize = atoi( argv[argNum++] );
    const int solveBlocksize = atoi( argv[argNum++] );
    const bool writeInfo = atoi( argv[argNum++] );
    if( writeInfo && argc == 10 )
    {
        if( commRank == 0 )
            Usage();
        clique::Finalize();
        return 0;
    }
    const char* basename = argv[argNum++];
    if( commRank == 0 )
        std::cout << "numPanels=" << numPanels << "\n"
                  << "(nx,ny,nz)=(" << nx << "," << ny << "," << nz << ")\n"
                  << "cutoff=" << cutoff << "\n"
                  << "factBlocksize=" << factBlocksize << "\n"
                  << "solveBlocksize=" << solveBlocksize << std::endl;

    try
    {
        std::ofstream infoFile;
        if( writeInfo )
        {
            std::ostringstream filename;
            filename << basename << "-" << commRank << ".dat";
            infoFile.open( filename.str().c_str() );
        }

        // Fill the distributed portion of the original structure
        clique::mpi::Barrier( comm );
        const double initStartTime = clique::mpi::Time();
        clique::symbolic::SymmOrig SOrig;
        FillOrigStruct
        ( nx, ny, nz, cutoff, comm, log2CommSize, SOrig );
        clique::mpi::Barrier( comm );
        const double initStopTime = clique::mpi::Time();
        if( commRank == 0 )
            std::cout << "Init time: " << initStopTime-initStartTime << " secs"
                      << std::endl;

        if( writeInfo )
        {
            const int numLocalSupernodes = SOrig.local.supernodes.size();
            const int numDistSupernodes = SOrig.dist.supernodes.size();
            infoFile << "Local original structure sizes:\n";
            for( int s=0; s<numLocalSupernodes; ++s )
                infoFile << "  " << s << ": supernode=" 
                         << SOrig.local.supernodes[s].size << ", "
                         << " lower=" 
                         << SOrig.local.supernodes[s].lowerStruct.size() 
                         << "\n";
            infoFile << "\n"
                     << "Dist original structure sizes:\n";
            for( int s=0; s<numDistSupernodes; ++s )
                infoFile << "  " << s << ": supernode="
                         << SOrig.dist.supernodes[s].size << ", "
                         << " lower="
                         << SOrig.dist.supernodes[s].lowerStruct.size() << "\n";
            infoFile << std::endl;
        }

        // Call the symbolic factorization routine
        const double symbFactStartTime = clique::mpi::Time();
        clique::symbolic::SymmFact S;
        clique::symbolic::SymmetricFactorization( SOrig, S, true );
        clique::mpi::Barrier( comm );
        const double symbFactStopTime = clique::mpi::Time();
        if( commRank == 0 )
            std::cout << "Symbolic factorization time: " 
                      << symbFactStopTime-symbFactStartTime << " secs"
                      << std::endl;

        if( writeInfo )
        {
            const int numLocalSupernodes = S.local.supernodes.size();
            const int numDistSupernodes = S.dist.supernodes.size();
            infoFile << "Local factor structure sizes:\n";
            for( int s=0; s<numLocalSupernodes; ++s )
                infoFile << "  " << s << ": supernode=" 
                         << S.local.supernodes[s].size << ", "
                         << " lower=" 
                         << S.local.supernodes[s].lowerStruct.size() << "\n";
            infoFile << "\n"
                     << "Dist factor structure sizes:\n";
            for( int s=0; s<numDistSupernodes; ++s )
                infoFile << "  " << s << ": supernode="
                         << S.dist.supernodes[s].size << ", "
                         << " lower="
                         << S.dist.supernodes[s].lowerStruct.size() << "\n";
            infoFile << std::endl;

            infoFile << "\n"
                     << "Local height 1d: " 
                     << S.dist.supernodes.back().localOffset1d + 
                        S.dist.supernodes.back().localSize1d << std::endl;
        }

        std::vector<clique::numeric::SymmFrontTree<F>*> symmFrontTrees;
        for( int panel=0; panel<numPanels; ++panel )
        {
            if( commRank == 0 )
                std::cout << "Panel " << panel << " of " << numPanels 
                          << std::endl;
            if( writeInfo )
                infoFile << "Panel " << panel << " of " << numPanels 
                         << std::endl;
            // Directly initialize the frontal matrices with the original 
            // sparse matrix (for now, use an original matrix equal to identity)
            const double fillStartTime = clique::mpi::Time();
            symmFrontTrees.push_back( new clique::numeric::SymmFrontTree<F> );
            clique::numeric::SymmFrontTree<F>& L = *symmFrontTrees.back();
            L.local.fronts.resize( S.local.supernodes.size() );
            for( unsigned s=0; s<S.local.supernodes.size(); ++s )
            {
                const clique::symbolic::LocalSymmFactSupernode& sn =
                    S.local.supernodes[s];
                elemental::Matrix<F>& frontL = L.local.fronts[s].frontL;

                const int frontSize = sn.size+sn.lowerStruct.size();
                frontL.ResizeTo( frontSize, sn.size );
                frontL.SetToRandom();
                if( writeInfo )
                    frontL.Print( infoFile, "frontL local" );
            }
            L.dist.mode = clique::MANY_RHS;
            L.dist.fronts.resize( log2CommSize+1 );
            clique::numeric::InitializeDistLeaf( S, L );
            for( unsigned s=1; s<log2CommSize+1; ++s )
            {
                const clique::symbolic::DistSymmFactSupernode& sn = 
                    S.dist.supernodes[s];
                elemental::DistMatrix<F,elemental::MC,elemental::MR>& front2dL =
                    L.dist.fronts[s].front2dL;

                front2dL.SetGrid( *sn.grid );
                const int frontSize = sn.size+sn.lowerStruct.size();
                front2dL.Align( 0, 0 );
                front2dL.ResizeTo( frontSize, sn.size );
                front2dL.SetToRandom();
                if( writeInfo )
                    front2dL.Print( infoFile, "frontL dist" );
            }
            clique::mpi::Barrier( comm );
            const double fillStopTime = clique::mpi::Time();
            if( commRank == 0 )
                std::cout << "Fill time: " << fillStopTime-fillStartTime
                          << " secs" << std::endl;

            // Generate a random RHS, multiply it by our matrix, and then 
            // compare the original RHS against the solution.
            const int NUM_RHS = 1;
            const double makeRhsStartTime = clique::mpi::Time();
            const int localHeight1d = 
                S.dist.supernodes.back().localOffset1d + 
                S.dist.supernodes.back().localSize1d;
            elemental::Matrix<F> localX( localHeight1d, NUM_RHS );
            localX.SetToRandom();
            elemental::Matrix<F> localYLower = localX;
            clique::numeric::SetSolveMode( L, clique::FEW_RHS );
            clique::numeric::LowerMultiply
            ( elemental::NORMAL, elemental::NON_UNIT, -1, S, L, localYLower );
            elemental::Matrix<F> localY = localX;
            clique::numeric::LowerMultiply
            ( elemental::TRANSPOSE, elemental::NON_UNIT, 0, S, L, localY );
            elemental::Axpy( (F)1, localYLower, localY );
            localYLower.Empty();
            clique::numeric::SetSolveMode( L, clique::MANY_RHS );
            clique::mpi::Barrier( comm );
            const double makeRhsStopTime = clique::mpi::Time();
            if( commRank == 0 )
                std::cout << "Make RHS time: " 
                          << makeRhsStopTime-makeRhsStartTime
                          << " secs" << std::endl;
            if( writeInfo )
            {
                localX.Print( infoFile, "localX" );
                localY.Print( infoFile, "localY" );
            }
            const double myYNorm = elemental::Norm( localY );
            double YNorm;
            clique::mpi::Reduce
            ( &myYNorm, &YNorm, 1, clique::mpi::SUM, 0, comm );

            // Call the numerical factorization routine
            elemental::SetBlocksize( factBlocksize );
            const double factStartTime = clique::mpi::Time();
            clique::numeric::LDL( elemental::TRANSPOSE, S, L );
            clique::mpi::Barrier( comm );
            const double factStopTime = clique::mpi::Time();
            if( commRank == 0 )
                std::cout << "Factorization time: " 
                          << factStopTime-factStartTime
                          << " secs" << std::endl;

            // Redistribute to 1d for fast solves with few right-hand sides
            const double redistStartTime = clique::mpi::Time();
            clique::numeric::SetSolveMode( L, clique::FEW_RHS_FAST_LDL );
            clique::mpi::Barrier( comm );
            const double redistStopTime = clique::mpi::Time();
            if( commRank == 0 )
                std::cout << "Redistribution time: " 
                          << redistStopTime-redistStartTime 
                          << " secs" << std::endl;

            // Solve
            elemental::SetBlocksize( solveBlocksize );
            clique::mpi::Barrier( comm );
            const double solveStartTime = clique::mpi::Time();
            clique::numeric::LDLSolve( clique::TRANSPOSE, S, L, localY );
            clique::mpi::Barrier( comm );
            const double solveStopTime = clique::mpi::Time();
            if( commRank == 0 )
                std::cout << "Solve time: " << solveStopTime-solveStartTime
                          << " secs" << std::endl;

            if( writeInfo )
                localY.Print( infoFile, "localY final" );
            elemental::Axpy( (F)-1, localX, localY );
            const double myErrorNorm = elemental::Norm( localY );
            double errorNorm;
            clique::mpi::Reduce
            ( &myErrorNorm, &errorNorm, 1, clique::mpi::SUM, 0, comm );
            if( commRank == 0 )
            {
                std::cout << "||y||_2: " << YNorm << "\n"
                          << "||e||_2: " << errorNorm << "\n" 
                          << "||e||_2/||y||_2: " << errorNorm/YNorm << "\n"
                          << std::endl;
            }
            if( writeInfo )
                localY.Print( infoFile, "error final" );
        }
    }
    catch( std::exception& e )
    {
#ifndef RELEASE
        elemental::DumpCallStack();
        clique::DumpCallStack();
#endif
        std::ostringstream msg;
        msg << "Process " << commRank << " caught message:\n"
            << e.what() << "\n";
        std::cerr << msg.str() << std::endl;
    }

    clique::Finalize();
    return 0;
}
 
int ReorderedIndexRecursion
( int x, int y, int z, int nx, int ny, int nz,
  int stepsLeft, int cutoff, int offset )
{
    const int size = nx*ny*nz;
#ifndef RELEASE
    if( stepsLeft != 0 && size == 0 )
        throw std::logic_error("Null supernode in the upper tree");
#endif
    if( stepsLeft == 0 && size <= cutoff )
    {
        // We have satisfied the nested dissection constraints
        return offset + (x+y*nx+z*nx*ny);
    }
    else if( nx >= ny )
    {
        // Partition the X dimension
        const int middle = (nx-1)/2;
        if( x < middle )
        {
            return ReorderedIndexRecursion
            ( x, y, z, middle, ny, nz, 
              std::max(stepsLeft-1,0), cutoff, offset );
        }
        else if( x == middle )
        {
            return offset + std::max(nx-1,0)*ny*nz + (y+z*ny);
        }
        else // x > middle
        {
            return ReorderedIndexRecursion
            ( x-middle-1, y, z, std::max(nx-middle-1,0), ny, nz,
              std::max(stepsLeft-1,0), cutoff, offset+middle*ny*nz );
        }
    }
    else
    {
        // Partition the Y dimension
        const int middle = (ny-1)/2;
        if( y < middle )
        {
            return ReorderedIndexRecursion
            ( x, y, z, nx, middle, nz,
              std::max(stepsLeft-1,0), cutoff, offset );
        }
        else if( y == middle )
        {
            return offset + nx*std::max(ny-1,0)*nz + (x+z*nx);
        }
        else // y > middle 
        {
            return ReorderedIndexRecursion
            ( x, y-middle-1, z, nx, std::max(ny-middle-1,0), nz,
              std::max(stepsLeft-1,0), cutoff, offset+nx*middle*nz );
        }
    }
}

int ReorderedIndex
( int x, int y, int z, int nx, int ny, int nz,
  int minBalancedDepth, int cutoff )
{
#ifndef RELEASE    
    clique::PushCallStack("ReorderedIndex");
#endif
    int index = ReorderedIndexRecursion
    ( x, y, z, nx, ny, nz, minBalancedDepth, cutoff, 0 );
#ifndef RELEASE
    clique::PopCallStack();
#endif
    return index;
}

void FillOrigStruct
( int nx, int ny, int nz, int cutoff, clique::mpi::Comm comm, int log2CommSize,
  clique::symbolic::SymmOrig& SOrig )
{
#ifndef RELEASE
    clique::PushCallStack("FillOrigStruct");
#endif
    int nxSub=nx, nySub=ny, xOffset=0, yOffset=0;
    FillDistOrigStruct
    ( nx, ny, nz, nxSub, nySub, xOffset, yOffset, cutoff, 
      comm, log2CommSize, SOrig );
    FillLocalOrigStruct
    ( nx, ny, nz, nxSub, nySub, xOffset, yOffset, cutoff, 
      log2CommSize, SOrig );
#ifndef RELEASE
    clique::PopCallStack();
#endif
}
  
void FillDistOrigStruct
( int nx, int ny, int nz, int& nxSub, int& nySub, int& xOffset, int& yOffset, 
  int cutoff, clique::mpi::Comm comm, int log2CommSize, 
  clique::symbolic::SymmOrig& S )
{
#ifndef RELEASE
    clique::PushCallStack("FillDistOrigStruct");
#endif
    const int commRank = clique::mpi::CommRank( comm );
    S.dist.comm = comm;
    S.dist.supernodes.resize( log2CommSize+1 );
    // Fill the distributed nodes
    for( int s=log2CommSize; s>0; --s )
    {
        clique::symbolic::DistSymmOrigSupernode& sn = S.dist.supernodes[s];
        const int powerOfTwo = 1u<<(s-1);
        const bool onLeft = (commRank&powerOfTwo) == 0;
        if( nxSub >= nySub )
        {
            // Form the structure of a partition of the X dimension
            const int middle = (nxSub-1)/2;
            sn.size = nySub*nz;
            sn.offset = 
                ReorderedIndex
                ( xOffset+middle, yOffset, 0, nx, ny, nz, 
                  log2CommSize, cutoff );

            // Allocate space for the lower structure
            int numJoins = 0;
            if( yOffset-1 >= 0 )
                ++numJoins;
            if( yOffset+nySub < ny )
                ++numJoins;
            sn.lowerStruct.resize( numJoins*nz );

            // Fill the (unsorted) lower structure
            int joinOffset = 0;
            if( yOffset-1 >= 0 )
            {
                for( int i=0; i<nz; ++i )
                    sn.lowerStruct[i] = ReorderedIndex
                    ( xOffset+middle, yOffset-1, i, nx, ny, nz, 
                      log2CommSize, cutoff );
                joinOffset += nz;
            }
            if( yOffset+nySub < ny )
            {
                for( int i=0; i<nz; ++i )
                    sn.lowerStruct[joinOffset+i] = ReorderedIndex
                    ( xOffset+middle, yOffset+nySub, i, nx, ny, nz, 
                      log2CommSize, cutoff );
            }

            // Sort the lower structure
            std::sort( sn.lowerStruct.begin(), sn.lowerStruct.end() );
            
            // Pick the new offsets and sizes based upon our rank
            if( onLeft )
            {
                xOffset = xOffset;
                nxSub = middle;
            }
            else
            {
                xOffset = xOffset+middle+1;
                nxSub = std::max(nxSub-middle-1,0);
            }
        }
        else
        {
            // Form the structure of a partition of the Y dimension
            const int middle = (nySub-1)/2;
            sn.size = nxSub*nz;
            sn.offset = 
                ReorderedIndex
                ( xOffset, yOffset+middle, 0, nx, ny, nz, 
                  log2CommSize, cutoff );

            // Allocate space for the lower structure
            int numJoins = 0;
            if( xOffset-1 >= 0 )
                ++numJoins;
            if( xOffset+nxSub < nx )
                ++numJoins;
            sn.lowerStruct.resize( numJoins*nz );

            // Fill the (unsorted) lower structure
            int joinOffset = 0;
            if( xOffset-1 >= 0 )
            {
                for( int i=0; i<nz; ++i )
                    sn.lowerStruct[i] = ReorderedIndex
                    ( xOffset-1, yOffset+middle, i, nx, ny, nz, 
                      log2CommSize, cutoff );
                joinOffset += nz;
            }
            if( xOffset+nxSub < nx )
            {
                for( int i=0; i<nz; ++i )
                    sn.lowerStruct[joinOffset+i] = ReorderedIndex
                    ( xOffset+nxSub, yOffset+middle, i, nx, ny, nz,
                      log2CommSize, cutoff );
            }

            // Sort the lower structure
            std::sort( sn.lowerStruct.begin(), sn.lowerStruct.end() );

            // Pick the new offsets and sizes based upon our rank
            if( onLeft )
            {
                yOffset = yOffset;
                nySub = middle;
            }
            else
            {
                yOffset = yOffset+middle+1;
                nySub = std::max(nySub-middle-1,0);
            }
        }
    }

    // Fill the bottom node, which is only owned by a single process
    clique::symbolic::DistSymmOrigSupernode& sn = S.dist.supernodes[0];
    if( nxSub*nySub*nz <= cutoff )
    {
        sn.size = nxSub*nySub*nz;
        sn.offset = ReorderedIndex
            ( xOffset, yOffset, 0, nx, ny, nz, log2CommSize, cutoff );

        // Count, allocate, and fill the lower struct
        int joinSize = 0;
        if( xOffset-1 >= 0 )
            joinSize += nySub*nz;
        if( xOffset+nxSub < nx )
            joinSize += nySub*nz;
        if( yOffset-1 >= 0 )
            joinSize += nxSub*nz;
        if( yOffset+nySub < ny )
            joinSize += nxSub*nz;
        sn.lowerStruct.resize( joinSize );

        int joinOffset = 0;
        if( xOffset-1 >= 0 )
        {
            for( int i=0; i<nz; ++i )
                for( int j=0; j<nySub; ++j )
                    sn.lowerStruct[i*nySub+j] = ReorderedIndex
                    ( xOffset-1, yOffset+j, i, 
                      nx, ny, nz, log2CommSize, cutoff );
            joinOffset += nySub*nz;
        }
        if( xOffset+nxSub < nx )
        {
            for( int i=0; i<nz; ++i )
                for( int j=0; j<nySub; ++j )
                    sn.lowerStruct[joinOffset+i*nySub+j] = ReorderedIndex
                    ( xOffset+nxSub, yOffset+j, i, 
                      nx, ny, nz, log2CommSize, cutoff );
            joinOffset += nySub*nz;
        }
        if( yOffset-1 >= 0 )
        {
            for( int i=0; i<nz; ++i )
                for( int j=0; j<nxSub; ++j )
                    sn.lowerStruct[joinOffset+i*nxSub+j] = ReorderedIndex
                    ( xOffset+j, yOffset-1, i, 
                      nx, ny, nz, log2CommSize, cutoff );
            joinOffset += nxSub*nz;
        }
        if( yOffset+nySub < ny )
        {
            for( int i=0; i<nz; ++i )
                for( int j=0; j<nxSub; ++j )
                    sn.lowerStruct[joinOffset+i*nxSub+j] = ReorderedIndex
                    ( xOffset+j, yOffset+nySub, i, 
                      nx, ny, nz, log2CommSize, cutoff );
        }

        // Sort the lower structure
        std::sort( sn.lowerStruct.begin(), sn.lowerStruct.end() );
    }
    else if( nxSub >= nySub )
    {
        // Form the structure of a partition of the X dimension
        const int middle = (nxSub-1)/2;
        sn.size = nySub*nz;
        sn.offset = 
            ReorderedIndex
            ( xOffset+middle, yOffset, 0, nx, ny, nz, log2CommSize, cutoff );

        // Allocate space for the lower structure
        int numJoins = 0;
        if( yOffset-1 >= 0 )
            ++numJoins;
        if( yOffset+nySub < ny )
            ++numJoins;
        sn.lowerStruct.resize( numJoins*nz );

        // Fill the (unsorted) lower structure
        int joinOffset = 0;
        if( yOffset-1 >= 0 )
        {
            for( int i=0; i<nz; ++i )
                sn.lowerStruct[i] = ReorderedIndex
                ( xOffset+middle, yOffset-1, i, nx, ny, nz, 
                  log2CommSize, cutoff );
            joinOffset += nz;
        }
        if( yOffset+nySub < ny )
        {
            for( int i=0; i<nz; ++i )
                sn.lowerStruct[joinOffset+i] = ReorderedIndex
                ( xOffset+middle, yOffset+nySub, i, nx, ny, nz, 
                  log2CommSize, cutoff );
        }

        // Sort the lower structure
        std::sort( sn.lowerStruct.begin(), sn.lowerStruct.end() );
    }
    else
    {
        // Form the structure of a partition of the Y dimension
        const int middle = (nySub-1)/2;
        sn.size = nxSub*nz;
        sn.offset = 
            ReorderedIndex
            ( xOffset, yOffset+middle, 0, nx, ny, nz, log2CommSize, cutoff );

        // Allocate space for the lower structure
        int numJoins = 0;
        if( xOffset-1 >= 0 )
            ++numJoins;
        if( xOffset+nxSub < nx )
            ++numJoins;
        sn.lowerStruct.resize( numJoins*nz );

        // Fill the (unsorted) lower structure
        int joinOffset = 0;
        if( xOffset-1 >= 0 )
        {
            for( int i=0; i<nz; ++i )
                sn.lowerStruct[i] = ReorderedIndex
                ( xOffset-1, yOffset+middle, i, nx, ny, nz, 
                  log2CommSize, cutoff );
            joinOffset += nz;
        }
        if( xOffset+nxSub < nx )
        {
            for( int i=0; i<nz; ++i )
                sn.lowerStruct[joinOffset+i] = ReorderedIndex
                ( xOffset+nxSub, yOffset+middle, i, nx, ny, nz,
                  log2CommSize, cutoff );
        }

        // Sort the lower structure
        std::sort( sn.lowerStruct.begin(), sn.lowerStruct.end() );
    }
#ifndef RELEASE
    clique::PopCallStack();
#endif
}

// Even though this is only used within the following function, it must be
// declared outside of the function in order to be used as a template argument
// (for std::stack) since template arguments must have external linkage.
//
// This data holds the bounding box information.
struct Box
{
    int parentIndex, nx, ny, xOffset, yOffset;
    bool leftChild;
};

void FillLocalOrigStruct
( int nx, int ny, int nz, int nxSub, int nySub, int xOffset, int yOffset, 
  int cutoff, int log2CommSize, clique::symbolic::SymmOrig& S )
{
#ifndef RELEASE
    clique::PushCallStack("FillLocalOrigStruct");
#endif
    // First count the depth, resize, and then run the actual fill
    int numSupernodes = 0;
    CountLocalTreeSize( nxSub, nySub, nz, cutoff, numSupernodes );
    S.local.supernodes.resize( numSupernodes );
    
    // Initialize with the root's box
    std::stack<Box> boxStack;
    {
        Box box;
        box.parentIndex = -1;
        box.nx = nxSub;
        box.ny = nySub;
        box.xOffset = xOffset;
        box.yOffset = yOffset;
        box.leftChild = false;
        boxStack.push(box);
    }

    // Fill the local tree
    for( int s=numSupernodes-1; s>=0; --s )
    {
        Box box = boxStack.top();
        boxStack.pop();

        clique::symbolic::LocalSymmOrigSupernode& sn = S.local.supernodes[s];
        sn.parent = box.parentIndex;
        if( sn.parent != -1 )
        {
            if( box.leftChild )
                S.local.supernodes[sn.parent].children[0] = s;
            else
                S.local.supernodes[sn.parent].children[1] = s;
        }

        if( box.nx*box.ny*nz <= cutoff )
        {
            sn.size = box.nx*box.ny*nz;
            sn.offset = ReorderedIndex
                ( box.xOffset, box.yOffset, 0, nx, ny, nz, 
                  log2CommSize, cutoff );
            sn.children.clear();

            // Count, allocate, and fill the lower struct
            int joinSize = 0;
            if( box.xOffset-1 >= 0 )
                joinSize += box.ny*nz;
            if( box.xOffset+box.nx < nx )
                joinSize += box.ny*nz;
            if( box.yOffset-1 >= 0 )
                joinSize += box.nx*nz;
            if( box.yOffset+box.ny < ny )
                joinSize += box.nx*nz;
            sn.lowerStruct.resize( joinSize );

            int joinOffset = 0;
            if( box.xOffset-1 >= 0 )
            {
                for( int i=0; i<nz; ++i )
                    for( int j=0; j<box.ny; ++j )
                        sn.lowerStruct[i*box.ny+j] = ReorderedIndex
                        ( box.xOffset-1, box.yOffset+j, i, 
                          nx, ny, nz, log2CommSize, cutoff );
                joinOffset += box.ny*nz;
            }
            if( box.xOffset+box.nx < nx )
            {
                for( int i=0; i<nz; ++i )
                    for( int j=0; j<box.ny; ++j )
                        sn.lowerStruct[joinOffset+i*box.ny+j] = ReorderedIndex
                        ( box.xOffset+box.nx, box.yOffset+j, i, 
                          nx, ny, nz, log2CommSize, cutoff );
                joinOffset += box.ny*nz;
            }
            if( box.yOffset-1 >= 0 )
            {
                for( int i=0; i<nz; ++i )
                    for( int j=0; j<box.nx; ++j )
                        sn.lowerStruct[joinOffset+i*box.nx+j] = ReorderedIndex
                        ( box.xOffset+j, box.yOffset-1, i, 
                          nx, ny, nz, log2CommSize, cutoff );
                joinOffset += box.nx*nz;
            }
            if( box.yOffset+box.ny < ny )
            {
                for( int i=0; i<nz; ++i )
                    for( int j=0; j<box.nx; ++j )
                        sn.lowerStruct[joinOffset+i*box.nx+j] = ReorderedIndex
                        ( box.xOffset+j, box.yOffset+box.ny, i,
                          nx, ny, nz, log2CommSize, cutoff );
            }

            // Sort the lower structure
            std::sort( sn.lowerStruct.begin(), sn.lowerStruct.end() );
        }
        else
        {
            sn.children.resize(2);
            if( box.nx >= box.ny )
            {
                // Partition the X dimension (this is the separator)
                const int middle = (box.nx-1)/2;
                sn.size = box.ny*nz;
                sn.offset = ReorderedIndex
                    ( box.xOffset+middle, box.yOffset, 0, nx, ny, nz,
                      log2CommSize, cutoff );

                // Count, allocate, and fill the lower struct
                int numJoins = 0;
                if( box.yOffset-1 >= 0 )
                    ++numJoins;
                if( box.yOffset+box.ny < ny )
                    ++numJoins;
                sn.lowerStruct.resize( numJoins*nz );

                int joinOffset = 0;
                if( box.yOffset-1 >= 0 )
                {
                    for( int i=0; i<nz; ++i )
                        sn.lowerStruct[i] = ReorderedIndex
                        ( box.xOffset+middle, box.yOffset-1, i, nx, ny, nz,
                          log2CommSize, cutoff );
                    joinOffset += nz;
                }
                if( box.yOffset+box.ny < ny )
                {
                    for( int i=0; i<nz; ++i )
                        sn.lowerStruct[joinOffset+i] = ReorderedIndex
                        ( box.xOffset+middle, box.yOffset+box.ny, i, nx, ny, nz,
                          log2CommSize, cutoff );
                }

                // Sort the lower structure
                std::sort( sn.lowerStruct.begin(), sn.lowerStruct.end() );

                // Push the left child box onto the stack
                Box leftBox;
                leftBox.parentIndex = s;
                leftBox.nx = middle;
                leftBox.ny = box.ny;
                leftBox.xOffset = box.xOffset;
                leftBox.yOffset = box.yOffset;
                leftBox.leftChild = true;
                boxStack.push( leftBox );

                // Push the right child box onto the stack
                Box rightBox;
                rightBox.parentIndex = s;
                rightBox.nx = std::max(box.nx-middle-1,0);
                rightBox.ny = box.ny;
                rightBox.xOffset = box.xOffset+middle+1;
                rightBox.yOffset = box.yOffset;
                rightBox.leftChild = false;
                boxStack.push( rightBox );
            }
            else 
            {
                // Partition the Y dimension (this is the separator)
                const int middle = (box.ny-1)/2;
                sn.size = box.nx*nz;
                sn.offset = ReorderedIndex
                    ( box.xOffset, box.yOffset+middle, 0, nx, ny, nz,
                      log2CommSize, cutoff );

                // Count, allocate, and fill the lower struct
                int numJoins = 0;
                if( box.xOffset-1 >= 0 )
                    ++numJoins;
                if( box.xOffset+box.nx < nx )
                    ++numJoins;
                sn.lowerStruct.resize( numJoins*nz );

                int joinOffset = 0;
                if( box.xOffset-1 >= 0 )
                {
                    for( int i=0; i<nz; ++i )
                        sn.lowerStruct[i] = ReorderedIndex
                        ( box.xOffset-1, box.yOffset+middle, i, nx, ny, nz,
                          log2CommSize, cutoff );
                    joinOffset += nz;
                }
                if( box.xOffset+box.nx < nx )
                {
                    for( int i=0; i<nz; ++i )
                        sn.lowerStruct[joinOffset+i] = ReorderedIndex
                        ( box.xOffset+box.nx, box.yOffset+middle, i, nx, ny, nz,
                          log2CommSize, cutoff );
                }

                // Sort the lower structure
                std::sort( sn.lowerStruct.begin(), sn.lowerStruct.end() );

                // Push the left child box onto the stack
                Box leftBox;
                leftBox.parentIndex = s;
                leftBox.nx = box.nx;
                leftBox.ny = middle;
                leftBox.xOffset = box.xOffset;
                leftBox.yOffset = box.yOffset;
                leftBox.leftChild = true;
                boxStack.push( leftBox );
                
                // Push the right child box onto the stack
                Box rightBox;
                rightBox.parentIndex = s;
                rightBox.nx = box.nx;
                rightBox.ny = std::max(box.ny-middle-1,0);
                rightBox.xOffset = box.xOffset;
                rightBox.yOffset = box.yOffset+middle+1;
                rightBox.leftChild = false;
                boxStack.push( rightBox );
            }
        }
    }
#ifndef RELEASE
    clique::PopCallStack();
#endif
}

void CountLocalTreeSize
( int nx, int ny, int nz, int cutoff, int& numSupernodes )
{
    ++numSupernodes;
    const int size = nx*ny*nz;
    if( size <= cutoff )
    {
        // no-op
    }
    else if( nx >= ny )
    {
        // Partition the X dimension
        const int middle = (nx-1)/2;

        // Recurse on the left partition
        CountLocalTreeSize( middle, ny, nz, cutoff, numSupernodes );

        // Recurse on the right partition
        CountLocalTreeSize
        ( std::max(nx-middle-1,0), ny, nz, cutoff, numSupernodes );
    }
    else
    {
        // Partition the Y dimension
        const int middle = (ny-1)/2;

        // Recurse on the top partition
        CountLocalTreeSize( nx, middle, nz, cutoff, numSupernodes );

        // Recurse on the bottom partition
        CountLocalTreeSize
        ( nx, std::max(ny-middle-1,0), nz, cutoff, numSupernodes );
    }
}

