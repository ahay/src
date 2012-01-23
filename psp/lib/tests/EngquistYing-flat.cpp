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
#include "psp.hpp"
using namespace psp;

void Usage()
{
    std::cout << "EngquistYing-flat <N> <omega> <imagShift> <velocity> "
                 "<numPlanesPerPanel> <fact blocksize> <solve blocksize> "
                 "<accelerate?> <SQMR?> <viz?>\n"
              << "  <N>: Size of grid in each dimension\n"
              << "  <omega>: Frequency (in rad/sec) of problem\n"
              << "  <imagShift>: imaginary shift [2 pi is standard]\n"
              << "  <velocity>: Which velocity field to use, {1,2}\n"
              << "  <numPlanesPerPanel>: depth of sparse-direct solves\n"
              << "  <fact blocksize>: factorization algorithmic blocksize\n"
              << "  <solve blocksize>: solve algorithmic blocksize\n"
              << "  <accelerate?>: accelerate solves iff !=0\n"
              << "  <SQMR?>: GMRES iff 0, SQMR otherwise\n"
              << "  <full viz?>: Full visualization iff != 0\n"
              << "\n"
              << "Please see \"Sweeping preconditioner for the Helmholtz "
                 "equation: moving perfectly matched layers\" for details\n"
              << std::endl;
}

int
main( int argc, char* argv[] )
{
    clique::Initialize( argc, argv );
    clique::mpi::Comm comm = clique::mpi::COMM_WORLD;
    const int commSize = clique::mpi::CommSize( comm );
    const int commRank = clique::mpi::CommRank( comm );

    if( argc < 11 )
    {
        if( commRank == 0 )
            Usage();
        clique::Finalize();
        return 0;
    }
    int argNum=1;
    const int N = atoi( argv[argNum++] );
    const double omega = atof( argv[argNum++] );
    const double imagShift = atof( argv[argNum++] );
    const int velocityModel = atoi( argv[argNum++] );
    const int numPlanesPerPanel = atoi( argv[argNum++] );
    const int factBlocksize = atoi( argv[argNum++] );
    const int solveBlocksize = atoi( argv[argNum++] );
    const bool accelerate = atoi( argv[argNum++] );
    const bool useSQMR = atoi( argv[argNum++] );
    const bool fullVisualize = atoi( argv[argNum++] );

    if( velocityModel < 1 || velocityModel > 2 )
    {
        if( commRank == 0 )
            std::cout << "Invalid velocity model choice, must be in {1,2};\n"
                      << "Please see \"Sweeping preconditioner for the "
                      << "Helmholtz equation: moving perfectly matched layers\""
                      << " for more details." << std::endl;
        clique::Finalize();
        return 0;
    }

    if( commRank == 0 )
    {
        std::cout << "Running with N=" << N << ", omega=" << omega 
                  << ", and numPlanesPerPanel=" << numPlanesPerPanel
                  << " with velocity model " << velocityModel << std::endl;
    }
    
    FiniteDiffControl<double> control;
    control.stencil = SEVEN_POINT;
    control.nx = N;
    control.ny = N;
    control.nz = N/8;
    control.wx = 1;
    control.wy = 1;
    control.wz = 0.125;
    control.omega = omega;
    control.Cx = 1.5;
    control.Cy = 1.5;
    control.Cz = 1.5;
    control.bx = 4;
    control.by = 4;
    control.bz = 4;
    /*
    control.bx = 6;
    control.by = 6;
    control.bz = 6;
    */
    control.imagShift = imagShift;
    control.cutoff = 96;
    control.numPlanesPerPanel = numPlanesPerPanel;
    control.frontBC = PML;
    control.rightBC = PML;
    control.backBC = PML;
    control.leftBC = PML;
    control.topBC = PML;

    try 
    {
        DistHelmholtz<double> helmholtz( control, comm );

        const int cubeRoot = 
            std::max(1,(int)std::floor(pow((double)commSize,0.333)));
        int px = cubeRoot;
        while( commSize % px != 0 )
            ++px;
        const int reduced = commSize / px;
        const int sqrtReduced = 
            std::max(1,(int)std::floor(sqrt((double)reduced)));
        int py = sqrtReduced;
        while( reduced % py != 0 )
            ++py;
        const int pz = reduced / py;
        if( px*py*pz != commSize )
            throw std::logic_error("Nonsensical process grid");
        else if( commRank == 0 )
            std::cout << "px=" << px << ", py=" << py << ", pz=" << pz 
                      << std::endl;

        GridData<double> velocity( 1, N, N, N/8, XYZ, px, py, pz, comm );
        double* localVelocity = velocity.LocalBuffer();
        const int xLocalSize = velocity.XLocalSize();
        const int yLocalSize = velocity.YLocalSize();
        const int zLocalSize = velocity.ZLocalSize();
        const int xShift = velocity.XShift();
        const int yShift = velocity.YShift();
        const int zShift = velocity.ZShift();
        if( velocityModel == 1 )
        {
            // Converging lens
            for( int zLocal=0; zLocal<zLocalSize; ++zLocal )
            {
                const int z = zShift + zLocal*pz;
                const double Z = z / (N+1.0);
                const double argZ = (Z-0.5/8)*(Z-0.5/8);
                for( int yLocal=0; yLocal<yLocalSize; ++yLocal )
                {
                    const int y = yShift + yLocal*py;
                    const double Y = y / (N+1.0);
                    const double argY = (Y-0.5)*(Y-0.5);
                    for( int xLocal=0; xLocal<xLocalSize; ++xLocal )
                    {
                        const int x = xShift + xLocal*px;
                        const double X = x / (N+1.0);
                        const double argX = (X-0.5)*(X-0.5);
                        
                        const int localIndex = 
                            xLocal + yLocal*xLocalSize + 
                            zLocal*xLocalSize*yLocalSize;
                        const double speed =
                            1.0 - 0.4*std::exp(-32.*(argX+argY+argZ));
                        localVelocity[localIndex] = speed / 0.8;
                    }
                }
            }
        }
        else // velocityModel == 2
        {
            // Wave guide
            for( int zLocal=0; zLocal<zLocalSize; ++zLocal )
            {
                for( int yLocal=0; yLocal<yLocalSize; ++yLocal )
                {
                    const int y = yShift + yLocal*py;
                    const double Y = y / (N+1.0);
                    const double argY = (Y-0.5)*(Y-0.5);
                    for( int xLocal=0; xLocal<xLocalSize; ++xLocal )
                    {
                        const int x = xShift + xLocal*px;
                        const double X = x / (N+1.0);
                        const double argX = (X-0.5)*(X-0.5);
                        
                        const int localIndex = 
                            xLocal + yLocal*xLocalSize + 
                            zLocal*xLocalSize*yLocalSize;
                        const double speed =
                            1.0 - 0.4*std::exp(-32.*(argX+argY));
                        localVelocity[localIndex] = speed / 0.8;
                    }
                }
            }
        }

        velocity.WritePlane( XY, N/16, "velocity-middleXY" );
        velocity.WritePlane( XZ, N/2,  "velocity-middleXZ" );
        velocity.WritePlane( YZ, N/2,  "velocity-middleYZ" );
        if( fullVisualize )
        {
            if( commRank == 0 )
            {
                std::cout << "Writing full velocity data...";
                std::cout.flush();
            }
            velocity.WriteVolume("velocity");
            elemental::mpi::Barrier( comm );
            if( commRank == 0 )
                std::cout << "done" << std::endl;
        }

        elemental::SetBlocksize( factBlocksize );
        if( commRank == 0 )
            std::cout << "Beginning to initialize..." << std::endl;
        clique::mpi::Barrier( comm );
        const double initialStartTime = clique::mpi::Time(); 
        helmholtz.Initialize( velocity, accelerate );
        clique::mpi::Barrier( comm );
        const double initialStopTime = clique::mpi::Time();
        const double initialTime = initialStopTime - initialStartTime;
        if( commRank == 0 )
            std::cout << "Finished initialization: " << initialTime 
                      << " seconds." << std::endl;

        GridData<std::complex<double> > 
            B( 2, N, N, N/8, XYZ, px, py, pz, comm );
        std::complex<double>* localB = B.LocalBuffer();
        const double dX = 0;
        const double dY = sqrt(2.0)/2.0;
        const double dZ = sqrt(2.0)/2.0;
        const std::complex<double> imagOne( 0.0, 1.0 );
        for( int zLocal=0; zLocal<zLocalSize; ++zLocal )
        {
            const int z = zShift + zLocal*pz;
            const double Z = z / (N+1.0);
            const double argZ = (Z-0.25/8)*(Z-0.25/8);
            for( int yLocal=0; yLocal<yLocalSize; ++yLocal )
            {
                const int y = yShift + yLocal*py;
                const double Y = y / (N+1.0);
                const double argYOne = (Y-0.5)*(Y-0.5);
                const double argYTwo = (Y-0.25)*(Y-0.25);
                for( int xLocal=0; xLocal<xLocalSize; ++xLocal )
                {
                    const int x = xShift + xLocal*px;
                    const double X = x / (N+1.0);
                    const double argX = (X-0.5)*(X-0.5);
                    
                    const int localIndex = 
                        2*(xLocal + yLocal*xLocalSize + 
                           zLocal*xLocalSize*yLocalSize);
                    const std::complex<double> fOne = 
                        N*std::exp(-N*N*(argX+argYOne+argZ));
                    const std::complex<double> fTwo = 
                        N*std::exp(-2*omega*(argX+argYTwo+argZ))*
                        std::exp(omega*imagOne*(X*dX+Y*dY+Z*dZ));
                    localB[localIndex+0] = fOne;
                    localB[localIndex+1] = fTwo;
                }
            }
        }

        B.WritePlane( XY, N/16, "source-middleXY" );
        B.WritePlane( XZ, N/2,  "source-middleXZ" );
        B.WritePlane( YZ, N/2,  "source-middleYZ" );
        if( fullVisualize )
        {
            if( commRank == 0 )
            {
                std::cout << "Writing source data...";
                std::cout.flush();
            }
            B.WriteVolume("source");
            if( commRank == 0 )
                std::cout << "done" << std::endl;
        }

        elemental::SetBlocksize( solveBlocksize );
        if( commRank == 0 )
            std::cout << "Beginning solve..." << std::endl;
        clique::mpi::Barrier( comm );
        const double solveStartTime = clique::mpi::Time();
        if( useSQMR )
            helmholtz.SolveWithSQMR( B );
        else
            helmholtz.SolveWithGMRES( B );
        clique::mpi::Barrier( comm );
        const double solveStopTime = clique::mpi::Time();
        const double solveTime = solveStopTime - solveStartTime;
        if( commRank == 0 )
            std::cout << "Finished solve: " << solveTime << " seconds." 
                      << std::endl;

        B.WritePlane( XY, N/16, "solution-middleXY" );
        B.WritePlane( XZ, N/2,  "solution-middleXZ" );
        B.WritePlane( YZ, N/2,  "solution-middleYZ" );
        if( fullVisualize )
        {
            if( commRank == 0 )
            {
                std::cout << "Writing solution data...";
                std::cout.flush();
            }
            B.WriteVolume("solution");
            if( commRank == 0 )
                std::cout << "done" << std::endl;
        }

        helmholtz.Finalize();
    }
    catch( std::exception& e )
    {
        std::cerr << "Caught exception on process " << commRank << ":\n"
                  << e.what() << std::endl;
#ifndef RELEASE
        elemental::DumpCallStack();
        clique::DumpCallStack();
#endif
    }

    clique::Finalize();
    return 0;
}
