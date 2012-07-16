/*
   Parallel Sweeping Preconditioner (PSP): a distributed-memory implementation
   of a sweeping preconditioner for 3d Helmholtz equations.

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
#include "rsf.hh"
#include "psp.hpp"
using namespace psp;

void Usage()
{
    std::cout << "Madagascar <omega>" << std::endl;
}

int
main( int argc, char* argv[] )
{
    psp::Initialize( argc, argv );
    mpi::Comm comm = mpi::COMM_WORLD;
    const int commSize = mpi::CommSize( comm );
    const int commRank = mpi::CommRank( comm );

    if( argc < 2 )
    {
        if( commRank == 0 )
            Usage();
        psp::Finalize();
        return 0;
    }

    try 
    {
        sf_init( argc, argv );
        iRSF in("inFile.rsf"), par(0);

        int Nx, Ny, Nz;
        in.get( "n1", Nx );
        in.get( "n2", Ny );
        in.get( "n3", Nz );

        float dx, dy, dz;
        in.get( "d1", dx );
        in.get( "d2", dy );
        in.get( "d3", dz );
        const double wx = dx*(Nx-1);
        const double wy = dy*(Ny-1);
        const double wz = dz*(Nz-1);

        if( commRank == 0 )
        {
            std::cout << "Nx=" << Nx << "\n"
                      << "Ny=" << Ny << "\n"
                      << "Nz=" << Nz << "\n"
                      << "dx=" << dx << "\n"
                      << "dy=" << dy << "\n"
                      << "dz=" << dz << "\n"
                      << "wx=" << wx << "\n"
                      << "wy=" << wy << "\n"
                      << "wz=" << wz << std::endl;
        }

        float omega;
        par.get( "omega", omega );

        if( commRank == 0 )
            std::cout << "omega=" << omega << std::endl;

        const double damping = 7.;
        const int numPlanesPerPanel = 4;
        const bool vtkVisualize = false;
        const int factBlocksize = 96;
        const int solveBlocksize = 64;

    // This uses 5 grid points of PML by default
    Discretization<double> disc
    ( omega, Nx, Ny, Nz, wx, wy, wz, PML, PML, PML, PML, DIRICHLET );


        DistHelmholtz<double> helmholtz
        ( disc, comm, damping, numPlanesPerPanel );

        DistUniformGrid<double> velocity( 1, Nx, Ny, Nz, XYZ, comm );
        double* localVelocity = velocity.LocalBuffer();
        const int xLocalSize = velocity.XLocalSize();
        const int yLocalSize = velocity.YLocalSize();
        const int zLocalSize = velocity.ZLocalSize();
        const int xShift = velocity.XShift();
        const int yShift = velocity.YShift();
        const int zShift = velocity.ZShift();
        const int px = velocity.XStride();
        const int py = velocity.YStride();
        const int pz = velocity.ZStride();
        for( int zLocal=0; zLocal<zLocalSize; ++zLocal )
        {
            const int z = zShift + zLocal*pz;
            for( int yLocal=0; yLocal<yLocalSize; ++yLocal )
            {
                const int y = yShift + yLocal*py;
                for( int xLocal=0; xLocal<xLocalSize; ++xLocal )
                {
                    const int x = xShift + xLocal*px;
                    const int filePos = (x + y*Nx + z*Nx*Ny)*sizeof(float);
                    in.seek( filePos, SEEK_SET );
                    const int localIndex = 
                        xLocal + yLocal*xLocalSize + 
                        zLocal*xLocalSize*yLocalSize;
                    float speed;
                    in >> speed;
                    localVelocity[localIndex] = speed;
                }
            }
        }

        if( vtkVisualize )
        {
            if( commRank == 0 )
            {
                std::cout << "Writing velocity data...";
                std::cout.flush();
            }
            velocity.WritePlane( XY, Nz/2, "velocity-middleXY" );
            velocity.WritePlane( XZ, Ny/2, "velocity-middleXZ" );
            velocity.WritePlane( YZ, Nx/2, "velocity-middleYZ" );
            velocity.WriteVolume("velocity");
            mpi::Barrier( comm );
            if( commRank == 0 )
                std::cout << "done" << std::endl;
        }

        elem::SetBlocksize( factBlocksize );
        if( commRank == 0 )
            std::cout << "Beginning to initialize..." << std::endl;
        mpi::Barrier( comm );
        const double initialStartTime = mpi::Time();
        helmholtz.Initialize( velocity );
        mpi::Barrier( comm );
        const double initialStopTime = mpi::Time();
        const double initialTime = initialStopTime - initialStartTime;
        if( commRank == 0 )
            std::cout << "Finished initialization: " << initialTime
                      << " seconds." << std::endl;

        DistUniformGrid<Complex<double> > B( 1, Nx, Ny, Nz, XYZ, comm );
        Complex<double>* localB = B.LocalBuffer();
        const double center[] = { 0.5, 0.5, 0.25 };
        double arg[3];
        for( int zLocal=0; zLocal<zLocalSize; ++zLocal )
        {
            const int z = zShift + zLocal*pz;
            const double Z = z / (Nz+1.0);
            arg[2] = (Z-center[2])*(Z-center[2]);
            for( int yLocal=0; yLocal<yLocalSize; ++yLocal )
            {
                const int y = yShift + yLocal*py;
                const double Y = y / (Ny+1.0);
                arg[1] = (Y-center[1])*(Y-center[1]);
                for( int xLocal=0; xLocal<xLocalSize; ++xLocal )
                {
                    const int x = xShift + xLocal*px;
                    const double X = x / (Nx+1.0);
                    arg[0] = (X-center[0])*(X-center[0]);

                    const int localIndex =
                        xLocal + yLocal*xLocalSize +
                        zLocal*xLocalSize*yLocalSize;
                    localB[localIndex] =
                        Nz*Exp(-Nx*Ny*(arg[0]+arg[1]+arg[2]));
                }
            }
        }

        if( vtkVisualize )
        {
            if( commRank == 0 )
            {
                std::cout << "Writing source data...";
                std::cout.flush();
            }
            B.WritePlane( XY, Nz/2, "source-middleXY" );
            B.WritePlane( XZ, Ny/2, "source-middleXZ" );
            B.WritePlane( YZ, Nx/2, "source-middleYZ" );
            B.WriteVolume("source");
            if( commRank == 0 )
                std::cout << "done" << std::endl;
        }

        elem::SetBlocksize( solveBlocksize );
        if( commRank == 0 )
            std::cout << "Beginning solve..." << std::endl;
        mpi::Barrier( comm );
        const double solveStartTime = mpi::Time();
        helmholtz.SolveWithGMRES( B );
        mpi::Barrier( comm );
        const double solveStopTime = mpi::Time();
        const double solveTime = solveStopTime - solveStartTime;
        if( commRank == 0 )
            std::cout << "Finished solve: " << solveTime << " seconds."
                      << std::endl;

        if( vtkVisualize )
        {
            if( commRank == 0 )
            {
                std::cout << "Writing solution data...";
                std::cout.flush();
            }
            B.WritePlane( XY, Nz/2, "solution-middleXY" );
            B.WritePlane( XZ, Ny/2, "solution-middleXZ" );
            B.WritePlane( YZ, Nx/2, "solution-middleYZ" );
            B.WriteVolume("solution");
            if( commRank == 0 )
                std::cout << "done" << std::endl;
        }

        helmholtz.Finalize();

        if( commRank == 0 )
            std::cout << "Beginning to write output data...please be patient"
                      << std::endl;
        oRSF out("outFile.rsf");
        out.type( SF_COMPLEX );
        out.put( "n1", Nx );
        out.put( "n2", Ny );
        out.put( "n3", Nz );
        out.put( "d1", dx );
        out.put( "d2", dy );
        out.put( "d3", dz );
        const int xMaxLocalSize = elem::MaxLocalLength( Nx, px );
        const int yMaxLocalSize = elem::MaxLocalLength( Ny, py );
        const int zMaxLocalSize = elem::MaxLocalLength( Nz, pz );
        const int maxLocalSize = xMaxLocalSize*yMaxLocalSize*zMaxLocalSize;
        std::vector<Complex<double> > receiveData;
        if( commRank != 0 )
        {
            mpi::Gather
            ( B.LocalBuffer(), maxLocalSize,
              &receiveData[0], maxLocalSize, 0, comm );
        }
        else
        {
            receiveData.resize( maxLocalSize*commSize );
            mpi::Gather
            ( B.LocalBuffer(), maxLocalSize,
              &receiveData[0], maxLocalSize, 0, comm );

            for( int proc=0; proc<commSize; ++proc )
            {
                Complex<double>* procSolution = &receiveData[proc*maxLocalSize];
                const int xProc = proc % px;
                const int yProc = (proc/px) % py;
                const int zProc = proc/(px*py);
                const int xSize = LocalLength( Nx, xProc, px );
                const int ySize = LocalLength( Ny, yProc, py );
                const int zSize = LocalLength( Nz, zProc, pz );
                for( int zLocal=0; zLocal<zSize; ++zLocal )
                {
                    const int z = zProc + zLocal*pz;
                    for( int yLocal=0; yLocal<ySize; ++yLocal )
                    {
                        const int y = yProc + yLocal*py;
                        for( int xLocal=0; xLocal<xSize; ++xLocal )
                        {
                            const int x = xProc + xLocal*px;
                            const int pos = (x+y*Nx+z*Nx*Ny)*sizeof(sf_complex);
                            out.seek( pos, SEEK_SET );
                            const int procIndex = 
                                xLocal + yLocal*xSize + zLocal*xLocalSize*ySize;
                            const float real = procSolution[procIndex].real;
                            const float imag = procSolution[procIndex].imag;
                            sf_complex value = sf_cmplx(real,imag);
                            out << value;
                        }
                    }
                }
            }
        }
    }
    catch( std::exception& e )
    {
        std::cerr << "Caught exception on process " << commRank << ":\n"
                  << e.what() << std::endl;
#ifndef RELEASE
        elem::DumpCallStack();
        cliq::DumpCallStack();
        psp::DumpCallStack();
#endif
    }

    psp::Finalize();
    return 0;
}
