// Parallel Sweeping Preconditioner (PSP) for solving 3D Helmholtz equations.

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

int
main( int argc, char* argv[] )
{
    psp::Initialize( argc, argv );
    mpi::Comm comm = mpi::COMM_WORLD;
    const int commSize = mpi::CommSize( comm );
    const int commRank = mpi::CommRank( comm );

    atexit(psp::Finalize);

    try 
    {
        sf_init( argc, argv );
        iRSF in("velocity"); // velocity RSF file
        iRSF par(0); // command-line parameters

        int origNx, origNy, origNz;
        in.get( "n1", origNx );
        in.get( "n2", origNy );
        in.get( "n3", origNz );

        float origDx, origDy, origDz;
        in.get( "d1", origDx );
        in.get( "d2", origDy );
        in.get( "d3", origDz );
        const double wx = origDx*(origNx-1);
        const double wy = origDy*(origNy-1);
        const double wz = origDz*(origNz-1);
        
        int Nx, Ny, Nz;
        par.get( "n1", Nx, origNx);
        par.get( "n2", Ny, origNy);
        par.get( "n3", Nz, origNz);
        const float dx = wx/(Nx-1);
        const float dy = wy/(Ny-1);
        const float dz = wz/(Nz-1);

	    iRSF src("src");
	    int ncomp,nrec,nsrc;
	    src.get("n1",ncomp);
	    src.get("n2",nrec);
        src.get("n3",nsrc);

	    float srcx[nsrc],srcy[nsrc],srcz[nsrc];
        float recx[nrec][nsrc],recy[nrec][nsrc],recz[nrec][nsrc];
	    for (int isrc=0; isrc<nsrc; isrc++ ) {
            for (int irec=0; irec<nrec; irec++ ) { 
		        src >> srcx[isrc];
		        src >> srcy[isrc];
		        src >> srcz[isrc];
                src >> recx[irec][isrc];
                src >> recy[irec][isrc];
                src >> recz[irec][isrc];
            }
	    }

        int adjSimulation;
        par.get("adj",adjSimulation);

        int numSimulation;
        float adjreal[nrec][nsrc],adjimag[nrec][nsrc];

        if ( adjSimulation == 1 ) { 
            numSimulation=2;

            iRSF adj("adjsrc");
            adj.get("n1",ncomp);
            adj.get("n2",nrec);
            adj.get("n3",nsrc);

            for (int isrc=0; isrc<nsrc; isrc++ ) { 
                for (int irec=0; irec<nrec; irec++ ) { 
                    adj >> adjreal[irec][isrc];
                    adj >> adjimag[irec][isrc];
                }
            }
        }
        else if ( adjSimulation == 0) { 
            numSimulation=1;
        } else {
            sf_error("wrong adjSimulation: 1-adj, 0-forward");
        }


        if( commRank == 0 )
        {
            std::cout << "origNx=" << origNx << "\n"
                      << "origNy=" << origNy << "\n"
                      << "origNz=" << origNz << "\n"
                      << "origDx=" << origDx << "\n"
                      << "origDy=" << origDy << "\n"
                      << "origDz=" << origDz << "\n"
                      << "wx=" << wx << "\n"
                      << "wy=" << wy << "\n"
                      << "wz=" << wz << "\n"
                      << "Nx=" << Nx << "\n"
                      << "Ny=" << Ny << "\n"
                      << "Nz=" << Nz << "\n" 
                      << "dx=" << dx << "\n"
                      << "dy=" << dy << "\n"
                      << "dz=" << dz << "\n"
                      << std::endl;
        }

        float omega, freq;
        par.get( "freq", freq ); // frequency in HZ
	    omega=freq*2*SF_PI;

        float sigma;
        par.get( "sigma", sigma, 1.5 ); // magnitude of PML stretching

        int pmlSize;
        par.get( "pmlSize", pmlSize, 5 ); // number of grid points of PML

        if( commRank == 0 )
            std::cout << "omega=" << omega << "\n"
                      << "sigma=" << sigma << "\n"
                      << "pmlSize=" << pmlSize << std::endl;

        const double damping = 7.;
        const int numPlanesPerPanel = 4;
        const bool vtkVisualize = false;
        const int factBlocksize = 96;
        const int solveBlocksize = 64;

        // This uses 5 grid points of PML by default
        //Discretization<double> disc
        //( omega, Nx, Ny, Nz, wx, wy, wz, 
        //  PML, PML, PML, PML, DIRICHLET, pmlSize, sigma );

        Discretization<double> disc
        ( omega, Nx, Ny, Nz, wx, wy, wz, 
          PML, PML, PML, PML, PML, pmlSize, sigma );

        DistHelmholtz<double> helmholtz
        ( disc, comm, damping, numPlanesPerPanel );

        DistUniformGrid<double> velocity
	    ( origNx, origNy, origNz, comm );
        double* localVelocity = velocity.Buffer();
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
                    const int filePos = 
                        (x + y*origNx + z*origNx*origNy)*sizeof(float);
                    in.seek( filePos, SEEK_SET );
                    const int localIndex = 
                        xLocal + yLocal*xLocalSize + 
                        zLocal*xLocalSize*yLocalSize;
                    float speed;
                    in >> speed;
                    if( localIndex >= xLocalSize*yLocalSize*zLocalSize )
                        throw std::logic_error("local index out of bounds");
                    localVelocity[localIndex] = speed;
                }
            }
        }

        if( commRank == 0 )
            std::cout << "Interpolating to " 
                      << Nx << " x " << Ny << " x " << Nz << std::endl;
        velocity.InterpolateTo( Nx, Ny, Nz );
        if( commRank == 0 )
            std::cout << "done" << std::endl;

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
        mpi::Barrier( comm );
        if( commRank == 0 )
            std::cout << "Beginning to initialize..." << std::endl;
        const double initialStartTime = mpi::Time();
        helmholtz.Initialize( velocity );
        mpi::Barrier( comm );
        const double initialStopTime = mpi::Time();
        const double initialTime = initialStopTime - initialStartTime;
        if( commRank == 0 )
            std::cout << "Finished initialization: " << initialTime
                      << " seconds." << std::endl;

        DistUniformGrid<Complex<double> > B
        ( Nx, Ny, Nz, comm, nsrc*numSimulation);
        Complex<double>* localB = B.Buffer();
        double arg[3];
        for( int zLocal=0; zLocal<zLocalSize; ++zLocal )
        {
            const int z = zShift + zLocal*pz;
            const double Z = z / (Nz+1.0);
            for( int yLocal=0; yLocal<yLocalSize; ++yLocal )
            {
                const int y = yShift + yLocal*py;
                const double Y = y / (Ny+1.0);
                for( int xLocal=0; xLocal<xLocalSize; ++xLocal )
                {
                    const int x = xShift + xLocal*px;
                    const double X = x / (Nx+1.0);

                    const int localIndex =
                        nsrc*numSimulation*(xLocal + yLocal*xLocalSize +
                        zLocal*xLocalSize*yLocalSize);

		            for (int isrc=0;isrc<nsrc;isrc++ ) { 
                    	arg[0] = (X-srcx[isrc])*(X-srcx[isrc]);
                    	arg[1] = (Y-srcy[isrc])*(Y-srcy[isrc]);
                    	arg[2] = (Z-srcz[isrc])*(Z-srcz[isrc]);

                        localB[localIndex+isrc*numSimulation] = 
                              Nz*Exp(-Nx*Ny*(arg[0]+arg[1]+arg[2]));

                        if ( adjSimulation == 1) {
                            double freal=0.0;
                            double fimag=0.0;
                            double arg0[3];
                            for (int irec=0; irec<nrec; irec++ ) { 
                                arg0[0]=(X-recx[irec][isrc])*(X-recx[irec][isrc]);
                                arg0[1]=(Y-recy[irec][isrc])*(Y-recy[irec][isrc]);
                                arg0[2]=(Z-recz[irec][isrc])*(Z-recz[irec][isrc]);

                                freal=freal+adjreal[irec][isrc]*Nz*Exp(-Nx*Ny*(arg0[0]+arg0[1]+arg0[2]));
                                fimag=fimag+adjimag[irec][isrc]*Nz*Exp(-Nx*Ny*(arg0[0]+arg0[1]+arg0[2]));
                            }

                            localB[localIndex+isrc*numSimulation+1].real = freal;
                            localB[localIndex+isrc*numSimulation+1].imag = -fimag;
                        }

		            }
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
        helmholtz.Solve( B );
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
        {
            std::cout << "Beginning to write output data...please be patient"
            << std::endl;
        }

        const int xMaxLocalSize = elem::MaxLength( Nx, px );
        const int yMaxLocalSize = elem::MaxLength( Ny, py );
        const int zMaxLocalSize = elem::MaxLength( Nz, pz );
        const int maxLocalSize = nsrc*numSimulation*xMaxLocalSize*yMaxLocalSize*zMaxLocalSize;
        std::vector<Complex<double> > receiveData;

        if( commRank != 0 )
        {
            mpi::Gather
            ( B.Buffer(), maxLocalSize,
              &receiveData[0], maxLocalSize, 0, comm );
        }
        else
        {
            oRSF out("solution"); // output RSF file
            out.type( SF_COMPLEX );
            out.put( "n1", nsrc);
            out.put( "n2", Nx );
            out.put( "n3", Ny );
            out.put( "n4", Nz );
            out.put( "d1", 1);
            out.put( "d2", dx );
            out.put( "d3", dy );
            out.put( "d4", dz );

            receiveData.resize( maxLocalSize*commSize );
            mpi::Gather
            ( B.Buffer(), maxLocalSize,
              &receiveData[0], maxLocalSize, 0, comm );

            for( int proc=0; proc<commSize; ++proc )
            {
                Complex<double>* procSolution = &receiveData[proc*maxLocalSize];
                const int xProc = proc % px;
                const int yProc = (proc/px) % py;
                const int zProc = proc/(px*py);
                const int xSize = Length( Nx, xProc, px );
                const int ySize = Length( Ny, yProc, py );
                const int zSize = Length( Nz, zProc, pz );
                for( int zLocal=0; zLocal<zSize; ++zLocal )
                {
                    const int z = zProc + zLocal*pz;
                    for( int yLocal=0; yLocal<ySize; ++yLocal )
                    {
                        const int y = yProc + yLocal*py;
                        for( int xLocal=0; xLocal<xSize; ++xLocal )
                        {
                            const int x = xProc + xLocal*px;

                            const int procIndex = 
                                nsrc*numSimulation*(xLocal + yLocal*xSize + zLocal*xSize*ySize);

                            const int localIndex = 
                                xLocal + yLocal*xSize + zLocal*xSize*ySize;
                            float speed=localVelocity[localIndex];
                            float weight=(omega*omega)/(speed*speed*speed);

                            for (int isrc=0; isrc<nsrc; isrc++ ) { 

                                const int pos = (isrc+x*nsrc+y*nsrc*Nx+z*nsrc*Nx*Ny)*sizeof(sf_complex);
                                out.seek( pos, SEEK_SET );
                                                        
                                const float real = procSolution[procIndex+isrc*numSimulation].real;
                                const float imag = procSolution[procIndex+isrc*numSimulation].imag;

                                if ( adjSimulation == 1) { 
                                    const float real0=procSolution[procIndex+isrc*numSimulation+1].real;
                                    const float imag0=procSolution[procIndex+isrc*numSimulation+1].imag;

                                    float kernel_real=(real*real0-imag*imag0)*weight;
                                    float kernel_imag=(real*imag0+imag*real0)*weight;
                                    sf_complex value=sf_cmplx(kernel_real,kernel_imag);
                                    out << value;
                                } else { 
                                    sf_complex value = sf_cmplx(real,imag);
                                    out << value;
                                }
                            }
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
