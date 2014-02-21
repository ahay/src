#include "rsf.hh"
#include "psp.hpp"


void helm_solver(const DistUniformGrid<double>& velocity,
                       DistUniformGrid<complex>& solution,
                 const complex& freq,
                 const int Adjoint)
{

        Discretization<double> disc
        ( omega, Nx, Ny, Nz, wx, wy, wz, 
          PML, PML, PML, PML, PML, pmlSize, sigma );

        DistHelmholtz<double> helmholtz
        ( disc, comm, damping, numPlanesPerPanel );


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


        helmholtz.Finalize();

}

