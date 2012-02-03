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
#include <rsf.hh>

#include "psp.hpp"
using namespace psp;

void Usage()
{
    std::cout << "Overthrust <omega> <imagShift> <numPlanesPerPanel> "
                 "<fact blocksize> <solve blocksize> <accelerate?> <SQMR?> "
                 "<viz?>\n"
              << "  <omega>: Frequency (in rad/sec) of problem\n"
              << "  <imagShift>: imaginary shift [2 pi is standard]\n"
              << "  <numPlanesPerPanel>: depth of sparse-direct solves\n"
              << "  <fact blocksize>: factorization algorithmic blocksize\n"
              << "  <solve blocksize>: solve algorithmic blocksize\n"
              << "  <accelerate?>: accelerate solves iff !=0\n"
              << "  <SQMR?>: GMRES iff 0, SQMR otherwise\n"
              << "  <full viz?>: Full visualization iff != 0\n"
              << std::endl;
}

int
main( int argc, char* argv[] )
{
    clique::Initialize( argc, argv );
    clique::mpi::Comm comm = clique::mpi::COMM_WORLD;
    const int commSize = clique::mpi::CommSize( comm );
    const int commRank = clique::mpi::CommRank( comm );

    sf_init(argc,argv); // Initialize RSF

    iRSF par(0), vel;   // input parameter, file
    oRSF wave;          // output file

    float fomega;
    par.get("omega",fomega); // frequency
    const double omega = fomega;

    float fimagShift;
    par.get("damping",fimagShift,2*SF_PI); // damping (imaginary frequency shift)
    const double imagShift = fimagShift;

    int numPlanesPerPanel;
    par.get("nppp",numPlanesPerPanel,4); // depth of sparse-direct solves

    int factBlocksize;
    par.get("fbs",factBlocksize,96); // block size for factorization

    int solveBlocksize;
    par.get("sbs",solveBlocksize,64); // block size for triangle solves

    bool accelerate;
    par.get("accelerate",accelerate,true); // if accelerate by inverting diagonal blocks

    bool useSQMR;
    par.get("useSQMR",useSQMR,false); // if use low-memory alternative to GMRES

    const bool fullVisualize = false;

    int Nx, Ny, Nz;

    vel.get("n1",Nx);
    vel.get("n2",Ny);
    vel.get("n3",Nz);

    if( commRank == 0 )
    {
        std::cerr << "Running with omega=" << omega 
                  << ", numPlanesPerPanel=" << numPlanesPerPanel
                  << ", and imagShift=" << imagShift << std::endl;
    }
    
    FiniteDiffControl<double> control;
    control.stencil = SEVEN_POINT;
    control.nx = Nx;
    control.ny = Ny;
    control.nz = Nz;

    float Dx, Dy, Dz;

    vel.get("d1",Dx);
    vel.get("d2",Dy);
    vel.get("d3",Dz);

    control.wx = Dx*(Nx-1);
    control.wy = Dy*(Ny-1);
    control.wz = Dz*(Nz-1);

    control.omega = omega;
    // sInv functions should be dimensionless, so C should have units of length,
    // so we should scale it by the size of our domain (i.e., 20)

    double wmax = SF_MAX(SF_MAX(control.wx,control.wy), control.wz);

    float PMLscale;
    par.get("pml_scale",PMLscale,1.5f); // dimensionless scale factor for setting PML

    control.Cx = PMLscale*wmax;
    control.Cy = PMLscale*wmax;
    control.Cz = PMLscale*wmax;

    int PMLdepth;
    par.get("pml_depth",PMLdepth,5); // depth of PML in grid points 

    control.bx = PMLdepth;
    control.by = PMLdepth;
    control.bz = PMLdepth;

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
            std::cerr << "px=" << px << ", py=" << py << ", pz=" << pz 
                      << std::endl;

        GridData<double> velocity( 1, Nx, Ny, Nz, XYZ, px, py, pz, comm );
        double* localVelocity = velocity.LocalBuffer();
        const int xLocalSize = velocity.XLocalSize();
        const int yLocalSize = velocity.YLocalSize();
        const int zLocalSize = velocity.ZLocalSize();
        const int xShift = velocity.XShift();
        const int yShift = velocity.YShift();
        const int zShift = velocity.ZShift();
        std::ostringstream os;
        os << "overthrust_" << commRank << ".dat";

        std::ifstream velocityFile;
        velocityFile.open( os.str().c_str(), std::ios::in|std::ios::binary );
        velocityFile.read
        ( (char*)localVelocity, 
          xLocalSize*yLocalSize*zLocalSize*sizeof(double) );
        velocityFile.close();

// visualize velocity
//        velocity.WriteVTKPlane( XY, Nz/2, "velocity-middleXY" );
//        velocity.WriteVTKPlane( XZ, Ny/2, "velocity-middleXZ" );
//        velocity.WriteVTKPlane( YZ, Nx/2, "velocity-middleYZ" );
        if( fullVisualize )
        {
            if( commRank == 0 )
            {
                std::cout << "Writing full velocity data...";
                std::cout.flush();
            }
            velocity.WriteVTKVolume("velocity");
            elemental::mpi::Barrier( comm );
            if( commRank == 0 )
                std::cout << "done" << std::endl;
        }

        elemental::SetBlocksize( factBlocksize );
        if( commRank == 0 )
            std::cerr << "Beginning to initialize..." << std::endl;
        clique::mpi::Barrier( comm );
        const double initialStartTime = clique::mpi::Time(); 
        helmholtz.Initialize( velocity, accelerate );
        clique::mpi::Barrier( comm );
        const double initialStopTime = clique::mpi::Time();
        const double initialTime = initialStopTime - initialStartTime;
        if( commRank == 0 )
            std::cerr << "Finished initialization: " << initialTime 
                      << " seconds." << std::endl;

	int nshots;
	par.get("nshots",nshots,1); // number of shots

	float Sx, Sy, Sz;
	par.get("sx",Sx,0.5);  // x shot location in relative units
	par.get("sy",Sy,0.5);  // y shot location in relative units
	par.get("sy",Sz,0.1);  // z shot location in relative units

	float gaussianDecay;
	par.get("decay",gaussianDecay,10.0f*Nx); // Gaussian decay factor 

        GridData<std::complex<double> > 
            B( nshots, Nx, Ny, Nz, XYZ, px, py, pz, comm );
        std::complex<double>* localB = B.LocalBuffer();

        for( int zLocal=0; zLocal<zLocalSize; ++zLocal )
        {
            const int z = zShift + zLocal*pz;
            const double zSqueeze = control.wz / wmax;
            const double Z = zSqueeze * z / (Nz+1.0);
            const double argZ = (Z-Sz*zSqueeze)*(Z-Sz*zSqueeze);

            for( int yLocal=0; yLocal<yLocalSize; ++yLocal )
            {
                const int y = yShift + yLocal*py;
		const double ySqueeze = control.wy / wmax;
                const double Y = ySqueeze * y / (Ny+1.0);
                const double argY = (Y-Sy*ySqueeze)*(Y-Sy*ySqueeze);
//                const double argYTwo = (Y-0.25)*(Y-0.25);
//                const double argYThree = (Y-0.75)*(Y-0.75);
                for( int xLocal=0; xLocal<xLocalSize; ++xLocal )
                {
                    const int x = xShift + xLocal*px;
		    const double xSqueeze = control.wx / wmax;
		    const double X = xSqueeze * x / (Nx+1.0);
		    const double argX = (X-Sx*xSqueeze)*(X-Sx*xSqueeze);
//                    const double argXTwo = (X-0.25)*(X-0.25);
//                    const double argXThree = (X-0.75)*(X-0.75);
                    
                    const int localIndex = 
                        nshots*(xLocal + yLocal*xLocalSize + 
				zLocal*xLocalSize*yLocalSize);
                    const std::complex<double> fOne = 
                        Nx*std::exp(-gaussianDecay*(argX+argY+argZ));
//                    const std::complex<double> fTwo = 
//                        Nx*std::exp(-10*Nx*(argXTwo+argYTwo+argZ));
//                    const std::complex<double> fThree = 
//                        Nx*std::exp(-10*Nx*(argXThree+argYThree+argZ));
                    localB[localIndex+0] = fOne;
//                    localB[localIndex+1] = fTwo;
//                    localB[localIndex+2] = fThree;
                }
            }
        }

	// B.WriteVTKPlane( XY, Nz/2, "source-middleXY" );
        // B.WriteVTKPlane( XZ, Ny/2, "source-middleXZ" );
        // B.WriteVTKPlane( YZ, Nx/2, "source-middleYZ" );
        if( fullVisualize )
        {
            if( commRank == 0 )
            {
                std::cout << "Writing source data...";
                std::cout.flush();
            }
            B.WriteVTKVolume("source");
            if( commRank == 0 )
                std::cout << "done" << std::endl;
        }

        elemental::SetBlocksize( solveBlocksize );
        if( commRank == 0 )
            std::cerr << "Beginning solve..." << std::endl;
        clique::mpi::Barrier( comm );
        const double solveStartTime = clique::mpi::Time();

	int GMRESrestart;
	par.get("gmres_restart",GMRESrestart,20); // GMRES subspace size
	
	float GMREStol;
	par.get("gmres_tol",GMREStol,1e-4); // GMRES relative tolerance

        if( useSQMR )
            helmholtz.SolveWithSQMR( B );
        else
            helmholtz.SolveWithGMRES( B, GMRESrestart, GMREStol);

        clique::mpi::Barrier( comm );
        const double solveStopTime = clique::mpi::Time();
        const double solveTime = solveStopTime - solveStartTime;
        if( commRank == 0 )
            std::cerr << "Finished solve: " << solveTime << " seconds." 
		      << std::endl;

//        B.WriteVTKPlane( XY, Nz/2, "solution-middleXY" );
//        B.WriteVTKPlane( XZ, Ny/2, "solution-middleXZ" );
//        B.WriteVTKPlane( YZ, Nx/2, "solution-middleYZ" );
        if( fullVisualize )
        {
            if( commRank == 0 )
            {
                std::cout << "Writing solution data...";
                std::cout.flush();
            }
            B.WriteVTKVolume("solution");
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

    exit(0);
}
