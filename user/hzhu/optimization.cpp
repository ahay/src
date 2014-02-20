#include <iostream>
#include "rsf.hh" 

double gen_adjsrc
( const DistUniformGrid<Complex<double>>& S, 
  const DistUniformGrid<Complex<double>>& O, 
  const Matrix<double>& xRecv,
  const Matrix<double>& yRecv,
  const Matrix<double>& zRecv 
        Matrix<Complex<double>>& A )
{
    // TODO: Compute local misfits and then sum the result with an
    //       MPI_Allreduce summation
    // NOTE: This routine should already work but may not be as fast
    //       as possible.
    double l2MisfitSquared=0.0; 
    const int nsrc = xRecv.Height();
    const int nrec = xRecv.Width();
    A.Resize( nsrc, nrec );
    for (int isrc=0;isrc<nsrc;isrc++ ) { 
        for (int irec=0;irec<nrec;irec++ ) { 
            int x=round(xRecv[irec][isrc]*S.XSize());
            int y=round(yRecv[irec][isrc]*S.YSize());
            int z=round(zRecv[irec][isrc]*S.ZSize());

	        Complex<double> obs = O.Get(x,y,z);
	        Complex<double> syn = S.Get(x,y,z);
	        A.Set( irec, isrc, obs-syn );

            l2MisfitSquared += Abs(adj.Get(irec,isrc));
        }
    }
    return l2MisfitSquared;
}

void smooth_kernel
( const DistUniformGrid<double>& K,
        DistUniformGrid<double>& KSmooth, double sigma2)
{
    const int nx = K.XSize();
    const int ny = K.YSize();
    const int nz = K.ZSize();
    // NOTE: Assuming that KSmooth is already nx x ny x nz since the
    //       DistUniformGrid class does not yet support resizing.
    
    // TODO: A distributed version of the following operations. Since
    //       K and KSmooth are distributed, and the smoothing requires
    //       knowledge of a neighborhood of points, there will need to
    //       be a communication step which gathers the neighbor data.
    // (STEP 6)
    /*
    for ( int k=0;k<nz;k++ ) { 
        for (int j=0;j<ny;j++ ) { 
            for (int i=0;i<nx;i++ ) { 
                double kernel_sum=0.0;
                double gauss_sum=0.0;

                int kstart=k-2*sigma2;
                if ( kstart < 0 ) { 
                    kstart=0 ;
                }
                int kend=k+2*sigma2;
                if ( kend > nz ) { 
                    kend=nz;
                }
                int jstart=j-2*sigma2;
                if (jstart < 0 ) {
                    jstart =0 ;
                }
                int jend=j+2*sigma2;
                if (jend > ny ) { 
                    jend=ny;
                }
                int istart = i-2*sigma2;
                if ( istart < 0 ) { 
                    istart =0 ; 
                }
                int iend = i+2*sigma2;
                if ( iend > nx ) { 
                    iend = nx;
                }

                for (int ik=kstart;ik<kend;ik++){
                    for (int ij=jstart;ij<jend;ij++) { 
                        for (int ii=istart;ii<iend;ii++) { 
                            double r2=pow(ii-i,2)+pow(ij-j,2)+pow(ik-k,2);
                            double gauss=exp(-r2/(2*pow(sigma2,2)));
                            kernel_sum=kernel_sum+kernel[ii][ij][ik]*gauss;
                            gauss_sum=gauss_sum+gauss;
                        }
                    }
                }
                double kernel_smooth[i][j][k]=kernel_sum/gauss_sum;
            }
        }
    }
    */
}

void calc_gradient_sd
( const DistUniformGrid<double>& KCur, DistUniformGrid<double>& KDir )
{
    const int nx = KCur.XSize();
    const int ny = KCur.YSize();
    const int nz = KCur.ZSize();
    // NOTE: Assuming that KDir is already nx x ny x nz since the
    //       DistUniformGrid class does not yet support resizing.
    
    // TODO: Perform local updates (no communication should be required)
    // (STEP 1)
    
    double* localKCur=KCur.Buffer();
    double* localKDir=KDir.Buffer();

    const int xLocalSize=KCur.XLocalSize();
    const int yLocalSize=KCur.YLocalSize();
    const int zLocalSize=KCur.ZLocalSize();

    for ( int zLocal=0; zLoca<zLocalSize; ++zLocal) 
    {
        for ( int yLocal=0; yLocal<yLocalSize; ++yLocal) 
        {
            for (int xLocal=0; xLocal<xLocalSize; ++xLocal) 
            {
                const int localIndex=xLocal+
                            yLocal*xLocalSize+
                            zlocal*xLocalSize*yLocalSize;
                localKDir[localIndex]=-localKCur[localIndex];
            }
        }
    }
}
    
    
    

void calc_gradient_cg
( const DistUniformGrid<double>& KCur, 
  const DistUniformGrid<double>& KPrev, 
  const DistUniformGrid<double>& DirPrev,
        DistUniformGrid<dobule>& DirCur,
        mpi::Comm comm)
{
    const int nx = KCur.XSize();
    const int ny = KCur.YSize();
    const int nz = KCur.ZSize();

    double beta_top=0.0;
    double beta_bot=0.0;

    // TODO: Local summation followed by MPI_Allreduce summation
    // (STEP 2)
    
    double* localKCur=KCur.Buffer();
    double* localKPrev=KPrev.Buffer();
    double* localDirCur=DirCur.Buffer();
    double* localDirPrev=DirPrev.Buffer();

    const int xLocalSize=KCur.XLocalSize();
    const int yLocalSize=KCur.YLocalSize();
    const int zLocalSize=KCur.ZLocalSize();

    for ( int zLocal=0; zLoca<zLocalSize; ++zLocal) 
    {
        for ( int yLocal=0; yLocal<yLocalSize; ++yLocal) 
        {
            for (int xLocal=0; xLocal<xLocalSize; ++xLocal) 
            {
                const int localIndex=xLocal+
                            yLocal*xLocalSize+
                            zlocal*xLocalSize*yLocalSize;

                double diff=localKCur[localIndex]-localKPrev[localIndex];
                beta_top=beta_top+localKCur[localIndex]*diff;
                beta_bot=beta_bot+localKPrev[localIndex]*localKPrev[localIndex];
            }
        }
    }
    

    double betalocal=beta_top/beta_bot;
    double betaglob;

    mpi::AllReduce(&betalocal, &betaglob,1, mpi::SUM ,comm ); 

    if (betaglob < 0.0) { 
        betaglob=0.0;
    }

    // TODO: Local subtraction
    // (STEP 3)
    for ( int zLocal=0; zLoca<zLocalSize; ++zLocal) 
    {
        for ( int yLocal=0; yLocal<yLocalSize; ++yLocal) 
        {
            for (int xLocal=0; xLocal<xLocalSize; ++xLocal) 
            {
                const int localIndex=xLocal+
                            yLocal*xLocalSize+
                            zlocal*xLocalSize*yLocalSize;
                localDirCur[localIndex]=-localKCur[localIndex]+betaglob*localDirPrev[localIndex];
            }
        }
    }
}

void update_model
( const DistUniformGrid<double>& MOld, 
  const DistUniformGrid<double>& Grad,
        DistUniformGrid<double>& MNew, double alpha, mpi::Comm comm )
{
    int nx = MOld.XSize();
    int ny = MOld.YSize();
    int nz = MOld.ZSize();

    // TODO: Compute infinity norm of vector in parallel via a local
    //       infinity norm calculation followed by a call to 
    //       MPI_Allreduce with MPI_MAX
    // (STEP 4)
    double grdmaxlocal=0.0;

    double* localGrad=Grad.Buffer();
    double* localMOld=MOld.Buffer();
    double* localMNew=MNew.Buffer();

    const int xLocalSize=Grad.XLocalSize();
    const int yLocalSize=Grad.YLocalSize();
    const int zLocalSize=Grad.ZLocalSize();

    for ( int zLocal=0; zLoca<zLocalSize; ++zLocal) 
    {
        for ( int yLocal=0; yLocal<yLocalSize; ++yLocal) 
        {
            for (int xLocal=0; xLocal<xLocalSize; ++xLocal) 
            {
                const int localIndex=xLocal+
                            yLocal*xLocalSize+
                            zlocal*xLocalSize*yLocalSize;
                if (fabs(localGrad[localIndex]) > gradmaxlocal) {
                    gradmaxlocal=fabs(localGrad[localIndex]);
                }
            }
        }
    }
    

    double grdmaxglob;
    mpi::AllReduce(&grdmaxlocal, &grdmaxglob, 1, mpi::MAX, comm);

    double ratio=alpha/grdmaxglob;
    
    // TODO: Perform a local update
    // (STEP 5)
    
    for ( int zLocal=0; zLoca<zLocalSize; ++zLocal) 
    {
        for ( int yLocal=0; yLocal<yLocalSize; ++yLocal) 
        {
            for (int xLocal=0; xLocal<xLocalSize; ++xLocal) 
            {
                const int localIndex=xLocal+
                            yLocal*xLocalSize+
                            zlocal*xLocalSize*yLocalSize;
                
                localMNew[localIndex]=localMOld[localIndex]+ratio*localGrad[localIndex]; 
            }
        }
    }
    
}    
    
    

void gen_acquisition
(double srcx, double srcy, double srcz, 
 Matrix<double>& xRecv, 
 Matrix<double>& yRecv,
 Matrix<double>& zRecv )
{
    double srcx=0.5;
    double srcy=0.5;
    double srcz=0.03;

    int nxr=1;
    int nyr=1;
    int nzr=1;

    int nrec=nxr*nyr*nzr;

    const double oxr=0.1, oyr=0.1, ozr=0.03; 
    const double dxr=0.02, dyr=0.02, dzr=0.0; 

    xRecv.Resize( nrec, nsrc );
    yRecv.Resize( nrec, nsrc );
    zRecv.Resize( nrec, nsrc );
    for (int isrc=0; isrc<nsrc; isrc++) { 
        for ( int k=0; k<nzr; k++ ) { 
            for (int j=0; j<nyr; j++) {
                for (int i=0; i<nxr; i++) {
		            xRecv.Set( irec, isrc, oxr+i*dxr );
		            yRecv.Set( irec, isrc, oyr+j*dyr );
		            zRecv.Set( irec, isrc, ozr+k*dzr );
                }
            }
        }
    }  
}
