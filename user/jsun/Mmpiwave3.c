/* Simple 3-D lowrank onestep wave propagation */
/*
  Copyright (C) 2009 University of Texas at Austin
  
  This program is free software; you can redistribute it and/or modify
  it under the terms of the GNU General Public License as published by
  the Free Software Foundation; either version 2 of the License, or
  (at your option) any later version.
  
  This program is distributed in the hope that it will be useful,
  but WITHOUT ANY WARRANTY; without even the implied warranty of
  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
  GNU General Public License for more details.
  
  You should have received a copy of the GNU General Public License
  along with this program; if not, write to the Free Software
  Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA
*/
#include <rsf.h>
#include <omp.h>
#include <fftw3-mpi.h>

int threads_ok;

static int n1, n2, n3, nk;
static float wt;
static sf_complex *cc=NULL;
static fftwf_plan cfg=NULL, icfg=NULL;
static ptrdiff_t alloc_local, local_n0, local_0_start;

int mcfft3_init(int pad1           /* padding on the first axis */,
	       int nx,   int ny,  int nz   /* input data size */, 
	       int *nx2, int *ny2, int *nz2 /* padded data size */,
               int *n_local, int *o_local /* local size & start */)
/*< initialize >*/
{
  int cpuid;
  MPI_Comm_rank(MPI_COMM_WORLD, &cpuid);

  if (threads_ok) threads_ok = fftwf_init_threads();

  fftwf_mpi_init();

  if (false)
    sf_warning("Using threaded FFTW3! \n");
  if (threads_ok)
    fftwf_plan_with_nthreads(omp_get_max_threads());

  /* axis 1 */
  nk = n1 = kiss_fft_next_fast_size(nx*pad1);
  /* axis 2 */
  n2 = kiss_fft_next_fast_size(ny);
  /* axis 3 */
  n3 = kiss_fft_next_fast_size(nz);

  alloc_local = fftwf_mpi_local_size_3d(n3, n2, n1, MPI_COMM_WORLD, &local_n0, &local_0_start);

  //cc = sf_complexalloc3(n1,n2,n3);
  cc = sf_complexalloc(alloc_local);

  cfg  = fftwf_mpi_plan_dft_3d(n3,n2,n1,
                               (fftwf_complex *) cc,
                               (fftwf_complex *) cc,
                               MPI_COMM_WORLD,
                               FFTW_FORWARD, FFTW_MEASURE);

  icfg = fftwf_mpi_plan_dft_3d(n3,n2,n1,
                               (fftwf_complex *) cc, 
                               (fftwf_complex *) cc,
                               MPI_COMM_WORLD,
                               FFTW_BACKWARD, FFTW_MEASURE);

  if (NULL == cfg || NULL == icfg) sf_error("FFTW failure.");

  *nx2 = n1;
  *ny2 = n2;
  *nz2 = n3;
  *n_local = (int) local_n0;
  *o_local = (int) local_0_start;
	
  wt =  1.0/(n3*n2*n1);

  return (nk*n2*n3);
}

void mcfft3(sf_complex *inp /* [n1*n2*n3] */, 
	   sf_complex *out /* [nk*n2*n3] */)
/*< 3-D FFT >*/
{
  int i1, i2, i3;

  /* FFT centering */
#pragma omp parallel for private(i3,i2,i1) default(shared)
  for (i3=0; i3<local_n0; i3++) {
    for (i2=0; i2<n2; i2++) {
      for (i1=0; i1<n1; i1++) {
#ifdef SF_HAS_COMPLEX_H
        cc[(i3*n2+i2)*n1+i1] = ((((i3+local_0_start)%2==0)==(i2%2==0))==(i1%2==0))? inp[(i3*n2+i2)*n1+i1]:-inp[(i3*n2+i2)*n1+i1];
#else
        cc[(i3*n2+i2)*n1+i1] = ((((i3+local_0_start)%2==0)==(i2%2==0))==(i1%2==0))? inp[(i3*n2+i2)*n1+i1]:sf_cneg(inp[(i3*n2+i2)*n1+i1]);
#endif
      }
    }
  }

  fftwf_execute(cfg);

#pragma omp parallel for private(i3,i2,i1) default(shared)
  for (i3=0; i3<local_n0; i3++)
    for (i2=0; i2<n2; i2++)
      for (i1=0; i1<nk; i1++)
        out[(i3*n2+i2)*n1+i1]=cc[(i3*n2+i2)*n1+i1];

}

void imcfft3(sf_complex *out /* [n1*n2*n3] */, 
	    sf_complex *inp /* [nk*n2*n3] */)
/*< 3-D inverse FFT >*/
{
  int i1, i2, i3;

#pragma omp parallel for private(i3,i2,i1) default(shared)
  for (i3=0; i3<local_n0; i3++) {
    for (i2=0; i2<n2; i2++) {
      for (i1=0; i1<nk; i1++) {
        cc[(i3*n2+i2)*n1+i1]=inp[(i3*n2+i2)*n1+i1];
      }
    }
  }
  fftwf_execute(icfg);
    
  /* FFT centering and normalization*/
#pragma omp parallel for private(i3,i2,i1) default(shared)
  for (i3=0; i3<local_n0; i3++) {
    for (i2=0; i2<n2; i2++) {
      for (i1=0; i1<n1; i1++) {
#ifdef SF_HAS_COMPLEX_H
        out[(i3*n2+i2)*n1+i1] = (((((i3+local_0_start)%2==0)==(i2%2==0))==(i1%2==0))? wt:-wt)*cc[(i3*n2+i2)*n1+i1];
#else
	out[(i3*n2+i2)*n1+i1] = sf_crmul(cc[(i3*n2+i2)*n1+i1],(((((i3+local_0_start)%2==0)==(i2%2==0))==(i1%2==0))? wt:-wt));
#endif
      }
    }
  }
}

void mcfft3_finalize()
/*< clean up fftw >*/
{
  /* make sure everything is back to its pristine state */

  fftwf_destroy_plan(cfg);
  fftwf_destroy_plan(icfg);
  fftwf_mpi_cleanup();
  //fftwf_cleanup_threads();
  //fftwf_cleanup();
  cfg=NULL;
  icfg=NULL;
  free(cc);
}

int main(int argc, char* argv[])
{
  bool verb;        
  int it,iz,im,ikz,ikx,iky,ix,iy,i,j,snap;     /* index variables */
  int nt,nz,nx,ny, m2, nk, nzx, nz2, nx2, ny2, nzx2, n2, pad1;
  float dt;
  sf_complex c;

  float  *rr;      /* I/O arrays*/
  sf_complex *cwave, *cwavem, *ww;
  sf_complex **wave, *curr;
  float *rcurr, *rcurr_all;

  sf_file Fw,Fr,Fo;    /* I/O files */
  sf_axis at,az,ax,ay;    /* cube axes */

  sf_complex **lt, **rt;
  sf_file left, right, snaps;

  /*MPI related*/
  int cpuid,numprocs;
  int provided;
  int n_local, o_local;
  int ozx2;
  float *sendbuf, *recvbuf;
  int *rcounts, *displs;

  MPI_Init_thread(&argc,&argv,MPI_THREAD_FUNNELED,&provided);
  threads_ok = provided >= MPI_THREAD_FUNNELED;

  sf_init(argc,argv);

  MPI_Comm_rank(MPI_COMM_WORLD, &cpuid);
  MPI_Comm_size(MPI_COMM_WORLD, &numprocs);

  if(!sf_getbool("verb",&verb)) verb=true; /* verbosity */

  /* setup I/O files */
  Fw = sf_input ("--input" );
  Fo = sf_output("--output");
  Fr = sf_input ("ref");

  /* Read/Write axes */
  at = sf_iaxa(Fw,1); nt = sf_n(at); dt = sf_d(at); 
  az = sf_iaxa(Fr,1); nz = sf_n(az); 
  ax = sf_iaxa(Fr,2); nx = sf_n(ax); 
  ay = sf_iaxa(Fr,3); ny = sf_n(ay); 

  if (!sf_getint("pad1",&pad1)) pad1=1; /* padding factor on the first axis */

  if (!sf_getint("snap",&snap)) snap=0;
  /* interval for snapshots */
    
  if (cpuid==0) {

    sf_oaxa(Fo,az,1); 
    sf_oaxa(Fo,ax,2);
    sf_oaxa(Fo,ay,3);
    
    sf_settype(Fo,SF_FLOAT);

    if (snap > 0) {
      snaps = sf_output("snaps");
      /* (optional) snapshot file */
	
      sf_oaxa(snaps,az,1); 
      sf_oaxa(snaps,ax,2);
      sf_oaxa(snaps,ay,3);
      sf_oaxa(snaps,at,4);
      sf_settype(snaps,SF_FLOAT);
      sf_putint(snaps,"n4",nt/snap);
      sf_putfloat(snaps,"d4",dt*snap);
      sf_putfloat(snaps,"o4",0.);
    } else {
      snaps = NULL;
    }

  }

  //nk = cfft3_init(pad1,nz,nx,ny,&nz2,&nx2,&ny2);
  //n_local = ny2;
  //o_local = 0;
  nk = mcfft3_init(pad1,nz,nx,ny,&nz2,&nx2,&ny2,&n_local,&o_local);
  sf_warning("Cpuid=%d,n2=%d,n1=%d,n0=%d,local_n0=%d,local_0_start=%d",cpuid,nz2,nx2,ny2,n_local,o_local);

  nzx = nz*nx*ny;
  //nzx2 = nz2*nx2*ny2;
  nzx2 = n_local*nz2*nx2;
  ozx2 = o_local*nz2*nx2;

  /* propagator matrices */
  left = sf_input("left");
  right = sf_input("right");

  if (!sf_histint(left,"n1",&n2) || n2 != nzx) sf_error("Need n1=%d in left",nzx);
  if (!sf_histint(left,"n2",&m2))  sf_error("Need n2=%d in left",m2);
    
  if (!sf_histint(right,"n1",&n2) || n2 != m2) sf_error("Need n1=%d in right",m2);
  if (!sf_histint(right,"n2",&n2) || n2 != nk) sf_error("Need n2=%d in right",nk);
 
  lt = sf_complexalloc2(nzx,m2);
  rt = sf_complexalloc2(m2,nk);

  sf_complexread(lt[0],nzx*m2,left);
  sf_complexread(rt[0],m2*nk,right);

  /* read wavelet & reflectivity */
  ww=sf_complexalloc(nt);  sf_complexread(ww,nt ,Fw);
  rr=sf_floatalloc(nzx); sf_floatread(rr,nzx,Fr);

  curr = sf_complexalloc(nzx2);
  rcurr= sf_floatalloc(nzx2);

  cwave  = sf_complexalloc(nzx2);
  cwavem = sf_complexalloc(nzx2);
  wave = sf_complexalloc2(nzx2,m2);

  //icfft3_allocate(cwavem);

  for (iz=0; iz < nzx2; iz++) {
    curr[iz]=sf_cmplx(0.,0.);
    rcurr[iz]=0.;
  }

  sendbuf = rcurr;
  if (cpuid==0) {
    rcurr_all = sf_floatalloc(nz2*nx2*ny2);
    recvbuf = rcurr_all;
    rcounts = sf_intalloc(numprocs);
    displs  = sf_intalloc(numprocs);
  } else {
    rcurr_all = NULL;
    recvbuf = NULL;
    rcounts = NULL;
    displs = NULL;
  }

  MPI_Gather(&nzx2, 1, MPI_INT, rcounts, 1, MPI_INT, 0, MPI_COMM_WORLD);
  MPI_Gather(&ozx2, 1, MPI_INT, displs, 1, MPI_INT, 0, MPI_COMM_WORLD);

  /* MAIN LOOP */
  for (it=0; it<nt; it++) {
    if(verb) sf_warning("it=%d;",it);

    /* matrix multiplication */
    mcfft3(curr,cwave);

    for (im = 0; im < m2; im++) {
      for (iky = 0; iky < n_local; iky++) {
        for (ikx = 0; ikx < nx2; ikx++) {
          for (ikz = 0; ikz < nz2; ikz++) {
            i = ikz + ikx*nz2 + (o_local+iky)*nx2*nz2;
            j = ikz + ikx*nz2 + iky*nx2*nz2;
#ifdef SF_HAS_COMPLEX_H
            cwavem[j] = cwave[j]*rt[i][im];
#else
            cwavem[j] = sf_cmul(cwave[j],rt[i][im]);
#endif
          }
        }
      }
      imcfft3(wave[im],cwavem);
    }

    for (iy = 0; iy < n_local && (iy+o_local)<ny; iy++) {
      for (ix = 0; ix < nx; ix++) {
        for (iz=0; iz < nz; iz++) {
          i = iz + ix*nz + (o_local+iy)*nx*nz;  /* original grid */
          j = iz + ix*nz2+ iy*nx2*nz2; /* padded grid */
#ifdef SF_HAS_COMPLEX_H		
          c = ww[it] * rr[i];
#else
          c = sf_crmul(ww[it],rr[i]);
#endif

          for (im = 0; im < m2; im++) {
#ifdef SF_HAS_COMPLEX_H
            c += lt[im][i]*wave[im][j];
#else
            c += sf_cmul(lt[im][i],wave[im][j]);
#endif
          }
		    
          curr[j] = c;
          rcurr[j]= crealf(c);
        }
      }
    }

    /* output movie */
    if (NULL != snaps && 0 == it%snap) {
      MPI_Gatherv(sendbuf, nzx2, MPI_FLOAT, recvbuf, rcounts, displs, MPI_FLOAT, 0, MPI_COMM_WORLD);

      if (cpuid==0) {
        for (iy = 0; iy < ny; iy++)
          for (ix = 0; ix < nx; ix++)
            sf_floatwrite(rcurr_all+nz2*(ix+nx2*iy),nz,snaps);
      }
    }

  }
  if(verb) sf_warning(".");    
	    	
  /* write wavefield to output */
  MPI_Gatherv(sendbuf, nzx2, MPI_FLOAT, recvbuf, rcounts, displs, MPI_FLOAT, 0, MPI_COMM_WORLD);
  if (cpuid==0) {
    for (iy = 0; iy < ny; iy++)
      for (ix = 0; ix < nx; ix++)
        sf_floatwrite(rcurr_all+nz2*(ix+nx2*iy),nz,Fo);
  }
    
  mcfft3_finalize();

  MPI_Finalize();
  exit (0);
}
