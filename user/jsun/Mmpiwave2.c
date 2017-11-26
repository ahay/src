/* Complex 2-D wave propagation (with multi-threaded FFTW3)*/
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
#include <fftw3-mpi.h>
#include <omp.h>

static int n1, n2, nk;
static float wt;
static sf_complex *cc,*dd;
static fftwf_plan cfg=NULL, icfg=NULL;
static ptrdiff_t alloc_local, local_n0, local_0_start;

int threads_ok;

int cfft2_init(int pad1           /* padding on the first axis */,
	       int nx,   int ny   /* input data size */, 
	       int *nx2, int *ny2 /* padded data size */,
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

  nk = n1 = kiss_fft_next_fast_size(nx*pad1);
  n2 = kiss_fft_next_fast_size(ny);

  alloc_local = fftwf_mpi_local_size_2d(n2, n1, MPI_COMM_WORLD, &local_n0, &local_0_start);

  //cc = sf_complexalloc2(n1,n2);
  //dd = sf_complexalloc2(nk,n2);
  cc = sf_complexalloc(alloc_local);
  dd = sf_complexalloc(alloc_local);

  cfg = fftwf_mpi_plan_dft_2d(n2,n1,
                              (fftwf_complex *) cc,
                              (fftwf_complex *) dd,
                              MPI_COMM_WORLD,
                              FFTW_FORWARD, FFTW_MEASURE);

  icfg = fftwf_mpi_plan_dft_2d(n2,n1,
                               (fftwf_complex *) dd, 
                               (fftwf_complex *) cc,
                               MPI_COMM_WORLD,
                               FFTW_BACKWARD, FFTW_MEASURE);

  if (NULL == cfg || NULL == icfg) sf_error("FFTW failure.");

  *nx2 = n1;
  *ny2 = n2;
  *n_local = (int) local_n0;
  *o_local = (int) local_0_start;
	
  wt =  1.0/(n1*n2);
	
  return (nk*n2);
}

void cfft2(sf_complex *inp /* [n1*n2] */, 
	   sf_complex *out /* [nk*n2] */)
/*< 2-D FFT >*/
{
  int i1, i2;

  /* FFT centering */
#pragma omp parallel for private(i2,i1) default(shared)
  for (i2=0; i2<local_n0; i2++) {
    for (i1=0; i1<n1; i1++) {
#ifdef SF_HAS_COMPLEX_H
      cc[i2*n1+i1] = (((i2+local_0_start)%2==0)==(i1%2==0))? inp[i2*n1+i1]:-inp[i2*n1+i1];
#else
      cc[i2*n1+i1] = (((i2+local_0_start)%2==0)==(i1%2==0))? inp[i2*n1+i1]:sf_cneg(inp[i2*n1+i1]);
#endif
    }
  }

  fftwf_execute(cfg);
#pragma omp parallel for private(i2,i1) default(shared)
  for (i2=0; i2<local_n0; i2++) {
    for (i1=0; i1<nk; i1++) {
      out[i2*nk+i1]=dd[i2*nk+i1];
    }
  }
}

void icfft2(sf_complex *out /* [n1*n2] */, 
	    sf_complex *inp /* [nk*n2] */)
/*< 2-D inverse FFT >*/
{
  int i1, i2;

#pragma omp parallel for private(i2,i1) default(shared)
  for (i2=0; i2<local_n0; i2++) {
    for (i1=0; i1<nk; i1++) {
      dd[i2*nk+i1]=inp[i2*nk+i1];
    }
  }
  fftwf_execute(icfg);
    
  /* FFT centering and normalization*/
#pragma omp parallel for private(i2,i1) default(shared)
  for (i2=0; i2<local_n0; i2++) {
    for (i1=0; i1<n1; i1++) {
#ifdef SF_HAS_COMPLEX_H
      out[i2*n1+i1] = ((((i2+local_0_start)%2==0)==(i1%2==0))? wt:-wt) * cc[i2*n1+i1];
#else
      out[i2*n1+i1] = sf_crmul(cc[i2*n1+i1],((((i2+local_0_start)%2==0)==(i1%2==0))? wt:-wt));
#endif
    }
  }
}

void cfft2_finalize()
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
  free(dd);
  free(cc);
}

int main(int argc, char* argv[])
{
  bool verb,sub,os;
  int it,iz,im,ikx,ikz,ix,i,j;     /* index variables */
  int nt,nz,nx, m2, nk, nzx, nz2, nx2, nzx2, n2, pad1,nth;
  sf_complex c,old;

  /* I/O arrays*/
  sf_complex *ww,*curr,*prev,*cwave,*cwavem,**wave,**lt, **rt;
  sf_complex *snap;
  float *rr;
    
  sf_file Fw,Fr,Fo;    /* I/O files */
  sf_axis at,az,ax;    /* cube axes */
  sf_file left, right;

  /*MPI related*/
  int cpuid,numprocs;
  int provided;
  int n_local, o_local;
  int ozx2;
  sf_complex *sendbuf, *recvbuf;
  int *rcounts, *displs;

  //MPI_Init(&argc,&argv);
  MPI_Init_thread(&argc,&argv,MPI_THREAD_FUNNELED,&provided);
  threads_ok = provided >= MPI_THREAD_FUNNELED;

  sf_init(argc,argv);

  MPI_Comm_rank(MPI_COMM_WORLD, &cpuid);
  MPI_Comm_size(MPI_COMM_WORLD, &numprocs);

  if(!sf_getbool("verb",&verb)) verb=false; /* verbosity */
  if(!sf_getbool("os",&os)) os=true; /* one-step flag */
  if (os) {
    sf_warning("One-step wave extrapolation");
    if(!sf_getbool("sub",&sub)) sub=false; /* subtraction flag */
  } else {
    sf_warning("Two-step wave extrapolation");
    if(!sf_getbool("sub",&sub)) sub=true; /* subtraction flag */
  }

  /* setup I/O files */
  Fw = sf_input ("--input" );
  Fo = sf_output("--output");
  Fr = sf_input ("ref");

  if (SF_COMPLEX != sf_gettype(Fw)) sf_error("Need complex input");
  if (SF_FLOAT != sf_gettype(Fr)) sf_error("Need float ref");

  sf_settype(Fo,SF_COMPLEX);

  /* Read/Write axes */
  at = sf_iaxa(Fw,1); nt = sf_n(at); 
  az = sf_iaxa(Fr,1); nz = sf_n(az); 
  ax = sf_iaxa(Fr,2); nx = sf_n(ax); 

  if (cpuid==0) {
    sf_oaxa(Fo,az,1); 
    sf_oaxa(Fo,ax,2); 
    //sf_setn(at,1);
    sf_oaxa(Fo,at,3);
  }
    
  if (!sf_getint("pad1",&pad1)) pad1=1; /* padding factor on the first axis */

#pragma omp parallel
  {
    nth = omp_get_num_threads();
  }
  if (verb) sf_warning(">>>> Using %d threads <<<<<", nth);

  nk = cfft2_init(pad1,nz,nx,&nz2,&nx2,&n_local,&o_local);
  sf_warning("Cpuid=%d,n0=%d,n1=%d,local_n0=%d,local_0_start=%d",cpuid,nz2,nx2,n_local,o_local);

  nzx = nz*nx;
  //  nzx2 = nz2*nx2;
  nzx2 = n_local*nz2;
  ozx2 = o_local*nz2;

  /* propagator matrices */
  left = sf_input("left");
  right = sf_input("right");

  if (!sf_histint(left,"n1",&n2) || n2 != nzx) sf_error("Need n1=%d in left",nzx);
  if (!sf_histint(left,"n2",&m2))  sf_error("Need n2= in left");
    
  if (!sf_histint(right,"n1",&n2) || n2 != m2) sf_error("Need n1=%d in right",m2);
  if (!sf_histint(right,"n2",&n2) || n2 != nk) sf_error("Need n2=%d in right",nk);
  
  lt = sf_complexalloc2(nzx,m2);
  rt = sf_complexalloc2(m2,nk);

  sf_complexread(lt[0],nzx*m2,left);
  sf_complexread(rt[0],m2*nk,right);

  sf_fileclose(left);
  sf_fileclose(right);

  /* read wavelet & reflectivity */
  ww=sf_complexalloc(nt);
  sf_complexread(ww,nt ,Fw);

  rr=sf_floatalloc(nzx); 
  sf_floatread(rr,nzx,Fr);

  curr   = sf_complexalloc(nzx2);
  if (!os) prev = sf_complexalloc(nzx2);
  else prev = NULL;

  cwave  = sf_complexalloc(nzx2);
  cwavem = sf_complexalloc(nzx2);
  wave   = sf_complexalloc2(nzx2,m2);

  for (iz=0; iz < nzx2; iz++) {
    curr[iz] = sf_cmplx(0.,0.);
    if (!os) prev[iz] = sf_cmplx(0.,0.);
  }

  sendbuf = curr;
  if (cpuid==0) {
    snap = sf_complexalloc(nz2*nx2);
    recvbuf = snap;
    rcounts = sf_intalloc(numprocs);
    displs  = sf_intalloc(numprocs);
  } else {
    snap = NULL;
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
    cfft2(curr,cwave);

    for (im = 0; im < m2; im++) {
      //for (ik = 0; ik < nk; ik++) {
      for (ikx = 0; ikx < n_local; ikx++) {
        for (ikz = 0; ikz < nz2; ikz++) {
          i = ikz + (o_local+ikx)*nz2;
          j = ikz + ikx*nz2;
#ifdef SF_HAS_COMPLEX_H
          cwavem[j] = cwave[j]*rt[i][im];
#else
          cwavem[j] = sf_cmul(cwave[j],rt[i][im]);
#endif
        }
      }
      icfft2(wave[im],cwavem);
    }

    for (ix = 0; ix < n_local && (ix+o_local)<nx ; ix++) {
      for (iz=0; iz < nz; iz++) {
        i = iz + (o_local+ix)*nz;  /* original grid */
        j = iz + ix*nz2; /* padded grid */
#ifdef SF_HAS_COMPLEX_H
        c = ww[it] * rr[i]; // source term
#else
        c = sf_crmul(ww[it], rr[i]); // source term
#endif
        if (sub) c += curr[j];
        if (!os) {
          old = curr[j];
#ifdef SF_HAS_COMPLEX_H
          c += sub? (old-prev[j]) : -prev[j];
#else
          c = sf_cadd(c,sub? sf_csub(old,prev[j]) : sf_cneg(prev[j]));
#endif
          prev[j] = old;
        }
        for (im = 0; im < m2; im++) {
#ifdef SF_HAS_COMPLEX_H
          c += lt[im][i]*wave[im][j];
#else
          c += sf_cmul(lt[im][i], wave[im][j]);
#endif
        }

        curr[j] = c;
      }
    }

    /* write wavefield to output */
    MPI_Gatherv(sendbuf, nzx2, MPI_COMPLEX, recvbuf, rcounts, displs, MPI_COMPLEX, 0, MPI_COMM_WORLD);

    if (cpuid==0) {
      for (ix = 0; ix < nx; ix++)
        sf_complexwrite(snap+ix*nz2,nz,Fo);
    }

  }

  if(verb) sf_warning("."); 
  cfft2_finalize();

  MPI_Finalize();
  exit (0);
}

