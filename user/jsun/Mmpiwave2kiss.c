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
static sf_complex **cc, *tmp2;

static kiss_fft_cfg *cfg1, *icfg1, *cfg2, *icfg2;
static kiss_fft_cpx *tmp, **ctrace2;

static fftwf_plan cfg=NULL, icfg=NULL;
static ptrdiff_t alloc_local, local_n0, local_0_start, local_n1, local_1_start;

int threads_ok;

int cfft2_init(int pad1           /* padding on the first axis */,
	       int nx,   int ny   /* input data size */, 
	       int *nx2, int *ny2 /* padded data size */,
               int *n_local, int *o_local /* local size & start */)
/*< initialize >*/
{
  int i, nth=1;
  int cpuid;
  ptrdiff_t n[2];
  
  MPI_Comm_rank(MPI_COMM_WORLD, &cpuid);

  fftwf_mpi_init();

  nk = n1 = kiss_fft_next_fast_size(nx*pad1);
  n2 = kiss_fft_next_fast_size(ny);

  n[0]=n2; n[1]=n1;
  //alloc_local = fftwf_mpi_local_size_many_transposed(2, n, 2, FFTW_MPI_DEFAULT_BLOCK, FFTW_MPI_DEFAULT_BLOCK, MPI_COMM_WORLD, &local_n0, &local_0_start, &local_n1, &local_1_start);
  alloc_local = fftwf_mpi_local_size_2d_transposed(n2, n1, MPI_COMM_WORLD, &local_n0, &local_0_start, &local_n1, &local_1_start);

  cc = sf_complexalloc2(n1,local_n0);
  //cc = sf_complexalloc(alloc_local);

  /* kiss-fft */

#ifdef _OPENMP
#pragma omp parallel
  {nth = omp_get_num_threads();}
#endif

  cfg1  = (kiss_fft_cfg *) sf_alloc(nth,sizeof(kiss_fft_cfg));
  icfg1 = (kiss_fft_cfg *) sf_alloc(nth,sizeof(kiss_fft_cfg));
  cfg2  = (kiss_fft_cfg *) sf_alloc(nth,sizeof(kiss_fft_cfg));
  icfg2 = (kiss_fft_cfg *) sf_alloc(nth,sizeof(kiss_fft_cfg));

  for (i=0; i < nth; i++) {
    cfg1[i] = kiss_fft_alloc(n1,0,NULL,NULL);
    icfg1[i]= kiss_fft_alloc(n1,1,NULL,NULL);
    cfg2[i] = kiss_fft_alloc(n2,0,NULL,NULL);
    icfg2[i]= kiss_fft_alloc(n2,1,NULL,NULL);
  }

  ctrace2= (kiss_fft_cpx **) sf_complexalloc2(n2,nth);

  tmp =    (kiss_fft_cpx *) sf_alloc(alloc_local,sizeof(kiss_fft_cpx));
  tmp2= (sf_complex *) tmp;

  /* fftw for transpose */

  cfg = fftwf_mpi_plan_many_transpose(n2,n1,2,
                              FFTW_MPI_DEFAULT_BLOCK,FFTW_MPI_DEFAULT_BLOCK,
                              (float *) tmp,
                              (float *) tmp,
                              MPI_COMM_WORLD,
                              FFTW_MEASURE);

  icfg= fftwf_mpi_plan_many_transpose(n1,n2,2,
                              FFTW_MPI_DEFAULT_BLOCK,FFTW_MPI_DEFAULT_BLOCK,
                              (float *) tmp,
                              (float *) tmp,
                              MPI_COMM_WORLD,
                              FFTW_MEASURE);

  if (NULL == cfg || NULL == icfg) sf_error("FFTW failure.");

  *nx2 = n1;
  *ny2 = n2;
  *n_local = (int) local_n0;
  *o_local = (int) local_0_start;
	
  wt =  1.0/(n1*n2);
	
  return (nk*n2);
}

void cfft2(sf_complex *inp /* [n1*local_n0] */, 
	   sf_complex *out /* [nk*local_n0] */)
/*< 2-D FFT >*/
{
  int i1, i2;
  int ith=0;

  /* FFT centering */
#pragma omp parallel for private(i2,i1) default(shared)
  for (i2=0; i2<local_n0; i2++) {
    for (i1=0; i1<n1; i1++) {
#ifdef SF_HAS_COMPLEX_H
      cc[i2][i1] = (((i2+local_0_start)%2==0)==(i1%2==0))? inp[i2*n1+i1]:-inp[i2*n1+i1];
#else
      cc[i2][i1] = (((i2+local_0_start)%2==0)==(i1%2==0))? inp[i2*n1+i1]:sf_cneg(inp[i2*n1+i1]);
#endif
    }
  }

#ifdef _OPENMP
#pragma omp parallel for private(i2,ith) default(shared)
#endif
  for (i2=0; i2 < local_n0; i2++) {
#ifdef _OPENMP
    ith = omp_get_thread_num();
#endif
    kiss_fft_stride(cfg1[ith],(kiss_fft_cpx *) cc[i2],tmp+i2*nk,1);
  }

  fftwf_execute(cfg);

#ifdef _OPENMP
#pragma omp parallel for private(i1,i2,ith) default(shared)
#endif
  for (i1=0; i1 < local_n1; i1++) {
#ifdef _OPENMP
    ith = omp_get_thread_num();
#endif
    kiss_fft_stride(cfg2[ith],tmp+i1*n2,ctrace2[ith],1);
    for (i2=0; i2<n2; i2++) {
      tmp[i1*n2+i2] = ctrace2[ith][i2];
    }
  }

  fftwf_execute(icfg);

#ifdef _OPENMP
#pragma omp parallel for private(i1,i2) default(shared)
#endif
  for (i2=0; i2<local_n0; i2++) {
    for (i1=0; i1<nk; i1++) {
      out[i2*nk+i1] = tmp2[i2*nk+i1];
    }
  }

}

void icfft2(sf_complex *out /* [n1*local_n0] */, 
	    sf_complex *inp /* [nk*local_n0] */)
/*< 2-D inverse FFT >*/
{
  int i1, i2, ith=0;

#ifdef _OPENMP
#pragma omp parallel for private(i2,ith) default(shared)
#endif
  for (i2=0; i2 < local_n0; i2++) {
#ifdef _OPENMP
    ith = omp_get_thread_num();
#endif
    kiss_fft_stride(icfg1[ith],(kiss_fft_cpx *) inp+i2*nk,tmp+i2*nk,1);
  }

  fftwf_execute(cfg);

#ifdef _OPENMP
#pragma omp parallel for private(i1,i2,ith) default(shared)
#endif
  for (i1=0; i1 < local_n1; i1++) {
#ifdef _OPENMP
    ith = omp_get_thread_num();
#endif
    kiss_fft_stride(icfg2[ith],tmp+i1*n2,ctrace2[ith],1);
    for (i2=0; i2<n2; i2++) {
      tmp[i1*n2+i2] = ctrace2[ith][i2];
    }
  }

  fftwf_execute(icfg);

  /* FFT centering and normalization*/
#pragma omp parallel for private(i2,i1) default(shared)
  for (i2=0; i2<local_n0; i2++) {
    for (i1=0; i1<n1; i1++) {
#ifdef SF_HAS_COMPLEX_H
      out[i2*n1+i1] = ((((i2+local_0_start)%2==0)==(i1%2==0))? wt:-wt) * tmp2[i2*n1+i1];
#else
      out[i2*n1+i1] = sf_crmul(tmp2[i2*n1+i1],((((i2+local_0_start)%2==0)==(i1%2==0))? wt:-wt));
#endif
    }
  }
}

void cfft2_finalize()
/*< clean up fftw >*/
{
  /* make sure everything is back to its pristine state */
  int nth=1,i;

  fftwf_destroy_plan(cfg); cfg=NULL;
  fftwf_destroy_plan(icfg); icfg=NULL;
  fftwf_mpi_cleanup();

#ifdef _OPENMP
#pragma omp parallel
  {nth = omp_get_num_threads();}
#endif
  for (i=0; i < nth; i++) {
    free(cfg1[i]);   cfg1[i]=NULL;
    free(icfg1[i]); icfg1[i]=NULL;
    free(cfg2[i]);   cfg2[i]=NULL;
    free(icfg2[i]); icfg2[i]=NULL;
  }
  free(*cc); free(cc);
  free(tmp);
  free(*ctrace2); free(ctrace2);
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
  Fw = sf_input ("input" );
  Fo = sf_output("output");
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

