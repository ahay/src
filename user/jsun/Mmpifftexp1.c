/* 2-D FFT-based prestack exploding reflector modeling/migration  */
/*
  Copyright (C) 2010 University of Texas at Austin
  
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

void mcfft3(float *inp /* [n1*n2*n3] */, 
            sf_complex *out /* [nk*n2*n3] */)
/*< 3-D FFT >*/
{
  int i1, i2, i3;

  /* FFT centering */
#pragma omp parallel for private(i3,i2,i1) default(shared)
  for (i3=0; i3<local_n0; i3++) {
    for (i2=0; i2<n2; i2++) {
      for (i1=0; i1<n1; i1++) {
        cc[(i3*n2+i2)*n1+i1] = sf_cmplx(((((i3+local_0_start)%2==0)==(i2%2==0))==(i1%2==0))? inp[(i3*n2+i2)*n1+i1]:-inp[(i3*n2+i2)*n1+i1],0.);
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

void imcfft3(float *out /* [n1*n2*n3] */, 
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
        out[(i3*n2+i2)*n1+i1] = (((((i3+local_0_start)%2==0)==(i2%2==0))==(i1%2==0))? wt:-wt)*crealf(cc[(i3*n2+i2)*n1+i1]);
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
    bool mig, sub;
    int it, nt, ix, nx, iz, nz, nx2, nz2, nzx, nzx2, ih, nh, nh2;
    int im, i, j, m2, it1, it2, its, ikz, ikx, ikh, n2, nk, snap;
    float dt, dx, dz, c, old, dh;
    float *curr, *prev, **img, **dat, **lft, **rht, **wave;
    sf_complex *cwave, *cwavem;
    sf_file data, image, left, right, snaps;

    /*MPI related*/
    int cpuid,numprocs;
    int provided;
    int n_local, o_local;
    int ozx2;
    float *sendbuf, *recvbuf, *wave_all;
    int *rcounts, *displs;

    MPI_Init_thread(&argc,&argv,MPI_THREAD_FUNNELED,&provided);
    threads_ok = provided >= MPI_THREAD_FUNNELED;

    sf_init(argc,argv);

    MPI_Comm_rank(MPI_COMM_WORLD, &cpuid);
    MPI_Comm_size(MPI_COMM_WORLD, &numprocs);

    if (!sf_getbool("mig",&mig)) mig=false;
    /* if n, modeling; if y, migration */

    if (!sf_getint("snap",&snap)) snap=0;
    /* interval for snapshots */

    snaps = (snap > 0)? sf_output("snaps"): NULL;
    /* (optional) snapshot file */

    if (mig) { /* migration */
	data = sf_input("input");
	image = sf_output("output");

	if (!sf_histint(data,"n1",&nh)) sf_error("No n1=");
	if (!sf_histfloat(data,"d1",&dh)) sf_error("No d1=");

	if (!sf_histint(data,"n2",&nx)) sf_error("No n2=");
	if (!sf_histfloat(data,"d2",&dx)) sf_error("No d2=");

	if (!sf_histint(data,"n3",&nt)) sf_error("No n3=");
	if (!sf_histfloat(data,"d3",&dt)) sf_error("No d3=");

	if (!sf_getint("nz",&nz)) sf_error("Need nz=");
	/* time samples (if migration) */
	if (!sf_getfloat("dz",&dz)) sf_error("Need dz=");
	/* time sampling (if migration) */
        
        if (cpuid==0) {
	sf_putint(image,"o1",0.);
	sf_putint(image,"n1",nz);
	sf_putfloat(image,"d1",dz);
	sf_putstring(image,"label1","Depth");
	sf_putint(image,"o2",0.);
	sf_putint(image,"n2",nx);
	sf_putfloat(image,"d2",dx);
	sf_putstring(image,"label2","Midpoint");
	sf_putint(image,"n3",1); /* stack for now */
        }
    } else { /* modeling */
	image = sf_input("input");
	data = sf_output("output");

	if (!sf_histint(image,"n1",&nz)) sf_error("No n1=");
	if (!sf_histfloat(image,"d1",&dz)) sf_error("No d1=");

	if (!sf_histint(image,"n2",&nx)) sf_error("No n2=");
	if (!sf_histfloat(image,"d2",&dx)) sf_error("No d2=");

	if (!sf_getint("nt",&nt)) sf_error("Need nt=");
	/* time samples (if modeling) */
	if (!sf_getfloat("dt",&dt)) sf_error("Need dt=");
	/* time sampling (if modeling) */

	if (!sf_getint("nh",&nh)) sf_error("Need nh=");
        /* offset samples (if modeling) */
	if (!sf_getfloat("dh",&dh)) sf_error("Need dh=");
	/* offset sampling (if modeling) */

        if (cpuid==0) {
	sf_putint(data,"n1",nh);
	sf_putfloat(data,"d1",dh);
	sf_putstring(data,"label1","Half-Offset");
	sf_putint(data,"o2",0.);
	sf_putint(data,"n2",nx);
	sf_putfloat(data,"d2",dx);
	sf_putstring(data,"label2","Midpoint");
	sf_putint(data,"n3",nt);
	sf_putfloat(data,"d3",dt);
	sf_putstring(data,"label3","Time");
	sf_putstring(data,"unit3","s");
        }
    }

    if (cpuid==0) {
    if (NULL != snaps) {
      sf_putint(snaps,"n1",nh);
      sf_putfloat(snaps,"d1",dh);
      sf_putstring(snaps,"label1","Half-Offset");

      sf_putint(snaps,"n2",nx);
      sf_putfloat(snaps,"d2",dx);
      sf_putstring(snaps,"label2","Midpoint");

      sf_putint(snaps,"n3",nz);
      sf_putfloat(snaps,"d3",dz);
      sf_putstring(snaps,"label3","Depth");

      sf_putint(snaps,"n4",nt/snap);
      sf_putfloat(snaps,"d4",dt*snap);
      if (mig) {
        sf_putfloat(snaps,"o4",(nt-1)*dt);
      } else {
        sf_putfloat(snaps,"o4",0.);
      }
      sf_putstring(snaps,"label4","Time");
    }
    }

    nk = mcfft3_init(1,nh,nx,nz,&nh2,&nx2,&nz2,&n_local,&o_local);
    sf_warning("Cpuid=%d,n2=%d,n1=%d,n0=%d,local_n0=%d,local_0_start=%d",cpuid,nh2,nx2,nz2,n_local,o_local);
    if (cpuid==0)
      if (o_local!=0) sf_error("Cpuid and o_local inconsistant!");

    nzx = nz*nx*nh;
    //nzx2 = nz2*nx2*nh2;
    nzx2 = n_local*nx2*nh2;
    ozx2 = o_local*nx2*nh2;

    img = sf_floatalloc2(nz,nx);
    dat = sf_floatalloc2(nh,nx);

    /* propagator matrices */
    left = sf_input("left");
    right = sf_input("right");

    if (!sf_histbool(left,"sub",&sub) && !sf_getbool("sub",&sub)) sub=true;
    /* if -1 is included in the matrix */

    if (!sf_histint(left,"n1",&n2) || n2 != nzx) sf_error("Need n1=%d in left",nzx);
    if (!sf_histint(left,"n2",&m2))  sf_error("No n2= in left");
    
    if (!sf_histint(right,"n1",&n2) || n2 != m2) sf_error("Need n1=%d in right",m2);
    if (!sf_histint(right,"n2",&n2) || n2 != nk) sf_error("Need n2=%d in right",nk);
 
    lft = sf_floatalloc2(nzx,m2);
    rht = sf_floatalloc2(m2,nk);

    sf_floatread(lft[0],nzx*m2,left);
    sf_floatread(rht[0],m2*nk,right);

    curr = sf_floatalloc(nzx2);
    prev = sf_floatalloc(nzx2);

    cwave  = sf_complexalloc(nk);
    cwavem = sf_complexalloc(nk);
    wave = sf_floatalloc2(nzx2,m2);

    for (iz=0; iz < nzx2; iz++) {
	curr[iz]=0.;
	prev[iz]=0.;
    }

    sendbuf = prev;
    if (cpuid==0) {
      wave_all = sf_floatalloc(nh2*nx2*nz2);
      recvbuf = wave_all;
      rcounts = sf_intalloc(numprocs);
      displs  = sf_intalloc(numprocs);
    } else {
      wave_all = NULL;
      recvbuf = NULL;
      rcounts = NULL;
      displs = NULL;
    }

    MPI_Gather(&nzx2, 1, MPI_INT, rcounts, 1, MPI_INT, 0, MPI_COMM_WORLD);
    MPI_Gather(&ozx2, 1, MPI_INT, displs, 1, MPI_INT, 0, MPI_COMM_WORLD);

    if (mig) { /* migration */
	/* step backward in time */
	it1 = nt-1;
	it2 = -1;
	its = -1;	
    } else { /* modeling */
	sf_floatread(img[0],nz*nx,image);

	/* transpose and initialize at zero offset */
	for (iz=0; iz < n_local && (iz+o_local)<nz; iz++) {
	    for (ix=0; ix < nx; ix++) {
		curr[nh2*(ix+iz*nx2)]=img[ix][iz+o_local];
	    }
	}
	
	/* step forward in time */
	it1 = 0;
	it2 = nt;
	its = +1;
    }

    /* time stepping */
    for (it=it1; it != it2; it += its) {
	sf_warning("it=%d;",it);

	if (mig) { /* migration <- read data */
	    sf_floatread(dat[0],nx*nh,data);
	} else {
	    for (ix=0; ix < nx; ix++) {
		for (ih=0; ih < nh; ih++) {
		    dat[ix][ih] = 0.;
		}
	    }
	}
	
	if (NULL != snaps && 0 == it%snap) {
          MPI_Gatherv(sendbuf, nzx2, MPI_FLOAT, recvbuf, rcounts, displs, MPI_FLOAT, 0, MPI_COMM_WORLD);
          if (cpuid==0) {
            for (iz = 0; iz < nz; iz++)
              for (ix = 0; ix < nx; ix++)
                sf_floatwrite(wave_all+nh2*(ix+nx2*iz),nh,snaps);
          }
        }

	/* at z=0 */
        if (cpuid==0) {
	for (ix=0; ix < nx; ix++) {
	    for (ih=0; ih < nh; ih++) {
		if (mig) {
		    curr[ix*nh2+ih] += dat[ix][ih];
		} else {
		    dat[ix][ih] = curr[ix*nh2+ih];
		}
	    }
	}
        }

	/* matrix multiplication */
	mcfft3(curr,cwave);

	for (im = 0; im < m2; im++) {
          //for (ik = 0; ik < nk; ik++) {
          for (ikz = 0; ikz < n_local; ikz++) {
            for (ikx = 0; ikx < nx2; ikx++) {
              for (ikh = 0; ikh < nh2; ikh++) {
                i = ikh + ikx*nh2 + (o_local+ikz)*nx2*nh2;
                j = ikh + ikx*nh2 + ikz*nx2*nh2;
#ifdef SF_HAS_COMPLEX_H
		cwavem[j] = cwave[j]*rht[i][im];
#else
		cwavem[j] = sf_crmul(cwave[j],rht[i][im]);
#endif
              }
            }
          }
          imcfft3(wave[im],cwavem);
	}

        //#ifdef _OPENMP
        //#pragma omp parallel for private(ix,iz,ih,i,j,im,old,c) shared(curr,prev,lft,wave)
        //#endif
        for (iz=0; iz < n_local && (iz+o_local)<nz; iz++) {
	    for (ix = 0; ix < nx; ix++) {
		for (ih=0; ih < nh; ih++) {	
                    i = ih + ix*nh + (o_local+iz)*nx*nh;  /* original grid */
                    j = ih + ix*nh2+ iz*nx2*nh2; /* padded grid */
		
		    old = curr[j];

		    c = sub? 2*old: 0.0f;

		    c -= prev[j];

		    prev[j] = old;

		    for (im = 0; im < m2; im++) {
			c += lft[im][i]*wave[im][j];
		    }
		    
		    curr[j] = c;
		}
	    }
	}
	
	if (!mig) { /* modeling -> write out data */
          if (cpuid==0)
	    sf_floatwrite(dat[0],nx*nh,data);
	}
    }
    sf_warning(".");

    if (mig) {
      sendbuf = curr;
      MPI_Gatherv(sendbuf, nzx2, MPI_FLOAT, recvbuf, rcounts, displs, MPI_FLOAT, 0, MPI_COMM_WORLD);

      if (cpuid==0) {
        /* transpose */
        for (iz=0; iz < nz; iz++) {
          for (ix=0; ix < nx; ix++) {
            img[ix][iz] = wave_all[nh2*(ix+iz*nx2)];
          }
        }
	sf_floatwrite(img[0],nz*nx,image);
      }

    }

    mcfft3_finalize();

    MPI_Finalize();

    exit(0);
}
