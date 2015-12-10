/* 2-D two-components wavefield modeling using original elastic displacement wave equation in TTI media by lowrank method. */
/*
  Copyright (C) 2014 UT Austin
     
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
#ifdef _OPENMP
#include <omp.h>
#endif

/* head files aumatically produced from *.c */
#include "ricker1.h"
#include "puthead.h"
#include "fft2w.h"
#include "absorb.h"

/*****************************************************************************************/
int main(int argc, char* argv[])
{
  int nx,nz;
  int iz, ix, ik, im, it, i, j;
  int nt, src, isx, isz; /* model/data dimension and source location */
  int nzx,nzx2,nz2,nx2,nk,nkz,nkx,pad1,n2,m2p,m2s; /*fft related*/

  bool cmplx,sep;
  int nth;

  bool abc;
  int nbt,nbb,nbl,nbr; /* abc width */
  float ct,cb,cl,cr; /* decay parameter */
  
  float f0, t, t0, dx, dz, dt, A; /* time and source related */
  float kz0, kx0, dkz, dkx;
  float *kx,*kz,kmod;
  float **plt, **prt, **slt, **srt; /* lowrank matrices */
  float **wavep, **waves, *uz, *ux, *upz1, *upx1, *usz1, *usx1, *upz2, *upx2, *usz2, *usx2; /* wavefield array */
  float c, old;
  sf_complex *cwavepz, *cwavepx, *cwavesz, *cwavesx, *cwavex, *cwavez, *cwavem;

  float x0, z0; /* starting location */

  sf_file Fvp0, Fpl, Fpr, Fsl, Fsr; /* input files, lowrank matrices */
  sf_file Fo1, Fo2; /* output files */
  sf_file Fpz, Fpx, Fsz, Fsx;
  sf_axis az, ax; /* axes */

  sf_init(argc,argv);

  /* setup I/O files */
  Fvp0 = sf_input ("in");  /* vp0 using standard input, just to get model dimension */
  Fpl = sf_input ("pleft");   /* left matrix for P */
  Fpr = sf_input ("pright");  /* right matrix for P */
  Fsl = sf_input ("sleft");   /* left matrix for Sv */
  Fsr = sf_input ("sright");  /* right matrix for Sv */

  /* Read/Write axes */
  az = sf_iaxa(Fvp0,1); nz = sf_n(az); dz = sf_d(az);
  ax = sf_iaxa(Fvp0,2); nx = sf_n(ax); dx = sf_d(ax);
  x0=sf_o(ax);
  z0=sf_o(az);

  /* paramters for the commandline */

  /* time sampling */
  if (!sf_getint("nt",&nt)) nt=301;
  if (!sf_getfloat("dt",&dt)) dt=0.001;
  /* absorbing boundary condition*/
  if(!sf_getbool("abc",&abc)) abc=false; /* absorbing flag */
  if (abc) {
    if(!sf_getint("nbt",&nbt)) sf_error("Need nbt!");
    if(!sf_getint("nbb",&nbb)) nbb = nbt;
    if(!sf_getint("nbl",&nbl)) nbl = nbt;
    if(!sf_getint("nbr",&nbr)) nbr = nbt;
    if(!sf_getfloat("ct",&ct)) sf_error("Need ct!");
    if(!sf_getfloat("cb",&cb)) cb = ct;
    if(!sf_getfloat("cl",&cl)) cl = ct;
    if(!sf_getfloat("cr",&cr)) cr = ct;
  }
  /* source definition */
  if (!sf_getint("isx",&isx)) isx=nx/2;
  if (!sf_getint("isz",&isz)) isz=nz/2;
  if (!sf_getfloat("t0",&t0)) t0=0.04;   /*  wavelet time lag */
  if (!sf_getfloat("f0",&f0)) f0=30.0;   /*  wavelet peak freq */
  if (!sf_getfloat("A",&A)) A=1.0;       /*  wavelet amplitude */
  /* flags and misc */
  if (!sf_getbool("cmplx",&cmplx)) cmplx=false; /* use complex FFT */
  if (!sf_getbool("sep",&sep)) sep=false; /* output separated wavefields */
  if (!sf_getint("pad1",&pad1)) pad1=1; /* padding factor on the first axis */
  if (!sf_getint("src",&src)) src=1; /* source mode: 1 - exploding force; 2 - equil-energy force */

  sf_warning("nt=%d dt=%f",nt,dt);

#ifdef _OPENMP
#pragma omp parallel
    {
      nth = omp_get_num_threads();
      /* omp_set_num_threads(nth); */
    }
    sf_warning(">>>> Using %d threads <<<<<", nth);
#endif

  /* wave modeling space */
  sf_warning("x0=%f z0=%f dx=%f dz=%f",x0,z0,dx,dz);
  sf_warning("nx=%d nz=%d", nx,nz);

  nk = fft2_init(cmplx,pad1,nz,nx,&nz2,&nx2);
  nzx = nz*nx;
  nzx2 = nz2*nx2;
  kz0 = cmplx ? -SF_PI/dz : 0;
  kx0 = -SF_PI/dx;
  dkz = 2*SF_PI/(dz*nz2);
  dkx = 2*SF_PI/(dx*nx2);
  nkz = cmplx ? nz2 : 1+nz2/2;
  nkx = nx2;

  kx=sf_floatalloc(nk);
  kz=sf_floatalloc(nk);
#ifdef _OPENMP
#pragma omp parallel for default(shared) private(ix,iz,ik,kmod)
#endif
  for (ix=0; ix < nkx; ix++) {
    for (iz=0; iz < nkz; iz++) {
      ik = iz + ix*nkz;
      kz[ik] = kz0+iz*dkz;
      kx[ik] = kx0+ix*dkx;
      kmod = sqrtf(kz[ik]*kz[ik]+kx[ik]*kx[ik]);
      kz[ik] /= kmod;
      kx[ik] /= kmod;
    }
  }

  /* propagator matrices */
  if (!sf_histint(Fpl,"n1",&n2) || n2 != nzx) sf_error("Need n1=%d in left for P wave",nzx);
  if (!sf_histint(Fpl,"n2",&m2p))  sf_error("Need n2= in left for P wave");
  if (!sf_histint(Fpr,"n1",&n2) || n2 != m2p) sf_error("Need n1=%d in right for P wave",m2p);
  if (!sf_histint(Fpr,"n2",&n2) || n2 != nk) sf_error("Need n2=%d in right for P wave",nk);

  if (!sf_histint(Fsl,"n1",&n2) || n2 != nzx) sf_error("Need n1=%d in left for P wave",nzx);
  if (!sf_histint(Fsl,"n2",&m2s))  sf_error("Need n2= in left for P wave");
  if (!sf_histint(Fsr,"n1",&n2) || n2 != m2s) sf_error("Need n1=%d in right for P wave",m2s);
  if (!sf_histint(Fsr,"n2",&n2) || n2 != nk) sf_error("Need n2=%d in right for P wave",nk);

  plt = sf_floatalloc2(nzx,m2p);
  prt = sf_floatalloc2(m2p,nk);
  slt = sf_floatalloc2(nzx,m2s);
  srt = sf_floatalloc2(m2s,nk);  

  sf_floatread(plt[0],nzx*m2p,Fpl);
  sf_floatread(prt[0],m2p*nk,Fpr);
  sf_floatread(slt[0],nzx*m2s,Fsl);
  sf_floatread(srt[0],m2s*nk,Fsr);

  /* setup I/O files */
  Fo1 = sf_output("out"); /* Elastic-wave x-component */
  Fo2 = sf_output("Elasticz"); /* Elastic-wave z-component */
  
  puthead3(Fo1, nz, nx, nt, dz, dx, dt, z0, x0, 0.0); /* should be m instead*/
  puthead3(Fo2, nz, nx, nt, dz, dx, dt, z0, x0, 0.0); /* should be m instead*/

  if (sep) {
    if (NULL != sf_getstring("Pwavez")) {
      Fpz = sf_output("Pwavez");
      puthead3(Fpz, nz, nx, nt, dz, dx, dt, z0, x0, 0.0);
    } else Fpz = NULL;
    if (NULL != sf_getstring("Pwavex")) {
      Fpx = sf_output("Pwavex");
      puthead3(Fpx, nz, nx, nt, dz, dx, dt, z0, x0, 0.0);
    } else Fpx = NULL;
    if (NULL != sf_getstring("Swavez")) {
      Fsz = sf_output("Swavez");
      puthead3(Fsz, nz, nx, nt, dz, dx, dt, z0, x0, 0.0);
    } else Fsz = NULL;
    if (NULL != sf_getstring("Swavex")) {
      Fsx = sf_output("Swavex");
      puthead3(Fsx, nz, nx, nt, dz, dx, dt, z0, x0, 0.0);
    } else Fsx = NULL;
  } else {
    Fpz = NULL;
    Fpx = NULL;
    Fsz = NULL;
    Fsx = NULL;
  }

  /**************** wavefield propagation ****************/
  sf_warning("==================================================");
  sf_warning("==      Propagation Using Elastic Wave Eq.      ==");
  sf_warning("==================================================");

  uz   = sf_floatalloc(nzx2);
  ux   = sf_floatalloc(nzx2);
  upz1 = sf_floatalloc(nzx2);
  upx1 = sf_floatalloc(nzx2);
  upz2 = sf_floatalloc(nzx2);
  upx2 = sf_floatalloc(nzx2);
  usz1 = sf_floatalloc(nzx2);
  usx1 = sf_floatalloc(nzx2);
  usz2 = sf_floatalloc(nzx2);
  usx2 = sf_floatalloc(nzx2);

  cwavepz= sf_complexalloc(nk);
  cwavepx= sf_complexalloc(nk);
  cwavesz= sf_complexalloc(nk);
  cwavesx= sf_complexalloc(nk);
  cwavez = sf_complexalloc(nk);
  cwavex = sf_complexalloc(nk);
  cwavem = sf_complexalloc(nk);
  wavep  = sf_floatalloc2(nzx2,m2p);
  waves  = sf_floatalloc2(nzx2,m2s);

#ifdef _OPENMP
#pragma omp parallel for default(shared) private(iz)
#endif
  for (iz=0; iz < nzx2; iz++) {
    uz[iz]  =0.0f;
    ux[iz]  =0.0f;
    upz1[iz]=0.0f;
    upx1[iz]=0.0f;
    upz2[iz]=0.0f;
    upx2[iz]=0.0f;
    usz1[iz]=0.0f;
    usx1[iz]=0.0f;
    usz2[iz]=0.0f;
    usx2[iz]=0.0f;
  }

  /* initialize abc */
  if (abc)
    abc_init(nz,nx,nz2,nx2,nbt,nbb,nbl,nbr,ct,cb,cl,cr);

  //  ifft2_allocate(cwavem);

  for(it=0;it<nt;it++) {
    t=it*dt;

    /* forward-propagating using two-step k-space solution in VTI media */

    if (src==1) {
      // 2D exploding force source
#ifdef _OPENMP
#pragma omp parallel for default(shared) private(ix,iz,j)
#endif
      for(ix=-1;ix<=1;ix++) {
	for(iz=-1;iz<=1;iz++) {
	  if(fabs(ix)+fabs(iz)==2) {
	    j = isz+iz+nz2*(isx+ix);
	    uz[j]+=iz*Ricker(t, f0, t0, A);
	    ux[j]+=ix*Ricker(t, f0, t0, A);
	  }
	}
      }
    } else {
      // 2D equil-energy force source
#ifdef _OPENMP
#pragma omp parallel for default(shared) private(ix,iz,j)
#endif
      for(ix=-1;ix<=1;ix++) {
        for(iz=-1;iz<=1;iz++) {
	  if(fabs(ix)+fabs(iz)==2) {
	    j = isz+iz+nz2*(isx+ix);
	    if(ix==-1&&iz==1)  
	      uz[j]+=sqrt(2.0)*Ricker(t, f0, t0, A);
	    if(ix==-1&&iz==-1) 
	      ux[j]+=-sqrt(2.0)*Ricker(t, f0, t0, A);
	    if(ix==1&&iz==1)  
	      ux[j]+=sqrt(2.0)*Ricker(t, f0, t0, A);
	    if(ix==1&&iz==-1) 
	      uz[j]+=-sqrt(2.0)*Ricker(t, f0, t0, A);
	  }
        }
      }
    }

    /* Going to wavenumber domain */
    fft2(uz,cwavez);
    fft2(ux,cwavex);

    /* First decompose the wavefield into P and S wave components */
#ifdef _OPENMP
#pragma omp parallel for default(shared) private(ix,iz,ik)
#endif
    for (ix=0; ix < nkx; ix++) {
      for (iz=0; iz < nkz; iz++) {
	ik = iz + ix*nkz;
	cwavepz[ik] = cwavez[ik]*kz[ik]*kz[ik] + cwavex[ik]*kx[ik]*kz[ik]; /* upz(k) */
	cwavesz[ik] = cwavez[ik] - cwavepz[ik]; /* usz(k) */
	cwavepx[ik] = cwavez[ik]*kz[ik]*kx[ik] + cwavex[ik]*kx[ik]*kx[ik]; /* upx(k) */
	cwavesx[ik] = cwavex[ik] - cwavepx[ik]; /* usx(k) */
      }
    }

    /****************************************************/
    /* Computing the Z component of P wave */

    for (im = 0; im < m2p; im++) {
#ifdef _OPENMP
#pragma omp parallel for default(shared) private(ik)
#endif
      for (ik = 0; ik < nk; ik++) {
#ifdef SF_HAS_COMPLEX_H
	cwavem[ik] = cwavepz[ik]*prt[ik][im];
#else
	cwavem[ik] = sf_crmul(cwavepz[ik],prt[ik][im]);
#endif
      }
      ifft2(wavep[im],cwavem);
    }

#ifdef _OPENMP
#pragma omp parallel for default(shared) private(ix,iz,i,j,old,c,im)
#endif
    for (ix = 0; ix < nx; ix++) {
      for (iz=0; iz < nz; iz++) {
	i = iz+ix*nz;  /* original grid */
	j = iz+ix*nz2; /* padded grid */

	old = c = upz2[j];
	c += c - upz1[j];
	upz1[j] = old;

	for (im = 0; im < m2p; im++) {
	  c += plt[im][i]*wavep[im][j];
	}

	upz2[j] = c;
      }
    }

    /* Computing the Z component of S wave */

    for (im = 0; im < m2s; im++) {
#ifdef _OPENMP
#pragma omp parallel for default(shared) private(ik)
#endif
      for (ik = 0; ik < nk; ik++) {
#ifdef SF_HAS_COMPLEX_H
	cwavem[ik] = cwavesz[ik]*srt[ik][im];
#else
	cwavem[ik] = sf_crmul(cwavesz[ik],srt[ik][im]);
#endif
      }
      ifft2(waves[im],cwavem);
    }

#ifdef _OPENMP
#pragma omp parallel for default(shared) private(ix,iz,i,j,old,c,im)
#endif
    for (ix = 0; ix < nx; ix++) {
      for (iz=0; iz < nz; iz++) {
	i = iz+ix*nz;  /* original grid */
	j = iz+ix*nz2; /* padded grid */

	old = c = usz2[j];
	c += c - usz1[j];
	usz1[j] = old;

	for (im = 0; im < m2s; im++) {
	  c += slt[im][i]*waves[im][j];
	}

	usz2[j] = c;
      }
    }

    /* apply abc */
    if (abc) {
      abc_apply(upz1);
      abc_apply(upz2);
      abc_apply(usz1);
      abc_apply(usz2);
    }

    /* accumulating the Z component*/
#ifdef _OPENMP
#pragma omp parallel for default(shared) private(ix,iz,j)
#endif
    for (ix = 0; ix < nx; ix++) {
      for (iz=0; iz < nz; iz++) {
	j = iz+ix*nz2; /* padded grid */
        uz[j] = upz2[j] + usz2[j];
      }
    }

    /****************************************************/
    /* Computing the X component of P wave */

    for (im = 0; im < m2p; im++) {
#ifdef _OPENMP
#pragma omp parallel for default(shared) private(ik)
#endif
      for (ik = 0; ik < nk; ik++) {
#ifdef SF_HAS_COMPLEX_H
	cwavem[ik] = cwavepx[ik]*prt[ik][im];
#else
	cwavem[ik] = sf_crmul(cwavepx[ik],prt[ik][im]);
#endif
      }
      ifft2(wavep[im],cwavem);
    }

#ifdef _OPENMP
#pragma omp parallel for default(shared) private(ix,iz,i,j,old,c,im)
#endif
    for (ix = 0; ix < nx; ix++) {
      for (iz=0; iz < nz; iz++) {
	i = iz+ix*nz;  /* original grid */
	j = iz+ix*nz2; /* padded grid */

	old = c = upx2[j];
	c += c - upx1[j];
	upx1[j] = old;

	for (im = 0; im < m2p; im++) {
	  c += plt[im][i]*wavep[im][j];
	}

	upx2[j] = c;
      }
    }

    /* Computing the X component of S wave */

    for (im = 0; im < m2s; im++) {
#ifdef _OPENMP
#pragma omp parallel for default(shared) private(ik)
#endif
      for (ik = 0; ik < nk; ik++) {
#ifdef SF_HAS_COMPLEX_H
	cwavem[ik] = cwavesx[ik]*srt[ik][im];
#else
	cwavem[ik] = sf_crmul(cwavesx[ik],srt[ik][im]);
#endif
      }
      ifft2(waves[im],cwavem);
    }

#ifdef _OPENMP
#pragma omp parallel for default(shared) private(ix,iz,i,j,old,c,im)
#endif
    for (ix = 0; ix < nx; ix++) {
      for (iz=0; iz < nz; iz++) {
	i = iz+ix*nz;  /* original grid */
	j = iz+ix*nz2; /* padded grid */

	old = c = usx2[j];
	c += c - usx1[j];
	usx1[j] = old;

	for (im = 0; im < m2s; im++) {
	  c += slt[im][i]*waves[im][j];
	}

	usx2[j] = c;
      }
    }

    /* apply abc */
    if (abc) {
      abc_apply(upx1);
      abc_apply(upx2);
      abc_apply(usx1);
      abc_apply(usx2);
    }

    /* accumulating the X component*/
#ifdef _OPENMP
#pragma omp parallel for default(shared) private(ix,iz,j)
#endif
    for (ix = 0; ix < nx; ix++) {
      for (iz=0; iz < nz; iz++) {
	j = iz+ix*nz2; /* padded grid */
        ux[j] = upx2[j] + usx2[j];
      }
    }

    /******* output wavefields: component and divergence *******/
    for(ix=0;ix<nx;ix++) {
      sf_floatwrite(ux+ix*nz2,nz,Fo1);
      sf_floatwrite(uz+ix*nz2,nz,Fo2);
    }

    if (sep) {
      if (NULL != sf_getstring("Pwavez")) {
	for(ix=0;ix<nx;ix++) sf_floatwrite(upz2+ix*nz2,nz,Fpz);
      }
      if (NULL != sf_getstring("Pwavex")) {
	for(ix=0;ix<nx;ix++) sf_floatwrite(upx2+ix*nz2,nz,Fpx);
      }
      if (NULL != sf_getstring("Swavez")) {
	for(ix=0;ix<nx;ix++) sf_floatwrite(usz2+ix*nz2,nz,Fsz);
      }
      if (NULL != sf_getstring("Swavex")) {
	for(ix=0;ix<nx;ix++) sf_floatwrite(usx2+ix*nz2,nz,Fsx);
      }
    }

    if(it%100==0)
      sf_warning("Elastic: it= %d",it);
  }/* it loop */

  if (abc)
    abc_close();

  fft2_finalize();
  exit(0);
}
