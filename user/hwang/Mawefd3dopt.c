/* 3D constant density acoustic wave equation time domain modeling 
with optimized fd scheme option and hybrid one-way ABC option */
/* adj flag only injects source backwards, exact adjoint propagator is not implemented */

/*
  Copyright (C) 2014 Colorado School of Mines
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
#include "fdutil.h"

static float* compute_fdcoef(int nop, float dz2, float dx2, float dy2, bool is_optimized);
static int factorial(int n);
static float* private_fdcoef(int N);

static void expand_domain(float*** vtmp, float*** v, int nz, int nx, int ny, int nbd);
static void step_forward(float*** u0, float*** u1, float*** vel, float* fdcoef,
                  int nop, int nzpad, int nxpad, int nypad);
static float* damp_make(int ntransit);
static void apply_abc(float*** uu2, float*** uu1, float*** vel,
              float dz, float dx, float dy,
              int nzpad, int nxpad, int nypad, int nbd, int nop, float* damp);

int
main(int argc, char** argv)
{
  sf_init(argc,argv);

  sf_file file_wav=NULL;
  sf_file file_vel=NULL;
  sf_file file_wfl=NULL;
  sf_file file_dat=NULL;
  sf_file file_src=NULL;
  sf_file file_rec=NULL;
  bool verb;  /* verbose flag */
  bool adj;  /* backward source injection */
  bool hybrid;  /* hybrid ABC if hybrid=y */
  bool snap;  /* output snapshot if snap=y */
  bool optfd; /* optimized fd stencil if optfd=y */
  bool expl; /* Multiple sources, one wvlt */
  int jsnap;  /* output snapshot every jsnap timestep */
  int nbd;  /* ABC boundary size */
  int fdorder;  /* finite difference spatial accuracy order */
  int nzpad,nxpad,nypad; /* boundary padded model size */
  int ix,iy,iz,it,nx,ny,nz,nt,ns,nr;
  float dx,dy,dz,dt,x0,y0,z0,t0;
  float dx2,dy2,dz2,dt2;
  float* damp=NULL; /* damping profile for hybrid bc */
  float** ws;  /* wavelet */
  float*** vel=NULL;  /* velocity */
  float*** u0=NULL;  /* wavefield array um (up) */
  float*** u1=NULL;  /* wavefield array uo */
  float* u_dat=NULL; /* output data */
  float*** ptr_tmp=NULL;   
  pt3d* src3d=NULL;  /* source position */
  pt3d* rec3d=NULL;  /*receiver position*/
  scoef3d cssinc = NULL, crsinc = NULL;  /* sinc interpolation */
  fdm3d fdm = NULL;

  if (!sf_getbool("verb",&verb))  verb=false;
  if (!sf_getbool("adj",&adj))  adj=false;
  if (!sf_getbool("optfd",&optfd))  optfd=true;
  if (!sf_getint("fdorder",&fdorder))  fdorder=8;
  if (!sf_getint("nb",&nbd))  nbd=10;
  if (!sf_getbool("hybridbc",&hybrid))  hybrid=true;
  if (!sf_getbool("snap",&snap))  snap=false;
  if(! sf_getbool("expl",&expl)) expl=false;
  if (snap) {
    if (!sf_getint("jsnap",&jsnap))  jsnap=1;
  }
  
  /* Initialize variables */
  file_wav = sf_input("in");
  file_dat = sf_output("out");
  file_vel = sf_input("vel");
  file_src = sf_input("sou");
  file_rec = sf_input("rec");

  if (snap)  file_wfl = sf_output("wfl");

  if (!sf_histint(file_wav,"n2",&nt))  sf_error("Need n1= in input");
  if (!sf_histfloat(file_wav,"d2",&dt))  sf_error("Need d1= in input");
  if (!sf_histfloat(file_wav,"o2",&t0))  sf_error("Need o1= in input");
  if (!sf_histint(file_vel,"n1",&nz))  sf_error("Need n1= in vel");
  if (!sf_histint(file_vel,"n2",&nx))  sf_error("Need n2= in vel");
  if (!sf_histint(file_vel,"n3",&ny))  sf_error("Need n3= in vel");
  if (!sf_histfloat(file_vel,"d1",&dz))  sf_error("Need d1= in vel");
  if (!sf_histfloat(file_vel,"d2",&dx))  sf_error("Need d2= in vel");
  if (!sf_histfloat(file_vel,"d3",&dy))  sf_error("Need d3= in vel");
  if (!sf_histfloat(file_vel,"o1",&z0))  sf_error("Need o1= in vel");
  if (!sf_histfloat(file_vel,"o2",&x0))  sf_error("Need o2= in vel");
  if (!sf_histfloat(file_vel,"o3",&y0))  sf_error("Need o3= in vel");
  if (!sf_histint(file_src,"n2",&ns)) sf_error("Need n2= in sou");
  if (!sf_histint(file_rec,"n2",&nr))  sf_error("Need n2= in rec");

  /* Precompute coefficients */
  dz2 = dz*dz;  dx2 = dx*dx;  dy2 = dy*dy;  dt2 = dt*dt;
  nzpad = nz+2*nbd;  nxpad = nx+2*nbd;  nypad = ny+2*nbd;

  /* set up output header files */
  if (snap) {
    sf_putint(file_wfl,"n1",nz);
    sf_putint(file_wfl,"n2",nx);
    sf_putint(file_wfl,"n3",ny);
    sf_putint(file_wfl,"n4",(nt-1)/jsnap+1);
    sf_putfloat(file_wfl,"d1",dz);
    sf_putfloat(file_wfl,"d2",dx);
    sf_putfloat(file_wfl,"d3",dy);
    sf_putfloat(file_wfl,"d4",dt*jsnap);
    sf_putfloat(file_wfl,"o1",z0);
    sf_putfloat(file_wfl,"o2",x0);
    sf_putfloat(file_wfl,"o3",y0);
    sf_putfloat(file_wfl,"o4",t0);
    sf_putstring(file_wfl,"label1","z"); 
    sf_putstring(file_wfl,"label2","x");
    sf_putstring(file_wfl,"label3","y");
    sf_putstring(file_wfl,"label4","t");
  }
  sf_putint(file_dat,"n1",nr);
  sf_putint(file_dat,"n2",nt);
  sf_putint(file_dat,"n3",1);
  sf_putfloat(file_dat,"d1",1);
  sf_putfloat(file_dat,"d2",dt);
  sf_putfloat(file_dat,"o2",t0);
  sf_putstring(file_dat,"label2","t");
  sf_putstring(file_dat,"unit2","s");
  sf_putstring(file_dat,"label1","r");
  sf_putstring(file_dat,"unit1","");

  sf_axis az,ax,ay;
  az = sf_maxa(nz,z0,dz);
  ax = sf_maxa(nx,x0,dx);
  ay = sf_maxa(ny,y0,dy);
 
  fdm = fdutil3d_init(verb,false,az,ax,ay,nbd,1);

  /* 2-2N finite difference coefficient */
  int nop = fdorder/2; /* fd half-length stencil */
  if (nbd<nop) nbd = nop;

  float* fdcoef = compute_fdcoef(nop,dz2,dx2,dy2,optfd);

  /* Allocate memories */
  if (expl) ws = sf_floatalloc2(1,nt);
  else ws = sf_floatalloc2(ns,nt);
  vel = sf_floatalloc3(nzpad,nxpad,nypad);
  u_dat = sf_floatalloc(nr);
  src3d = pt3dalloc1(ns);
  rec3d = pt3dalloc1(nr);

  /* temperary array */
  float*** tmp_array;
  tmp_array = sf_floatalloc3(nz,nx,ny);

  /* source and receiver position */
  pt3dread1(file_src,src3d,ns,3);  /* read format: (x,y,z) */
  cssinc = sinc3d_make(ns,src3d,fdm);

  pt3dread1(file_rec,rec3d,nr,3);  /* read format: (x,y,z) */
  crsinc = sinc3d_make(nr,rec3d,fdm);

  /* read velocity and pad */
  sf_floatread(tmp_array[0][0],nz*nx*ny,file_vel);
  expand_domain(tmp_array,vel,nz,nx,ny,nbd);

  free(**tmp_array);  free(*tmp_array);  free(tmp_array);
  /* read wavelet */                 
  if (expl) sf_floatread(ws[0],nt,file_wav);
  else sf_floatread(ws[0],ns*nt,file_wav);
  if (expl) {
    for (int it=0; it<nt; it++)
      ws[it][0] *= dt2;
  } else {
    for (int it=0; it<nt; it++)
      for (int is=0; is<ns; is++)
        ws[it][is] *= dt2;
  }

  if (hybrid) {
    damp = damp_make(nbd-nop); /* compute damping profiles for hybrid bc */
    if (nbd==nop) nbd = 2*nop;
  }

  /* allocate memory for wavefield variables */
  float** oslice = sf_floatalloc2(nz,nx); /* output 3D wavefield slice-by-slice */
  u0 = sf_floatalloc3(nzpad,nxpad,nypad);
  u1 = sf_floatalloc3(nzpad,nxpad,nypad);
  /* initialize variables */
  memset(u0[0][0],0,sizeof(float)*nzpad*nxpad*nypad);
  memset(u1[0][0],0,sizeof(float)*nzpad*nxpad*nypad);
  memset(u_dat,0,sizeof(float)*nr);

  /* v = (v*dt)^2 */
  for (ix=0;ix<nzpad*nxpad*nypad;ix++)
    *(vel[0][0]+ix) *= *(vel[0][0]+ix)*dt2;

  for (it=0;it<nt;it++) {
    if (verb)  sf_warning("it=%d;",it+1);

    step_forward(u0,u1,vel,fdcoef,nop,nzpad,nxpad,nypad);

    /* backward inject source wavelet */
    if (adj) {
      if (expl) sinc3d_inject1(u0,ws[nt-1-it][0],cssinc);
      else sinc3d_inject(u0,ws[nt-1-it],cssinc);
    } else {
    /* forward inject source wavelet */
      if (expl) sinc3d_inject1(u0,ws[it][0],cssinc);
      else sinc3d_inject(u0,ws[it],cssinc);
    }

    /* apply abc */
    apply_abc(u0,u1,vel,dz,dx,dy,nzpad,nxpad,nypad,nbd,nop,damp);

    /* loop over pointers */
    ptr_tmp = u0;  u0 = u1;  u1 = ptr_tmp;

    /* extract snapshot */
    if (snap && it%jsnap==0) {
      for (iy=nbd; iy<nypad-nbd; iy++) {
        for (ix=nbd; ix<nxpad-nbd; ix++)
          for (iz=nbd; iz<nzpad-nbd; iz++)
            oslice[ix-nbd][iz-nbd] = u0[iy][ix][iz];
        sf_floatwrite(oslice[0],nx*nz,file_wfl);
      }
    }

    /* extract receiver data */
    sinc3d_extract(u0,u_dat,crsinc);

    sf_floatwrite(u_dat,nr,file_dat);

  }
  if (verb)  fprintf(stderr,"\n");

  free(**u0); free(*u0); free(u0);
  free(**u1); free(*u1); free(u1);
  free(**vel); free(*vel); free(vel);
  free(u_dat);
  free(*ws); free(ws);
  free(fdcoef);
  if (hybrid) free(damp);

  return 0;
}


static void
expand_domain(float*** vtmp, float*** v, int nz, int nx, int ny, int nbd)
{
  int iz,ix,iy;
  int nzpad=nz+2*nbd;
  int nxpad=nx+2*nbd;
  int nypad=ny+2*nbd;

  for (iy=0;iy<ny;iy++) {
    for (ix=0;ix<nx;ix++) {
      for (iz=0;iz<nz;iz++) {
        v[iy+nbd][ix+nbd][iz+nbd] = vtmp[iy][ix][iz];
      }
    }
  }
  
  /* Top and bottom */
  for (iy=0;iy<nypad;iy++) {
    for (ix=0;ix<nxpad;ix++) {
      for (iz=0;iz<nbd;iz++) {
        v[iy][ix][iz] = v[iy][ix][nbd];
        v[iy][ix][iz+nz+nbd] = v[iy][ix][nz+nbd-1];
      }
    }
  }

  /* front and back */
  for (ix=0;ix<nxpad;ix++) {
    for (iz=0;iz<nzpad;iz++) {
      for (iy=0;iy<nbd;iy++) {
        v[iy][ix][iz] = v[nbd][ix][iz];
        v[iy+nbd+ny][ix][iz] = v[ny+nbd-1][ix][iz];
      }
    }
  }

  /* left and right */
  for (iy=0;iy<nypad;iy++) {
    for (iz=0;iz<nzpad;iz++) {
      for (ix=0;ix<nbd;ix++) {
        v[iy][ix][iz] = v[iy][nbd][iz];
        v[iy][ix+nbd+nx][iz] = v[iy][nx+nbd-1][iz];
      }
    }
  }
  return;
}


static void 
step_forward(float*** u0, float*** u1, float*** vel, float* fdcoef, int nop,
            int nzpad, int nxpad, int nypad)
{
#ifdef _OPENMP
#pragma omp parallel for \
  schedule(static,1) \
  shared(nxpad,nypad,nzpad,u0,u1,vel,fdcoef)
#endif
  for (int iy=nop; iy<nypad-nop; iy++) {
    for (int ix=nop; ix<nxpad-nop; ix++) {
      for (int iz=nop; iz<nzpad-nop; iz++) {
        float lap = u1[iy][ix][iz]*fdcoef[0];
        for (int iop=1; iop<=nop; iop++) {
          lap += (u1[iy][ix][iz-iop] + u1[iy][ix][iz+iop]) * fdcoef[iop]
               + (u1[iy][ix-iop][iz] + u1[iy][ix+iop][iz]) * fdcoef[iop+nop]
               + (u1[iy-iop][ix][iz] + u1[iy+iop][ix][iz]) * fdcoef[iop+nop+nop];
        }
        u0[iy][ix][iz] = 2.*u1[iy][ix][iz] - u0[iy][ix][iz] + vel[iy][ix][iz]*lap;
      }
    }
  }
  return;
}

static float*
damp_make(int ntransit)
{

  float* damp = NULL;
  if (ntransit>0) damp = sf_floatalloc(ntransit);
  float sb = 4.0*ntransit;
  for(int ib=0; ib<ntransit; ib++) {
    float fb = ib/(sqrt(2.0)*sb);
    damp[ib] = exp(-fb*fb);
  }
  return damp;
}


static void 
apply_abc(float*** uu2, float*** uu1, float*** vel,
          float dz, float dx, float dy,
          int nzpad, int nxpad, int nypad, int nbd, int nop, float* damp)
{ /* absorbing boundary condition (A1 type) Forward-Eular scheme */
#ifdef _OPENMP
#pragma omp parallel for schedule(static,1)
#endif
  /* up and bottom boundary */
  for (int iy=0;iy<nypad;iy++) {
    for (int ix=0;ix<nxpad;ix++) {
      for (int iz=0;iz<nop;iz++) { /* one-way zones */
        uu2[iy][ix][iz] =
          uu1[iy][ix][iz] + (uu1[iy][ix][iz+1] - uu1[iy][ix][iz])*sqrtf(vel[iy][ix][iz])/dz;

        uu2[iy][ix][nzpad-iz-1] = 
          uu1[iy][ix][nzpad-iz-1]
          - (uu1[iy][ix][nzpad-iz-1]-uu1[iy][ix][nzpad-iz-2])
          * sqrtf(vel[iy][ix][nzpad-iz-1])/dz;
      }
      if (damp != NULL) {
        for (int iz=nop; iz<nbd; iz++) { /* transition zones */
          float u2_bc =
            uu1[iy][ix][iz] + (uu1[iy][ix][iz+1] - uu1[iy][ix][iz])*sqrtf(vel[iy][ix][iz])/dz;
          uu2[iy][ix][iz] = u2_bc*(1.f - damp[nbd-iz-1]) + uu2[iy][ix][iz]*damp[nbd-iz-1];

          u2_bc =
            uu1[iy][ix][nzpad-iz-1]
            - (uu1[iy][ix][nzpad-iz-1]-uu1[iy][ix][nzpad-iz-2])
            * sqrtf(vel[iy][ix][nzpad-iz-1])/dz;
          uu2[iy][ix][nzpad-iz-1] = 
            u2_bc*(1.f - damp[nbd-iz-1]) + damp[nbd-iz-1]*uu2[iy][ix][nzpad-iz-1];
        }
      }
    }
  }
#ifdef _OPENMP
#pragma omp parallel for schedule(static,1)
#endif
  /* left and right boundary */
  for (int iy=0;iy<nypad;iy++) {
    for (int iz=0;iz<nzpad;iz++) {
      for (int ix=0;ix<nop;ix++) { /* one-way zones */
        uu2[iy][ix][iz] =
          uu1[iy][ix][iz] + (uu1[iy][ix+1][iz] - uu1[iy][ix][iz])*sqrtf(vel[iy][ix][iz])/dx;

        uu2[iy][nxpad-ix-1][iz] = 
          uu1[iy][nxpad-ix-1][iz]
          - (uu1[iy][nxpad-ix-1][iz]-uu1[iy][nxpad-ix-2][iz])
          * sqrtf(vel[iy][nxpad-ix-1][iz])/dx;
      }
      if (damp != NULL) {
        for (int ix=nop; ix<nbd; ix++) { /* transition zones */
          float u2_bc =
            uu1[iy][ix][iz] 
            + (uu1[iy][ix+1][iz] - uu1[iy][ix][iz])*sqrtf(vel[iy][ix][iz])/dx;
          uu2[iy][ix][iz] = u2_bc*(1.f - damp[nbd-ix-1]) + uu2[iy][ix][iz]*damp[nbd-ix-1];

          u2_bc =
            uu1[iy][nxpad-ix-1][iz]
            - (uu1[iy][nxpad-ix-1][iz]-uu1[iy][nxpad-ix-2][iz])
            * sqrtf(vel[iy][nxpad-ix-1][iz])/dx;
          uu2[iy][nxpad-ix-1][iz] =
            u2_bc*(1.f - damp[nbd-ix-1]) + damp[nbd-ix-1]*uu2[iy][nxpad-ix-1][iz];
        }
      }
    }
  }
#ifdef _OPENMP
#pragma omp parallel for schedule(static,1)
#endif
  /* front and back boundary */
  for (int ix=0;ix<nxpad;ix++) {
    for (int iz=0;iz<nzpad;iz++) {
      for (int iy=0;iy<nop;iy++) { /* one-way zones */
        uu2[iy][ix][iz] =
          uu1[iy][ix][iz] + (uu1[iy+1][ix][iz] - uu1[iy][ix][iz])*sqrtf(vel[iy][ix][iz])/dy;

        uu2[nypad-iy-1][ix][iz] =
          uu1[nypad-iy-1][ix][iz]
          - (uu1[nypad-iy-1][ix][iz] - uu1[nypad-iy-2][ix][iz])
          * sqrtf(vel[nypad-iy-1][ix][iz])/dy;
      }
      if (damp != NULL) {
        for (int iy=nop; iy<nbd; iy++) { /* transition zones */
          float u2_bc = 
            uu1[iy][ix][iz] + (uu1[iy+1][ix][iz] - uu1[iy][ix][iz])*sqrtf(vel[iy][ix][iz])/dy;
          uu2[iy][ix][iz] = u2_bc*(1.f - damp[nbd-iy-1]) + uu2[iy][ix][iz]*damp[nbd-iy-1];

          u2_bc =
            uu1[nypad-iy-1][ix][iz]
            - (uu1[nypad-iy-1][ix][iz] - uu1[nypad-iy-2][ix][iz])
            * sqrtf(vel[nypad-iy-1][ix][iz])/dy;
          uu2[nypad-iy-1][ix][iz] = 
            u2_bc*(1.f - damp[nbd-iy-1]) + damp[nbd-iy-1]*uu2[nypad-iy-1][ix][iz];
        }
      }

    }
  }
  return;
}

static float*
compute_fdcoef(int nop, float dz2, float dx2, float dy2, bool is_optimized)
/* Optimized fd coeffifients from
1. Yang Liu. "Globally optimal finite-difference schemes based on least squares" Geophysics 78.4 (2013)
  Conventional fd coefficients formula from:
2. Chunlei Chu, Paul Stoffa "Determination of finite-difference weights using scaled binomial windows" Geophysics 77.3 (2012)
*/
{
  int ndim = 3;
  float* fdcoef = sf_floatalloc(ndim*nop+1);
  float d2[ndim];
  d2[0] = dz2; d2[1] = dx2; d2[2] = dy2;
  if (!is_optimized) {
    float *ccc = private_fdcoef(nop*2);
    for (int idim=0; idim<ndim; idim++) {
      for (int ii=1; ii<=nop; ii++)
        fdcoef[idim*nop+ii] = ccc[ii]/d2[idim];
    }
  } else { // optimized
    if (nop==1) { /* 2nd-order (unoptimized) */
      for (int idim=0; idim<ndim; idim++)
        fdcoef[idim*nop+1] = 1.f / d2[idim];
    } else if (nop==2) {  /* 4th-order */
      for (int idim=0; idim<ndim; idim++) {
        fdcoef[idim*nop+1] = 1.369074 / d2[idim];
        fdcoef[idim*nop+2] = - 0.09266816 / d2[idim];
      }
    } else if (nop==3) { /* 6th-order */
      for (int idim=0; idim<ndim; idim++) {
        fdcoef[idim*nop+1] = 1.573661f / d2[idim];
        fdcoef[idim*nop+2] = - 0.1820268f / d2[idim];
        fdcoef[idim*nop+3] = 0.01728053f / d2[idim];
      }
    } else if (nop==4) { /* 8th-order */
      for (int idim=0; idim<ndim; idim++) {
        fdcoef[idim*nop+1] = 1.700010f / d2[idim];
        fdcoef[idim*nop+2] = - 0.2554615f / d2[idim];
        fdcoef[idim*nop+3] = 0.04445392f / d2[idim];
        fdcoef[idim*nop+4] = - 0.004946851f / d2[idim];
      }
    } else if (nop==5) { /* 10th-order */
      for (int idim=0; idim<ndim; idim++) {
        fdcoef[idim*nop+1] = 1.782836f / d2[idim];
        fdcoef[idim*nop+2] = - 0.3124513f / d2[idim];
        fdcoef[idim*nop+3] = 0.07379487f / d2[idim];
        fdcoef[idim*nop+4] = - 0.01532122f / d2[idim];
        fdcoef[idim*nop+5] = 0.001954439f / d2[idim];
      }
    } else if (nop==6) { /* 12th-order */
      for (int idim=0; idim<ndim; idim++) {
        fdcoef[idim*nop+1] = 1.837023f / d2[idim];
        fdcoef[idim*nop+2] = - 0.3538895f / d2[idim];
        fdcoef[idim*nop+3] = 0.09978343f / d2[idim];
        fdcoef[idim*nop+4] = - 0.02815486f / d2[idim];
        fdcoef[idim*nop+5] = 0.006556587f / d2[idim];
        fdcoef[idim*nop+6] = - 0.0009405699f / d2[idim];
      }
    } else if (nop==7) { /* 14th-order */
      for (int idim=0; idim<ndim; idim++) {
        fdcoef[idim*nop+1] = 1.874503f / d2[idim];
        fdcoef[idim*nop+2] = - 0.3845794f / d2[idim];
        fdcoef[idim*nop+3] = 0.1215162f / d2[idim];
        fdcoef[idim*nop+4] = - 0.04121749f / d2[idim];
        fdcoef[idim*nop+5] = 0.01295522f / d2[idim];
        fdcoef[idim*nop+6] = - 0.003313813f / d2[idim];
        fdcoef[idim*nop+7] = 0.0005310053f / d2[idim];
      }
    } else if (nop==8) { /* 16th-order */
      for (int idim=0; idim<ndim; idim++) {
        fdcoef[idim*nop+1] = 1.901160f / d2[idim];
        fdcoef[idim*nop+2] = - 0.4074304f / d2[idim];
        fdcoef[idim*nop+3] = 0.1390909f / d2[idim];
        fdcoef[idim*nop+4] = - 0.05318775f / d2[idim];
        fdcoef[idim*nop+5] = 0.02004823f / d2[idim];
        fdcoef[idim*nop+6] = - 0.006828249f / d2[idim];
        fdcoef[idim*nop+7] = 0.001895771f / d2[idim];
        fdcoef[idim*nop+8] = - 0.0003369052f / d2[idim];
      }
    } else if (nop==9) {  /* 18-th order */
      for (int idim=0; idim<ndim; idim++) {
        fdcoef[idim*nop+1] = 1.919909f / d2[idim];
        fdcoef[idim*nop+2] = - 0.4240446f / d2[idim];
        fdcoef[idim*nop+3] = 0.1526043f / d2[idim];
        fdcoef[idim*nop+4] = - 0.06322328f / d2[idim];
        fdcoef[idim*nop+5] = 0.02676005f / d2[idim];
        fdcoef[idim*nop+6] = - 0.01080739f / d2[idim];
        fdcoef[idim*nop+7] = 0.003907747f / d2[idim];
        fdcoef[idim*nop+8] = - 0.001158024f / d2[idim];
        fdcoef[idim*nop+9] = 0.0002240247f / d2[idim];
      }
    } else if (nop==10) {  /* 20-th order */
      for (int idim=0; idim<ndim; idim++) {
        fdcoef[idim*nop+1] = 1.918204f / d2[idim];
        fdcoef[idim*nop+2] = - 0.4225858f / d2[idim];
        fdcoef[idim*nop+3] = 0.1514992f / d2[idim];
        fdcoef[idim*nop+4] = - 0.06249474f / d2[idim];
        fdcoef[idim*nop+5] = 0.02637196f / d2[idim];
        fdcoef[idim*nop+6] = - 0.01066631f / d2[idim];
        fdcoef[idim*nop+7] = 0.003915625f / d2[idim];
        fdcoef[idim*nop+8] = - 0.001219872f / d2[idim];
        fdcoef[idim*nop+9] = 0.0002863976f / d2[idim];
        fdcoef[idim*nop+10] = - 0.00003744830f / d2[idim];
      }
    }
  }
  fdcoef[0] = 0.f;
  for (int ii=1; ii<ndim*nop+1; ii++)
    fdcoef[0] += fdcoef[ii];
  fdcoef[0] *= - 2.0f;
  return fdcoef;
}

static int
factorial(int n)
{
  int result = 1;
  for (int i=1; i<=n; ++i)
    result *= i;
  return result;
}

static float*
private_fdcoef(int N)
{
  float *cc = calloc(N/2+1,sizeof(float));
  int halfN = N/2;
  int halfN_fact = factorial(halfN);
  halfN_fact *= halfN_fact;
  for (int n=1; n<=N/2+1; n++) {
    cc[n] = - 2.f / (n*n) * cos(n*M_PI) * halfN_fact/factorial(halfN+n)/factorial(halfN-n); 
  }
  return cc;
}


