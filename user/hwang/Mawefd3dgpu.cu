/* 3D acoustic wave equation finite difference time domain modeling */
/* Sigle GPU version; 2nd order in time, 8th order in space */
/*
   Copyright (C) 2013 King Abduallah University of Science and Techonology

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
  
  /*
   * :WARNING: 
   *  1.Only work for sigle GPU card therefore only small models feasible
   *  2.Might need computing capability >=2.x
   * ========================================================================== 
   * I/O File                             Description
   * --------------------------------------------------------------------------
   *  in        :    [nt]                 wavelet 
   *  vel       :    [ny][nx][nz]         velocity 
   *  den       :    [ny][nx][nz]         density (if cden=n) 
   *  sou       :    [4][ns]              shot position 
   *  rec       :    [3][nr][ns]          receiver position 
   *  wfl       :    [ns][nt][ny][nx][nz] wavefield snapshot (if snap=y) 
   *  out       :    [ns][nt][nr]         shot gathers 
   * ========================================================================== 
   */

extern "C" {
#include <rsf.h>
}
#include <cuda.h>

typedef struct Grid3d grid3d; 
typedef struct LinearWeit3d weit3d; 
typedef struct Abc3d abc3d;

struct Grid3d { /* mapping info of src/rec onto regular grids*/
  int jx,jy,jz;
  float fx,fy,fz;
};

struct LinearWeit3d { /* linear interpolation weights for src/rec mapping */
  float slice0[2][2];
  float slice1[2][2];
};

__constant__ float d_fdcoef[13];
__constant__ float cz1,cz2,cz3,cz4,cx1,cx2,cx3,cx4,cy1,cy2,cy3,cy4;
texture<float,1> tex_v;
texture<float,1> tex_d;

void expand_domain(float*** vtmp, float*** v,
                  int nz, int nx, int ny, int nbd);
void compute_interpolator_weights(pt3d* v3d, grid3d* geo_v, weit3d* weit_v,
                                  int npt, int nbd, int nz, int nx, int ny,
                                  float dz, float dx, float dy,
                                  float z0, float x0, float y0);
__global__ void make_abc3d_z_kernel(float* d_bzl, float* d_bzh, float dz,
                                    int nbd, int nzpad, int nxpad, int nypad);
__global__ void make_abc3d_x_kernel(float* d_bxl, float* d_bxh, float dx,
                                    int nbd, int nzpad, int nxpad, int nypad);
__global__ void make_abc3d_y_kernel(float* d_byl, float* d_byh, float dy,
                                    int nbd, int nzpad, int nxpad, int nypad);
__global__ void apply_abc3d_z_kernel(float* d_u1, float* d_u0, float* bzl, float* bzh,
                                    int nzpad, int nxpad, int nypad, int nbd);
__global__ void apply_abc3d_x_kernel(float* d_u1, float* d_u0, float* bxl, float* bxh,
                                    int nzpad, int nxpad, int nypad, int nbd);
__global__ void apply_abc3d_y_kernel(float* d_u1, float* d_u0, float* byl, float* byh,
                                    int nzpad, int nxpad, int nypad, int nbd);
__global__ void extract_data_kernel(float* d_dat, float* d_u0,
                                    weit3d* d_weit_r, grid3d* d_geo_r,
                                    int nr, int n1, int n2, int n3);
__global__ void inject_source_kernel(float* d_u0, float wt, weit3d d_weit_s,
                                    int s1, int s2, int s3,
                                    int n1, int n2, int n3);
__global__ void step_forward_kernel(float* d_u0, float* d_u1, 
                                    int nzpad, int nxpad, int nypad, bool cden);

#define gpu_errchk(ans) { gpu_assert((ans),__FILE__,__LINE__); }
inline void gpu_assert(cudaError_t err, char* file, int line, bool abort=true) {
  if (err != cudaSuccess) {
    fprintf(stderr,"CUDA error: %s %s %d\n",cudaGetErrorString(err),file,line);
    if (abort)  exit(err);
  }
}
int
main(int argc, char** argv)
{
  sf_init(argc,argv);

  sf_file file_wav=NULL;
  sf_file file_vel=NULL;
  sf_file file_den=NULL;
  sf_file file_wfl=NULL;
  sf_file file_dat=NULL;
  sf_file file_src=NULL;
  sf_file file_rec=NULL;
  bool verb;  /* verbose flag */
  bool cden;  /* constant if cden=y */
  bool snap;  /* output snapshot if snap=y */
  int jsnap;  /* output snapshot every jsnap timestep */
  int nbd;  /* ABC boundary size */
  int nzpad,nxpad,nypad; /* boundary padded model size */
  int ix,is,it,nx,ny,nz,nt,ns,nr;
  float dx,dy,dz,dt,x0,y0,z0,t0;
  float dx2,dy2,dz2,dt2;
  float fdcoef[13];
  float h_cz1,h_cz2,h_cz3,h_cz4,h_cx1,h_cx2,h_cx3,h_cx4,h_cy1,h_cy2,h_cy3,h_cy4;
  float* wav;  /* wavelet */
  float*** den=NULL;
  float*** vel=NULL;  /* velocity */
  float*** u1=NULL;  /* wavefield array */
  float* u_dat=NULL; /* output data */
  pt3d* src3d=NULL;  /* source position */
  pt3d** rec3d=NULL;  /* receiver position */
  grid3d* geo_s=NULL; /* source position onto regular grids */
  grid3d** geo_r=NULL; /* receiver position onto regular grids */
  weit3d* weit_s=NULL;  /* source position linear interpolation weight */
  weit3d** weit_r=NULL;  /* receiver position linear interpolation weight */

  if (!sf_getbool("verb",&verb))  verb=false;
  if (!sf_getint("nbd",&nbd))  nbd=20;
  if (!sf_getbool("snap",&snap))  snap=false;
  if (!sf_getbool("cden",&cden))  cden=true;
  if (snap) {
    if (!sf_getint("jsnap",&jsnap))  jsnap=1;
  }
  
  /* Initialize variables */
  file_wav = sf_input("in");
  file_dat = sf_output("out");
  file_vel = sf_input("vel");
  file_src = sf_input("sou");
  file_rec = sf_input("rec");

  if (!cden) file_den = sf_input("den");
  if (snap)  file_wfl = sf_output("wfl");

  if (!sf_histint(file_wav,"n1",&nt))  sf_error("Need n1= in input");
  if (!sf_histfloat(file_wav,"d1",&dt))  sf_error("Need d1= in input");
  if (!sf_histfloat(file_wav,"o1",&t0))  sf_error("Need o1= in input");
  if (!sf_histint(file_vel,"n1",&nz))  sf_error("Need n1= in vel");
  if (!sf_histint(file_vel,"n2",&nx))  sf_error("Need n2= in vel");
  if (!sf_histint(file_vel,"n3",&ny))  sf_error("Need n3= in vel");
  if (!sf_histfloat(file_vel,"d1",&dz))  sf_error("Need d1= in vel");
  if (!sf_histfloat(file_vel,"d2",&dx))  sf_error("Need d2= in vel");
  if (!sf_histfloat(file_vel,"d3",&dy))  sf_error("Need d3= in vel");
  if (!sf_histfloat(file_vel,"o1",&z0))  sf_error("Need o1= in vel");
  if (!sf_histfloat(file_vel,"o2",&x0))  sf_error("Need o2= in vel");
  if (!sf_histfloat(file_vel,"o3",&y0))  sf_error("Need o3= in vel");
  if (!sf_histint(file_src,"n2",&ns)) {
    sf_error("Need n2= in sou");
  } else {
    int nstmp;
    if (!sf_histint(file_rec,"n3",&nstmp))  sf_error("Need n3= in rec");
    if (nstmp!=ns) sf_error("Shot number inconsistency");
  }
  if (!sf_histint(file_rec,"n2",&nr))  sf_error("Need n2= in rec");

  /* Precompute coefficients */
  dz2 = dz*dz;  dx2 = dx*dx;  dy2 = dy*dy;  dt2 = dt*dt;
  nzpad = nz+2*nbd;  nxpad = nx+2*nbd;  nypad = ny+2*nbd;

  /* set up output header files */
  if (snap) {
    sf_putint(file_wfl,"n1",nzpad);
    sf_putint(file_wfl,"n2",nxpad);
    sf_putint(file_wfl,"n3",nypad);
    sf_putint(file_wfl,"n4",nt/jsnap);
    sf_putint(file_wfl,"n5",ns);
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
    sf_putstring(file_wfl,"label5","shot");
  }
  sf_putint(file_dat,"n1",nr);
  sf_putint(file_dat,"n2",nt);
  sf_putint(file_dat,"n3",ns);
  sf_putfloat(file_dat,"d2",dt);
  sf_putfloat(file_dat,"o2",t0);
  sf_putstring(file_dat,"label2","t");
  sf_putstring(file_dat,"unit2","s");

  /* 2-8 FD coefficients */
  fdcoef[1] = 1.6/dz2;
  fdcoef[2] = -0.2/dz2;
  fdcoef[3] = 8./(315*dz2);
  fdcoef[4] = -1./(560*dz2);
  fdcoef[5] = 1.6/dx2;
  fdcoef[6] = -0.2/dx2;
  fdcoef[7] = 8./(315*dx2);
  fdcoef[8] = -1./(560*dx2);
  fdcoef[9] = 1.6/dy2;
  fdcoef[10] = -0.2/dy2;
  fdcoef[11] = 8./(315*dy2);
  fdcoef[12] = -1./(560*dy2);
  fdcoef[0] = 0.f;
  /* FD coefficients for gradient */
  h_cz1 = 0.8f/dz;  h_cz2 = -0.2f/dz;  h_cz3 = 4./(105*dz); h_cz4 = -1.f/(280*dz);
  h_cx1 = 0.8f/dx;  h_cx2 = -0.2f/dx;  h_cx3 = 4./(105*dx); h_cx4 = -1.f/(280*dx);
  h_cy1 = 0.8f/dy;  h_cy2 = -0.2f/dy;  h_cy3 = 4./(105*dy); h_cy4 = -1.f/(280*dy);
  
  for (ix=1;ix<13;ix++)
    fdcoef[0] += fdcoef[ix]; 
  fdcoef[0] *= - 2.0;


  /* Allocate memory */
  wav = sf_floatalloc(nt);
  if (!cden)  den = sf_floatalloc3(nzpad,nxpad,nypad);
  vel = sf_floatalloc3(nzpad,nxpad,nypad);
  u1 = sf_floatalloc3(nzpad,nxpad,nypad);
  u_dat = sf_floatalloc(nr);
  src3d = pt3dalloc1(ns);
  rec3d = pt3dalloc2(nr,ns);
  geo_s = (grid3d* ) malloc(ns*sizeof(*geo_s));
  weit_s = (weit3d* ) malloc(ns*sizeof(*weit_s));
  geo_r = (grid3d** ) malloc(ns*sizeof(*geo_r));
  weit_r = (weit3d** ) malloc(ns*sizeof(*weit_r));
  geo_r[0] = (grid3d* ) malloc(ns*nr*sizeof(*geo_r[0]));
  weit_r[0] = (weit3d* ) malloc(ns*nr*sizeof(*weit_r[0]));
  for (is=1;is<ns;is++) {
    geo_r[is] = geo_r[0] + is*nr;
    weit_r[is] = weit_r[0] + is*nr;
  }

  /* temperary array */
  float*** tmp_array;
  tmp_array = sf_floatalloc3(nz,nx,ny);

  /* Read wavelet, source and receiver position */
  sf_floatread(wav,nt,file_wav);
  pt3dread1(file_src,src3d,ns,4);  /* read format: (x,y,z,v) */
  pt3dread2(file_rec,rec3d,nr,ns,3);
  
  /* read velocity and pad */
  sf_floatread(tmp_array[0][0],nz*nx*ny,file_vel);
  expand_domain(tmp_array,vel,nz,nx,ny,nbd);

  /* read density and pad if cden=n */
  if (!cden) {
    sf_floatread(tmp_array[0][0],nz*nx*ny,file_den);
    expand_domain(tmp_array,den,nz,nx,ny,nbd);
  }
  free(**tmp_array);  free(*tmp_array);  free(tmp_array);

  /* compute linear interpolation weights */
  compute_interpolator_weights(src3d,geo_s,weit_s,ns,nbd,nz,nx,ny,dz,dx,dy,z0,x0,y0);
  for (is=0;is<ns;is++) {
    compute_interpolator_weights(rec3d[is],geo_r[is],weit_r[is],nr,nbd,nz,nx,ny,dz,dx,dy,z0,x0,y0);
  }

  for (ix=0;ix<nzpad*nxpad*nypad;ix++)
    *(vel[0][0]+ix) *= *(vel[0][0]+ix)*dt2; /* v=(v*dt)^2 */

  /* device side */
  float *d_bzl, *d_bzh, *d_bxl, *d_bxh, *d_byl, *d_byh; /* ABC coefficients */
  float *d_u0,*d_u1,*d_vel,*d_den;
  float *d_dat;
  weit3d *d_weit_r;
  grid3d *d_geo_r;
  int memsize = sizeof(float)*nzpad*nxpad*nypad;
  gpu_errchk( cudaMalloc(&d_u0,memsize) );
  gpu_errchk( cudaMalloc(&d_u1,memsize) );
  gpu_errchk( cudaMalloc(&d_vel,memsize) );
  if (!cden)  gpu_errchk( cudaMalloc(&d_den,memsize) );
  gpu_errchk( cudaMalloc(&d_dat,sizeof(float)*nr) );
  gpu_errchk( cudaMalloc(&d_weit_r,sizeof(weit3d)*nr*ns) );
  gpu_errchk( cudaMalloc(&d_geo_r,sizeof(grid3d)*nr*ns) );

  gpu_errchk( cudaMemcpy(d_weit_r,weit_r[0],sizeof(weit3d)*nr*ns,cudaMemcpyHostToDevice) ); 
  gpu_errchk( cudaMemcpy(d_geo_r,geo_r[0],sizeof(grid3d)*nr*ns,cudaMemcpyHostToDevice) ); 
  gpu_errchk( cudaMemcpyToSymbol(d_fdcoef,fdcoef,sizeof(float)*13,0,cudaMemcpyHostToDevice) );
  gpu_errchk( cudaMemcpyToSymbol(cz1,&h_cz1,sizeof(float),0,cudaMemcpyHostToDevice) );
  gpu_errchk( cudaMemcpyToSymbol(cz2,&h_cz2,sizeof(float),0,cudaMemcpyHostToDevice) );
  gpu_errchk( cudaMemcpyToSymbol(cz3,&h_cz3,sizeof(float),0,cudaMemcpyHostToDevice) );
  gpu_errchk( cudaMemcpyToSymbol(cz4,&h_cz4,sizeof(float),0,cudaMemcpyHostToDevice) );
  gpu_errchk( cudaMemcpyToSymbol(cx1,&h_cx1,sizeof(float),0,cudaMemcpyHostToDevice) );
  gpu_errchk( cudaMemcpyToSymbol(cx2,&h_cx2,sizeof(float),0,cudaMemcpyHostToDevice) );
  gpu_errchk( cudaMemcpyToSymbol(cx3,&h_cx3,sizeof(float),0,cudaMemcpyHostToDevice) );
  gpu_errchk( cudaMemcpyToSymbol(cx4,&h_cx4,sizeof(float),0,cudaMemcpyHostToDevice) );
  gpu_errchk( cudaMemcpyToSymbol(cy1,&h_cy1,sizeof(float),0,cudaMemcpyHostToDevice) );
  gpu_errchk( cudaMemcpyToSymbol(cy2,&h_cy2,sizeof(float),0,cudaMemcpyHostToDevice) );
  gpu_errchk( cudaMemcpyToSymbol(cy3,&h_cy3,sizeof(float),0,cudaMemcpyHostToDevice) );
  gpu_errchk( cudaMemcpyToSymbol(cy4,&h_cy4,sizeof(float),0,cudaMemcpyHostToDevice) );
  gpu_errchk( cudaMemcpy(d_vel,vel[0][0],memsize,cudaMemcpyHostToDevice) );
  gpu_errchk( cudaBindTexture(0,tex_v,d_vel,memsize) );
  if (!cden) {
    gpu_errchk( cudaMemcpy(d_den,den[0][0],memsize,cudaMemcpyHostToDevice) );
    gpu_errchk( cudaBindTexture(0,tex_d,d_den,memsize) );
  }
  int n1_threads = 32;
  int n2_threads = 16;
  int n3_threads = 16;
  int nz_blocks = nzpad/n1_threads + (nzpad%n1_threads==0 ? 0:1);
  int nx_blocks = nxpad/n2_threads + (nxpad%n2_threads==0 ? 0:1);
  int ny_blocks = nypad/n3_threads + (nypad%n3_threads==0 ? 0:1);

  dim3 block_dim(n1_threads,n2_threads);
  dim3 grid_dim(nz_blocks,nx_blocks);
  dim3 block_dim_src(n1_threads,n2_threads);
  dim3 grid_dim_src(nz_blocks,nx_blocks,nypad);
  int block_rec = nr<1024 ? nr:1024;
  int grid_rec = ceilf(nr/1024.f);
  int smemsize = sizeof(float)*(n1_threads+8)*(n2_threads+8);
  dim3 abc_z_block(n2_threads,n3_threads);
  dim3 abc_z_grid(nx_blocks,ny_blocks);
  dim3 abc_x_block(n1_threads,n3_threads);
  dim3 abc_x_grid(nz_blocks,ny_blocks);
  dim3 abc_y_block(n1_threads,n2_threads);
  dim3 abc_y_grid(nz_blocks,nx_blocks);

  /* allocate ABC coefficients memory and initialize */
  gpu_errchk( cudaMalloc(&d_bzl,sizeof(float)*nxpad*nypad) );
  gpu_errchk( cudaMalloc(&d_bzh,sizeof(float)*nxpad*nypad) );
  gpu_errchk( cudaMalloc(&d_bxl,sizeof(float)*nzpad*nypad) );
  gpu_errchk( cudaMalloc(&d_bxh,sizeof(float)*nzpad*nypad) );
  gpu_errchk( cudaMalloc(&d_byl,sizeof(float)*nzpad*nxpad) );
  gpu_errchk( cudaMalloc(&d_byh,sizeof(float)*nzpad*nxpad) );
  make_abc3d_z_kernel<<<abc_z_grid,abc_z_block>>>(d_bzl,d_bzh,dz,nbd,nzpad,nxpad,nypad);
  gpu_errchk( cudaGetLastError() );
  make_abc3d_x_kernel<<<abc_x_grid,abc_x_block>>>(d_bxl,d_bxh,dx,nbd,nzpad,nxpad,nypad);
  gpu_errchk( cudaGetLastError() );
  make_abc3d_y_kernel<<<abc_y_grid,abc_y_block>>>(d_byl,d_byh,dy,nbd,nzpad,nxpad,nypad);
  gpu_errchk( cudaGetLastError() );

  abc_z_grid.z = abc_x_grid.z = abc_y_grid.z = nbd;

  cudaMemset(d_dat,0,sizeof(float)*nr);
  /* ---------------------------------------------------------------------- */
  for (is=0;is<ns;is++) { /* SHOT LOOP */
    if (verb)  fprintf(stderr,"\nShot:<%d of %d> Location: x=%f y=%f z=%f \n",is+1,ns,src3d[is].x,src3d[is].y,src3d[is].z);
    gpu_errchk( cudaMemset(d_u0,0,memsize) );
    gpu_errchk( cudaMemset(d_u1,0,memsize) ); 

    for (it=0;it<nt;it++) { /* TIME LOOP */
      if (verb)  sf_warning("it=%d;",it);

      /* fd forward 1 step*/
      step_forward_kernel<<<grid_dim,block_dim,smemsize>>>(d_u0,d_u1,nzpad,nxpad,nypad,cden);
      gpu_errchk( cudaGetLastError() );

      /* inject source */
      inject_source_kernel<<<grid_dim_src,block_dim_src>>>(d_u0,wav[it],weit_s[is],geo_s[is].jz,geo_s[is].jx,geo_s[is].jy,nzpad,nxpad,nypad);
      gpu_errchk( cudaGetLastError() );

      /* circulate pointers */
      float* d_ptr_tmp;
      d_ptr_tmp = d_u0;
      d_u0 = d_u1;
      d_u1 = d_ptr_tmp;

      /* apply ABC  */
      apply_abc3d_z_kernel<<<abc_z_grid,abc_z_block>>>(d_u1,d_u0,d_bzl,d_bzh,nzpad,nxpad,nypad,nbd);
      gpu_errchk( cudaGetLastError() );
      apply_abc3d_x_kernel<<<abc_x_grid,abc_x_block>>>(d_u1,d_u0,d_bxl,d_bxh,nzpad,nxpad,nypad,nbd);
      gpu_errchk( cudaGetLastError() );
      apply_abc3d_y_kernel<<<abc_y_grid,abc_y_block>>>(d_u1,d_u0,d_byl,d_byh,nzpad,nxpad,nypad,nbd);
      gpu_errchk( cudaGetLastError() );

      extract_data_kernel<<<grid_rec,block_rec>>>(d_dat,d_u1,&d_weit_r[is*nr],&d_geo_r[is*nr],nr,nzpad,nxpad,nypad);
      gpu_errchk( cudaGetLastError() );

      /* extract shot gather */
      gpu_errchk( cudaMemcpyAsync(u_dat,d_dat,sizeof(float)*nr,cudaMemcpyDeviceToHost) );
      sf_floatwrite(u_dat,nr,file_dat);

      /* extract snapshot */
      if (snap && it%jsnap==0) {
        sf_floatwrite(u1[0][0],nzpad*nxpad*nypad,file_wfl);
        gpu_errchk( cudaMemcpyAsync(u1[0][0],d_u1,memsize,cudaMemcpyDeviceToHost) );
      }
    } /* END TIME LOOP */
  } /* END SHOT LOOP */

  return 0;
}

void
expand_domain(float*** vtmp, float*** v, int nz, int nx, int ny, int nbd)
{
  int iz,ix,iy;
  int nzpad=nz+2*nbd;
  int nxpad=nx+2*nbd;
  int nypad=ny+2*nbd;

#ifdef _OPENMP
#pragma omp for \
  schedule(dynamic,1) \
  private(ix,iy,iz)
#endif
  for (iy=0;iy<ny;iy++) {
    for (ix=0;ix<nx;ix++) {
      for (iz=0;iz<nz;iz++) {
        v[iy+nbd][ix+nbd][iz+nbd] = vtmp[iy][ix][iz];
      }
    }
  }
  
#ifdef _OPENMP
#pragma omp for \
  schedule(dynamic,1) \
  private(ix,iy,iz)
#endif
  /* Top and bottom */
  for (iy=0;iy<nypad;iy++) {
    for (ix=0;ix<nxpad;ix++) {
      for (iz=0;iz<nbd;iz++) {
        v[iy][ix][iz] = v[iy][ix][nbd];
        v[iy][ix][iz+nz+nbd] = v[iy][ix][nz+nbd-1];
      }
    }
  }

#ifdef _OPENMP
#pragma omp for \
  schedule(dynamic,1) \
  private(ix,iy,iz)
#endif
  /* front and back */
  for (ix=0;ix<nxpad;ix++) {
    for (iz=0;iz<nzpad;iz++) {
      for (iy=0;iy<nbd;iy++) {
        v[iy][ix][iz] = v[nbd][ix][iz];
        v[iy+nbd+ny][ix][iz] = v[ny+nbd-1][ix][iz];
      }
    }
  }

#ifdef _OPENMP
#pragma omp for \
  schedule(dynamic,1) \
  private(ix,iy,iz)
#endif
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

void 
compute_interpolator_weights(pt3d* v3d, grid3d* geo_v, weit3d* weit_v,
                            int npt, int nbd, int nz, int nx, int ny, 
                            float dz, float dx, float dy,
                            float z0, float x0, float y0)
{
  int ii;
  float zmax = z0 + (nz-1)*dz;
  float xmax = x0 + (nx-1)*dx;
  float ymax = y0 + (ny-1)*dy;
  for (ii=0;ii<npt;ii++) {
    if (v3d[ii].z>=z0 && v3d[ii].z<=zmax &&
        v3d[ii].x>=x0 && v3d[ii].x<=xmax &&
        v3d[ii].y>=y0 && v3d[ii].y<=ymax) 
    {
      geo_v[ii].jz = (int) ((v3d[ii].z - z0)/dz);
      geo_v[ii].jx = (int) ((v3d[ii].x - x0)/dx);
      geo_v[ii].jy = (int) ((v3d[ii].y - y0)/dy);
      geo_v[ii].fz = (v3d[ii].z - z0)/dz - geo_v[ii].jz;
      geo_v[ii].fx = (v3d[ii].x - x0)/dx - geo_v[ii].jx;
      geo_v[ii].fy = (v3d[ii].y - y0)/dy - geo_v[ii].jy;
      geo_v[ii].jx += nbd;
      geo_v[ii].jy += nbd;
      geo_v[ii].jz += nbd;
    } else {
      geo_v[ii].jz = geo_v[ii].jx = geo_v[ii].jy = nbd;
      geo_v[ii].fz = geo_v[ii].fx = geo_v[ii].fy = 0.0;
    }
    weit_v[ii].slice0[0][0] = (1-geo_v[ii].fz)*(1-geo_v[ii].fx)*(1-geo_v[ii].fy); 
    weit_v[ii].slice0[0][1] = (  geo_v[ii].fz)*(1-geo_v[ii].fx)*(1-geo_v[ii].fy); 
    weit_v[ii].slice0[1][0] = (1-geo_v[ii].fz)*(  geo_v[ii].fx)*(1-geo_v[ii].fy); 
    weit_v[ii].slice0[1][1] = (  geo_v[ii].fz)*(  geo_v[ii].fx)*(1-geo_v[ii].fy); 
    weit_v[ii].slice1[0][0] = (1-geo_v[ii].fz)*(1-geo_v[ii].fx)*(  geo_v[ii].fy); 
    weit_v[ii].slice1[0][1] = (  geo_v[ii].fz)*(1-geo_v[ii].fx)*(  geo_v[ii].fy); 
    weit_v[ii].slice1[1][0] = (1-geo_v[ii].fz)*(  geo_v[ii].fx)*(  geo_v[ii].fy); 
    weit_v[ii].slice1[1][1] = (  geo_v[ii].fz)*(  geo_v[ii].fx)*(  geo_v[ii].fy); 
  }
  return;
}


__global__ void
step_forward_kernel(float* d_u0, float* d_u1, int n1, int n2, int n3, bool cden)
{
/* Modified from P.Micikevicius, 3D Finite Difference Computation on GPUs using CUDA. 2009 */
#define radius 4 /* half-width for halos */
  extern __shared__ float s_u[];

  int i1 = blockIdx.x*blockDim.x + threadIdx.x;
  int i2 = blockIdx.y*blockDim.y + threadIdx.y;
  int in_idx = i2*n1 + i1;  /* index in global memory array */
  int out_idx;
  int stride = n1*n2;  /* sliding stride */

  int t_i1 = threadIdx.x + radius; /* z index in shared memory tile */
  int t_i2 = threadIdx.y + radius; /* x index in shared memory tile */
  int t_n1 = blockDim.x + radius + radius;

  /* fill the "in-front" and "behind" data */
  float behind4, behind3, behind2, behind1;
  float infront4, infront3, infront2, infront1;
  float current;

  behind3 = d_u1[in_idx];  in_idx += stride;
  behind2 = d_u1[in_idx];  in_idx += stride;
  behind1 = d_u1[in_idx];  in_idx += stride;  out_idx = in_idx;
  current = d_u1[in_idx];  in_idx += stride;
  infront1 = d_u1[in_idx];  in_idx += stride;
  infront2 = d_u1[in_idx];  in_idx += stride;
  infront3 = d_u1[in_idx];  in_idx += stride;
  infront4 = d_u1[in_idx];  in_idx += stride;

  for (int i=radius; i<n3-radius; i++) {
    behind4 = behind3;
    behind3 = behind2;
    behind2 = behind1;
    behind1 = current;
    current = infront1;
    infront1 = infront2;
    infront2 = infront3;
    infront3 = infront4;
    infront4 = d_u1[in_idx];

    out_idx += stride;
    in_idx += stride;
    __syncthreads();  /* ensure all threads in block are in same y-slide */
    
    /* load elements into tile from global memory array */
    int t_idx = t_i2*t_n1+t_i1;
    s_u[t_idx] = current;  /* interior tile elements load */

    /* Near boundary threads need additional job */
    if (threadIdx.x < radius) { /* first-axis halos (left & right) */
      s_u[t_i2*t_n1+threadIdx.x] = d_u1[out_idx - radius];
      s_u[t_i2*t_n1+threadIdx.x+blockDim.x+radius] = d_u1[out_idx + blockDim.x];
    }

    if (threadIdx.y < radius) { /* second-axis halos (top & bottom) */
      s_u[threadIdx.y*t_n1+t_i1] = d_u1[out_idx - radius*n1];
      s_u[(threadIdx.y+blockDim.y+radius)*t_n1+t_i1] = d_u1[out_idx + blockDim.y*n1];
    }
    __syncthreads();  /*wait for all threads finish loading*/

    /*fd evolve one time step*/
    float vv = tex1Dfetch(tex_v,out_idx);
    d_u0[out_idx] = 2.f*current - d_u0[out_idx] +
              vv*(current*d_fdcoef[0] + 
              (s_u[t_idx+1] + s_u[t_idx-1])*d_fdcoef[1] +
              (s_u[t_idx+2] + s_u[t_idx-2])*d_fdcoef[2] +
              (s_u[t_idx+3] + s_u[t_idx-3])*d_fdcoef[3] +
              (s_u[t_idx+4] + s_u[t_idx-4])*d_fdcoef[4] +
              (s_u[t_idx+  t_n1] + s_u[t_idx-  t_n1])*d_fdcoef[5] +
              (s_u[t_idx+2*t_n1] + s_u[t_idx-2*t_n1])*d_fdcoef[6] +
              (s_u[t_idx+3*t_n1] + s_u[t_idx-3*t_n1])*d_fdcoef[7] +
              (s_u[t_idx+4*t_n1] + s_u[t_idx-4*t_n1])*d_fdcoef[8] +
              (behind1 + infront1)*d_fdcoef[9 ] +
              (behind2 + infront2)*d_fdcoef[10] +
              (behind3 + infront3)*d_fdcoef[11] +
              (behind4 + infront4)*d_fdcoef[12] );
    if (!cden) {  /* variable density term grad(\rho) \cdot grad(u)*/
      float grad_d_y = ( tex1Dfetch(tex_d,out_idx+4*stride) - tex1Dfetch(tex_d,out_idx-4*stride) )*cy4 +
                       ( tex1Dfetch(tex_d,out_idx+3*stride) - tex1Dfetch(tex_d,out_idx-3*stride) )*cy3 +
                       ( tex1Dfetch(tex_d,out_idx+2*stride) - tex1Dfetch(tex_d,out_idx-2*stride) )*cy2 +
                       ( tex1Dfetch(tex_d,out_idx+  stride) - tex1Dfetch(tex_d,out_idx-1*stride) )*cy1 ;
      float grad_d_x = ( tex1Dfetch(tex_d,out_idx+4*n1) - tex1Dfetch(tex_d,out_idx-4*n1) )*cx4 +
                       ( tex1Dfetch(tex_d,out_idx+3*n1) - tex1Dfetch(tex_d,out_idx-3*n1) )*cx3 +
                       ( tex1Dfetch(tex_d,out_idx+2*n1) - tex1Dfetch(tex_d,out_idx-2*n1) )*cx2 +
                       ( tex1Dfetch(tex_d,out_idx+  n1) - tex1Dfetch(tex_d,out_idx-  n1) )*cx1 ;
      float grad_d_z = ( tex1Dfetch(tex_d,out_idx+4) - tex1Dfetch(tex_d,out_idx-4) )*cz4 +
                       ( tex1Dfetch(tex_d,out_idx+3) - tex1Dfetch(tex_d,out_idx-3) )*cz3 +
                       ( tex1Dfetch(tex_d,out_idx+2) - tex1Dfetch(tex_d,out_idx-2) )*cz2 +
                       ( tex1Dfetch(tex_d,out_idx+1) - tex1Dfetch(tex_d,out_idx-1) )*cz1 ;
      float grad_u_z = (s_u[t_idx+1] - s_u[t_idx-1])*cz1 +
                       (s_u[t_idx+2] - s_u[t_idx-2])*cz2 +
                       (s_u[t_idx+3] - s_u[t_idx-3])*cz3 +
                       (s_u[t_idx+4] - s_u[t_idx-4])*cz4 ;
      float grad_u_x = (s_u[t_idx+  t_n1] - s_u[t_idx-  t_n1])*cz1 +
                       (s_u[t_idx+2*t_n1] - s_u[t_idx-2*t_n1])*cz2 +
                       (s_u[t_idx+3*t_n1] - s_u[t_idx-3*t_n1])*cz3 +
                       (s_u[t_idx+4*t_n1] - s_u[t_idx-4*t_n1])*cz4 ;
      float grad_u_y = (behind1 - infront1)*cy1 +
                       (behind2 - infront2)*cy2 +
                       (behind3 - infront3)*cy2 +
                       (behind4 - infront4)*cy4 ;
      d_u0[out_idx] -= vv*(grad_u_z*grad_d_z + 
                        grad_u_x*grad_d_x + 
                        grad_u_y*grad_d_y)/tex1Dfetch(tex_d,out_idx);
    }
  }
  return;
}

__global__ void
inject_source_kernel(float* d_u0, float wt, weit3d d_weit_s,
                    int s1, int s2, int s3, int n1, int n2, int n3)
{
  int i1 = blockIdx.x*blockDim.x + threadIdx.x;
  int i2 = blockIdx.y*blockDim.y + threadIdx.y;
  int i3 = blockIdx.z;
  int idx = (i3*n2+i2)*n1+i1;

  wt *= tex1Dfetch(tex_v,idx);

  float w[9];
  w[0] = d_weit_s.slice0[0][0]*wt;
  w[1] = d_weit_s.slice0[0][1]*wt;
  w[2] = d_weit_s.slice0[1][0]*wt;
  w[3] = d_weit_s.slice0[1][1]*wt;
  w[4] = d_weit_s.slice1[0][0]*wt;
  w[5] = d_weit_s.slice1[0][1]*wt;
  w[6] = d_weit_s.slice1[1][0]*wt;
  w[7] = d_weit_s.slice1[1][1]*wt;
  w[8] = 0.f;

  int zflag = (i1==s1)?0:(i1==(s1+1)?1:-99);
  int xflag = (i2==s2)?0:(i2==(s2+1)?2:-99);
  int yflag = (i3==s3)?0:(i3==(s3+1)?4:-99);
  int flag = xflag+yflag+zflag;
  if (flag<-1)  flag = 8;
  d_u0[idx] -= w[flag];
}

__global__ void
extract_data_kernel(float* d_dat, float* d_u, weit3d* d_weit_r, grid3d* d_geo_r,
                    int nr, int n1, int n2, int n3)
{ /* parallel along receiver axis  */
  int ir = blockIdx.x*blockDim.x + threadIdx.x;
  grid3d d_geo_current = d_geo_r[ir];
  weit3d d_weit_current = d_weit_r[ir];
  int idx = (d_geo_current.jy*n2+d_geo_current.jx)*n1+d_geo_current.jz;
  int stride = n1*n2;
  d_dat[ir] = 
    d_u[idx]*d_weit_current.slice0[0][0]+
    d_u[idx+1]*d_weit_current.slice0[0][1]+
    d_u[idx+n1]*d_weit_current.slice0[1][0]+
    d_u[idx+n1+1]*d_weit_current.slice0[1][1]+
    d_u[idx+stride]*d_weit_current.slice1[0][0]+
    d_u[idx+stride+1]*d_weit_current.slice1[0][1]+
    d_u[idx+stride+n1]*d_weit_current.slice1[1][0]+
    d_u[idx+stride+n1+1]*d_weit_current.slice1[1][1];
}

__global__ void
make_abc3d_z_kernel(float* d_bzl, float* d_bzh, float dz,
                    int nbd, int nzpad, int nxpad, int nypad)
{
  int i1 = blockIdx.x*blockDim.x + threadIdx.x; /* x axis */
  int i2 = blockIdx.y*blockDim.y + threadIdx.y; /* y axis */
  /* int iop = blockIdx.z; [> z axis: (0,nbd)<] */

  int bz_idx = i2*nxpad+i1;
  int v_idx = (i2*nxpad+i1)*nzpad+nbd;
  float d = tex1Dfetch(tex_v,v_idx);
  d = sqrtf(d)/dz;
  d_bzl[bz_idx] = (1.f-d)/(1.f+d);

  v_idx = (i2*nxpad+i1)*nzpad+nzpad-1-nbd;
  d = tex1Dfetch(tex_v,v_idx);
  d = sqrtf(d)/dz;
  d_bzh[bz_idx] = (1.f-d)/(1.f+d);
}

__global__ void
make_abc3d_x_kernel(float* d_bxl, float* d_bxh, float dx,
                    int nbd, int nzpad, int nxpad, int nypad)
{
  int i1 = blockIdx.x*blockDim.x + threadIdx.x; /* z axis */
  int i2 = blockIdx.y*blockDim.y + threadIdx.y; /* y axis */

  int bx_idx = i2*nzpad+i1;
  int v_idx = (i2*nxpad+nbd)*nzpad+i1; 
  float d = tex1Dfetch(tex_v,v_idx);
  d = sqrtf(d)/dx;
  d_bxl[bx_idx] = (1.f-d)/(1.f+d);

  v_idx = (i2*nxpad+nxpad-1-nbd)*nzpad+i1;
  d = tex1Dfetch(tex_v,v_idx);
  d = sqrtf(d)/dx;
  d_bxh[bx_idx] = (1.f-d)/(1.f+d);
}

__global__ void
make_abc3d_y_kernel(float* d_byl, float* d_byh, float dy,
                    int nbd, int nzpad, int nxpad, int nypad)
{
  int i1 = blockIdx.x*blockDim.x + threadIdx.x; /* z axis */
  int i2 = blockIdx.y*blockDim.y + threadIdx.y; /* x axis */

  int by_idx = i2*nzpad+i1;
  int v_idx = (nbd*nxpad+i2)*nzpad+i1;
  float d = tex1Dfetch(tex_v,v_idx);
  d = sqrtf(d)/dy;
  d_byl[by_idx] = (1.f-d)/(1.f+d);

  v_idx = ((nypad-1-nbd)*nxpad+i2)*nzpad+i1;
  d = tex1Dfetch(tex_v,v_idx);
  d = sqrtf(d)/dy;
  d_byh[by_idx] =(1.f-d)/(1.f+d);
}

__global__ void 
apply_abc3d_z_kernel(float* d_u1, float* d_u0, float* bzl, float* bzh,
                    int nzpad, int nxpad, int nypad, int nbd)
{
  int i1 = blockIdx.x*blockDim.x + threadIdx.x; /* x axis */
  int i2 = blockIdx.y*blockDim.y + threadIdx.y; /* y axis */
  int iop = blockIdx.z; /* z axis: (0,nbd)*/
  int bz_idx = i2*nxpad+i1;
  int iz = nbd - iop;
  int idx = (i2*nxpad+i1)*nzpad+iz;
  d_u1[idx] = d_u0[idx+1] + (d_u0[idx] - d_u1[idx+1])*bzl[bz_idx];
  iz = nzpad-1-nbd+iop;
  idx = (i2*nxpad+i1)*nzpad+iz;
  d_u1[idx] = d_u0[idx-1] + (d_u0[idx] - d_u1[idx-1])*bzh[bz_idx];
}

__global__ void
apply_abc3d_x_kernel(float* d_u1, float* d_u0, float* bxl, float* bxh,
                    int nzpad, int nxpad, int nypad, int nbd)
{
  int i1 = blockIdx.x*blockDim.x + threadIdx.x; /* z axis */
  int i2 = blockIdx.y*blockDim.y + threadIdx.y; /* y axis */
  int iop = blockIdx.z; /* x axis: (0,nbd)*/
  int bx_idx = i2*nzpad+i1;
  int ix = nbd - iop;
  int idx = (i2*nxpad+ix)*nzpad+i1;
  d_u1[idx] = d_u0[idx+nzpad] + (d_u0[idx]- d_u1[idx+nzpad])*bxl[bx_idx];
  ix = nxpad-1-nbd+iop;
  idx = (i2*nxpad+ix)*nzpad+i1; 
  d_u1[idx] = d_u0[idx-nzpad] + (d_u0[idx]- d_u1[idx-nzpad])*bxh[bx_idx];
}

__global__ void
apply_abc3d_y_kernel(float* d_u1, float* d_u0, float* byl, float* byh,
                    int nzpad, int nxpad, int nypad, int nbd)
{
  int i1 = blockIdx.x*blockDim.x + threadIdx.x; /* z axis */
  int i2 = blockIdx.y*blockDim.y + threadIdx.y; /* x axis */
  int iop = blockIdx.z; /* y axis: (0,nbd)*/
  int stride = nzpad*nxpad;
  int by_idx = i2*nzpad+i1;
  int iy = nbd - iop;
  int idx = (iy*nxpad+i2)*nzpad+i1;
  d_u1[idx] = d_u0[idx+stride] + (d_u0[idx] - d_u1[idx+stride])*byl[by_idx];
  iy = nypad-1-nbd+iop;
  idx = (iy*nxpad+i2)*nzpad+i1;
  d_u1[idx] = d_u0[idx-stride] + (d_u0[idx] - d_u1[idx-stride])*byh[by_idx];
}
