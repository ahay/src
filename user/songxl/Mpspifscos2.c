/* 1-D finite-difference wave extrapolation */
/*
  Copyright (C) 2008 University of Texas at Austin
  
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
#include <math.h>
#include <limits.h>
#ifdef _OPENMP
#include <omp.h>
#endif
int main(int argc, char* argv[]) 
{
    int nx, nt, nkx, nkz,  ix, it, ikx, ikz, iv, nv, nz, iz, isx, isz;
    float dt, dx, dkx, kx,  dz, dkz, kz, dv, tmpdt, pi=SF_PI, vmax, vmin, wsum;
    float **nxt,  **old,  **cur,  **nxtc, *wav, **uk, **uktmp;
    //sf_complex  **uk, **uktmp, **curcmp, *ctracex, *ctracez; 
    //kiss_fft_cfg cfgx, cfgxi, cfgz, cfgzi;
    float  **v, *vc, ***weight; 
    sf_file out, vel, source;
    bool opt;    /* optimal padding */
   // #ifdef _OPENMP
   // int nth;
   // #endif
     

    sf_init(argc,argv);
    //inp = sf_input("in");
    out = sf_output("out");
    vel = sf_input("vel");   /* velocity */
    //source = sf_input("source");   /* source wavlet*/
    source = sf_input("in");   /* source wavlet*/

//    if (SF_FLOAT != sf_gettype(inp)) sf_error("Need float input");
    if (SF_FLOAT != sf_gettype(vel)) sf_error("Need float input");
    if (SF_FLOAT != sf_gettype(source)) sf_error("Need float input");
    if (!sf_histint(vel,"n1",&nx)) sf_error("No n1= in input");
    if (!sf_histfloat(vel,"d1",&dx)) sf_error("No d1= in input");
    if (!sf_histint(vel,"n2",&nz)) sf_error("No n2= in input");
    if (!sf_histfloat(vel,"d2",&dz)) sf_error("No d2= in input");
  //  if (!sf_histint(inp,"n2",&nt)) sf_error("No n2= in input");
  //  if (!sf_histfloat(inp,"d2",&dt)) sf_error("No d2= in input");
    if (!sf_getbool("opt",&opt)) opt=true;
    /* if y, determine optimal size for efficiency */
    if (!sf_getfloat("dt",&dt)) sf_error("Need dt input");
    if (!sf_getint("nt",&nt)) sf_error("Need nt input");
    if (!sf_getint("nv",&nv)) sf_error("Need nv input");
    if (!sf_getint("isx",&isx)) sf_error("Need isx input");
    if (!sf_getint("isz",&isz)) sf_error("Need isz input");

    sf_putint(out,"n1",nx);
    sf_putfloat(out,"d1",dx);
//    sf_putfloat(out,"o1",x0);
    sf_putint(out,"n2",nz);
    sf_putfloat(out,"d2",dz);
    sf_putint(out,"n3",nt);
    sf_putfloat(out,"d3",dt);
    sf_putfloat(out,"o3",0.0); 

    //nkx = opt? kiss_fft_next_fast_size(nx): nx;
    //nkz = opt? kiss_fft_next_fast_size(nz): nz;
    nkx = nx;
    nkz = nz;
    dkx = 1./(2.0*kiss_fft_next_fast_size(nx-1)*dx);
    dkz = 1./(2.0*kiss_fft_next_fast_size(nz-1)*dz);



    old    =  sf_floatalloc2(nx,nz);
    cur    =  sf_floatalloc2(nx,nz);
    nxt    =  sf_floatalloc2(nx,nz);
    uk     =  sf_floatalloc2(nx,nz);
    uktmp  =  sf_floatalloc2(nx,nz);
    nxtc   =  sf_floatalloc2(nx,nz);
    vc     =  sf_floatalloc(nv);
    wav    =  sf_floatalloc(nt);
    weight =  sf_floatalloc3(nx,nz,nv);
    
    v = sf_floatalloc2(nx,nz);


    sf_floatread(v[0],nx*nz,vel);
    sf_floatread(wav,nt,source);
//    sf_floatread(vx,nx,grad);
    vmax = -FLT_MAX;
    vmin = +FLT_MAX;
    for (iz=0; iz < nz; iz++) {
        for (ix=0; ix < nx; ix++) {
            if (v[iz][ix] > vmax) vmax = v[iz][ix];
            if (v[iz][ix] < vmin) vmin = v[iz][ix];
          }
         }
    dv = (vmax-vmin)/(nv-1);
    for (iv=0; iv < nv; iv++) vc[iv] = vmin + dv * iv;
 
    for (iz=0; iz < nz; iz++){
         for (ix=0; ix < nx; ix++) {
             wsum = 0.0;
             for (iv=0; iv < nv; iv++) {
                weight[iv][iz][ix] = 1.0/((v[iz][ix]-vc[iv])*(v[iz][ix]-vc[iv])/10000.0+1.0); // Weight Function
                wsum += weight[iv][iz][ix];
             }
         for (iv=0; iv < nv; iv++) weight[iv][iz][ix] /= wsum;
        }
      }


    for (iz=0; iz < nz; iz++) {
        for (ix=0; ix < nx; ix++) {
            old[iz][ix] =  0.0; 
            cur[iz][ix] =  0.0; 
            nxtc[iz][ix] =  0.0; 
           }
         }
    cur[isz][isx] = wav[0];
    sf_floatwrite(cur[0],nx*nz,out);
/*
    #ifdef _OPENMP
    #pragma omp parallel
   {nth = omp_get_num_threads();}
    sf_warning("using %d threads",nth);
    #endif
*/
    /* propagation in time */
    for (it=1; it < nt; it++) {

         for (iz=0; iz < nz; iz++){
             for (ix=0; ix < nx; ix++){ 
                  nxt[iz][ix] = 0.0; 
                  uk[iz][ix] = cur[iz][ix]; 
                }
           }  


/* compute u(kx,kz) */
         sf_cosft_init(nx);
         for (iz=0; iz < nz; iz++){
             /* Fourier transform x to kx */
             sf_cosft_frw(uk[iz],0,1);
             }
         sf_cosft_close();

         sf_cosft_init(nz);
         for (ikx=0; ikx < nkx; ikx++){
             /* Fourier transform z to kz */
             sf_cosft_frw(uk[0],ikx,nx);
             }
         sf_cosft_close();

/*    #ifdef _OPENMP
    #pragma omp parallel for private(ik,ix,x,k,tmp,tmpex,tmpdt) 
    #endif
*/
         for (iv=0; iv < nv; iv++){

              for (ikz=0; ikz < nkz; ikz++) {
                   kz = ikz* dkz*2.0*pi;

                   for (ikx=0; ikx < nkx; ikx++) {
                        kx = ikx * dkx*2.0*pi;
                        tmpdt = 2.0*(cosf(vc[iv]*sqrtf(kx*kx+kz*kz)*dt)-1.0)/(vc[iv]*vc[iv]*dt*dt);
                        uktmp[ikz][ikx] = uk[ikz][ikx]*tmpdt;
                   }

              }   
/* Inverse FFT*/
              sf_cosft_init(nz);
              for (ikx=0; ikx < nkx; ikx++){
             /* Inverse Fourier transform kz to z */
                  sf_cosft_inv(uktmp[0],ikx,nx);
              }
              sf_cosft_close();
              sf_cosft_init(nx);
              for (iz=0; iz < nz; iz++){
             /* Inverse Fourier transform kx to x */
              sf_cosft_inv(uktmp[iz],0,1);
              }
              sf_cosft_close();

	      for (iz=0; iz < nz; iz++) {  
	          for (ix=0; ix < nx; ix++) {  
                      nxtc[iz][ix]  = uktmp[iz][ix];
                   //   nxtc[iz][ix] /= nkx * nkz; 
		      nxtc[iz][ix] *= weight[iv][iz][ix];
                      nxt[iz][ix] += nxtc[iz][ix];
                  }
              }  
         }
	 for (iz=0; iz < nz; iz++) {  
             for(ix=0; ix<nx; ix++){
                nxt[iz][ix] *= v[iz][ix]*v[iz][ix]*dt*dt;
                nxt[iz][ix] += 2.0*cur[iz][ix]-old[iz][ix];
             }
         }
         //nxt[isz][isx] += wav[it];
         nxt[isz][isx] += wav[it];
         //nxt[isz][isx] += wav[it]*v[isz][isx]*v[isz][isx]*dt*dt;
	 for (iz=0; iz < nz; iz++) {  
             for(ix=0; ix<nx; ix++){
	        old[iz][ix] = cur[iz][ix]; 
	        //old[iz][ix] = creal(curcmp[iz][ix]);
	        cur[iz][ix] = nxt[iz][ix]; 
	       // curcmp[iz][ix] = sf_cmplx(nxt[iz][ix],0.0);
             }
         }
         sf_floatwrite(nxt[0],nz*nx,out);
    }

    free(v);     
    free(vc);     
    free(nxt);     
    free(nxtc);     
    free(old);     
    free(uk);     
    free(uktmp);     
    free(weight);
 //   sf_fileclose(vel);
 //   sf_fileclose(inp);
 //   sf_fileclose(out);
 
    exit(0); 
}           
           
