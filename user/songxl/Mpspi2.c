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
    int nx, nt, nkx, nkz,  ix, it, ikx, ikz, iv, nv, nz, iz ;
    float dt, dx, dkx, kx, kx0, dz, dkz, kz, kz0, dv, tmpdt, pi=SF_PI, vmax, vmin, wsum;
    float **sig,  **new,  **old,  **newc;
    sf_complex  **uk, **uktmp, **curcmp, *ctracex, *ctracez; 
    kiss_fft_cfg cfgx, cfgxi, cfgz, cfgzi;
    float  **v, *vc, ***weight; 
    sf_file inp, out, vel;
    bool opt;    /* optimal padding */
   // #ifdef _OPENMP
   // int nth;
   // #endif
     

    sf_init(argc,argv);
    inp = sf_input("in");
    out = sf_output("out");
    vel = sf_input("vel");   /* velocity */
  //  grad = sf_input("grad"); /* velocity gradient */

    if (SF_FLOAT != sf_gettype(inp)) sf_error("Need float input");
    if (SF_FLOAT != sf_gettype(vel)) sf_error("Need float input");
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

    sf_putint(out,"n1",nx);
    sf_putfloat(out,"d1",dx);
//    sf_putfloat(out,"o1",x0);
    sf_putint(out,"n2",nz);
    sf_putfloat(out,"d2",dz);
    sf_putint(out,"n3",nt);
    sf_putfloat(out,"d3",dt);
    sf_putfloat(out,"o3",0.0); 

    nkx = opt? kiss_fft_next_fast_size(nx): nx;
    nkz = opt? kiss_fft_next_fast_size(nz): nz;
    if (nkx != nx) sf_warning("nkx padded to %d",nkx);
    if (nkz != nz) sf_warning("nkz padded to %d",nkz);
    dkx = 1./(nkx*dx);
    kx0 = -0.5/dx;
    dkz = 1./(nkz*dz);
    kz0 = -0.5/dz;
    cfgx = kiss_fft_alloc(nkx,0,NULL,NULL);
    cfgxi = kiss_fft_alloc(nkx,1,NULL,NULL);
    cfgz = kiss_fft_alloc(nkz,0,NULL,NULL);
    cfgzi = kiss_fft_alloc(nkz,1,NULL,NULL);


    uk = sf_complexalloc2(nkx,nkz);
    uktmp = sf_complexalloc2(nkx,nkz);
    curcmp = sf_complexalloc2(nkx,nkz);
    ctracex = sf_complexalloc(nkx);
    ctracez = sf_complexalloc(nkz);

    sig    =  sf_floatalloc2(nx,nz);
    old    =  sf_floatalloc2(nx,nz);
    new    =  sf_floatalloc2(nx,nz);
    newc   =  sf_floatalloc2(nx,nz);
    vc     =  sf_floatalloc(nv);
    weight =  sf_floatalloc3(nx,nz,nv);
    
    v = sf_floatalloc2(nx,nz);

    sf_floatread(v[0],nx*nz,vel);
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

    sf_floatread(sig[0],nx*nz,inp);		
    sf_floatwrite(sig[0],nx*nz,out);

    for (iz=0; iz < nz; iz++) {
        for (ix=0; ix < nx; ix++) {
            curcmp[iz][ix] =  sf_cmplx(sig[iz][ix],0.0);
            old[iz][ix] =  0.0; 
           }
         }
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
                  new[iz][ix] = 0.0; 
                }
           }  
/* compute u(kx,kz) */
         for (iz=0; iz < nz; iz++){
             /* Fourier transform x to kx */
                kiss_fft_stride(cfgx,(kiss_fft_cpx *)curcmp[iz],(kiss_fft_cpx *)ctracex,1); 
                for (ikx=0; ikx<nkx; ikx++) uk[iz][ikx] = ctracex[ikx]; 
             }
         for (ikx=0; ikx < nkx; ikx++){
             /* Fourier transform z to kz */
                kiss_fft_stride(cfgz,(kiss_fft_cpx *)uk[0]+ikx,(kiss_fft_cpx *)ctracez,nkx); 
                for (ikz=0; ikz<nkz; ikz++) uk[ikz][ikx] = ctracez[ikz]; 
             }

/*    #ifdef _OPENMP
    #pragma omp parallel for private(ik,ix,x,k,tmp,tmpex,tmpdt) 
    #endif
*/
        for (iv=0; iv < nv; iv++){

              for (ikz=0; ikz < nkz; ikz++) {
                   kz = (kz0 +  ikz * dkz)*2.0*pi;

#ifdef SF_HAS_COMPLEX_H
                   for (ikx=0; ikx < nkx; ikx++) {
                        kx = (kx0 +  ikx * dkx)*2.0*pi;
                        tmpdt = vc[iv]*sqrt(kx*kx+kz*kz)*dt;
                        uktmp[ikz][ikx] = uk[ikz][ikx] *2.0*cosf(tmpdt);
                  }

#else
                   for (ikx=0; ikx < nkx; ikx++) {
                        kx = kx0 +  ikx * dkx*2.0*pi;
                        tmpdt = vc[iv]*sqrt(kx*kx+kz*kz)*dt;
                        uktmp[ikz][ikx] = sf_crmul(uk[ikz][ikx],2.0*cosf(tmpdt));
                  }
#endif
               }   
/* Inverse FFT*/
         for (ikx=0; ikx < nkx; ikx++){
             /* Inverse Fourier transform kz to z */
                kiss_fft_stride(cfgzi,(kiss_fft_cpx *)uktmp[0]+ikx,(kiss_fft_cpx *)ctracez,nkx); 
                for (ikz=0; ikz < nkz; ikz++) uktmp[ikz][ikx] = ctracez[ikz]; 
             }
         for (ikz=0; ikz < nkz; ikz++){
             /* Inverse Fourier transform kx to x */
                kiss_fft_stride(cfgxi,(kiss_fft_cpx *)uktmp[ikz],(kiss_fft_cpx *)ctracex,1); 
                for (ikx=0; ikx < nkx; ikx++) uktmp[ikz][ikx] = ctracex[ikx]; 
             }
	      for (iz=0; iz < nz; iz++) {  
	           for (ix=0; ix < nx; ix++) {  
                         newc[iz][ix]  = crealf(uktmp[iz][ix]);
                         newc[iz][ix] /= nkx * nkz; 
		         newc[iz][ix] -= old[iz][ix];
		         newc[iz][ix] *= weight[iv][iz][ix];
                         new[iz][ix] += newc[iz][ix];
                   }
               }  
      }
	      for (iz=0; iz < nz; iz++) {  
                  for(ix=0; ix<nx; ix++){
	             old[iz][ix] = crealf(curcmp[iz][ix]);
	             curcmp[iz][ix] = sf_cmplx(new[iz][ix],0.0);
                   }
              }
         sf_floatwrite(new[0],nz*nx,out);
  }

   free(v);     
   free(vc);     
   free(new);     
   free(newc);     
   free(curcmp);     
   free(old);     
   free(uk);     
   free(uktmp);     
   free(sig);
   free(weight);

   sf_close();
   exit(0); 
}           
           
