/* 2-D Fourier finite-difference wave extrapolation */
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
#include <math.h>
#include <limits.h>
#include "abcpass.h"
#include "srcsm.h"
float dehf(int k /*current frequency*/,
          int kn /*highest frequency*/,
          float a /*suppress factor*/,
          float factor /*propotion*/);
int main(int argc, char* argv[]) 
{
    int nx, nt, nkx, nkz, ix, it, ikx, ikz, nz, iz, isx, isz;
    float dt, dx, dkx, kx, dz, dkz, kz, tmpdt, pi=SF_PI, o1, o2, kx0, kz0, v0;
    float **nxt,  **old,  **cur,  **ukr, **v, *wav;
    float k, kmax;
    float **a, **b1, **b2, **c1, **c2, **d1, **d2; 
    kiss_fft_cpx **uk, *ctracex, *ctracez;
    kiss_fft_cfg cfgx, cfgxi, cfgz, cfgzi;
    sf_file out, vel, source, G;
    float ax, az, factor;
     

    sf_init(argc,argv);
    out = sf_output("out");
    vel = sf_input("vel");   /* velocity */
    G = sf_input("G");   /* velocity */
    source = sf_input("in");   /* source wavlet*/

//    if (SF_FLOAT != sf_gettype(inp)) sf_error("Need float input");
    if (SF_FLOAT != sf_gettype(vel)) sf_error("Need float input");
    if (SF_FLOAT != sf_gettype(source)) sf_error("Need float input");
    if (!sf_histint(vel,"n1",&nz)) sf_error("No n1= in input");
    if (!sf_histfloat(vel,"d1",&dz)) sf_error("No d1= in input");
    if (!sf_histint(vel,"n2",&nx)) sf_error("No n2= in input");
    if (!sf_histfloat(vel,"d2",&dx)) sf_error("No d2= in input");
    if (!sf_histfloat(vel,"o1",&o1)) o1=0.0;
    if (!sf_histfloat(vel,"o2",&o2)) o2=0.0;
    if (!sf_getfloat("dt",&dt)) sf_error("Need dt input");
    if (!sf_getint("nt",&nt)) sf_error("Need nt input");
    if (!sf_getint("isx",&isx)) sf_error("Need isx input");
    if (!sf_getint("isz",&isz)) sf_error("Need isz input");
    if (!sf_getfloat("ax",&ax)) ax= 5.0; /*suppress HF parameter*/
    if (!sf_getfloat("az",&az)) az= 5.0; /*suppress HF parameter*/
    if (!sf_getfloat("factor",&factor)) factor= 5.0/6.0; /*suppress HF parameter*/


    sf_putint(out,"n1",nz);
    sf_putfloat(out,"d1",dz);
    sf_putint(out,"n2",nx);
    sf_putfloat(out,"d2",dx);
    sf_putint(out,"n3",nt);
    sf_putfloat(out,"d3",dt);
    sf_putfloat(out,"o1",o1); 
    sf_putfloat(out,"o2",o2); 
    sf_putfloat(out,"o3",0.0); 

    nkx = nx;
    nkz = nz;
    dkx = 1./(nkx*dx);
    kx0 = -0.5/dx;
    dkz = 1./(nkz*dz);
    kz0 = -0.5/dz;

    kmax = (0.5/dx > 0.5/dz)?0.5/dx*2*pi : 0.5/dz*2*pi;
    cfgx = kiss_fft_alloc(nkx,0,NULL,NULL);
    cfgxi = kiss_fft_alloc(nkx,1,NULL,NULL);
    cfgz = kiss_fft_alloc(nkz,0,NULL,NULL);
    cfgzi = kiss_fft_alloc(nkz,1,NULL,NULL);
    //dkz = 1./(2.0*kiss_fft_next_fast_size(nzb-1)*dz);
    //dkz = 1./(2.0*kiss_fft_next_fast_size(nzb-1)*dz);
    uk = (kiss_fft_cpx **) sf_complexalloc2(nkz,nkx);
    ctracex = (kiss_fft_cpx *) sf_complexalloc(nkx);
    ctracez = (kiss_fft_cpx *) sf_complexalloc(nkz);



    wav    =  sf_floatalloc(nt);
    sf_floatread(wav,nt,source);

    old    =  sf_floatalloc2(nz,nx);
    cur    =  sf_floatalloc2(nz,nx);
    nxt    =  sf_floatalloc2(nz,nx);
    ukr     =  sf_floatalloc2(nz,nx);
    a     =  sf_floatalloc2(nz,nx);
    b1    =  sf_floatalloc2(nz,nx);
    b2    =  sf_floatalloc2(nz,nx);
    c1    =  sf_floatalloc2(nz,nx);
    c2    =  sf_floatalloc2(nz,nx);
    d1    =  sf_floatalloc2(nz,nx);
    d2    =  sf_floatalloc2(nz,nx);
    
    sf_floatread(a[0],nz*nx,G);
    sf_floatread(b1[0],nz*nx,G);
    sf_floatread(b2[0],nz*nx,G);
    sf_floatread(c1[0],nz*nx,G);
    sf_floatread(c2[0],nz*nx,G);
    sf_floatread(d1[0],nz*nx,G);
    sf_floatread(d2[0],nz*nx,G);
    

    v = sf_floatalloc2(nz,nx);
    sf_floatread(v[0],nx*nz,vel);
    v0 =0.0;
    for (ix=0; ix < nx; ix++) {
        for (iz=0; iz < nz; iz++) {
            v0 += v[ix][iz]*v[ix][iz];
         }
    }
    v0 = sqrtf(v0/(nx*nz));
 

    for (ix=0; ix < nx; ix++) {
        for (iz=0; iz < nz; iz++) {
            cur[ix][iz] = 0.0;
            old[ix][iz] =  0.0; 
        }
    }
    cur[isx][isz] = wav[0];
    srcsm_init(dx,dz);
//    source_smooth(cur,isx,isz,wav[0]);
    /* propagation in time */
    for (it=0; it < nt; it++) {
        sf_floatwrite(cur[0],nx*nz,out);

         for (ix=0; ix < nx; ix++){
             for (iz=0; iz < nz; iz++){ 
                  nxt[ix][iz] = 0.0; 
                  uk[ix][iz].r = cur[ix][iz];
                  uk[ix][iz].i = 0.0; 
                }
         }  

         for (ix=0; ix < nx;  ix++){
             /* Fourier transform x to kx */
                for (iz=1; iz < nz; iz+=2){
                    uk[ix][iz] = sf_cneg(uk[ix][iz]);
                    }
                kiss_fft_stride(cfgz,uk[ix],ctracez,1); 
                for (ikz=0; ikz<nkz; ikz++) uk[ix][ikz] = ctracez[ikz]; 
             }
         for (ikz=0; ikz < nkz; ikz++){
             /* Fourier transform z to kz */
                for (ikx=1; ikx<nkx; ikx+=2){
                    uk[ikx][ikz] = sf_cneg(uk[ikx][ikz]);
                    }
                kiss_fft_stride(cfgx,uk[0]+ikz,ctracex,nkz); 
                for (ikx=0; ikx<nkx; ikx++) uk[ikx][ikz] = ctracex[ikx]; 
             }

         for (ikx=0; ikx < nkx; ikx++) {
             kx = (kx0+ikx*dkx)*2.0*pi;

             for (ikz=0; ikz < nkz; ikz++) {
                 kz = (kz0+ikz*dkz)*2.0*pi;
                 k = sqrtf(kx*kx+kz*kz);
                 //tmpdt = 2.0*(cosf(v0*k*dt)-1.0);
                 tmpdt = 2.0*(cosf(v0*k*dt)-1);
                 tmpdt *= dehf(ikx,nkx,ax,factor)*dehf(ikz,nkx,az,factor);
                 uk[ikx][ikz] = sf_crmul(uk[ikx][ikz],tmpdt);
             }

         }   
/* Inverse FFT*/
         for (ikz=0; ikz < nkz; ikz++){
         /* Inverse Fourier transform kz to z */
             kiss_fft_stride(cfgxi,(kiss_fft_cpx *)uk[0]+ikz,(kiss_fft_cpx *)ctracex,nkz); 
             for (ikx=0; ikx < nkx; ikx++) uk[ikx][ikz] = sf_crmul(ctracex[ikx],ikx%2?-1.0:1.0); 
              }
             for (ikx=0; ikx < nkx; ikx++){
             /* Inverse Fourier transform kx to x */
                 kiss_fft_stride(cfgzi,(kiss_fft_cpx *)uk[ikx],(kiss_fft_cpx *)ctracez,1); 
                 for (ikz=0; ikz < nkz; ikz++) uk[ikx][ikz] = sf_crmul(ctracez[ikz],ikz%2?-1.0:1.0); 
             }

         for (ix=0; ix < nx; ix++){
             for (iz=0; iz < nz; iz++){ 
                  ukr[ix][iz] = sf_crealf(uk[ix][iz]); 
                  ukr[ix][iz] /= (nkx*nkz); 
                }
         }
         for (ix=2; ix < nx-2; ix++) {  
	     for (iz=2; iz < nz-2; iz++) {  
                 nxt[ix][iz]  = ukr[ix][iz]*a[ix][iz]
                              + 0.5*(ukr[ix][iz-1]+ukr[ix][iz+1])*b1[ix][iz]
                              + 0.5*(ukr[ix][iz-2]+ukr[ix][iz+2])*b2[ix][iz]
                              + 0.5*(ukr[ix-1][iz]+ukr[ix+1][iz])*c1[ix][iz]
                              + 0.5*(ukr[ix-2][iz]+ukr[ix+2][iz])*c2[ix][iz]
                              + 0.5*(ukr[ix+1][iz+1]+ukr[ix-1][iz-1])*d1[ix][iz]
                              + 0.5*(ukr[ix+1][iz-1]+ukr[ix-1][iz+1])*d2[ix][iz]
                              + 2.0*cur[ix][iz]
                              - old[ix][iz];
             }
         }  
         nxt[isx][isz] += wav[it];
      //   source_smooth(nxt,isx,isz,wav[it]);

                 
	 for (ix=0; ix < nx; ix++) {  
             for(iz=0; iz < nz; iz++) {
	        old[ix][iz] = cur[ix][iz]; 
	        cur[ix][iz] = nxt[ix][iz]; 
             }
         }
    }
    bd_close();
    srcsm_close();
    free(*v);     
    free(*nxt);     
    free(*cur);     
    free(*old);     
    free(*uk);     
    free(v);     
    free(nxt);     
    free(cur);     
    free(old);     
    free(uk);     
}           

float dehf(int k /*current frequency*/,
          int kn /*highest frequency*/,
          float a /*suppress factor*/,
          float factor /*propotion*/)
/*< high frequency depressing>*/
{
    int kmax;
    float depress;
    kmax = (int) (kn*factor);
    if (k < kmax) {
       depress = 1.0;
       }
    else
       //depress =cosf(((float)(k-kmax))/((float)(kn-kmax))*pi/2.0);
       //depress = exp(-a*(float)((k-kmax)*(k-kmax))/((float)((kn-kmax)*(kn-kmax))));
       depress = exp(-a*(float)((k-kmax)*(k-kmax)*(k-kmax))/((float)((kn-kmax)*(kn-kmax)*(kn-kmax))));
    return(depress);
}
           
