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
float dehf(float k /*current frequency*/,
          float kn /*highest frequency*/,
          float a /*suppress factor*/,
          float factor /*propotion*/);
int main(int argc, char* argv[]) 
{
    int nx, nt, nkx, nkz, ix, it, ikx, ikz, nz, iz, isx, isz;
    float dt, dx, dkx, kx, dz, dkz, kz, tmpdt, tmpk, pi=SF_PI, o1, o2, kx0, kz0;
    float **new,  **old,  **cur,  **ukr, *wav;
    float vx0, vx02,  vz0, vz02,  yi0, se0;
    float cosg0, sing0;
    float vk, vk2, tmpvk, kx1, kz1;
    float ***B;  
    int len;
    int *s1, *s2, is;
    float *fs1, *fs2;
    sf_file out, source, G, files1, files2, paras;
    sf_file velz;
    kiss_fft_cpx **uk, *ctracex, *ctracez;
    kiss_fft_cfg cfgx, cfgxi, cfgz, cfgzi;
    float ax, az, factor;
    float **fcos,err,*p0;
    bool de;
     

    sf_init(argc,argv);
    out = sf_output("out");
    source = sf_input("in");   /* source wavlet */
    velz = sf_input("velz");   /* vertical velocity */
    G = sf_input("G");   /* coefficients */
    paras = sf_input("paras");   /* reference values*/
    files1 = sf_input("s1");   /* Z-stencil */
    files2 = sf_input("s2");   /* X-stencil */
    

    if (SF_FLOAT != sf_gettype(velz)) sf_error("Need float input");
    if (SF_FLOAT != sf_gettype(source)) sf_error("Need float input");
    if (SF_FLOAT != sf_gettype(G)) sf_error("Need float input");
    if (SF_FLOAT != sf_gettype(paras)) sf_error("Need float input");
    if (!sf_histint(velz,"n1",&nz)) sf_error("No n1= in input");
    if (!sf_histfloat(velz,"d1",&dz)) sf_error("No d1= in input");
    if (!sf_histint(velz,"n2",&nx)) sf_error("No n2= in input");
    if (!sf_histfloat(velz,"d2",&dx)) sf_error("No d2= in input");
    if (!sf_histfloat(velz,"o1",&o1)) o1=0.0;
    if (!sf_histfloat(velz,"o2",&o2)) o2=0.0;
    if (!sf_getfloat("dt",&dt)) sf_error("Need dt input");
    if (!sf_getint("nt",&nt)) sf_error("Need nt input");
    if (!sf_getint("isx",&isx)) sf_error("Need isx input");
    if (!sf_getint("isz",&isz)) sf_error("Need isz input");
    if (!sf_getfloat("ax",&ax)) ax= 2.0; /*suppress HF parameter*/
    if (!sf_getfloat("az",&az)) az= 2.0; /*suppress HF parameter*/
    if (!sf_getfloat("factor",&factor)) factor= 1.0; /*suppress HF parameter*/
    if (!sf_getfloat("err",&err)) err= 0.0001; /*suppress HF parameter*/
    if (!sf_getbool("de",&de)) de=true;
    if (!sf_histint(files1,"n1",&len)) sf_error("No n1= in input");
    sf_warning("len=%d",len);


    sf_putint(out,"n1",nz);
    sf_putfloat(out,"d1",dz);
    sf_putint(out,"n2",nx);
    sf_putfloat(out,"d2",dx);
    sf_putint(out,"n3",nt/20);
    sf_putfloat(out,"d3",dt*20);
    sf_putfloat(out,"o1",o1); 
    sf_putfloat(out,"o2",o2); 
    sf_putfloat(out,"o3",0.0); 

    nkx = nx;
    nkz = nz;
    dkx = 1./(nkx*dx);
    kx0 = -0.5/dx;
    dkz = 1./(nkz*dz);
    kz0 = -0.5/dz;

    /* kmax = (0.5/dx > 0.5/dz)?0.5/dx*2*pi : 0.5/dz*2*pi; */
    cfgx = kiss_fft_alloc(nkx,0,NULL,NULL);
    cfgxi = kiss_fft_alloc(nkx,1,NULL,NULL);
    cfgz = kiss_fft_alloc(nkz,0,NULL,NULL);
    cfgzi = kiss_fft_alloc(nkz,1,NULL,NULL);
    uk = (kiss_fft_cpx **) sf_complexalloc2(nkz,nkx);
    ctracex = (kiss_fft_cpx *) sf_complexalloc(nkx);
    ctracez = (kiss_fft_cpx *) sf_complexalloc(nkz);



    wav    =  sf_floatalloc(nt);
    sf_floatread(wav,nt,source);
    p0    =  sf_floatalloc(4);
    sf_floatread(p0,4,paras);
    vz0 = p0[0];
    vx0 = p0[1];
    yi0 = p0[2];
    se0 = p0[3];

    old    =  sf_floatalloc2(nz,nx);
    cur    =  sf_floatalloc2(nz,nx);
    new    =  sf_floatalloc2(nz,nx);
    ukr     =  sf_floatalloc2(nz,nx);
    fcos     =  sf_floatalloc2(nkz,nkx);

    B   =  sf_floatalloc3(nz,nx,len+1);
    sf_floatread(B[0][0],nz*nx*(len+1),G);
    
    fs1    =  sf_floatalloc(len);
    fs2    =  sf_floatalloc(len);
    s1    =  sf_intalloc(len);
    s2    =  sf_intalloc(len);
    sf_floatread(fs1,len,files1);
    sf_floatread(fs2,len,files2);
    for (ix=0; ix < len; ix++) {s1[ix] = (int) fs1[ix];}
    for (ix=0; ix < len; ix++) {s2[ix] = (int) fs2[ix];}
    free(fs1); free(fs2);
    

    sf_warning("vz0=%g,vx0=%g,yi0=%g,se0=%g",vz0,vx0,yi0,se0);
    cosg0 = cosf(se0);
    /* cosg02 = cosg0*cosg0; */
    sing0 = sinf(se0);
    /* sing02 = sing0*sing0; */ 
    vx02=vx0*vx0; 
    vz02=vz0*vz0; 
    for (ikx=0; ikx < nkx; ikx++) {
        kx1 = (kx0+ikx*dkx)*2.0*pi;
        for (ikz=0; ikz < nkz; ikz++) {
            kz1 = (kz0+ikz*dkz)*2.0*pi;
            kx = kx1*cosg0+kz1*sing0;
            kz = kz1*cosg0-kx1*sing0;
            tmpvk = (vx02*kx*kx+vz02*kz*kz);
            /* k2 = kx1*kx1+kz1*kz1; */
            vk2 = 0.5*tmpvk+0.5*sqrtf(tmpvk*tmpvk-8.0*yi0/(1.0+2.0*yi0)*vx02*vz02*kx*kx*kz*kz);
            vk = sqrtf(vk2);
            tmpk = vk *dt;
/*            tmp = vz0*dt; */
/*
            if(k2==0 || tmpk < 0.000001) 
		tmpdt = -(tmp)*(tmp);
            else
		tmpdt = 2.0*(cosf(tmpk)-1.0)/k2;
*/
            tmpdt = 2.0*(cosf(tmpk));
            fcos[ikx][ikz] = tmpdt;
         }
    }

    for (ix=0; ix < nx; ix++) {
        for (iz=0; iz < nz; iz++) {
            cur[ix][iz] = 0.0;
            old[ix][iz] =  0.0; 
        }
    }
    srcsm_init(dx,dz);
//    source_smooth(cur,isx,isz,wav[0]);
    /* propagation in time */
    for (it=0; it < nt; it++) {
        cur[isx][isz] += wav[it];
        if(!(it%20)) sf_floatwrite(cur[0],nx*nz,out);

         for (ix=0; ix < nx; ix++){
             for (iz=0; iz < nz; iz++){ 
                  new[ix][iz] = 0.0; 
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
             for (ikz=0; ikz < nkz; ikz++) {
                 tmpdt = fcos[ikx][ikz];
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
	// for (ix=4; ix < nx-4; ix++) {  
	  //   for (iz=4; iz < nz-4; iz++) {  
	 for (ix=4; ix < nx-4; ix++) {  
	     for (iz=4; iz < nz-4; iz++) {  
                 for (is=0; is < len; is++) {
                     new[ix][iz]  += 0.5*(ukr[ix+s2[is]][iz+s1[is]]+ukr[ix-s2[is]][iz-s1[is]])*B[is][ix][iz];
                 }
                     new[ix][iz]  += 0.5*(ukr[ix+1][iz+1]+ukr[ix-1][iz-1]+ukr[ix+1][iz-1]+ukr[ix-1][iz+1])*B[len][ix][iz];
                 //new[ix][iz] += 2.0*cur[ix][iz]-old[ix][iz];
             }
         }
                 
	 for (ix=0; ix < nx; ix++) {  
             for(iz=0; iz < nz; iz++) {
                new[ix][iz] += -old[ix][iz];
               // new[ix][iz] += 2.0*cur[ix][iz]-old[ix][iz];
	        old[ix][iz] = cur[ix][iz]; 
	        cur[ix][iz] = new[ix][iz]; 
             }
         }
    }
    bd_close();
    srcsm_close();
    free(*new);     
    free(*cur);     
    free(*old);     
    free(*uk);     
    free(new);     
    free(cur);     
    free(old);     
    free(uk);     
}           

float dehf(float k /*current frequency*/,
          float kn /*highest frequency*/,
          float a /*suppress factor*/,
          float factor /*propotion*/)
/*< high frequency depressing>*/
{
    float kmax;
    float depress;
    kmax = (kn*factor);
    if (k < kmax || k == kmax) {
       depress = 1.0;
       }
    else
         {
       //depress =cosf(((float)(k-kmax))/((float)(kn-kmax))*pi/2.0);
       //depress = exp(-a*(float)((k-kmax)*(k-kmax))/((float)((kn-kmax)*(kn-kmax))));
       depress = exp(-a*((k-kmax)*(k-kmax))/((kn-kmax)*(kn-kmax)));
    }
    //sf_warning("%f",depress);
    return(depress);
}
           
