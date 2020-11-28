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
#include "abcpass2.h"
#include "srcsm.h"
#ifdef _OPENMP
#include <omp.h>
#endif


float dehf(float k /*current frequency*/,
          float kn /*highest frequency*/,
          float a /*suppress factor*/,
          float factor /*propotion*/);


int main(int argc, char* argv[]) 
{
    int nx, nt, nkx, nkz, ix, it, ikx, ikz, nz, iz, isx, isz, tl, sht;
    float dt, dx, dkx, kx, dz, dkz, kz, tmpdt, pi=SF_PI, o1, o2, kx0, kz0, v0;
    float **new_p, **cur_p, **new_uz, **cur_uz,  **new_ux, **cur_ux,  **ukr, **v, **d, *wav;
    float k, tmp, *tmparray, dt2, dz2, dx2, dth, dxh, dzh, v02;
    float **a, **b1, **b2; 
    kiss_fft_cpx **uk, **ukx, **ukz, **ctracex, **ctracez, tmpc;
    //kiss_fft_cfg cfgx, cfgxi, cfgz, cfgzi;
    sf_file out, vel, source, den, snaps=NULL;
    float ax, az, factor;
    int nbl, nbr, nbt, nbb;
    float ct, cb, cl, cr;
    float **fsinc,*fsinx,*fsinz,*fcosx,*fcosz;
    int nth=1, ith=0;
    int i,jm, irz;
    bool snap;
    //float **buffer;

    sf_init(argc,argv);
    out = sf_output("out");
    vel = sf_input("vel");   /* velocity */
    den = sf_input("den");   /* density */
    source = sf_input("in");   /* source wavlet*/

    if (SF_FLOAT != sf_gettype(source)) sf_error("Need float input");
    if (SF_FLOAT != sf_gettype(vel)) sf_error("Need float input");
    if (SF_FLOAT != sf_gettype(den)) sf_error("Need float input");

    if (!sf_histint(source,"n1",&tl)) sf_error("No n1= in input");
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
    if (!sf_getint("irz",&irz)) irz=isz;
    if (!sf_getint("jm",&jm)) jm=20;
    if (!sf_getfloat("ax",&ax)) ax= 2.0; /*suppress HF parameter*/
    if (!sf_getfloat("az",&az)) az= 2.0; /*suppress HF parameter*/
    if (!sf_getfloat("factor",&factor)) factor= 2.0/3.0; /*suppress HF parameter*/
    if (!sf_getbool("snap",&snap)) snap=true; /*Output snapshots*/

    if (!sf_getint("nbt",&nbt)) nbt=0;
    if (!sf_getint("nbb",&nbb)) nbb=0;
    if (!sf_getint("nbl",&nbl)) nbl=0;
    if (!sf_getint("nbr",&nbr)) nbr=0;
    if (!sf_getint("sht",&sht)) sht=0;

    if (!sf_getfloat("ct",&ct)) ct = 0.01; /*decaying parameter*/
    if (!sf_getfloat("cb",&cb)) cb = 0.01; /*decaying parameter*/
    if (!sf_getfloat("cl",&cl)) cl = 0.01; /*decaying parameter*/
    if (!sf_getfloat("cr",&cr)) cr = 0.01; /*decaying parameter*/

    if (snap) { 
       snaps = sf_output("snaps");
       sf_putint(snaps,"n1",nz-nbt-nbb);
       sf_putfloat(snaps,"d1",dz);
       sf_putint(snaps,"n2",(nx-nbl-nbr));
       sf_putfloat(snaps,"d2",dx);
       sf_putint(snaps,"n3",(int)(1+(nt-1)/jm));
       sf_putfloat(snaps,"d3",(float)jm*dt);
       sf_putfloat(snaps,"o1",0.0); 
       sf_putfloat(snaps,"o2",0.0); 
       sf_putfloat(snaps,"o3",0.0); 
     }
 
    sf_putint(out,"n1",nx-nbl-nbr);
    sf_putfloat(out,"d1",dx);
    sf_putint(out,"n2",nt);
    sf_putfloat(out,"d2",dt);
    sf_putfloat(out,"o1",0.0); 
    sf_putfloat(out,"o2",0.0); 

    nkx = nx;
    nkz = nz;
    dkx = 1./(nkx*dx);
    kx0 = -0.5/dx;
    dkz = 1./(nkz*dz);
    kz0 = -0.5/dz;

    #ifdef _OPENMP
    #pragma omp parallel
   {nth = omp_get_num_threads();}
    sf_warning("using %d threads",nth);
    #endif

    kiss_fft_cfg cfgx[nth], cfgxi[nth], cfgz[nth], cfgzi[nth];
    for (i=0; i < nth; i++) {
        cfgx[i] = kiss_fft_alloc(nkx,0,NULL,NULL);
        cfgxi[i] = kiss_fft_alloc(nkx,1,NULL,NULL);
        cfgz[i] = kiss_fft_alloc(nkz,0,NULL,NULL);
        cfgzi[i]= kiss_fft_alloc(nkz,1,NULL,NULL);
    }


    uk  = (kiss_fft_cpx **) sf_complexalloc2(nkz,nkx);
    ukx = (kiss_fft_cpx **) sf_complexalloc2(nkz,nkx);
    ukz = (kiss_fft_cpx **) sf_complexalloc2(nkz,nkx);
    ctracex = (kiss_fft_cpx **) sf_complexalloc2(nkx,nth);
    ctracez = (kiss_fft_cpx **) sf_complexalloc2(nkz,nth);

    wav    =  sf_floatalloc(nt);
    tmparray    =  sf_floatalloc(nx-nbl-nbr);
    sf_floatread(wav,tl,source);

    cur_p    =  sf_floatalloc2(nz,nx);
    new_p   =  sf_floatalloc2(nz,nx);
    cur_ux    =  sf_floatalloc2(nz,nx);
    new_ux   =  sf_floatalloc2(nz,nx);
    cur_uz    =  sf_floatalloc2(nz,nx);
    new_uz   =  sf_floatalloc2(nz,nx);
    ukr     =  sf_floatalloc2(nz,nx);
    //buffer     =  sf_floatalloc2(nz*nx,5);

    fsinc    =  sf_floatalloc2(nkz,nkx);
    fsinx    =  sf_floatalloc(nkx);
    fcosx    =  sf_floatalloc(nkx);
    fsinz    =  sf_floatalloc(nkz);
    fcosz    =  sf_floatalloc(nkz);

    a     =  sf_floatalloc2(nz,nx);
    b1    =  sf_floatalloc2(nz,nx);
    b2    =  sf_floatalloc2(nz,nx);

    v = sf_floatalloc2(nz,nx);
    d = sf_floatalloc2(nz,nx);
    sf_floatread(v[0],nx*nz,vel);
    sf_floatread(d[0],nx*nz,den);

    v0 =0.0;
    for (ix=0; ix < nx; ix++) {
        for (iz=0; iz < nz; iz++) {
            v0 += v[ix][iz];
            //v0 += v[ix][iz]*v[ix][iz];
         }
    }
    //v0 = sqrtf(v0/(nx*nz));
    v0 = (v0/(nx*nz));
    
    v02 = v0*v0;
    dt2 = dt*dt;
    dz2 = dz*dz;
    dx2 = dx*dx; 
    dth = dt*0.5;
    dzh = dz*0.5;
    dxh = dx*0.5;
    sf_warning("v0=%f",v0);
    sf_warning("dt=%f",dt);
    sf_warning("dx=%f",dx);
    sf_warning("dz=%f",dz);
    sf_warning("dkx=%f",dkx);
    sf_warning("dkz=%f",dkz);


    for (ikx=0; ikx < nkx; ikx++) {
        kx = (kx0+ikx*dkx)*2.0*pi;
        for (ikz=0; ikz < nkz; ikz++) {
            kz = (kz0+ikz*dkz)*2.0*pi;
            k = sqrtf(kx*kx+kz*kz);
            tmpdt = sinf(v0*k*dth)/(v0*k*dth);
            if (k<0.0000001) tmpdt=1.0;
            //tmpdt *= dehf(ikx,nkx,ax,factor)*dehf(ikz,nkz,az,factor);
          //  tmpdt *= dehf(fabs(kx),knx,ax,factor)*dehf(fabs(kz),knz,az,factor);
            fsinc[ikx][ikz] = tmpdt;
         }
    }

    for (ikx=0; ikx < nkx; ikx++) {
        kx = (kx0+ikx*dkx)*2.0*pi;
        fsinx[ikx] = kx*sinf(kx*dxh);
        fcosx[ikx] = kx*cosf(kx*dxh);
    }
    for (ikz=0; ikz < nkz; ikz++) {
        kz = (kz0+ikz*dkz)*2.0*pi;
        fsinz[ikz] = kz*sinf(kz*dzh);
        fcosz[ikz] = kz*cosf(kz*dzh);
    }

    for (ix=0; ix < nx; ix++) {
        for (iz=0; iz < nz; iz++) {
            tmp = dt2*(v[ix][iz]*v[ix][iz]-v02)/24.0;
            b1[ix][iz] = tmp/dz2;
            b2[ix][iz] = tmp/dx2;
             a[ix][iz] = 1.0-2.0*b1[ix][iz]-2.0*b2[ix][iz];
        }
     }
 

    for (ix=0; ix < nx; ix++) {
        for (iz=0; iz < nz; iz++) {
            cur_p[ix][iz] = 0.0;
            cur_ux[ix][iz] = 0.0;
            cur_uz[ix][iz] = 0.0;
        }
    }
//    srcsm_init(dx,dz);
//    source_smooth(cur,isx,isz,wav[0]);
    bd_init(nx-nbl-nbr,nz-nbt-nbb,nbt,nbb,nbl,nbr,ct,cb,cl,cr);
    /* propagation in time */
    for (it=0; it < nt; it++) {
        sf_warning("it=%d",it);
        if(it<tl) { cur_p[isx+nbl][isz+nbt] += wav[it]; }
        if((!(it%jm)) && snap) { 
           for (ix=0; ix<nx-nbl-nbr; ix++) sf_floatwrite(cur_p[nbl+ix]+nbt,nz-nbt-nbb,snaps);
               
        }
        for(i=0; i < nx-nbl-nbr; i++) tmparray[i] = cur_p[nbl+i][irz+nbt];
        sf_floatwrite(tmparray,nx-nbl-nbr,out); 
        #ifdef _OPENMP
        #pragma omp parallel for private(iz,ix) 
        #endif
         for (ix=0; ix < nx; ix++){
             for (iz=0; iz < nz; iz++){ 
                  new_p[ix][iz]  = 0.0; 
                  new_uz[ix][iz] = 0.0; 
                  new_ux[ix][iz] = 0.0; 
                  uk[ix][iz].r = cur_p[ix][iz];
                  uk[ix][iz].i = 0.0; 
                }
         }  

        #ifdef _OPENMP
        #pragma omp parallel for private(iz,ix,ikz,ith) 
        #endif
         for (ix=0; ix < nx;  ix++){
             /* Fourier transform z to kz */
                for (iz=1; iz < nz; iz+=2){
                    uk[ix][iz] = sf_cneg(uk[ix][iz]);
                    }
        #ifdef _OPENMP
                ith = omp_get_thread_num();
        #endif
                kiss_fft_stride(cfgz[ith],uk[ix],ctracez[ith],1); 
                for (ikz=0; ikz<nkz; ikz++) uk[ix][ikz] = ctracez[ith][ikz]; 
             }
        #ifdef _OPENMP
        #pragma omp parallel for private(ikx,ikz,ith) 
        #endif
         for (ikz=0; ikz < nkz; ikz++){
             /* Fourier transform x to kx */
                for (ikx=1; ikx<nkx; ikx+=2){
                    uk[ikx][ikz] = sf_cneg(uk[ikx][ikz]);
                    }
	#ifdef _OPENMP
                ith = omp_get_thread_num();
        #endif
                kiss_fft_stride(cfgx[ith],uk[0]+ikz,ctracex[ith],nkz); 
                for (ikx=0; ikx<nkx; ikx++) uk[ikx][ikz] = ctracex[ith][ikx]; 
             }

         #ifdef _OPENMP
         #pragma omp parallel for private(ikz,ikx,tmpdt,tmpc) 
         #endif
         for (ikx=0; ikx < nkx; ikx++) {
             for (ikz=0; ikz < nkz; ikz++) {

                 tmpdt = fsinc[ikx][ikz];
                 tmpc.r = -fsinz[ikz]*tmpdt; 
                 tmpc.i = fcosz[ikz]*tmpdt; 
                 ukz[ikx][ikz] = sf_cmul(uk[ikx][ikz],tmpc);
                 tmpc.r = -fsinx[ikx]*tmpdt; 
                 tmpc.i = fcosx[ikx]*tmpdt; 
                 ukx[ikx][ikz] = sf_cmul(uk[ikx][ikz],tmpc);

             }

         }   
/* Inverse FFT*/
        #ifdef _OPENMP
        #pragma omp parallel for private(ikz,ikx,ith) 
        #endif 
         for (ikz=0; ikz < nkz; ikz++){
         /* Inverse Fourier transform kx to x */
#ifdef _OPENMP
             ith = omp_get_thread_num();
#endif
             kiss_fft_stride(cfgxi[ith],(kiss_fft_cpx *)ukz[0]+ikz,ctracex[ith],nkz); 
             for (ikx=0; ikx < nkx; ikx++) ukz[ikx][ikz] = sf_crmul(ctracex[ith][ikx],ikx%2?-1.0:1.0); 
              }
        #ifdef _OPENMP
        #pragma omp parallel for private(ikz,ikx,ith) 
        #endif
          for (ikx=0; ikx < nkx; ikx++){
         /* Inverse Fourier transform kz to z */
#ifdef _OPENMP
              ith = omp_get_thread_num();
#endif
              kiss_fft_stride(cfgzi[ith],(kiss_fft_cpx *)ukz[ikx],ctracez[ith],1); 
              for (ikz=0; ikz < nkz; ikz++) ukz[ikx][ikz] = sf_crmul(ctracez[ith][ikz],ikz%2?-1.0:1.0); 
             }
        #ifdef _OPENMP
        #pragma omp parallel for private(ikz,ikx,ith) 
        #endif
         for (ikz=0; ikz < nkz; ikz++){
         /* Inverse Fourier transform kx to x */
#ifdef _OPENMP
             ith = omp_get_thread_num();
#endif
             kiss_fft_stride(cfgxi[ith],(kiss_fft_cpx *)ukx[0]+ikz,ctracex[ith],nkz); 
             for (ikx=0; ikx < nkx; ikx++) ukx[ikx][ikz] = sf_crmul(ctracex[ith][ikx],ikx%2?-1.0:1.0); 
              }
        #ifdef _OPENMP
        #pragma omp parallel for private(ikz,ikx,ith) 
        #endif
         for (ikx=0; ikx < nkx; ikx++){
         /* Inverse Fourier transform kz to z */
#ifdef _OPENMP
             ith = omp_get_thread_num();
#endif
             kiss_fft_stride(cfgzi[ith],(kiss_fft_cpx *)ukx[ikx],ctracez[ith],1);
             for (ikz=0; ikz < nkz; ikz++) ukx[ikx][ikz] = sf_crmul(ctracez[ith][ikz],ikz%2?-1.0:1.0); 
             }

        #ifdef _OPENMP
        #pragma omp parallel for private(iz,ix) 
        #endif
         for (ix=0; ix < nx; ix++){
             for (iz=0; iz < nz; iz++){ 
                  ukr[ix][iz] = sf_crealf(ukz[ix][iz]); 
                  ukr[ix][iz] /= (nkx*nkz); 
                }
         }
        #ifdef _OPENMP
        #pragma omp parallel for private(iz,ix) 
        #endif
         for (ix=1; ix < nx-1; ix++) {  
	     for (iz=1; iz < nz-1; iz++) {  
                 new_uz[ix][iz]  = -(ukr[ix][iz]*a[ix][iz]
                              + (ukr[ix][iz-1]+ukr[ix][iz+1])*b1[ix][iz]
                              + (ukr[ix-1][iz]+ukr[ix+1][iz])*b2[ix][iz])*2.0*dt/(d[ix][iz]+d[ix][iz+1])
                              + cur_uz[ix][iz];
             }
         }  
        #ifdef _OPENMP
        #pragma omp parallel for private(iz,ix) 
        #endif
         for (ix=0; ix < nx; ix++){
             for (iz=0; iz < nz; iz++){ 
                  ukr[ix][iz] = sf_crealf(ukx[ix][iz]); 
                  ukr[ix][iz] /= (nkx*nkz); 
                }
         }
        #ifdef _OPENMP
        #pragma omp parallel for private(iz,ix) 
        #endif
         for (ix=1; ix < nx-1; ix++) {  
	     for (iz=1; iz < nz-1; iz++) {  
                 new_ux[ix][iz]  = -(ukr[ix][iz]*a[ix][iz]
                              + (ukr[ix][iz-1]+ukr[ix][iz+1])*b1[ix][iz]
                              + (ukr[ix-1][iz]+ukr[ix+1][iz])*b2[ix][iz])*2.0*dt/(d[ix][iz]+d[ix+1][iz])
                              + cur_ux[ix][iz];
             }
         }  
         bd_decay(new_uz); 
         bd_decay(new_ux); 
        #ifdef _OPENMP
        #pragma omp parallel for private(iz,ix) 
        #endif
         for (ix=0; ix < nx; ix++){
             for (iz=0; iz < nz; iz++){ 
                 cur_uz[ix][iz] = new_uz[ix][iz];
                 cur_ux[ix][iz] = new_ux[ix][iz];
             }
         }  
        #ifdef _OPENMP
        #pragma omp parallel for private(iz,ix) 
        #endif
         for (ix=0; ix < nx; ix++){
             for (iz=0; iz < nz; iz++){ 
                  new_uz[ix][iz] = 0.0; 
                  uk[ix][iz].r = cur_uz[ix][iz];
                  uk[ix][iz].i = 0.0; 
                }
         }  
        #ifdef _OPENMP
        #pragma omp parallel for private(iz,ix,ikz,ith) 
        #endif
         for (ix=0; ix < nx;  ix++){
             /* Fourier transform z to kz */
                for (iz=1; iz < nz; iz+=2){
                    uk[ix][iz] = sf_cneg(uk[ix][iz]);
                    }
#ifdef _OPENMP
                ith = omp_get_thread_num();
#endif
                kiss_fft_stride(cfgz[ith],uk[ix],ctracez[ith],1); 
                for (ikz=0; ikz<nkz; ikz++) uk[ix][ikz] = ctracez[ith][ikz]; 
             }
        #ifdef _OPENMP
        #pragma omp parallel for private(ikz,ikx,ith) 
        #endif
         for (ikz=0; ikz < nkz; ikz++){
             /* Fourier transform x to kx */
                for (ikx=1; ikx<nkx; ikx+=2){
                    uk[ikx][ikz] = sf_cneg(uk[ikx][ikz]);
                    }
#ifdef _OPENMP
                ith = omp_get_thread_num();
#endif
                kiss_fft_stride(cfgx[ith],uk[0]+ikz,ctracex[ith],nkz); 
                for (ikx=0; ikx<nkx; ikx++) uk[ikx][ikz] = ctracex[ith][ikx]; 
             }
        #ifdef _OPENMP
        #pragma omp parallel for private(ikz,ikx,tmpdt,tmpc) 
        #endif
         for (ikx=0; ikx < nkx; ikx++) {
             for (ikz=0; ikz < nkz; ikz++) {
                 tmpdt = fsinc[ikx][ikz];
                 tmpc.r = fsinz[ikz]*tmpdt; 
                 tmpc.i = fcosz[ikz]*tmpdt; 
                 ukz[ikx][ikz] = sf_cmul(uk[ikx][ikz],tmpc);
             }

         }   
/* Inverse FFT*/
         /* Inverse Fourier transform kx to x */
        #ifdef _OPENMP
        #pragma omp parallel for private(ikz,ikx,ith) 
        #endif
         for (ikz=0; ikz < nkz; ikz++){
#ifdef _OPENMP
             ith = omp_get_thread_num();
#endif
             kiss_fft_stride(cfgxi[ith],(kiss_fft_cpx *)ukz[0]+ikz,ctracex[ith],nkz); 
             for (ikx=0; ikx < nkx; ikx++) ukz[ikx][ikz] = sf_crmul(ctracex[ith][ikx],ikx%2?-1.0:1.0); 
              }
        #ifdef _OPENMP
        #pragma omp parallel for private(ikz,ikx,ith) 
        #endif
             for (ikx=0; ikx < nkx; ikx++){
         /* Inverse Fourier transform kz to z */
#ifdef _OPENMP
                 ith = omp_get_thread_num();
#endif
                 kiss_fft_stride(cfgzi[ith],(kiss_fft_cpx *)ukz[ikx],ctracez[ith],1); 
                 for (ikz=0; ikz < nkz; ikz++) ukz[ikx][ikz] = sf_crmul(ctracez[ith][ikz],ikz%2?-1.0:1.0); 
             }
        #ifdef _OPENMP
        #pragma omp parallel for private(iz,ix) 
        #endif
         for (ix=0; ix < nx; ix++){
             for (iz=0; iz < nz; iz++){ 
                  ukr[ix][iz] = sf_crealf(ukz[ix][iz]); 
                  ukr[ix][iz] /= (nkx*nkz); 
                }
         }
        #ifdef _OPENMP
        #pragma omp parallel for private(iz,ix) 
        #endif
         for (ix=1; ix < nx-1; ix++) {  
	     for (iz=1; iz < nz-1; iz++) {  
                 new_uz[ix][iz]  = ukr[ix][iz]*a[ix][iz]
                              + (ukr[ix][iz-1]+ukr[ix][iz+1])*b1[ix][iz]
                              + (ukr[ix-1][iz]+ukr[ix+1][iz])*b2[ix][iz];
             }
         }  
         bd_decay(new_uz); 
        #ifdef _OPENMP
        #pragma omp parallel for private(iz,ix) 
        #endif
         for (ix=0; ix < nx; ix++){
             for (iz=0; iz < nz; iz++){ 
                  new_ux[ix][iz] = 0.0; 
                  uk[ix][iz].r = cur_ux[ix][iz];
                  uk[ix][iz].i = 0.0; 
                }
         }  
        #ifdef _OPENMP
        #pragma omp parallel for private(iz,ix,ikz,ith) 
        #endif
         for (ix=0; ix < nx;  ix++){
             /* Fourier transform z to kz */
                for (iz=1; iz < nz; iz+=2){
                    uk[ix][iz] = sf_cneg(uk[ix][iz]);
                    }
#ifdef _OPENMP
                ith = omp_get_thread_num();
#endif
                kiss_fft_stride(cfgz[ith],uk[ix],ctracez[ith],1); 
                for (ikz=0; ikz<nkz; ikz++) uk[ix][ikz] = ctracez[ith][ikz]; 
             }
             /* Fourier transform x to kx */
        #ifdef _OPENMP
        #pragma omp parallel for private(ith,ikz,ikx) 
        #endif
         for (ikz=0; ikz < nkz; ikz++){
                for (ikx=1; ikx<nkx; ikx+=2){
                    uk[ikx][ikz] = sf_cneg(uk[ikx][ikz]);
                    }
#ifdef _OPENMP
                ith = omp_get_thread_num();
#endif
                kiss_fft_stride(cfgx[ith],uk[0]+ikz,ctracex[ith],nkz); 
                for (ikx=0; ikx<nkx; ikx++) uk[ikx][ikz] = ctracex[ith][ikx]; 
             }
        #ifdef _OPENMP
        #pragma omp parallel for private(ikz,ikx,tmpdt,tmpc) 
        #endif
         for (ikx=0; ikx < nkx; ikx++) {
             for (ikz=0; ikz < nkz; ikz++) {
                 tmpdt = fsinc[ikx][ikz];
                 tmpc.r = fsinx[ikx]*tmpdt; 
                 tmpc.i = fcosx[ikx]*tmpdt; 
                 ukx[ikx][ikz] = sf_cmul(uk[ikx][ikz],tmpc);
             }

         }   
/* Inverse FFT*/
        #ifdef _OPENMP
        #pragma omp parallel for private(ith,ikz,ikx) 
        #endif
         for (ikz=0; ikz < nkz; ikz++){
         /* Inverse Fourier transform kx to x */
#ifdef _OPENMP
             ith = omp_get_thread_num();
#endif
             kiss_fft_stride(cfgxi[ith],(kiss_fft_cpx *)ukx[0]+ikz,ctracex[ith],nkz); 
             for (ikx=0; ikx < nkx; ikx++) ukx[ikx][ikz] = sf_crmul(ctracex[ith][ikx],ikx%2?-1.0:1.0); 
              }
         /* Inverse Fourier transform kz to z */
        #ifdef _OPENMP
        #pragma omp parallel for private(ith,ikz,ikx) 
        #endif
             for (ikx=0; ikx < nkx; ikx++){
#ifdef _OPENMP
                 ith = omp_get_thread_num();
#endif
                 kiss_fft_stride(cfgzi[ith],(kiss_fft_cpx *)ukx[ikx],ctracez[ith],1); 
                 for (ikz=0; ikz < nkz; ikz++) ukx[ikx][ikz] = sf_crmul(ctracez[ith][ikz],ikz%2?-1.0:1.0); 
             }
        #ifdef _OPENMP
        #pragma omp parallel for private(ix,iz) 
        #endif
         for (ix=0; ix < nx; ix++){
             for (iz=0; iz < nz; iz++){ 
                  ukr[ix][iz] = sf_crealf(ukx[ix][iz]); 
                  ukr[ix][iz] /= (nkx*nkz); 
                }
         }
        #ifdef _OPENMP
        #pragma omp parallel for private(ix,iz) 
        #endif
         for (ix=1; ix < nx-1; ix++) {  
	     for (iz=1; iz < nz-1; iz++) {  
                 new_ux[ix][iz]  = ukr[ix][iz]*a[ix][iz]
                              + (ukr[ix][iz-1]+ukr[ix][iz+1])*b1[ix][iz]
                              + (ukr[ix-1][iz]+ukr[ix+1][iz])*b2[ix][iz];
             }
         }  
         bd_decay(new_ux); 
        #ifdef _OPENMP
        #pragma omp parallel for private(ix,iz) 
        #endif
         for (ix=0; ix < nx; ix++){
             for (iz=0; iz < nz; iz++){ 
                 new_p[ix][iz] = -(new_uz[ix][iz] + new_ux[ix][iz])
                                  *d[ix][iz]*v[ix][iz]*v[ix][iz]*dt
                                 + cur_p[ix][iz];
             }
         }  
         bd_decay(new_p); 
        #ifdef _OPENMP
        #pragma omp parallel for private(ix,iz) 
        #endif
         for (ix=0; ix < nx; ix++){
             for (iz=0; iz < nz; iz++){ 
                 cur_p[ix][iz] = new_p[ix][iz];
             }
         }  
    }

//    srcsm_close();
    bd_close();
    free(*v);     
    free(*d);     
    free(*new_p);     
    free(*new_uz);     
    free(*new_ux);     
    free(*cur_p);     
    free(*cur_uz);     
    free(*cur_ux);     
    free(*uk);     
    free(*ukx);     
    free(*ukz);     
    free(v);     
    free(d);     
    free(new_p);     
    free(new_uz);     
    free(new_ux);     
    free(cur_p);     
    free(cur_uz);     
    free(cur_ux);     
    free(uk);     
    free(ukx);     
    free(ukz);     
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
    if (k < kmax) {
       depress = 1.0;
       }
    else
       //depress =cosf(((float)(k-kmax))/((float)(kn-kmax))*pi/2.0);
       //depress = exp(-a*(float)((k-kmax)*(k-kmax))/((float)((kn-kmax)*(kn-kmax))));
       depress = exp(-a*((k-kmax)*(k-kmax))/(((kn-kmax)*(kn-kmax))));
    return(depress);
}
           
