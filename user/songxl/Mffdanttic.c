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
#ifdef _OPENMP
#include <omp.h>
#endif
int main(int argc, char* argv[]) 
{
    int nx, nt, nkx, nkz, ix, it, ikx, ikz, nz, iz, nbt, nbb, nbl, nbr, nxb, nzb, isx, isz;
    float dt, dx, dkx, kx, dz, dkz, kz, tmpdt, pi=SF_PI, o1, o2, kx0, kz0;
    float **new,  **old,  **cur,  **dercur, **derold, *wav;
    float **vx, vx0, vx02, vx04, **vz, vz0, vz02, vz04, yi0, se0;
    float dx2, dz2, dx4, dz4, ct, cb, cl, cr; //top, bottom, left, right 
    float cosg0, cosg02, sing0, sing02;
    float vk, vk2, tmpvk, err, dt2, kx1, kz1;
//    sf_complex  *ctracex, *ctracez, **uk; 
    kiss_fft_cpx **uk, *ctracex, *ctracez;
    kiss_fft_cfg cfgx, cfgxi, cfgz, cfgzi;
    sf_file out, velx, velz, source;
    bool opt;    /* optimal padding */
     

    sf_init(argc,argv);
    out = sf_output("out");
    velx = sf_input("velx");   /* velocity */
    velz = sf_input("velz");   /* velocity */
    source = sf_input("in");   /* source wavlet*/

//    if (SF_FLOAT != sf_gettype(inp)) sf_error("Need float input");
    if (SF_FLOAT != sf_gettype(velx)) sf_error("Need float input");
    if (SF_FLOAT != sf_gettype(velz)) sf_error("Need float input");
    if (SF_FLOAT != sf_gettype(source)) sf_error("Need float input");

    if (!sf_histint(velx,"n1",&nx)) sf_error("No n1= in input");
    if (!sf_histfloat(velx,"d1",&dx)) sf_error("No d1= in input");
    if (!sf_histint(velx,"n2",&nz)) sf_error("No n2= in input");
    if (!sf_histfloat(velx,"d2",&dz)) sf_error("No d2= in input");
    if (!sf_histfloat(velx,"o1",&o1)) o1=0.0;
    if (!sf_histfloat(velx,"o2",&o2)) o2=0.0;
  //  if (!sf_histint(inp,"n2",&nt)) sf_error("No n2= in input");
  //  if (!sf_histfloat(inp,"d2",&dt)) sf_error("No d2= in input");
    if (!sf_getbool("opt",&opt)) opt=true;
    /* if y, determine optimal size for efficiency */
    if (!sf_getfloat("dt",&dt)) sf_error("Need dt input");
    if (!sf_getint("nt",&nt)) sf_error("Need nt input");
    if (!sf_getint("isx",&isx)) sf_error("Need isx input");
    if (!sf_getint("isz",&isz)) sf_error("Need isz input");
    if (!sf_getfloat("se0",&se0)) sf_error("Need se0 input");
    if (!sf_getfloat("yi0",&yi0)) sf_error("Need yi0 input");
    if (!sf_getfloat("err",&err)) err = 0.0001;

    if (!sf_getint("nbt",&nbt)) nbt=44;
    if (!sf_getint("nbb",&nbb)) nbb=44;
    if (!sf_getint("nbl",&nbl)) nbl=44;
    if (!sf_getint("nbr",&nbr)) nbr=44;

    if (!sf_getfloat("ct",&ct)) ct = 0.01; /*decaying parameter*/
    if (!sf_getfloat("cb",&cb)) cb = 0.01; /*decaying parameter*/
    if (!sf_getfloat("cl",&cl)) cl = 0.01; /*decaying parameter*/
    if (!sf_getfloat("cr",&cr)) cr = 0.01; /*decaying parameter*/

    

    sf_putint(out,"n1",nx);
    sf_putfloat(out,"d1",dx);
//    sf_putfloat(out,"o1",x0);
    sf_putint(out,"n2",nz);
    sf_putfloat(out,"d2",dz);
    sf_putint(out,"n3",nt);
    sf_putfloat(out,"d3",dt);
    sf_putfloat(out,"o1",o1); 
    sf_putfloat(out,"o2",o2); 
    sf_putfloat(out,"o3",0.0); 

    nxb = nx + nbl + nbr;
    nzb = nz + nbt + nbb;
    nkx = opt? kiss_fft_next_fast_size(nxb): nxb;
    nkz = opt? kiss_fft_next_fast_size(nzb): nzb;
    if (nkx != nxb) sf_warning("nkx padded to %d",nkx);
    if (nkz != nzb) sf_warning("nkz padded to %d",nkz);
    dkx = 1./(nkx*dx);
    kx0 = -0.5/dx;
    dkz = 1./(nkz*dz);
    kz0 = -0.5/dz;
    cfgx = kiss_fft_alloc(nkx,0,NULL,NULL);
    cfgxi = kiss_fft_alloc(nkx,1,NULL,NULL);
    cfgz = kiss_fft_alloc(nkz,0,NULL,NULL);
    cfgzi = kiss_fft_alloc(nkz,1,NULL,NULL);


    uk = (kiss_fft_cpx **) sf_complexalloc2(nkx,nkz);
    ctracex = (kiss_fft_cpx *) sf_complexalloc(nkx);
    ctracez = (kiss_fft_cpx *) sf_complexalloc(nkz);
//    uk =  sf_complexalloc2(nkx,nkz);
//    ctracex = sf_complexalloc(nkx);
//    ctracez = sf_complexalloc(nkz);



    wav    =  sf_floatalloc(nt);
    sf_floatread(wav,nt,source);

    old    =  sf_floatalloc2(nxb,nzb);
    cur    =  sf_floatalloc2(nxb,nzb);
    new    =  sf_floatalloc2(nxb,nzb);
    derold    =  sf_floatalloc2(nxb,nzb);
    dercur    =  sf_floatalloc2(nxb,nzb);
    
    bd_init(nx,nz,nbt,nbb,nbl,nbr,ct,cb,cl,cr);

    vx = sf_floatalloc2(nxb,nzb);
    vz = sf_floatalloc2(nxb,nzb);

    /*input & extend velocity model*/
    for (iz=nbt; iz<nz+nbt; iz++){
        sf_floatread(vx[iz]+nbl,nx,velx);
        sf_floatread(vz[iz]+nbl,nx,velz);
         for (ix=0; ix<nbl; ix++){
             vx[iz][ix] = vx[iz][nbl];
             vz[iz][ix] = vz[iz][nbl];
         }
         for (ix=0; ix<nbr; ix++){
             vx[iz][nx+nbl+ix] = vx[iz][nx+nbl-1];
             vz[iz][nx+nbl+ix] = vz[iz][nx+nbl-1];
         }     
    }
    for (iz=0; iz<nbt; iz++){
        for (ix=0; ix<nxb; ix++){
            vx[iz][ix] = vx[nbt][ix];
            vz[iz][ix] = vz[nbt][ix];
        }
    }
    for (iz=0; iz<nbb; iz++){
        for (ix=0; ix<nxb; ix++){
            vx[nz+nbt+iz][ix] = vx[nz+nbt-1][ix];
            vz[nz+nbt+iz][ix] = vz[nz+nbt-1][ix];
        }
    }

    vx0 =0.0;
    vz0 =0.0;
    for (iz=0; iz < nzb; iz++) {
        for (ix=0; ix < nxb; ix++) {
            vx0 += vx[iz][ix]*vx[iz][ix];
            vz0 += vz[iz][ix]*vz[iz][ix];
         }
    }
    vx0 = sqrtf(vx0/(nxb*nzb));
    vz0 = sqrtf(vz0/(nxb*nzb));
    
    vx02=vx0*vx0; 
    vz02=vz0*vz0; 
    vx04=vx02*vx02; 
    vz04=vz02*vz02; 

    se0 *= pi/180.0;


    cosg0 = cosf(se0);
    cosg02 = cosg0*cosg0;
    sing0 = sinf(se0);
    sing02 = sing0*sing0; 

    dt2 = dt*dt;
    dx2 = dx*dx;
    dx4 = dx2*dx2;
    dz2 = dz*dz;
    dz4 = dz2*dz2;


    for (iz=0; iz < nzb; iz++) {
        for (ix=0; ix < nxb; ix++) {
            cur[iz][ix] = 0.0;
        }
    }
    cur[isz+nbt][isx+nbl] = wav[0];
    for (iz=0; iz < nzb; iz++) {
        for (ix=0; ix < nxb; ix++) {
            old[iz][ix] =  0.0; 
            derold[iz][ix] =cur[iz][ix]/dt;
           }
         }
    for (iz=nbt; iz<nz+nbt; iz++){
        sf_floatwrite(cur[iz]+nbl,nx,out);
    }
    /* propagation in time */
    for (it=1; it < nt; it++) {

         for (iz=0; iz < nzb; iz++){
             for (ix=0; ix < nxb; ix++){ 
                  new[iz][ix] = 0.0; 
               //   uk[iz][ix] = sf_cmplx(cur[iz][ix],0.0);
                  uk[iz][ix].r = cur[iz][ix];
                  uk[iz][ix].i = 0.0;
                }
         }  


/* compute u(kx,kz) */
         for (iz=0; iz < nzb; iz++){
             /* Fourier transform x to kx */
                for (ix=1; ix < nxb; ix+=2){
                    uk[iz][ix] = sf_cneg(uk[iz][ix]);
                    }
                kiss_fft_stride(cfgx,uk[iz],ctracex,1); 
                for (ikx=0; ikx<nkx; ikx++) uk[iz][ikx] = ctracex[ikx]; 
             }
         for (ikx=0; ikx < nkx; ikx++){
             /* Fourier transform z to kz */
                for (ikz=1; ikz<nkz; ikz+=2){
                    uk[ikz][ikx] = sf_cneg(uk[ikz][ikx]);
                    }
                kiss_fft_stride(cfgz,uk[0]+ikx,ctracez,nkx); 
                for (ikz=0; ikz<nkz; ikz++) uk[ikz][ikx] = ctracez[ikz]; 
             }


         for (ikz=0; ikz < nkz; ikz++) {
             kz1 = (kz0 + ikz* dkz)*2.0*pi;

             for (ikx=0; ikx < nkx; ikx++) {
                 kx1 = (kx0 + ikx * dkx)*2.0*pi;
                 kx = kx1*cosg0+kz1*sing0;
                 kz = kz1*cosg0-kx1*sing0;
                 tmpvk = (vx02*kx*kx+vz02*kz*kz);
                 vk2 = 0.5*tmpvk+0.5*sqrtf(tmpvk*tmpvk-8.0*yi0/(1.0+2.0*yi0)*vx02*vz02*kx*kx*kz*kz);
                 vk = sqrtf(vk2); 
                 tmpdt = 2.0*(cosf(vk*dt)-1.0);
                 uk[ikz][ikx] = sf_crmul(uk[ikz][ikx],tmpdt);
             }

         }   
/* Inverse FFT*/
              for (ikx=0; ikx < nkx; ikx++){
             /* Inverse Fourier transform kz to z */
                  kiss_fft_stride(cfgzi,(kiss_fft_cpx *)uk[0]+ikx,(kiss_fft_cpx *)ctracez,nkx); 
                  for (ikz=0; ikz < nkz; ikz++) uk[ikz][ikx] = sf_crmul(ctracez[ikz],ikz%2?-1.0:1.0); 
              }
              for (ikz=0; ikz < nkz; ikz++){
             /* Inverse Fourier transform kx to x */
                  kiss_fft_stride(cfgxi,(kiss_fft_cpx *)uk[ikz],(kiss_fft_cpx *)ctracex,1); 
                  for (ikx=0; ikx < nkx; ikx++) uk[ikz][ikx] = sf_crmul(ctracex[ikx],ikx%2?-1.0:1.0); 
              }

         for (iz=0; iz < nzb; iz++){
             for (ix=0; ix < nxb; ix++){ 
                  new[iz][ix] = sf_crealf(uk[iz][ix]); 
                  new[iz][ix] /= (nkx*nkz); 
                }
         }  
         
         new[isz+nbt][isx+nbl] += wav[it];

	 for (iz=0; iz < nzb; iz++) {  
             for (ix=0; ix < nxb; ix++) {
                 new[iz][ix] += 2.0*cur[iz][ix] - old[iz][ix]; 
             }
         }
 
                 
	 for (iz=0; iz < nzb; iz++) {  
             for(ix=0; ix < nxb; ix++) {
	        old[iz][ix] = cur[iz][ix]; 
	        cur[iz][ix] = new[iz][ix]; 
             }
         }
         for (iz=nbt; iz<nz+nbt; iz++){
             sf_floatwrite(new[iz]+nbl,nx,out);
         }  
    }
    bd_close();
    free(*vx);     
    free(*vz);     
    free(*new);     
    free(*cur);     
    free(*old);     
    free(*dercur);     
    free(*derold);     
    free(*uk);     
    free(vx);     
    free(vz);     
    free(new);     
    free(cur);     
    free(old);     
    free(dercur);     
    free(derold);     
    free(uk);     
 //   sf_fileclose(vel);
 //   sf_fileclose(inp);
 //   sf_fileclose(out);
 
    exit(0); 
}           
           
