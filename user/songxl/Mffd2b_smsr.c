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
#ifdef _OPENMP
#include <omp.h>
#endif
int main(int argc, char* argv[]) 
{
    int nx, nt, nkx, nkz, ix, it, ikx, ikz, nz, iz, nbt, nbb, nbl, nbr, nxb, nzb, isx, isz;
    float dt, dx, dkx, kx, dz, dkz, kz, tmpdt, pi=SF_PI, o1, o2;
    float **new,  **old,  **cur,  **uk, **dercur, **derold, *wav;
    float  **v, v0, ***aa, w, g1, g2, ct, cb, cl, cr; //top, bottom, left, right 
    sf_file out, vel, source;
    bool opt;    /* optimal padding */
   // #ifdef _OPENMP
   // int nth;
   // #endif
     

    sf_init(argc,argv);
    out = sf_output("out");
    vel = sf_input("vel");   /* velocity */
    source = sf_input("in");   /* source wavlet*/

//    if (SF_FLOAT != sf_gettype(inp)) sf_error("Need float input");
    if (SF_FLOAT != sf_gettype(vel)) sf_error("Need float input");
    if (SF_FLOAT != sf_gettype(source)) sf_error("Need float input");
    if (!sf_histint(vel,"n1",&nx)) sf_error("No n1= in input");
    if (!sf_histfloat(vel,"d1",&dx)) sf_error("No d1= in input");
    if (!sf_histint(vel,"n2",&nz)) sf_error("No n2= in input");
    if (!sf_histfloat(vel,"d2",&dz)) sf_error("No d2= in input");
    if (!sf_histfloat(vel,"o1",&o1)) o1=0.0;
    if (!sf_histfloat(vel,"o2",&o2)) o1=0.0;
  //  if (!sf_histint(inp,"n2",&nt)) sf_error("No n2= in input");
  //  if (!sf_histfloat(inp,"d2",&dt)) sf_error("No d2= in input");
    if (!sf_getbool("opt",&opt)) opt=true;
    /* if y, determine optimal size for efficiency */
    if (!sf_getfloat("dt",&dt)) sf_error("Need dt input");
    if (!sf_getint("nt",&nt)) sf_error("Need nt input");
    if (!sf_getint("isx",&isx)) sf_error("Need isx input");
    if (!sf_getint("isz",&isz)) sf_error("Need isz input");

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
    nkx = nxb;
    nkz = nzb;
    dkx = 1./(2.0*kiss_fft_next_fast_size(nxb-1)*dx);
    dkz = 1./(2.0*kiss_fft_next_fast_size(nzb-1)*dz);



    wav    =  sf_floatalloc(nt);
    sf_floatread(wav,nt,source);

    old    =  sf_floatalloc2(nxb,nzb);
    cur    =  sf_floatalloc2(nxb,nzb);
    new    =  sf_floatalloc2(nxb,nzb);
    uk     =  sf_floatalloc2(nxb,nzb);
    derold    =  sf_floatalloc2(nxb,nzb);
    dercur    =  sf_floatalloc2(nxb,nzb);
    aa     =  sf_floatalloc3(3,nxb,nzb);
    
    bd_init(nx,nz,nbt,nbb,nbl,nbr,ct,cb,cl,cr);

    v = sf_floatalloc2(nxb,nzb);
    /*input & extend velocity model*/
    for (iz=nbt; iz<nz+nbt; iz++){
        sf_floatread(v[iz]+nbl,nx,vel);
         for (ix=0; ix<nbl; ix++){
             v[iz][ix] = v[iz][nbl];
         }
         for (ix=0; ix<nbr; ix++){
             v[iz][nx+nbl+ix] = v[iz][nx+nbl-1];
         }     
    }
    for (iz=0; iz<nbt; iz++){
        for (ix=0; ix<nxb; ix++){
            v[iz][ix] = v[nbt][ix];
        }
    }
    for (iz=0; iz<nbb; iz++){
        for (ix=0; ix<nxb; ix++){
            v[nz+nbt+iz][ix] = v[nz+nbt-1][ix];
        }
    }

    v0 =0.0;
    for (iz=0; iz < nzb; iz++) {
        for (ix=0; ix < nxb; ix++) {
            v0 += v[iz][ix]*v[iz][ix];
         }
    }
    v0 = sqrtf(v0/(nxb*nzb));
 
    for (iz=0; iz < nzb; iz++){
         for (ix=0; ix < nxb; ix++) {
         w = dt*dt*v[iz][ix]*v[iz][ix];
         g1 = dt*dt*(v[iz][ix]*v[iz][ix]-v0*v0)/(12.0*dx*dx);
         g2 = dt*dt*(v[iz][ix]*v[iz][ix]-v0*v0)/(12.0*dz*dz);
         aa[iz][ix][1] = w*g1;
         aa[iz][ix][2] = w*g2;
         aa[iz][ix][0] = w-2.0*aa[iz][ix][1]-2.0*aa[iz][ix][2];
        }
      }

    for (iz=0; iz < nzb; iz++) {
        for (ix=0; ix < nxb; ix++) {
            cur[iz][ix] = 0.0;
        }
    }
    cur[isz+nbt][isx+nbl] = wav[0];
    srcsm_init(dz,dx);
    source_smooth(cur,isz+nbt,isx+nbl,wav[0]);
    for (iz=0; iz < nzb; iz++) {
        for (ix=0; ix < nxb; ix++) {
            old[iz][ix] =  0.0; 
            derold[iz][ix] =cur[iz][ix]/dt;
           }
         }
    for (iz=nbt; iz<nz+nbt; iz++){
        sf_floatwrite(cur[iz]+nbl,nx,out);
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

         for (iz=0; iz < nzb; iz++){
             for (ix=0; ix < nxb; ix++){ 
                  new[iz][ix] = 0.0; 
                  uk[iz][ix] = cur[iz][ix]; 
                }
         }  


/* compute u(kx,kz) */
         sf_cosft_init(nxb);
         for (iz=0; iz < nzb; iz++){
             /* Fourier transform x to kx */
             sf_cosft_frw(uk[iz],0,1);
         }
         sf_cosft_close();

         sf_cosft_init(nzb);
         for (ikx=0; ikx < nkx; ikx++){
             /* Fourier transform z to kz */
             sf_cosft_frw(uk[0],ikx,nxb);
         }
         sf_cosft_close();

/*    #ifdef _OPENMP
    #pragma omp parallel for private(ik,ix,x,k,tmp,tmpex,tmpdt) 
    #endif
*/

         for (ikz=0; ikz < nkz; ikz++) {
             kz = ikz* dkz*2.0*pi;

             for (ikx=0; ikx < nkx; ikx++) {
                 kx = ikx * dkx*2.0*pi;
                 tmpdt = 2.0*(cosf(v0*sqrtf(kx*kx+kz*kz)*dt)-1.0)/(v0*v0*dt*dt);
                 uk[ikz][ikx] = uk[ikz][ikx]*tmpdt;
             }

         }   
/* Inverse FFT*/
         sf_cosft_init(nzb);
         for (ikx=0; ikx < nkx; ikx++){
             /* Inverse Fourier transform kz to z */
             sf_cosft_inv(uk[0],ikx,nxb);
         }
         sf_cosft_close();
         sf_cosft_init(nxb);
         for (iz=0; iz < nzb; iz++){
             /* Inverse Fourier transform kx to x */
              sf_cosft_inv(uk[iz],0,1);
         }
         sf_cosft_close();

	 for (iz=1; iz < nzb-1; iz++) {  
	     for (ix=1; ix < nxb-1; ix++) {  
                 new[iz][ix]  = uk[iz][ix]*aa[iz][ix][0]
                              + (uk[iz][ix-1]+uk[iz][ix+1])*aa[iz][ix][1]
                              + (uk[iz-1][ix]+uk[iz+1][ix])*aa[iz][ix][2];
             }
         }  
    
         new[0][0] = uk[0][0]*aa[0][0][0] + uk[0][1]*aa[0][0][1] + uk[1][0]*aa[0][0][2];
         new[0][nxb-1] = uk[0][nxb-1]*aa[0][nxb-1][0] + uk[0][nxb-2]*aa[0][nxb-1][1] + uk[1][nxb-1]*aa[0][nxb-1][2];
         new[nzb-1][0] = uk[nzb-1][0]*aa[nzb-1][0][0] + uk[nzb-1][1]*aa[nzb-1][0][1] + uk[nzb-2][0]*aa[nzb-1][0][2];
         new[nzb-1][nxb-1] = uk[nzb-1][nxb-1]*aa[nzb-1][nxb-1][0] + uk[nzb-1][nxb-2]*aa[nzb-1][nxb-1][1] + uk[nzb-2][nxb-1]*aa[nzb-1][nxb-1][2];
          
	 for (ix=1; ix < nxb-1; ix++) {  
             new[0][ix] = uk[0][ix]*aa[0][ix][0] + (uk[0][ix-1]+uk[0][ix+1])*aa[0][ix][1] + uk[1][ix]*aa[0][ix][2];
             new[nz-1][ix] = uk[nz-1][ix]*aa[nz-1][ix][0] + (uk[nz-1][ix-1]+uk[nz-1][ix+1])*aa[nz-1][ix][1] + uk[nz-2][ix]*aa[nz-1][ix][2];
         }
	 for (iz=1; iz < nzb-1; iz++) {  
             new[iz][0] = uk[iz][0]*aa[iz][0][0] + uk[iz][1]*aa[iz][0][1] + (uk[iz-1][0]+uk[iz+1][0])*aa[iz][0][2]; 
             new[iz][nx-1] = uk[iz][nx-1]*aa[iz][nx-1][0] + uk[iz][nx-2]*aa[iz][nx-1][1] + (uk[iz-1][nx-1]+uk[iz+1][nx-1])*aa[iz][nx-1][2]; 
         }
          
         new[isz+nbt][isx+nbl] += wav[it];
         source_smooth(new,isz+nbt,isx+nbl,wav[it]);

	 for (iz=0; iz < nzb; iz++) {  
             for (ix=0; ix < nxb; ix++) {
                 dercur[iz][ix]= derold[iz][ix] + new[iz][ix]/dt;
                 new[iz][ix] = cur[iz][ix] + dercur[iz][ix]*dt; 
             }
         }
 
    //     new[isz+nb][isx+nb] += wav[it];
         bd_decay(new); 
         bd_decay(dercur); 
                 
	 for (iz=0; iz < nzb; iz++) {  
             for(ix=0; ix < nxb; ix++) {
	        old[iz][ix] = cur[iz][ix]; 
	        cur[iz][ix] = new[iz][ix]; 
	        derold[iz][ix] = dercur[iz][ix]; 
             }
         }
         for (iz=nbt; iz<nz+nbt; iz++){
             sf_floatwrite(new[iz]+nbl,nx,out);
         }  
    }
    bd_close();
    srcsm_close();
    free(**aa);
    free(*aa);
    free(aa);
    free(*v);     
    free(*new);     
    free(*cur);     
    free(*old);     
    free(*dercur);     
    free(*derold);     
    free(*uk);     
    free(v);     
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
           
