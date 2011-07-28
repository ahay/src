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
    float dt, dx, dkx, kx, dz, dkz, kz, tmpdt, pi=SF_PI, o1, o2;
    float **nxt,  **old,  **cur,  **uk, **dercur, **derold, *wav;
    float **vx, vx2, vx4, vx0, vx02, vx04, **vz, vz2, vz4, vz0, vz02, vz04, **yi, yi0;
    float ***aa, w, w0, g11, g22, dx2, dz2, dx4, dz4, ct, cb, cl, cr; //top, bottom, left, right 
    float vk, vk2, tmpvk, k2, err, dt2;
    sf_file out, velx, velz, source, yita;
    bool opt;    /* optimal padding */
   // #ifdef _OPENMP
   // int nth;
   // #endif
     

    sf_init(argc,argv);
    out = sf_output("out");
    velx = sf_input("velx");   /* velocity */
    velz = sf_input("velz");   /* velocity */
    yita = sf_input("yita");   /* anistropic parameter*/
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

    if (!sf_getint("nbt",&nbt)) nbt=44;
    if (!sf_getint("nbb",&nbb)) nbb=44;
    if (!sf_getint("nbl",&nbl)) nbl=44;
    if (!sf_getint("nbr",&nbr)) nbr=44;

    if (!sf_getfloat("ct",&ct)) ct = 0.01; /*decaying parameter*/
    if (!sf_getfloat("cb",&cb)) cb = 0.01; /*decaying parameter*/
    if (!sf_getfloat("cl",&cl)) cl = 0.01; /*decaying parameter*/
    if (!sf_getfloat("cr",&cr)) cr = 0.01; /*decaying parameter*/

    err = 0.00000001;

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
    nxt    =  sf_floatalloc2(nxb,nzb);
    uk     =  sf_floatalloc2(nxb,nzb);
    derold    =  sf_floatalloc2(nxb,nzb);
    dercur    =  sf_floatalloc2(nxb,nzb);
    aa     =  sf_floatalloc3(6,nxb,nzb);
    
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
    /*input & extend anistropic model*/
    yi = sf_floatalloc2(nxb,nzb);
    for (iz=nbt; iz<nz+nbt; iz++){
        sf_floatread(yi[iz]+nbl,nx,yita);
         for (ix=0; ix<nbl; ix++){
             yi[iz][ix] = yi[iz][nbl];
         }
         for (ix=0; ix<nbr; ix++){
             yi[iz][nx+nbl+ix] = yi[iz][nx+nbl-1];
         }     
    }
    for (iz=0; iz<nbt; iz++){
        for (ix=0; ix<nxb; ix++){
            yi[iz][ix] = yi[nbt][ix];
        }
    }
    for (iz=0; iz<nbb; iz++){
        for (ix=0; ix<nxb; ix++){
            yi[nz+nbt+iz][ix] = yi[nz+nbt-1][ix];
        }
    }

    yi0 = 0.0;
    for (iz=0; iz < nzb; iz++) {
        for (ix=0; ix < nxb; ix++) {
            yi0+= yi[iz][ix]*yi[iz][ix];
         }
    }
    yi0 = sqrtf(yi0/(nxb*nzb));
    
    w0 = sqrtf(((1.0+2.0*yi0)*(vx04+vz04)+(2.0-4.0*yi0)*vx02*vz02)/(1.0+2.0*yi0));
    dt2 = dt*dt;
    dx2 = dx*dx;
    dx4 = dx2*dx2;
    dz2 = dz*dz;
    dz4 = dz2*dz2;
    for (iz=0; iz < nzb; iz++){
         for (ix=0; ix < nxb; ix++) {
         vx2 = vx[iz][ix]*vx[iz][ix];
         vx4 = vx2*vx2;
         vz2 = vz[iz][ix]*vz[iz][ix];
         vz4 = vz2*vz2;
         w = sqrtf(((1.0+2.0*yi[iz][ix])*(vx4+vz4)+(2.0-4.0*yi[iz][ix])*vx2*vz2)/(1.0+2.0*yi[iz][ix]));
         g11 = vx2/vx02;
         g22 = vz2/vz02;
         aa[iz][ix][4] = (-dt2*vx4+dt2*vx2*vx02+vx2*dx2)/(12.0*vx02*dx4);
         aa[iz][ix][5] = (-dt2*vz4+dt2*vz2*vz02+vz2*dz2)/(12.0*vz02*dz4);
         aa[iz][ix][3] = -(-dt2*(vx2+vz2+w)+dt2*(vx2+vz2+w)*(vx2+vz2+w)/(vx02+vz02+w0)-g11*dx2-g22*dz2+12.0*aa[iz][ix][4]*dx4+12.0*aa[iz][ix][5]*dz4)/(12.0*dx2*dz2);
         aa[iz][ix][1] = -g11/dx2-2.0*aa[iz][ix][3]-4.0*aa[iz][ix][4];
         aa[iz][ix][2] = -g22/dz2-2.0*aa[iz][ix][3]-4.0*aa[iz][ix][5];
         aa[iz][ix][0] = -2.0*aa[iz][ix][1]-2.0*aa[iz][ix][2]-4.0*aa[iz][ix][3]-2.0*aa[iz][ix][4]-2.0*aa[iz][ix][5];
        }
      }

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
                  nxt[iz][ix] = 0.0; 
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
                 tmpvk = (vx02*kx*kx+vz02*kz*kz);
                 k2 = kx*kx+kz*kz;
                 vk2 = 0.5*tmpvk+0.5*sqrtf(tmpvk*tmpvk-8.0*yi0/(1+2.0*yi0)*vx02*vz02*kx*kx*kz*kz);
                 vk = sqrtf(vk2);
                 tmpdt = 2.0*(cosf(vk*dt)-1.0);
                 if((!ikx)&&(!ikz)) 
                   tmpdt /=(k2+err);
                 else 
                   tmpdt /= (k2);
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

	 for (iz=2; iz < nzb-2; iz++) {  
	     for (ix=2; ix < nxb-2; ix++) {  
                 nxt[iz][ix]  = uk[iz][ix]*aa[iz][ix][0]
                              + (uk[iz][ix-1]+uk[iz][ix+1])*aa[iz][ix][1]
                              + (uk[iz-1][ix]+uk[iz+1][ix])*aa[iz][ix][2]
                              + (uk[iz-1][ix-1]+uk[iz-1][ix+1]+uk[iz+1][ix-1]+uk[iz+1][ix+1])*aa[iz][ix][3]
                              + (uk[iz][ix-2]+uk[iz][ix+2])*aa[iz][ix][4]
                              + (uk[iz-2][ix]+uk[iz+2][ix])*aa[iz][ix][5];
             }
         }  
         
/* 
         nxt[0][0] = uk[0][0]*aa[0][0][0] + uk[0][1]*aa[0][0][1] + uk[1][0]*aa[0][0][2];
         nxt[0][nxb-1] = uk[0][nxb-1]*aa[0][nxb-1][0] + uk[0][nxb-2]*aa[0][nxb-1][1] + uk[1][nxb-1]*aa[0][nxb-1][2];
         nxt[nzb-1][0] = uk[nzb-1][0]*aa[nzb-1][0][0] + uk[nzb-1][1]*aa[nzb-1][0][1] + uk[nzb-2][0]*aa[nzb-1][0][2];
         nxt[nzb-1][nxb-1] = uk[nzb-1][nxb-1]*aa[nzb-1][nxb-1][0] + uk[nzb-1][nxb-2]*aa[nzb-1][nxb-1][1] + uk[nzb-2][nxb-1]*aa[nzb-1][nxb-1][2];
          
	 for (ix=1; ix < nxb-1; ix++) {  
             nxt[0][ix] = uk[0][ix]*aa[0][ix][0] + (uk[0][ix-1]+uk[0][ix+1])*aa[0][ix][1] + uk[1][ix]*aa[0][ix][2];
             nxt[nz-1][ix] = uk[nz-1][ix]*aa[nz-1][ix][0] + (uk[nz-1][ix-1]+uk[nz-1][ix+1])*aa[nz-1][ix][1] + uk[nz-2][ix]*aa[nz-1][ix][2];
         }
	 for (iz=1; iz < nzb-1; iz++) {  
             nxt[iz][0] = uk[iz][0]*aa[iz][0][0] + uk[iz][1]*aa[iz][0][1] + (uk[iz-1][0]+uk[iz+1][0])*aa[iz][0][2]; 
             nxt[iz][nx-1] = uk[iz][nx-1]*aa[iz][nx-1][0] + uk[iz][nx-2]*aa[iz][nx-1][1] + (uk[iz-1][nx-1]+uk[iz+1][nx-1])*aa[iz][nx-1][2]; 
         }
 */         
         nxt[isz+nbt][isx+nbl] += wav[it];

	 for (iz=0; iz < nzb; iz++) {  
             for (ix=0; ix < nxb; ix++) {
            //     dercur[iz][ix]= derold[iz][ix] + nxt[iz][ix]/dt;
            //     nxt[iz][ix] = cur[iz][ix] + dercur[iz][ix]*dt; 
                 nxt[iz][ix] += 2.0*cur[iz][ix] -old[iz][ix]; 
             }
         }
 
    //     nxt[isz+nb][isx+nb] += wav[it];
    //     bd_decay(nxt); 
    //     bd_decay(dercur); 
                 
	 for (iz=0; iz < nzb; iz++) {  
             for(ix=0; ix < nxb; ix++) {
	        old[iz][ix] = cur[iz][ix]; 
	        cur[iz][ix] = nxt[iz][ix]; 
	        derold[iz][ix] = dercur[iz][ix]; 
             }
         }
         for (iz=nbt; iz<nz+nbt; iz++){
             sf_floatwrite(nxt[iz]+nbl,nx,out);
         }  
    }
    bd_close();
    free(**aa);
    free(*aa);
    free(aa);
    free(*vx);     
    free(*vz);     
    free(*yi);     
    free(*nxt);     
    free(*cur);     
    free(*old);     
    free(*dercur);     
    free(*derold);     
    free(*uk);     
    free(vx);     
    free(vz);     
    free(yi);     
    free(nxt);     
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
           
