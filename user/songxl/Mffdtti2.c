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
#include "fft2.h"
#ifdef _OPENMP
#include <omp.h>
#endif
int main(int argc, char* argv[]) 
{
    int nx, nt, nkx, nkz, ix, it, ikx, ikz, nz, iz, nbt, nbb, nbl, nbr, nxb, nzb, isx, isz;
    float dt, dx, dkx, kx, dz, dkz, kz, tmpdt, pi=SF_PI, o1, o2, kx0, kz0;
    //float **nxt,  **old,  **cur, **ukr, **dercur, **derold, *wav;
    float **nxt,  **old,  **cur,  **dercur, **derold, *wav, *tmpcur;
    float **vx, vx2, /* vx4, */ vx0, vx02, /* vx04, */ **vz, vz2, /* vz4, */ vz0, vz02, /* vz04, */ **yi, yi0, **se, se0;
    float ***aa, dx2, dz2, dx4, dz4, ct, cb, cl, cr; /* top, bottom, left, right */
    float w1, w10, w2, w20, w3, w30, h1, h10, h2, h20, h3, h30;
    float cosg, cosg0, cosg2, cosg02, sing, sing0, sing2, sing02;
    float vk, vk2, tmpvk, k2, err, dt2, kx1, kz1;
    float alpha; /* source smoothing */
    //kiss_fft_cpx **uk, *ctracex, *ctracez;
    sf_complex *cwave;
   // kiss_fft_cfg cfgx, cfgxi, cfgz, cfgzi;
    sf_file out, velx, velz, source, yita, seta;
    bool opt;    /* optimal padding */     
    int nk, pad1, n;
    bool cmplx;

    sf_init(argc,argv);
    out = sf_output("out");
    velx = sf_input("velx");   /* velocity */
    velz = sf_input("velz");   /* velocity */
    yita = sf_input("yita");   /* anistropic parameter*/
    source = sf_input("in");   /* source wavlet*/
    seta = sf_input("seta");   /* TTI angle*/

/*    if (SF_FLOAT != sf_gettype(inp)) sf_error("Need float input"); */
    if (SF_FLOAT != sf_gettype(velx)) sf_error("Need float input");
    if (SF_FLOAT != sf_gettype(velz)) sf_error("Need float input");
    if (SF_FLOAT != sf_gettype(source)) sf_error("Need float input");
    if (SF_FLOAT != sf_gettype(seta)) sf_error("Need float input");

    if (!sf_histint(velx,"n1",&nx)) sf_error("No n1= in input");
    if (!sf_histfloat(velx,"d1",&dx)) sf_error("No d1= in input");
    if (!sf_histint(velx,"n2",&nz)) sf_error("No n2= in input");
    if (!sf_histfloat(velx,"d2",&dz)) sf_error("No d2= in input");
    if (!sf_histfloat(velx,"o1",&o1)) o1=0.0;
    if (!sf_histfloat(velx,"o2",&o2)) o2=0.0;
    /*  if (!sf_histint(inp,"n2",&nt)) sf_error("No n2= in input"); */
    /*  if (!sf_histfloat(inp,"d2",&dt)) sf_error("No d2= in input"); */
    if (!sf_getbool("opt",&opt)) opt=true;
    /* if y, determine optimal size for efficiency */
    if (!sf_getfloat("dt",&dt)) sf_error("Need dt input");
    if (!sf_getint("nt",&nt)) sf_error("Need nt input");
    if (!sf_getint("isx",&isx)) sf_error("Need isx input");
    if (!sf_getint("isz",&isz)) isz=0;
    if (!sf_getfloat("err",&err)) err = 0.0001;
    if (!sf_getfloat("alpha",&alpha)) alpha=-0.7;

    if (!sf_getint("nbt",&nbt)) nbt=126;
    if (!sf_getint("nbb",&nbb)) nbb=126;
    if (!sf_getint("nbl",&nbl)) nbl=128;
    if (!sf_getint("nbr",&nbr)) nbr=128;

    if (!sf_getfloat("ct",&ct)) ct = 0.002; /*decaying parameter*/
    if (!sf_getfloat("cb",&cb)) cb = 0.002; /*decaying parameter*/
    if (!sf_getfloat("cl",&cl)) cl = 0.002; /*decaying parameter*/
    if (!sf_getfloat("cr",&cr)) cr = 0.002; /*decaying parameter*/

    
    if (!sf_getint("pad1",&pad1)) pad1=1; /* padding factor on the first axis */
    //if (!sf_getbool("cmplx",&cmplx)) cmplx=true; /* use complex FFT */
    //if (!sf_getbool("cmplx",&cmplx)) cmplx=true; /* use complex FFT */
    if (!sf_getbool("cmplx",&cmplx)) cmplx=false; /* use complex FFT */

    sf_putint(out,"n1",nx);
    sf_putfloat(out,"d1",dx);
/*    sf_putfloat(out,"o1",x0); */
    sf_putint(out,"n2",nz);
    sf_putfloat(out,"d2",dz);
    sf_putint(out,"n3",nt);
    sf_putfloat(out,"d3",dt);
    sf_putfloat(out,"o1",o1); 
    sf_putfloat(out,"o2",o2); 
    sf_putfloat(out,"o3",0.0); 

    nxb = nx + nbl + nbr;
    nzb = nz + nbt + nbb;

    nk = fft2_init(cmplx,pad1,nxb,nzb,&nkx,&nkz);

    sf_warning("nx=%d\n",nx);
    sf_warning("nz=%d\n",nz);
    sf_warning("nxb=%d\n",nxb);
    sf_warning("nzb=%d\n",nzb);
    sf_warning("nkx=%d\n",nkx);
    sf_warning("nkz=%d\n",nkz);

    dkx = 1./(nkx*dx);
    kx0 = -0.5/dx;
    dkz = 1./(nkz*dz);
    kz0 = -0.5/dz;
/*
    nkx = opt? kiss_fft_next_fast_size(nxb): nxb;
    nkz = opt? kiss_fft_next_fast_size(nzb): nzb;
    if (nkx != nxb) sf_warning("nkx padded to %d",nkx);
    if (nkz != nzb) sf_warning("nkz padded to %d",nkz);
    cfgx = kiss_fft_alloc(nkx,0,NULL,NULL);
    cfgxi = kiss_fft_alloc(nkx,1,NULL,NULL);
    cfgz = kiss_fft_alloc(nkz,0,NULL,NULL);
    cfgzi = kiss_fft_alloc(nkz,1,NULL,NULL);


    uk = (kiss_fft_cpx **) sf_complexalloc2(nkx,nkz);
    ctracex = (kiss_fft_cpx *) sf_complexalloc(nkx);
    ctracez = (kiss_fft_cpx *) sf_complexalloc(nkz);
*/
    wav    =  sf_floatalloc(nt);
    sf_floatread(wav,nt,source);

    old    =  sf_floatalloc2(nxb,nzb);
    cur    =  sf_floatalloc2(nxb,nzb);
    tmpcur    =  sf_floatalloc(nxb*nzb);
    cwave  = sf_complexalloc(nk);
    nxt    =  sf_floatalloc2(nxb,nzb);
    derold    =  sf_floatalloc2(nxb,nzb);
    dercur    =  sf_floatalloc2(nxb,nzb);
    aa     =  sf_floatalloc3(6,nxb,nzb);
    
    ifft2_allocate(cwave);
    
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
/*    vx04=vx02*vx02; 
      vz04=vz02*vz02; */

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
    
    se = sf_floatalloc2(nxb,nzb);
    for (iz=nbt; iz<nz+nbt; iz++){
        sf_floatread(se[iz]+nbl,nx,seta);
         for (ix=0; ix<nbl; ix++){
             se[iz][ix] = se[iz][nbl];
         }
         for (ix=0; ix<nbr; ix++){
             se[iz][nx+nbl+ix] = se[iz][nx+nbl-1];
         }     
    }
    for (iz=0; iz<nbt; iz++){
        for (ix=0; ix<nxb; ix++){
            se[iz][ix] = se[nbt][ix];
        }
    }
    for (iz=0; iz<nbb; iz++){
        for (ix=0; ix<nxb; ix++){
            se[nz+nbt+iz][ix] = se[nz+nbt-1][ix];
        }
    }

    se0 = 0.0;
    for (iz=0; iz < nzb; iz++) {
        for (ix=0; ix < nxb; ix++) {
            se0+= se[iz][ix];
         }
    }
    se0 /= (nxb*nzb);
    se0 *= pi/180.0;

    for (iz=0; iz < nzb; iz++) {
        for (ix=0; ix < nxb; ix++) {
            se[iz][ix] *= pi/180.0;
        }
    }

    cosg0 = cosf(se0);
    cosg02 = cosg0*cosg0;
    sing0 = sinf(se0);
    sing02 = sing0*sing0; 

    w10 = vx02*cosg02+vz02*sing02;
    w20 = vz02*cosg02+vx02*sing02;
    w30 = vx02+vz02+(vx02-vz02)*sinf(2.0*se0);
    h10 = sqrtf(-8.0*yi0*vx02*vz02*cosg02*sing02/(1.0+2.0*yi0)+w10*w10);
    h20 = sqrtf(-8.0*yi0*vx02*vz02*cosg02*sing02/(1.0+2.0*yi0)+w20*w20);
    h30 = sqrtf(-2.0*yi0*vx02*vz02*cosf(2.0*se0)*cosf(2.0*se0)/(1.0+2.0*yi0)+0.25*w30*w30);
    dt2 = dt*dt;
    dx2 = dx*dx;
    dx4 = dx2*dx2;
    dz2 = dz*dz;
    dz4 = dz2*dz2;

    for (iz=0; iz < nzb; iz++){
         for (ix=0; ix < nxb; ix++) {
         vx2 = vx[iz][ix]*vx[iz][ix];
         /* vx4 = vx2*vx2; */
         vz2 = vz[iz][ix]*vz[iz][ix];
         /* vz4 = vz2*vz2; */
         cosg = cosf(se[iz][ix]);
         sing = sinf(se[iz][ix]);
         cosg2 = cosg*cosg;
         sing2 = sing*sing;
         w1 = vx2*cosg2+vz2*sing2;
         w2 = vz2*cosg2+vx2*sing2;
         w3 = vx2+vz2+(vx2-vz2)*sinf(2.0*se[iz][ix]);
         h1 = sqrtf(-8.0*yi[iz][ix]*vx2*vz2*cosg2*sing2/(1.0+2.0*yi[iz][ix])+w1*w1);
         h2 = sqrtf(-8.0*yi[iz][ix]*vx2*vz2*cosg2*sing2/(1.0+2.0*yi[iz][ix])+w2*w2);
         h3 = sqrtf(-2.0*yi[iz][ix]*vx2*vz2*cosf(2.0*se[iz][ix])*cosf(2.0*se[iz][ix])/(1.0+2.0*yi[iz][ix])+0.25*w3*w3);
         aa[iz][ix][4] = (w1+h1)*(dt2+(2.0*dx2-dt2*(w1+h1))/(w10+h10))/(24.0*dx4);
         aa[iz][ix][5] = (w2+h2)*(dt2+(2.0*dz2-dt2*(w2+h2))/(w20+h20))/(24.0*dz4);
         aa[iz][ix][3] = -aa[iz][ix][4]*dx2/dz2-aa[iz][ix][5]*dz2/dx2+(dt2*(w3+2.0*h3)+dx2*(w1+h1)/(w10+h10)+dz2*(w2+h2)/(w20+h20)-dt2*(w3+2.0*h3)*(w3+2.0*h3)/(w30+2.0*h30))/(12.0*dx2*dz2);
         aa[iz][ix][1] = -2.0*aa[iz][ix][3]-4.0*aa[iz][ix][4]-(w1+h1)/(dx2*(w10+h10));
         aa[iz][ix][2] = -2.0*aa[iz][ix][3]-4.0*aa[iz][ix][5]-(w2+h2)/(dz2*(w20+h20));
         aa[iz][ix][0] = -2.0*aa[iz][ix][1]-2.0*aa[iz][ix][2]-4.0*aa[iz][ix][3]-2.0*aa[iz][ix][4]-2.0*aa[iz][ix][5];
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

    for (iz=0; iz < nzb; iz++) {
        n=nxb*iz;
        for (ix=0; ix < nxb; ix++) {
            tmpcur[ix+n] = cur[iz][ix];
        }
    }
         for (iz=0; iz < nzb; iz++){
             for (ix=0; ix < nxb; ix++){ 
                  nxt[iz][ix] = 0.0; 
            //      uk[iz][ix].r = cur[iz][ix];
            //      uk[iz][ix].i = 0.0; 
                }
         }  
         fft2(tmpcur,cwave);


/* compute u(kx,kz) */
 //        for (iz=0; iz < nzb; iz++){
             /* Fourier transform x to kx */
             /*
                for (ix=1; ix < nxb; ix+=2){
                    uk[iz][ix] = sf_cneg(uk[iz][ix]);
                    }
                kiss_fft_stride(cfgx,uk[iz],ctracex,1); 
                for (ikx=0; ikx<nkx; ikx++) uk[iz][ikx] = ctracex[ikx]; 
             }
*/
  //       for (ikx=0; ikx < nkx; ikx++){
             /* Fourier transform z to kz */
             /*
                for (ikz=1; ikz<nkz; ikz+=2){
                    uk[ikz][ikx] = sf_cneg(uk[ikz][ikx]);
                    }
                kiss_fft_stride(cfgz,uk[0]+ikx,ctracez,nkx); 
                for (ikz=0; ikz<nkz; ikz++) uk[ikz][ikx] = ctracez[ikz]; 
             }
*/

/*    #ifdef _OPENMP
    #pragma omp parallel for private(ik,ix,x,k,tmp,tmpex,tmpdt) 
    #endif
*/

         for (ikz=0; ikz < nkz; ikz++) {
             kz1 = (kz0+ikz*dkz)*2.0*pi;
                 n = ikz*nkx;

             for (ikx=0; ikx < nkx; ikx++) {
                 kx1 = (kx0+ikx*dkx)*2.0*pi;
                 kx = kx1*cosg0+kz1*sing0;
                 kz = kz1*cosg0-kx1*sing0;
                 tmpvk = (vx02*kx*kx+vz02*kz*kz);
                 k2 = kx1*kx1+kz1*kz1;
                 vk2 = 0.5*tmpvk+0.5*sqrtf(tmpvk*tmpvk-8.0*yi0/(1.0+2.0*yi0)*vx02*vz02*kx*kx*kz*kz);
                 vk = sqrtf(vk2);
                 tmpdt = 2.0*(cosf(vk*dt)-1.0);
                 if(!k2) 
                   tmpdt /=(k2+err);
                 else 
                   tmpdt /= (k2);
                 //uk[ikz][ikx] = sf_crmul(uk[ikz][ikx],tmpdt);
                 //n = ikz*nkx+ikx;
 #ifdef SF_HAS_COMPLEX_H
                 cwave[n+ikx] = cwave[n+ikx]*tmpdt;
#else
                 cwave[n+ikx] = sf_crmul(cwave[n+ikx],tmpdt);
#endif
             }

         }   

/* Inverse FFT*/
         ifft2(tmpcur,cwave);

//         for (ikx=0; ikx < nkx; ikx++){
         /* Inverse Fourier transform kz to z */
/*
             kiss_fft_stride(cfgzi,(kiss_fft_cpx *)uk[0]+ikx,(kiss_fft_cpx *)ctracez,nkx); 
             for (ikz=0; ikz < nkz; ikz++) uk[ikz][ikx] = sf_crmul(ctracez[ikz],ikz%2?-1.0:1.0); 
              }
             for (ikz=0; ikz < nkz; ikz++){
*/
             /* Inverse Fourier transform kx to x */
/*
                 kiss_fft_stride(cfgxi,(kiss_fft_cpx *)uk[ikz],(kiss_fft_cpx *)ctracex,1); 
                 for (ikx=0; ikx < nkx; ikx++) uk[ikz][ikx] = sf_crmul(ctracex[ikx],ikx%2?-1.0:1.0); 
             }

         for (iz=0; iz < nzb; iz++){
             for (ix=0; ix < nxb; ix++){ 
                  ukr[iz][ix] = sf_crealf(uk[iz][ix]); 
                  ukr[iz][ix] /= (nkx*nkz); 
                }
         }  
*/

	 for (iz=2; iz < nzb-2; iz++) {  
	     for (ix=2; ix < nxb-2; ix++) {  
                 nxt[iz][ix]  = tmpcur[iz*nxb+ix]*aa[iz][ix][0]
                              + (tmpcur[iz*nxb+ix-1]+tmpcur[iz*nxb+ix+1])*aa[iz][ix][1]
                              + (tmpcur[(iz-1)*nxb+ix]+tmpcur[(iz+1)*nxb+ix])*aa[iz][ix][2]
                              + (tmpcur[(iz-1)*nxb+ix-1]+tmpcur[(iz-1)*nxb+ix+1]+tmpcur[(iz+1)*nxb+ix-1]+tmpcur[(iz+1)*nxb+ix+1])*aa[iz][ix][3]
                              + (tmpcur[iz*nxb+ix-2]+tmpcur[iz*nxb+ix+2])*aa[iz][ix][4]
                              + (tmpcur[(iz-2)*nxb+ix]+tmpcur[(iz+2)*nxb+ix])*aa[iz][ix][5];
             }
         }  
         

	 for (iz=0; iz < nzb; iz++) {  
             for (ix=0; ix < nxb; ix++) {
                 dercur[iz][ix]= derold[iz][ix] + nxt[iz][ix]/dt;
                 nxt[iz][ix] = cur[iz][ix] + dercur[iz][ix]*dt; 
		 /*     nxt[iz][ix] += 2.0*cur[iz][ix] -old[iz][ix]; */
             }
         }
 
           nxt[isz+nbt][isx+nbl] += wav[it];
           
           source_smooth(nxt,isz+nbt,isx+nbl,wav[it]);
           bd_decay(nxt); 
           bd_decay(dercur); 
                 
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
    fft2_finalize();
    bd_close();
    srcsm_close();
    free(**aa);
    free(*aa);
    free(aa);
    free(*vx);     
    free(*vz);     
    free(*yi);     
    free(*se);     
    free(*nxt);     
    free(*cur);     
    free(*old);     
    free(*dercur);     
    free(*derold);     
    free(vx);     
    free(vz);     
    free(yi);     
    free(se);     
    free(nxt);     
    free(cur);     
    free(old);     
    free(dercur);     
    free(derold);     
    free(tmpcur);     

    exit(0); 
}           
           
