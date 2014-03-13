/* 2-D FD RTM and its adjoint */
/*
 Copyright (C) 2014 University of Texas at Austin
 
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
/*^*/

#include "prertm.h"

#ifdef _OPENMP
#include <omp.h>
#endif

static bool verb;
static int nz, nx, nt, nr, ns, nw, nsource, dsource, ndelay;
static int padx, padz, padnx, padnz, dr_v, ds_v, r0_v, s0_v, zr_v, zs_v;
static float **padvv, *ww;
static float c0, c11, c12, c21, c22;

void prertm_init(bool verb1, int nz1, int nx1, int nt1, int nr1, int ns1, int nw1, int nsource1,
                 int dsource1, int ndelay1, float dx, float dz, int padx1, int padz1, int padnx1,
                 int padnz1, int dr_v1, int ds_v1, int r0_v1, int s0_v1, int zr_v1, int zs_v1,
                 float **padvv1, float *ww1)
/*< initialize >*/
{

#ifdef _OPENMP
	omp_init();
#endif
    
    float idx2, idz2;
    
    idx2=1.0/(dx*dx);
    idz2=1.0/(dz*dz);
    
    c11=4.0*idz2/3.0;
    c12=-idz2/12.0;
    c21=4.0*idx2/3.0;
    c22=-idx2/12.0;
    c0=-2.0*(c11+c12+c21+c22);
    
    verb=verb1;
    nz=nz1;
    nx=nx1;
    nt=nt1;
    nr=nr1;
    ns=ns1;
    
    nw=nw1;
    nsource=nsource1;
    dsource=dsource1;
    ndelay=ndelay1;
    
    padx=padx1;
    padz=padz1;
    padnx=padnx1;
    padnz=padnz1;
    
    dr_v=dr_v1;
    r0_v=r0_v1;
    zr_v=zr_v1;
    ds_v=ds_v1;
    s0_v=s0_v1;
    zs_v=zs_v1;
    
    padvv=padvv1;
    ww=ww1;    
}

void prertm_close()
/*< free allocated storage >*/
{
    free(ww);
    free(*padvv);
    free(padvv);
}

void laplacian(bool adj, float **u0, float **u1, float **u2)
{
    int ix, iz;
    
    if(adj){
#ifdef _OPENMP
#pragma omp parallel for  default(none) \
private(ix, iz)       \
shared(padnx, padnz, u0, u1, u2, padvv, c0, c11, c12, c21, c22)
#endif
        for(ix=2; ix<padnx-2; ix++){
            for(iz=2; iz<padnz-2; iz++){
                u2[ix][iz]=
                (c11*(u1[ix][iz-1]+u1[ix][iz+1])+
                 c12*(u1[ix][iz-2]+u1[ix][iz+2])+
                 c0*u1[ix][iz]+
                 c21*(u1[ix-1][iz]+u1[ix+1][iz])+
                 c22*(u1[ix-2][iz]+u1[ix+2][iz]))
                *padvv[ix][iz]+2*u1[ix][iz]-u0[ix][iz];
            }
        }
    }else{
#ifdef _OPENMP
#pragma omp parallel for  default(none) \
private(ix, iz)       \
shared(padnx, padnz, u0, u1, u2, padvv, c0, c11, c12, c21, c22)
#endif
        for(ix=2; ix<padnx-2; ix++){
            for(iz=2; iz<padnz-2; iz++){
                u2[ix][iz]=
                (c11*(u1[ix][iz-1]*padvv[ix][iz-1]+u1[ix][iz+1]*padvv[ix][iz+1])+
                 c12*(u1[ix][iz-2]*padvv[ix][iz-2]+u1[ix][iz+2]*padvv[ix][iz+2])+
                 c0*u1[ix][iz]*padvv[ix][iz]+
                 c21*(u1[ix-1][iz]*padvv[ix-1][iz]+u1[ix+1][iz]*padvv[ix+1][iz])+
                 c22*(u1[ix-2][iz]*padvv[ix-2][iz]+u1[ix+2][iz]*padvv[ix+2][iz]))
                +2*u1[ix][iz]-u0[ix][iz];
            }
        }
    }
}


void prertm2_oper(bool adj, bool add, int nm , int nd, float *mm, float *dd)
/*< linear operator >*/
{
    int ix, iz, is, it, ir, sx, rx, nnt;
    float **u0, **u1, **u2, **temp, **sou2, **perm, ***wave;
    
    sf_adjnull(adj, add, nm, nd, mm, dd);
    
    nnt=1+(nt-1)/nw;
    
    u0=sf_floatalloc2(padnz, padnx);
    u1=sf_floatalloc2(padnz, padnx);
    u2=sf_floatalloc2(padnz, padnx);
    sou2=sf_floatalloc2(nz, nx);
    perm=sf_floatalloc2(nz, nx);
    wave=sf_floatalloc3(nz, nx, nnt);
    
    if(adj){/* migration */
        
        for(is=0; is<ns; is++){
            if(verb) sf_warning("ishot @@@ nshot:  %d  %d", is+1, ns);
            sx=is*ds_v+s0_v;
            
            memset(u0[0], 0, padnz*padnx*sizeof(float));
            memset(u1[0], 0, padnz*padnx*sizeof(float));
            memset(u2[0], 0, padnz*padnx*sizeof(float));
            memset(sou2[0], 0, nz*nx*sizeof(float));
            memset(perm[0], 0, nz*nx*sizeof(float));
            memset(wave[0][0], 0, nz*nx*nnt*sizeof(float));
            
            for(it=0; it<nt; it++){
                //sf_warning("RTM_Source: it=%d;",it);
                laplacian(true, u0, u1, u2);
                
                temp=u0;
                u0=u1;
                u1=u2;
                u2=temp;
                
                for(ix=0; ix<nsource; ix++){
                if(it>=ix*ndelay)
                u1[sx+ix*dsource][zs_v]+=ww[it-ix*ndelay];
                }
                
                if(it%nw==0){
                for(ix=0; ix<nx; ix++){
                    for(iz=0; iz<nz; iz++){
                        sou2[ix][iz]+=u1[ix+padx][iz+padz]*u1[ix+padx][iz+padz];
                        wave[it/nw][ix][iz]=u1[ix+padx][iz+padz];
                    }
                }
                }
            }// end of it
            
            memset(u0[0], 0, padnz*padnx*sizeof(float));
            memset(u1[0], 0, padnz*padnx*sizeof(float));
            memset(u2[0], 0, padnz*padnx*sizeof(float));
            
            for(it=nt-1; it>=0; it--){
                //sf_warning("RTM_Receiver: it=%d;",it);
                laplacian(false, u0, u1, u2);
                
                temp=u0;
                u0=u1;
                u1=u2;
                u2=temp;
                
                for(ir=0; ir<nr; ir++){
                    rx=ir*dr_v+r0_v;
                    u1[rx][zr_v]+=dd[is*nr*nt+ir*nt+it];
                }
                if(it%nw==0){
                for(ix=0; ix<nx; ix++){
                    for(iz=0; iz<nz; iz++){
                        perm[ix][iz]+=wave[it/nw][ix][iz]*u1[ix+padx][iz+padz];
                    }
                }
                }
            }// end of it
            for(ix=0; ix<nx; ix++)
                for(iz=0; iz<nz; iz++)
                    mm[ix*nz+iz]+=perm[ix][iz]/(sou2[ix][iz]+FLT_EPSILON);
        }// end of shot
    }else{/* modeling */
        
        for(is=0; is<ns; is++){
            if(verb) sf_warning("ishot @@@ nshot:  %d  %d", is+1, ns);
            sx=is*ds_v+s0_v;
            
            memset(u0[0], 0, padnz*padnx*sizeof(float));
            memset(u1[0], 0, padnz*padnx*sizeof(float));
            memset(u2[0], 0, padnz*padnx*sizeof(float));
            memset(sou2[0], 0, nz*nx*sizeof(float));
            memset(perm[0], 0, nz*nx*sizeof(float));
            memset(wave[0][0], 0, nz*nx*nnt*sizeof(float));
            
            for(it=0; it<nt; it++){
                //sf_warning("RTM_Source: it=%d;",it);
                laplacian(true, u0, u1, u2);
                
                temp=u0;
                u0=u1;
                u1=u2;
                u2=temp;
                
                for(ix=0; ix<nsource; ix++){
                if(it>=ix*ndelay)
                u1[sx+ix*dsource][zs_v]+=ww[it-ix*ndelay];
                }
                
                if(it%nw==0){
                for(ix=0; ix<nx; ix++){
                    for(iz=0; iz<nz; iz++){
                        sou2[ix][iz]+=u1[ix+padx][iz+padz]*u1[ix+padx][iz+padz];
                        wave[it/nw][ix][iz]=u1[ix+padx][iz+padz];
                    }
                }
                }
            }// end of it
            
            memset(u0[0], 0, padnz*padnx*sizeof(float));
            memset(u1[0], 0, padnz*padnx*sizeof(float));
            memset(u2[0], 0, padnz*padnx*sizeof(float));
            
            for(ix=0; ix<nx; ix++)
                for(iz=0; iz<nz; iz++)
                    perm[ix][iz]+=mm[ix*nz+iz]/(sou2[ix][iz]+FLT_EPSILON);
            
            for(it=0; it<nt; it++){
                //sf_warning("Modeling_Receiver: it=%d;",it);
                laplacian(true, u0, u1, u2);
                
                temp=u0;
                u0=u1;
                u1=u2;
                u2=temp;
                
                if(it%nw==0){
                for(ix=0; ix<nx; ix++){
                    for(iz=0; iz<nz; iz++){
                        u1[ix+padx][iz+padz]+=wave[it/nw][ix][iz]*perm[ix][iz];
                    }
                }
                }
                
                for(ir=0; ir<nr; ir++){
                    rx=ir*dr_v+r0_v;
                    dd[is*nr*nt+ir*nt+it]+=u1[rx][zr_v];
                }
            } //end of it
        }// end of shot
    }// end of if
}
