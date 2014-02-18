/* 2-D exploding-reflector RTM of incomplete data and its adjoint */
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

#include "postrtm.h"

#ifdef _OPENMP
#include <omp.h>
#endif

static int nx, nz, nt, n0, padx, padz, padnx, padnz, nxz;
static float c0, c11, c12, c21, c22;
static float **padvv;
static int *head;

void postrtm_init(int nx1, int nz1, int nt1, int n01, int padx1, int padz1, int padnx1,
                     int padnz1, float dx, float dz, int *head1, float **padvv1)
/*< initialization >*/
{
#ifdef _OPENMP
	omp_init();
#endif
    float idx2, idz2;
    
    nx=nx1;
    nz=nz1;
    nt=nt1;
    n0=n01;
    
    padx=padx1;
    padz=padz1;
    padnx=padnx1;
    padnz=padnz1;
    nxz=padnx1*padnz1;
    
    head=head1;
    padvv=padvv1;
    
    idx2=1.0/(dx*dx);
    idz2=1.0/(dz*dz);
    
    c11=4.0*idz2/3.0;
    c12=-idz2/12.0;
    c21=4.0*idx2/3.0;
    c22=-idx2/12.0;
    c0=-2.0*(c11+c12+c21+c22);
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

void postrtm_lop(bool adj, bool add, int nm, int nd, float *mm, float *dd)
/*< linear operator >*/
{
    int it, ix, iz;
    float **u0, **u1, **u2, **tmp;
    
    u0=sf_floatalloc2(padnz, padnx);
    u1=sf_floatalloc2(padnz, padnx);
    u2=sf_floatalloc2(padnz, padnx);
    
    memset(u0[0], 0.0, nxz*sizeof(float));
    memset(u1[0], 0.0, nxz*sizeof(float));
    memset(u2[0], 0.0, nxz*sizeof(float));
    
    sf_adjnull(adj, add, nm, nd, mm, dd);
    
    if(adj){ //migration
        
        /* mask operator */
        for(ix=0; ix<nx; ix++){
            if(head[ix]==0)
                for(it=0; it<nt; it++)
                    dd[ix*nt+it]=0.0;
        }
        
        for(it=nt-1; it>=0; it--){
            sf_warning("Migration: %d;", it);
            
            laplacian(adj, u0, u1, u2);
            tmp=u0; u0=u1; u1=u2; u2=tmp;

#ifdef _OPENMP
#pragma omp parallel for default(none) private(ix) shared(padx, nx, dd, nt, it, u1, n0)
#endif
            for(ix=padx; ix<padx+nx; ix++)
                /* inject data */
                u1[ix][n0]+=dd[it+(ix-padx)*nt];
        }
        sf_warning(".");
        
        /* output image */
        for(ix=0; ix<nx; ix++)
            for(iz=0; iz<nz; iz++)
                mm[ix*nz+iz]+=u1[padx+ix][padz+iz];
    
    }else{ //modeling
        
        /* read reflector */
        for(ix=0; ix<nx; ix++)
            for(iz=0; iz<nz; iz++)
                u1[padx+ix][padz+iz]+=mm[ix*nz+iz];
        
        for(it=0; it<nt; it++){
            sf_warning("Modeling: %d;", it);
            
#ifdef _OPENMP
#pragma omp parallel for default(none) private(ix) shared(padx, nx, dd, nt, it, u1, n0)
#endif
            for(ix=padx; ix<padx+nx; ix++)
                /* record data */
                dd[it+(ix-padx)*nt]+=u1[ix][n0];
            
            laplacian(adj, u0, u1, u2);
            tmp=u0; u0=u1; u1=u2;; u2=tmp;
        }
        sf_warning(".");
        
        /* mask operator */
        for(ix=0; ix<nx; ix++){
            if(!head[ix])
                for(it=0; it<nt; it++)
                    dd[ix*nt+it]=0.0;
        }
        
    }
}
