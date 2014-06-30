/* Correct time-shift gathers */
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
#include <omp.h>

static int padnx, padnz;
static float c0, c11, c21, c12, c22;
static float **padvv;

static void laplacian(float **uin, float **uout)
{
    int ix, iz;
    
#pragma omp parallel for default(none) private(ix, iz) \
shared(padnx, padnz, uin, uout, padvv, c0, c11, c12, c21, c22)
    for(ix=2; ix<padnx-2; ix++){
        for(iz=2; iz<padnz-2; iz++){
            uout[ix][iz]=
            (c11*(uin[ix][iz-1]+uin[ix][iz+1])+
            c12*(uin[ix][iz-2]+uin[ix][iz+2])+
            c21*(uin[ix-1][iz]+uin[ix+1][iz])+
            c22*(uin[ix-2][iz]+uin[ix+2][iz])+
            c0*uin[ix][iz])*padvv[ix][iz];
        }
    }
}

void zero2( float **data, int n1, int n2)
{
    int i1, i2;
#pragma omp parallel for private(i2, i1)
    for(i2=0; i2<n2; i2++){
        for(i1=0; i1<n1; i1++){
            data[i2][i1]=0.0;
        }
    }
}

int main(int argc, char *argv[])
{
    int ix, iz, it, itau;
    int nx, nz, ntau, nt, pad;
    float dt, dtau, dx, dz, dt2, idz2, idx2;
    float tau0, tau;
    
    float ***dd, ***mm, **vv, **v0;
    float **u0, **u1, **u2, **ud, **tmp;
    
    sf_axis ax, az, atau;
    sf_file tgather, cgather, vel;
    
    sf_init(argc, argv);
    
    tgather=sf_input("in");
    cgather=sf_output("out");
    vel=sf_input("velocity");
    
    az=sf_iaxa(tgather, 1);
    ax=sf_iaxa(tgather, 2);
    atau=sf_iaxa(tgather, 3);
    
    nz=sf_n(az); dz=sf_d(az);
    nx=sf_n(ax); dx=sf_d(ax);
    ntau=sf_n(atau); dtau=sf_d(atau); tau0=sf_o(atau);
    
    if(!sf_getfloat("dt", &dt)) dt=0.001;
    if(!sf_getint("pad", &pad)) pad=30;
    padnx=nx+2*pad;
    padnz=nz+2*pad;
    
    idz2=1.0/(dz*dz);
    idx2=1.0/(dx*dx);
    
    c11=4.0*idz2/3.0;
    c12=-idz2/12.0;
    c21=4.0*idx2/3.0;
    c22=-idx2/12.0;
    c0=-2.0*(c11+c12+c21+c22);
    
    dd=sf_floatalloc3(nz, nx, ntau);
    mm=sf_floatalloc3(nz, nx, ntau);
    vv=sf_floatalloc2(nz, nx);
    
    v0=sf_floatalloc2(padnz, padnx);
    padvv=sf_floatalloc2(padnz, padnx);
    
    u0=sf_floatalloc2(padnz, padnx);
    u1=sf_floatalloc2(padnz, padnx);
    u2=sf_floatalloc2(padnz, padnx);
    ud=sf_floatalloc2(padnz, padnx);
    
    sf_floatread(dd[0][0], ntau*nx*nz, tgather);
    sf_floatread(vv[0], nx*nz, vel);
    
    dt2=dt*dt;
    for(ix=0; ix<nx; ix++)
        for(iz=0; iz<nz; iz++)
            padvv[ix+pad][iz+pad]=vv[ix][iz]*vv[ix][iz]*dt2;
    for(iz=0; iz<pad; iz++)
        for(ix=pad; ix<nx+pad; ix++){
            padvv[ix][iz]=padvv[ix][pad];
            padvv[ix][pad+nz+iz]=padvv[ix][pad+nz-1];
        }
    
    for(ix=0; ix<pad; ix++)
        for(iz=0; iz<padnz; iz++){
            padvv[ix][iz]=padvv[pad][iz];
            padvv[ix+pad+nx][iz]=padvv[pad+nx-1][iz];
        }
    
    for(itau=0; itau<ntau; itau++){
		sf_warning("itau=%d/%d", itau+1, ntau);
        
        zero2(u0, padnz, padnx);
        zero2(u1, padnz, padnx);
        zero2(u2, padnz, padnx);
        zero2(ud, padnz, padnx);
        zero2(v0, padnz, padnx);
        
        tau=tau0+itau*dtau;

		// tau=0
        if(tau==0.){
            for(ix=0; ix<nx; ix++)
                for(iz=0; iz<nz; iz++)
                    mm[itau][ix][iz]=dd[itau][ix][iz];
            continue;
        }
        
        // calculate v0 (derivative with respect to tau)
        if(itau==0){
            for(ix=0; ix<nx; ix++){
                for(iz=0; iz<nz; iz++){
                    v0[ix+pad][iz+pad]=(dd[1][ix][iz]-dd[0][ix][iz])/dtau;
                }
            }
        } else if (itau==ntau-1){
            for(ix=0; ix<nx; ix++){
                for(iz=0; iz<nz; iz++){
                    v0[ix+pad][iz+pad]=(dd[ntau-1][ix][iz]-dd[ntau-2][ix][iz])/dtau;
                }
            }
        } else {
#pragma omp parallel for private(ix, iz)
            for(ix=0; ix<nx; ix++){
                for(iz=0; iz<nz; iz++){
                    v0[ix+pad][iz+pad]=(dd[itau+1][ix][iz]-dd[itau-1][ix][iz])/dtau/2.0;
                }
            }
        }
        
        // calculate u1
#pragma omp parallel for private(ix, iz)
        for(ix=0; ix<nx; ix++)
            for(iz=0; iz<nz; iz++)
                u1[ix+pad][iz+pad]=dd[itau][ix][iz];
        
		// tau>0
        if(tau>0.){
            
            laplacian(u1, ud);
#pragma omp parallel for private(ix, iz)
            for(ix=0; ix<padnx; ix++){
                for(iz=0; iz<padnz; iz++){
                    u0[ix][iz]=u1[ix][iz]+ud[ix][iz]/2.0-v0[ix][iz]*dt;
                }
            }
            
            nt=tau/dt;
            
            for(it=1; it<nt; it++){
				sf_warning("it=%d/%d;", it+1, nt);
                tmp=u2; u2=u1; u1=u0; u0=tmp;
                
                laplacian(u1, ud);
#pragma omp parallel for private(ix, iz)
                for(ix=0; ix<padnx; ix++){
                    for(iz=0; iz<padnz; iz++){
                        u0[ix][iz]=2*u1[ix][iz]-u2[ix][iz]+ud[ix][iz];
                    }
                }
            } //end of it
#pragma omp parallel for private(ix, iz)
            for(ix=0; ix<nx; ix++){
                for(iz=0; iz<nz; iz++){
                    mm[itau][ix][iz]=u0[ix+pad][iz+pad];
                }
            }
        }
        
		// tau<0
        if(tau<0.){
            
            laplacian(u1, ud);
#pragma omp parallel for private(ix, iz)
            for(ix=0; ix<padnx; ix++){
                for(iz=0; iz<padnz; iz++){
                    u2[ix][iz]=u1[ix][iz]+dt*v0[ix][iz]+ud[ix][iz]/2.0;
                }
            }
            
            nt=-tau/dt;
            
            for(it=1; it<nt; it++){
				sf_warning("it=%d/%d;", it+1, nt);
                tmp=u0; u0=u1; u1=u2; u2=tmp;
                
                laplacian(u1, ud);
#pragma omp parallel for private(ix, iz)
                for(ix=0; ix<padnx; ix++){
                    for(iz=0; iz<padnz; iz++){
                        u2[ix][iz]=2*u1[ix][iz]-u0[ix][iz]+ud[ix][iz];
                    }
                }
            }//end of it
#pragma omp parallel for private(ix, iz)
            for(ix=0; ix<nx; ix++){
                for(iz=0; iz<nz; iz++){
                    mm[itau][ix][iz]=u2[ix+pad][iz+pad];
                }
            }
        }
    } //end of itau
    
    sf_floatwrite(mm[0][0], ntau*nx*nz, cgather);
    
    exit(0);
}
