/* Get the derivative of time-shift gathers */
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
#ifdef _OPENMP
#include <omp.h>
#endif

static int padnx, padnz;
static float c0, c11, c21, c12, c22;
static float **padvv;

static void laplacian(float **uin, float **uout)
{
    int ix, iz;
    
#ifdef _OPENMP
#pragma omp parallel for default(none) private(ix, iz) \
shared(padnx, padnz, uin, uout, padvv, c0, c11, c12, c21, c22)
#endif
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
#ifdef _OPENMP
#pragma omp parallel for private(i2, i1)
#endif
    for(i2=0; i2<n2; i2++){
        for(i1=0; i1<n1; i1++){
            data[i2][i1]=0.0;
        }
    }
}

int main(int argc, char *argv[])
{
    int ix, iz, itau;
    int nx, nz, ntau, pad, htau;
    float dt, dtau, dx, dz, idz2, idx2;
    
    float ***dd, ***der, **vv, **v0;
    float **cur, **nxt, **ud;
    
    sf_axis ax, az, atau;
    sf_file Ftg, Fvel, Fder;
    
    sf_init(argc, argv);
    
    Ftg=sf_input("in");
    Fder=sf_output("out");
    Fvel=sf_input("velocity");
    
    az=sf_iaxa(Ftg, 1);
    ax=sf_iaxa(Ftg, 2);
    atau=sf_iaxa(Ftg, 3);
    
    nz=sf_n(az); dz=sf_d(az);
    nx=sf_n(ax); dx=sf_d(ax);
    ntau=sf_n(atau); dtau=sf_d(atau);
    
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
    
	/* allocate storage */
    dd=sf_floatalloc3(nz, nx, ntau);
	der=sf_floatalloc3(nz, nx, ntau);
    vv=sf_floatalloc2(nz, nx);
    
    v0=sf_floatalloc2(padnz, padnx);
    padvv=sf_floatalloc2(padnz, padnx);

    cur=sf_floatalloc2(padnz, padnx);
    nxt=sf_floatalloc2(padnz, padnx);
    ud=sf_floatalloc2(padnz, padnx);

	/* read data and velocity */
    sf_floatread(dd[0][0], ntau*nx*nz, Ftg);
    sf_floatread(vv[0], nx*nz, Fvel);
    
	/* pad model to avoid boundary reflection */
    for(ix=0; ix<nx; ix++)
        for(iz=0; iz<nz; iz++)
            padvv[ix+pad][iz+pad]=vv[ix][iz]*vv[ix][iz]*dt*dt;
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
	
	
	/* calculate derivative and save it into ***der1 */
    htau=(ntau-1)/2+0.5;
	
	  // when tau <0
	for(itau=0; itau<htau; itau++){	
		zero2(nxt, padnz, padnx);
		zero2(cur, padnz, padnx);
		zero2(ud,  padnz, padnx);
		zero2(v0,  padnz, padnx);
		
#ifdef _OPENMP
#pragma omp parallel for private(ix, iz)
#endif
		for(ix=0; ix<nx; ix++){
			for(iz=0; iz<nz; iz++){
				if(itau == 0){
					v0[ix+pad][iz+pad]=(dd[1][ix][iz]-dd[0][ix][iz])/dtau;
				}else{
					v0[ix+pad][iz+pad]=(dd[itau+1][ix][iz]-dd[itau-1][ix][iz])/dtau/2.0;
				}
			}
		}

		// calculate u(tau)
#ifdef _OPENMP
#pragma omp parallel for private(ix, iz)
#endif
		for(ix=0; ix<nx; ix++)
			for(iz=0; iz<nz; iz++)
				nxt[ix+pad][iz+pad]=dd[itau][ix][iz];

		// calculate u(tau-det(t)) and derivative
		laplacian(nxt, ud);

#ifdef _OPENMP
#pragma omp parallel for private(ix, iz)
#endif
		for(ix=0; ix<padnx; ix++){
			for(iz=0; iz<padnz; iz++){
				cur[ix][iz]=nxt[ix][iz]+ud[ix][iz]/2.0 - v0[ix][iz]*dt;
			}
		}

#ifdef _OPENMP
#pragma omp parallel for private(ix, iz)
#endif
		for(ix=0; ix<nx; ix++){
			for(iz=0; iz<nz; iz++){
				der[itau][ix][iz]=(nxt[ix+pad][iz+pad]-cur[ix+pad][iz+pad])/dt;
			}
		}
	} // end of itau

	   // when tau=0
#ifdef _OPENMP
#pragma omp parallel for private(ix, iz)
#endif
	for(ix=0; ix<nx; ix++){
		for(iz=0; iz<nz; iz++){
			der[htau][ix][iz]=0.;
		}
	}

	   // when tau >0
	for(itau=htau+1; itau<ntau; itau++){
		zero2(nxt, padnz, padnx);
		zero2(cur, padnz, padnx);
		zero2(ud,  padnz, padnx);
		zero2(v0,  padnz, padnx);

#ifdef _OPENMP
#pragma omp parallel for private(ix, iz)
#endif
		for(ix=0; ix<nx; ix++){
			for(iz=0; iz<nz; iz++){
				if(itau == ntau-1){
					v0[ix+pad][iz+pad]=(dd[ntau-1][ix][iz]-dd[ntau-2][ix][iz])/dtau;
				}else{
					v0[ix+pad][iz+pad]=(dd[itau+1][ix][iz]-dd[itau-1][ix][iz])/dtau/2.0;
				}
			}
		}

		// calculate u(tau)
#ifdef _OPENMP
#pragma omp parallel for private(ix, iz)
#endif
		for(ix=0; ix<nx; ix++)
			for(iz=0; iz<nz; iz++)
				nxt[ix+pad][iz+pad]=dd[itau][ix][iz];

		// calculate u(tau+det(t)) and derivative
		laplacian(nxt, ud);

#ifdef _OPENMP
#pragma omp parallel for private(ix, iz)
#endif
		for(ix=0; ix<padnx; ix++){
			for(iz=0; iz<padnz; iz++){
				cur[ix][iz]=nxt[ix][iz]+ud[ix][iz]/2.0 + v0[ix][iz]*dt;
			}
		}

#ifdef _OPENMP
#pragma omp parallel for private(ix, iz)
#endif
		for(ix=0; ix<nx; ix++){
			for(iz=0; iz<nz; iz++){
				der[itau][ix][iz]=(cur[ix+pad][iz+pad]-nxt[ix+pad][iz+pad])/dt;
			}
		}
	} // end of itau
	sf_floatwrite(der[0][0], ntau*nx*nz, Fder);
   
	free(**dd); free(*dd); free(dd);
	free(**der); free(*der); free(der);
	free(*vv); free(vv);
	free(*padvv); free(padvv);
	free(*v0); free(v0);
	free(*cur); free(cur);
	free(*nxt); free(nxt);
	free(*ud); free(ud);

    exit(0);
}
