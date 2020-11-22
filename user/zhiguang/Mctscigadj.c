/* Correcting time-shift gathers and its adjoint */
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
	bool adj;
    int ix, iz, it, itau;
    int nx, nz, ntau, nt, pad, htau;
    float dt, dtau, dx, dz, idz2, idx2;
    float tau0, tau, sign = 0.;
    
    float ***dd, ***mm, ***dertau, ***der0, **vv;
    float **cur, **nxt, **dercur, **dernxt, **tmp;
    
    sf_axis ax, az, atau;
    sf_file Ftg, Fcg, Fvel, Fdertau, Fder0;
    
    sf_init(argc, argv);
	
	if(!sf_getbool("adj", &adj)) adj=true;
    if(!sf_getfloat("dt", &dt)) dt=0.001;
    if(!sf_getint("pad", &pad)) pad=100;
    
	/* files */
	if(adj){ /* migration */
		Ftg=sf_input("in");
		Fdertau=sf_input("Fdertau");
		Fcg=sf_output("out");
		Fder0=sf_output("Fder0");
	}else{ /* modeling */
		Fcg=sf_input("in");
		Fder0=sf_input("Fder0");
		Ftg=sf_output("out");
		Fdertau=sf_output("Fdertau");
	}
    Fvel=sf_input("velocity");
    
	if(adj){
		az=sf_iaxa(Ftg, 1);
		ax=sf_iaxa(Ftg, 2);
		atau=sf_iaxa(Ftg, 3);
	}else{
		az=sf_iaxa(Fcg, 1);
		ax=sf_iaxa(Fcg, 2);
		atau=sf_iaxa(Fcg, 3);
	}
    
    nz=sf_n(az); dz=sf_d(az);
    nx=sf_n(ax); dx=sf_d(ax);
    ntau=sf_n(atau); dtau=sf_d(atau); tau0=sf_o(atau);
    
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
    mm=sf_floatalloc3(nz, nx, ntau);
	dertau=sf_floatalloc3(nz, nx, ntau);
	der0=sf_floatalloc3(nz, nx, ntau);
    vv=sf_floatalloc2(nz, nx);
    
    padvv=sf_floatalloc2(padnz, padnx);
    cur=sf_floatalloc2(padnz, padnx);
    nxt=sf_floatalloc2(padnz, padnx);
    dercur=sf_floatalloc2(padnz, padnx);
    dernxt=sf_floatalloc2(padnz, padnx);

	/* read data, derivative and velocity */
	if(adj){
		sf_floatread(dd[0][0], ntau*nx*nz, Ftg);
		sf_floatread(dertau[0][0], ntau*nx*nz, Fdertau);
	}else{
		sf_floatread(mm[0][0], ntau*nx*nz, Fcg);
		sf_floatread(der0[0][0], ntau*nx*nz, Fder0);
	}
    sf_floatread(vv[0], nx*nz, Fvel);
    
	/* pad model to avoid boundary reflection */
    for(ix=0; ix<nx; ix++)
        for(iz=0; iz<nz; iz++)
            padvv[ix+pad][iz+pad]=vv[ix][iz]*vv[ix][iz]*dt;
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
	

    htau=(ntau-1)/2;
	if(adj){ /* migration */
		
		for(itau=0; itau<ntau; itau++){
		
			 // when tau=0
			if(itau==htau){
				sf_warning("Migration: itau=%d/%d, nt=%d", htau+1, ntau, 0);
#ifdef _OPENMP
#pragma omp parallel for private(ix, iz)
#endif
				for(ix=0; ix<nx; ix++){
					for(iz=0; iz<nz; iz++){
						mm[htau][ix][iz]=dd[htau][ix][iz];
						der0[htau][ix][iz]=dertau[htau][ix][iz];
					}
				}
				continue;
			}

			if(itau<htau) sign=-1;
			if(itau>htau) sign=1;
		
			zero2(cur, padnz, padnx);
			zero2(nxt, padnz, padnx);
			zero2(dercur, padnz, padnx);
			zero2(dernxt, padnz, padnx);
			
			tau=tau0+itau*dtau;
			nt=sign*tau/dt+0.5;
			sf_warning("Migration: itau=%d/%d, nt=%d", itau+1, ntau, nt);
			
			// current wavefield and derivative
#ifdef _OPENMP
#pragma omp parallel for private(ix, iz)
#endif
			for(ix=0; ix<nx; ix++){
				for(iz=0; iz<nz; iz++){
					cur[ix+pad][iz+pad]=dd[itau][ix][iz];
					dercur[ix+pad][iz+pad]=dertau[itau][ix][iz];
				}
			}
			
			for(it=0; it<nt; it++){
				laplacian(cur, nxt);
#ifdef _OPENMP
#pragma omp parallel for private(ix, iz)
#endif
				for(ix=0; ix<padnx; ix++){
					for(iz=0; iz<padnz; iz++){
						dernxt[ix][iz]=dercur[ix][iz]-sign*nxt[ix][iz];
						nxt[ix][iz]=cur[ix][iz]-sign*dt*dernxt[ix][iz];
					}
				}
				
				tmp=cur; cur=nxt; nxt=tmp;
				tmp=dercur; dercur=dernxt; dernxt=tmp;
			} //end of it

#ifdef _OPENMP
#pragma omp parallel for private(ix, iz)
#endif
			for(ix=0; ix<nx; ix++){
				for(iz=0; iz<nz; iz++){
					mm[itau][ix][iz]=cur[ix+pad][iz+pad];
					der0[itau][ix][iz]=dercur[ix+pad][iz+pad];
				}
			}
		} // end of itau
	}else{ /* modeling */

		for (itau=0; itau<ntau; itau++){

			 // when tau=0
			if(itau==htau){
				sf_warning("Modeling: itau=%d/%d, nt=%d", htau+1, ntau, 0);
#ifdef _OPENMP
#pragma omp parallel for private(ix, iz)
#endif
				for(ix=0; ix<nx; ix++){
					for(iz=0; iz<nz; iz++){
						dd[htau][ix][iz]=mm[htau][ix][iz];
						dertau[htau][ix][iz]=der0[htau][ix][iz];
					}
				}
				continue;
			}

			if(itau<htau) sign=-1;
			if(itau>htau) sign=1;
		
			zero2(cur, padnz, padnx);
			zero2(nxt, padnz, padnx);
			zero2(dercur, padnz, padnx);
			zero2(dernxt, padnz, padnx);
			
			tau=tau0+itau*dtau;
			nt=sign*tau/dt+0.5;
			sf_warning("Modeling: itau=%d/%d, nt=%d", itau+1, ntau, nt);
			
			// current wavefield and derivative
#ifdef _OPENMP
#pragma omp parallel for private(ix, iz)
#endif
			for(ix=0; ix<nx; ix++){
				for(iz=0; iz<nz; iz++){
					cur[ix+pad][iz+pad]=mm[itau][ix][iz];
					dercur[ix+pad][iz+pad]=der0[itau][ix][iz];
				}
			}
			
			for(it=0; it<nt; it++){
#ifdef _OPENMP
#pragma omp parallel for private(ix, iz)
#endif
				for(ix=0; ix<padnx; ix++){
					for(iz=0; iz<padnz; iz++){
						nxt[ix][iz]=cur[ix][iz]+sign*dt*dercur[ix][iz];
					}
				}

				laplacian(nxt, cur);

#ifdef _OPENMP
#pragma omp parallel for private(ix, iz)
#endif
				for(ix=0; ix<padnx; ix++){
					for(iz=0; iz<padnz; iz++){
						dernxt[ix][iz]=dercur[ix][iz]+sign*cur[ix][iz];
					}
				}
				
				tmp=cur; cur=nxt; nxt=tmp;
				tmp=dercur; dercur=dernxt; dernxt=tmp;
			} // end of it

#ifdef _OPENMP
#pragma omp parallel for private(ix, iz)
#endif
			for(ix=0; ix<nx; ix++){
				for(iz=0; iz<nz; iz++){
					dd[itau][ix][iz]=cur[ix+pad][iz+pad];
					dertau[itau][ix][iz]=dercur[ix+pad][iz+pad];
				}
			}
		} // end of itau	
	} // end of modeling
    
	if(adj){
		sf_floatwrite(mm[0][0], ntau*nx*nz, Fcg);
		sf_floatwrite(der0[0][0], ntau*nx*nz, Fder0);
	}else{
		sf_floatwrite(dd[0][0], ntau*nx*nz, Ftg);
		sf_floatwrite(dertau[0][0], ntau*nx*nz, Fdertau);
	}
    
    exit(0);
}
