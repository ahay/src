/* 2D 8-th order elastic wave propagation using sponge ABC
*/
/*
  Copyright (C) 2014  Xi'an Jiaotong University (Pengliang Yang)

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

static int nb, nz, nx, nt, nzpad, nxpad;
static float dz, dx, _dz, _dx, dt, fm;

void expand2d(float** b, float** a)
/*< expand domain of 'a' to 'b': source(a)-->destination(b) >*/
{
    int iz,ix;

#ifdef _OPENMP
#pragma omp parallel for default(none)	\
	private(ix,iz)			\
	shared(b,a,nb,nz,nx)
#endif
    for     (ix=0;ix<nx;ix++) {
	for (iz=0;iz<nz;iz++) {
	    b[nb+ix][nb+iz] = a[ix][iz];
	}
    }

    for     (ix=0; ix<nxpad; ix++) {
	for (iz=0; iz<nb;    iz++) {
	    b[ix][      iz  ] = b[ix][nb  ];
	    b[ix][nzpad-iz-1] = b[ix][nzpad-nb-1];
	}
    }

    for     (ix=0; ix<nb;    ix++) {
	for (iz=0; iz<nzpad; iz++) {
	    b[ix 	 ][iz] = b[nb  		][iz];
	    b[nxpad-ix-1 ][iz] = b[nxpad-nb-1	][iz];
	}
    }
}


void window2d(float **a, float **b)
/*< window 'b' to 'a': source(b)-->destination(a) >*/
{
    int iz,ix;

#ifdef _OPENMP
#pragma omp parallel for default(none)	\
	private(ix,iz)			\
	shared(b,a,nb,nz,nx)
#endif
    for     (ix=0;ix<nx;ix++) {
	for (iz=0;iz<nz;iz++) {
	    a[ix][iz]=b[nb+ix][nb+iz] ;
	}
    }
}


void forward_uvx_uvz(float **uvx, float **uvz, float **txx, float **tzz, float **txz, float **rho)
/*< forward step: update uvx, uvz >*/
{
	int i1, i2;
	float diff1, diff2, diff3, diff4;

#ifdef _OPENMP
#pragma omp parallel for default(none) 	\
	private(i1,i2,diff1,diff2,diff3,diff4)		\
	shared(uvx,uvz,txx,tzz,txz,rho,nxpad,nzpad,dt,_dx,_dz)
#endif
	for(i2=4; i2<nxpad-3; i2++)
	for(i1=4; i1<nzpad-3; i1++)
	{
		diff1 = 1.1962890625000f*(tzz[i2][i1]-tzz[i2][i1-1])
			-0.0797526041667f*(tzz[i2][i1+1]-tzz[i2][i1-2])
			+0.0095703125000f*(tzz[i2][i1+2]-tzz[i2][i1-3])
			-0.0006975446429f*(tzz[i2][i1+3]-tzz[i2][i1-4]);
		diff2 = 1.1962890625000f*(txz[i2][i1]-txz[i2-1][i1])
			-0.0797526041667f*(txz[i2+1][i1]-txz[i2-2][i1])
			+0.0095703125000f*(txz[i2+2][i1]-txz[i2-3][i1])
			-0.0006975446429f*(txz[i2+3][i1]-txz[i2-4][i1]);
		diff3 = 1.1962890625000f*(txx[i2][i1]-txx[i2-1][i1])
			-0.0797526041667f*(txx[i2+1][i1]-txx[i2-2][i1])
			+0.0095703125000f*(txx[i2+2][i1]-txx[i2-3][i1])
			-0.0006975446429f*(txx[i2+3][i1]-txx[i2-4][i1]);
		diff4 = 1.1962890625000f*(txz[i2][i1]-txz[i2][i1-1])
			-0.0797526041667f*(txz[i2][i1+1]-txz[i2][i1-2])
			+0.0095703125000f*(txz[i2][i1+2]-txz[i2][i1-3])
			-0.0006975446429f*(txz[i2][i1+3]-txz[i2][i1-4]);
		uvz[i2][i1]+=dt*rho[i2][i1]*(_dz*diff1+_dx*diff2);
		uvx[i2][i1]+=dt*rho[i2][i1]*(_dx*diff3+_dz*diff4);
	}
}


void forward_txx_tzz_txz(float **uvx, float **uvz, float **txx, float **tzz, float **txz, float **vp, float **vs)
/*< forward step: update txx, tzz, txz >*/
{
	int i1, i2;
	float diff1, diff2, diff3, diff4;

#ifdef _OPENMP
#pragma omp parallel for default(none) 	\
	private(i1,i2,diff1,diff2,diff3,diff4)		\
	shared(uvx,uvz,txx,tzz,txz,vp,vs,nxpad,nzpad,dt,_dx,_dz)
#endif
	for(i2=3; i2<nxpad-4; i2++)
	for(i1=3; i1<nzpad-4; i1++)
	{
		diff1 = 1.1962890625000f*(uvz[i2][i1+1]-uvz[i2][i1])
			-0.0797526041667f*(uvz[i2][i1+2]-uvz[i2][i1-1])
			+0.0095703125000f*(uvz[i2][i1+3]-uvz[i2][i1-2])
			-0.0006975446429f*(uvz[i2][i1+4]-uvz[i2][i1-3]);
		diff2 = 1.1962890625000f*(uvx[i2][i1+1]-uvx[i2][i1])
			-0.0797526041667f*(uvx[i2][i1+2]-uvx[i2][i1-1])
			+0.0095703125000f*(uvx[i2][i1+3]-uvx[i2][i1-2])
			-0.0006975446429f*(uvx[i2][i1+4]-uvx[i2][i1-3]);
		diff3 = 1.1962890625000f*(uvz[i2+1][i1]-uvz[i2][i1])
			-0.0797526041667f*(uvz[i2+2][i1]-uvz[i2-1][i1])
			+0.0095703125000f*(uvz[i2+3][i1]-uvz[i2-2][i1])
			-0.0006975446429f*(uvz[i2+4][i1]-uvz[i2-3][i1]);
		diff4 = 1.1962890625000f*(uvx[i2+1][i1]-uvx[i2][i1])
			-0.0797526041667f*(uvx[i2+2][i1]-uvx[i2-1][i1])
			+0.0095703125000f*(uvx[i2+3][i1]-uvx[i2-2][i1])
			-0.0006975446429f*(uvx[i2+4][i1]-uvx[i2-3][i1]);
		txx[i2][i1]+=dt*(vp[i2][i1]*_dx*diff4+(vp[i2][i1]-2*vs[i2][i1])*_dz*diff1);
		tzz[i2][i1]+=dt*(vp[i2][i1]*_dz*diff1+(vp[i2][i1]-2*vs[i2][i1])*_dx*diff4);
		txz[i2][i1]+=dt*vs[i2][i1]*(_dz*diff2+_dx*diff3);
	}
}


void apply_sponge(float **u, float *bndr)
/*< apply absorbing boundary condition >*/
{
	int ix,iz;

#ifdef _OPENMP
#pragma omp parallel for	    \
    private(ix,iz)		    \
    shared(bndr,u)
#endif
	for(ix=0; ix<nxpad; ix++)
	{
		for(iz=0;iz<nb;iz++){	// top ABC			
			u[ix][iz]=bndr[iz]*u[ix][iz];
		}
		for(iz=nz+nb;iz<nzpad;iz++){// bottom ABC			
			u[ix][iz]=bndr[nzpad-iz-1]*u[ix][iz];
		}
	}

#ifdef _OPENMP
#pragma omp parallel for	    \
    private(ix,iz)		    \
    shared(bndr,u)
#endif
	for(iz=0; iz<nzpad; iz++)
	{
		for(ix=0;ix<nb;ix++){	// left ABC			
			u[ix][iz]=bndr[ix]*u[ix][iz];
		}	
		for(ix=nx+nb;ix<nxpad;ix++){// right ABC			
			u[ix][iz]=bndr[nxpad-ix-1]*u[ix][iz];
		}	
	}
}

int main(int argc, char* argv[])
{
	bool verb;
	int jt, ft, kt, it, ib, sx, sz, ix, iz;
	float a, *wlt, *bndr;
	float **vp0, **vs0, **rho0, **vp, **vs, **rho, **uvx, **uvz, **txx, **tzz, **txz;

	sf_file Fvp, Fvs, Frho, Fwavx, Fwavz;
    
    	sf_init(argc,argv);
#ifdef _OPENMP
    	omp_init();
#endif

	Fvp = sf_input("in");/* p-wave veloctiy */
	Fvs = sf_input("vs");/* s-wave veloctiy */
	Frho = sf_input("rho");/* density */
	Fwavz = sf_output("out");/* z-component of wavefield */
	Fwavx = sf_output("wavx");/* x-component of wavefield */

    	if(!sf_getbool("verb",&verb)) verb=false;    /* verbosity */
    	if (!sf_histint(Fvp,"n1",&nz)) sf_error("No n1= in input");/* veloctiy model: nz */
    	if (!sf_histint(Fvp,"n2",&nx)) sf_error("No n2= in input");/* veloctiy model: nx */
    	if (!sf_histfloat(Fvp,"d1",&dz)) sf_error("No d1= in input");/* veloctiy model: dz */
    	if (!sf_histfloat(Fvp,"d2",&dx)) sf_error("No d2= in input");/* veloctiy model: dx */
    	if (!sf_getint("nb",&nb)) nb=30; /* thickness of PML boundary */
    	if (!sf_getint("nt",&nt)) sf_error("nt required");/* number of time steps */
    	if (!sf_getint("kt",&kt)) sf_error("kt required");/* record wavefield at time kt */
	if (kt>nt) sf_error("make sure kt<=nt");
    	if (!sf_getfloat("dt",&dt)) sf_error("dt required");/* time sampling interval */
    	if (!sf_getfloat("fm",&fm)) fm=20.0; /*dominant freq of Ricker wavelet */
   	if (!sf_getint("ft",&ft)) ft=0; /* first recorded time */
    	if (!sf_getint("jt",&jt)) jt=1;	/* time interval */

	nzpad=nz+2*nb;
	nxpad=nx+2*nb;	
	_dz=1.0/dz;
	_dx=1.0/dx;
	sx=nxpad/2;
	sz=nzpad/2;

	/* allocate memory for variables */
	wlt=sf_floatalloc(nt);
	bndr=sf_floatalloc(nb);
	vp0=sf_floatalloc2(nz,nx); 
	vs0=sf_floatalloc2(nz,nx); 
	rho0=sf_floatalloc2(nz,nx);
	vp=sf_floatalloc2(nzpad, nxpad);
	vs=sf_floatalloc2(nzpad, nxpad);
	rho=sf_floatalloc2(nzpad, nxpad);
	uvx=sf_floatalloc2(nzpad, nxpad);
	uvz=sf_floatalloc2(nzpad, nxpad);
	txx=sf_floatalloc2(nzpad, nxpad);
	tzz=sf_floatalloc2(nzpad, nxpad);
	txz=sf_floatalloc2(nzpad, nxpad);

	/* initialization */
	for(it=0;it<nt;it++)
	{
		a=SF_PI*fm*(it*dt-1.0/fm);a*=a;
		wlt[it]=(1.0-2.0*a)*expf(-a);
	}
	for(ib=0;ib<nb;ib++)
	{
		a=0.015*(nb-ib);
		bndr[ib]=expf(-a*a);
	}
	sf_floatread(vp0[0], nz*nx, Fvp);
	sf_floatread(vs0[0], nz*nx, Fvs);
	sf_floatread(rho0[0], nz*nx, Frho);
	for(ix=0; ix<nx; ix++)
	for(iz=0; iz<nz; iz++)
	{
		vp0[ix][iz]=rho0[ix][iz]*vp0[ix][iz]*vp0[ix][iz];
		vs0[ix][iz]=rho0[ix][iz]*vs0[ix][iz]*vs0[ix][iz];
		rho0[ix][iz]=1.0/rho0[ix][iz];
	}
	expand2d(vp, vp0);
	expand2d(vs, vs0);
	expand2d(rho, rho0);
	memset(uvx[0],0,nzpad*nxpad*sizeof(float));
	memset(uvz[0],0,nzpad*nxpad*sizeof(float));
	memset(txx[0],0,nzpad*nxpad*sizeof(float));
	memset(tzz[0],0,nzpad*nxpad*sizeof(float));
	memset(txz[0],0,nzpad*nxpad*sizeof(float));


	for(it=0; it<nt; it++)
	{
		txx[sx][sz]+=wlt[it];
		tzz[sx][sz]+=wlt[it];

		forward_uvx_uvz(uvx, uvz, txx, tzz, txz, rho);
		forward_txx_tzz_txz(uvx, uvz, txx, tzz, txz, vp, vs);

		apply_sponge(uvz, bndr);
		apply_sponge(uvx, bndr);
		apply_sponge(tzz, bndr);
		apply_sponge(txx, bndr);
		apply_sponge(txz, bndr);

		if (it==kt)
		{
			window2d(vp0, uvx);
			sf_floatwrite(vp0[0], nz*nx, Fwavx);
			window2d(vs0, uvz);
			sf_floatwrite(vs0[0], nz*nx, Fwavz);
		}
		if (verb) sf_warning("%d of %d;", it, nt);
	}

	free(wlt);
	free(bndr);
	free(*vp0); free(vp0);
	free(*vs0); free(vs0);
	free(*rho0); free(rho0);
	free(*vp); free(vp);
	free(*vs); free(vs);
	free(*rho); free(rho);
	free(*uvx); free(uvx);
	free(*uvz); free(uvz);
	free(*txx); free(txx);
	free(*tzz); free(tzz);
	free(*txz); free(txz);

    	exit(0);
}

