/* Forward modeling with Maxwell attenuation model and sponge ABC
*/
/*
  Copyright (C) 2015  University Joseph Fourier, Grenoble (Pengliang Yang)

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

static int nb, nz, nx, nzpad, nxpad, nt;
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


void step_forward(float **p,float **vz, float **vx, float **vv, float **rho)
{
	int i1, i2;
	float tmp, diff1, diff2;

	for(i2=3; i2<nxpad-4; i2++)
	for(i1=3; i1<nzpad-4; i1++)
	{
		diff1=	 1.196289062500000*(p[i2][i1+1]-p[i2][i1])
			-0.079752604166667*(p[i2][i1+2]-p[i2][i1-1])
			+0.009570312500000*(p[i2][i1+3]-p[i2][i1-2])
			-0.000697544642857*(p[i2][i1+4]-p[i2][i1-3]);
		diff2=	 1.196289062500000*(p[i2+1][i1]-p[i2][i1])
			-0.079752604166667*(p[i2+2][i1]-p[i2-1][i1])
			+0.009570312500000*(p[i2+3][i1]-p[i2-2][i1])
			-0.000697544642857*(p[i2+4][i1]-p[i2-3][i1]);
		vz[i2][i1]=vz[i2][i1]-dt*_dz*diff1/rho[i2][i1];
		vx[i2][i1]=vx[i2][i1]-dt*_dx*diff2/rho[i2][i1];
	}

	for(i2=4; i2<nxpad-3; i2++)
	for(i1=4; i1<nzpad-3; i1++)
	{
		tmp=vv[i2][i1]; tmp=rho[i2][i1]*tmp*tmp;
		diff1=	 1.196289062500000*(vz[i2][i1]-vz[i2][i1-1])
			-0.079752604166667*(vz[i2][i1+1]-vz[i2][i1-2])
			+0.009570312500000*(vz[i2][i1+2]-vz[i2][i1-3])
			-0.000697544642857*(vz[i2][i1+3]-vz[i2][i1-4]);
		diff2=	 1.196289062500000*(vx[i2][i1]-vx[i2-1][i1])
			-0.079752604166667*(vx[i2+1][i1]-vx[i2-2][i1])
			+0.009570312500000*(vx[i2+2][i1]-vx[i2-3][i1])
			-0.000697544642857*(vx[i2+3][i1]-vx[i2-4][i1]);
		p[i2][i1]=p[i2][i1]-dt*tmp*(_dz*diff1+_dx*diff2);
	}
}

void apply_attenuation(float **p, float **eta, float **rho, float **vv, float dt)
{
	int i1, i2;
	float a,tau;

	for(i2=0; i2<nxpad; i2++)
	for(i1=0; i1<nzpad; i1++)
	{
		a=rho[i2][i1]*vv[i2][i1]*vv[i2][i1];
		tau=eta[i2][i1]/a; 
		a=expf(-dt/tau);
		p[i2][i1]=a*p[i2][i1];
	}
}

void add_sources(float **p, float dt, float wlt, int sz, int sx)
{
	p[sx][sz]+=dt*wlt;
}


int main(int argc, char* argv[])
{
	bool order1;
	int it, ib, sx, sz;
	float tmp, *wlt, *bndr;
	float **v0, **vv, **rho, **eta, **p, **vz, **vx;
	sf_file Fv, Fw, Frho, Feta;

    	sf_init(argc,argv);
#ifdef _OPENMP
    	omp_init();
#endif

	Fv = sf_input("in");/* veloctiy model */
	Fw = sf_output("out");/* wavefield snaps */
	Frho=sf_input("rho"); /* density */
	Feta=sf_input("eta");/* eta */

    	if (!sf_histint(Fv,"n1",&nz)) sf_error("No n1= in input");/* veloctiy model: nz */
    	if (!sf_histint(Fv,"n2",&nx)) sf_error("No n2= in input");/* veloctiy model: nx */
    	if (!sf_histfloat(Fv,"d1",&dz)) sf_error("No d1= in input");/* veloctiy model: dz */
    	if (!sf_histfloat(Fv,"d2",&dx)) sf_error("No d2= in input");/* veloctiy model: dx */
    	if (!sf_getint("nb",&nb)) nb=30; /* thickness of PML ABC */
    	if (!sf_getint("nt",&nt)) sf_error("nt required");/* number of time steps */
    	if (!sf_getfloat("dt",&dt)) sf_error("dt required");/* time sampling interval */
    	if (!sf_getfloat("fm",&fm)) fm=20.0; /*dominant freq of Ricker wavelet */
    	if(!sf_getbool("order1",&order1)) order1=true;/* (1st order) or (2nd order) accuracy */

	sf_putint(Fw,"n1",nz);
	sf_putint(Fw,"n2",nx);
    	sf_putint(Fw,"n3",nt);
    	sf_putfloat(Fw,"d3",dt);
    	sf_putfloat(Fw,"o3",0);

	_dx=1./dx;
	_dz=1./dz;
	nzpad=nz+2*nb;
	nxpad=nx+2*nb;
	sx=nxpad/2;
	sz=nzpad/2;

	wlt=sf_floatalloc(nt);
	bndr=sf_floatalloc(nb);
	v0=sf_floatalloc2(nz,nx); 	
	vv=sf_floatalloc2(nzpad, nxpad);
	rho=sf_floatalloc2(nzpad, nxpad);
	eta=sf_floatalloc2(nzpad, nxpad);
	p =sf_floatalloc2(nzpad, nxpad);
	vz=sf_floatalloc2(nzpad, nxpad);
	vx=sf_floatalloc2(nzpad, nxpad);

	for(it=0;it<nt;it++){
		tmp=SF_PI*fm*(it*dt-1.0/fm);tmp*=tmp;
		wlt[it]=(1.0-2.0*tmp)*expf(-tmp);
	}
	for(ib=0;ib<nb;ib++){
		tmp=0.015*(nb-ib);
		bndr[ib]=expf(-tmp*tmp);
	}
	sf_floatread(v0[0],nz*nx,Fv);
	expand2d(vv, v0);
	sf_floatread(v0[0],nz*nx,Frho);
	expand2d(rho, v0);
	sf_floatread(v0[0],nz*nx,Feta);
	expand2d(eta, v0);
	memset(p [0],0,nzpad*nxpad*sizeof(float));
	memset(vx[0],0,nzpad*nxpad*sizeof(float));
	memset(vz[0],0,nzpad*nxpad*sizeof(float));

	for(it=0; it<nt; it++)
	{
		window2d(v0,p);
		sf_floatwrite(v0[0],nz*nx,Fw);

		if(order1){// scheme 1, 1st order accuracy, default
			step_forward(p, vz, vx, vv, rho);
			apply_attenuation(p, eta, rho, vv, dt);			
		}else{// scheme 2, 2nd order accuracy
			apply_attenuation(p, eta, rho, vv, 0.5*dt);
			step_forward(p, vz, vx, vv, rho);
			apply_attenuation(p, eta, rho, vv, 0.5*dt);
		}
		add_sources(p, dt, wlt[it], sz, sx);

		// apply sponge/Gaussian taper boundary condition
		apply_sponge(p, bndr);
		apply_sponge(vx, bndr);
		apply_sponge(vz, bndr);
	}

	free(wlt);
	free(bndr);
	free(*v0); free(v0);
	free(*vv); free(vv);
	free(*rho); free(rho);
	free(*eta); free(eta);
	free(*p); free(p);
	free(*vx); free(vx);
	free(*vz); free(vz);

    	exit(0);
}

