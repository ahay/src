/* RTM and angle gather (ADCIG) extraction using poynting vector
NB: SPML boundary condition combined with 4-th order finite difference,
effective boundary saving strategy used!
*/
/*
  Copyright (C) 2012  Xi'an Jiaotong University (Pengliang Yang)

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

void  pmlcoeff_init(float *d1z, float *d2x, float vmax)
/*< initialize PML abosorbing coefficients >*/
{
	int ix, iz;
 	float Rc=1.e-5;
    	float x, z, L=nb* SF_MAX(dx,dz);
    	float d0=-3.*vmax*logf(Rc)/(2.*L*L*L);

	for(ix=0; ix<nxpad; ix++)
	{
		x=0;
	    	if (ix>=0 && ix<nb)   x=(ix-nb)*dx;
	    	else if(ix>=nxpad-nb && ix<nxpad) x=(ix-(nxpad-nb-1))*dx;
	    	d2x[ix] = d0*x*x;
	}
	for(iz=0; iz<nzpad; iz++)
	{
		z=0;
	    	if (iz>=0 && iz<nb)   z=(iz-nb)*dz;
	    	else if(iz>=nzpad-nb && iz<nzpad) z=(iz-(nzpad-nb-1))*dz; 
	    	d1z[iz] = d0*z*z;  
	}
}



void step_forward(float **p, float **pz, float **px, float **vz, float **vx, float **vv, float *d1z, float *d2x)
{
	int i1, i2;
	float tmp, diff1, diff2;

#ifdef _OPENMP
#pragma omp parallel for default(none)			\
	private(i1,i2,diff1,diff2)			\
	shared(nzpad, nxpad, p, vz, vx, dt, _dz, _dx, d1z, d2x)
#endif
	for(i2=3; i2<nxpad-4; i2++)
	for(i1=3; i1<nzpad-4; i1++)
	{
		diff1=1.125*(p[i2][i1+1]-p[i2][i1])-0.041666666666667*(p[i2][i1+2]-p[i2][i1-1]);
		diff2=1.125*(p[i2+1][i1]-p[i2][i1])-0.041666666666667f*(p[i2+2][i1]-p[i2-1][i1]);
		vz[i2][i1]=((1.-0.5*dt*d1z[i1])*vz[i2][i1]+dt*_dz*diff1)/(1.+0.5*dt*d1z[i1]);
		vx[i2][i1]=((1.-0.5*dt*d2x[i2])*vx[i2][i1]+dt*_dx*diff2)/(1.+0.5*dt*d2x[i2]);
	}

#ifdef _OPENMP
#pragma omp parallel for default(none)			\
	private(i1,i2,diff1, diff2, tmp)		\
	shared(nzpad, nxpad, vv, p, pz, px, vz, vx, dt, _dz, _dx, d1z, d2x)
#endif
	for(i2=4; i2<nxpad-3; i2++)
	for(i1=4; i1<nzpad-3; i1++)
	{
		tmp=vv[i2][i1]; tmp=tmp*tmp;
		diff2=vx[i2][i1]-vx[i2-1][i1];
		diff1=vz[i2][i1]-vz[i2][i1-1];

		diff1=1.125f*(vz[i2][i1]-vz[i2][i1-1])-0.041666666666667f*(vz[i2][i1+1]-vz[i2][i1-2]);
		diff2=1.125f*(vx[i2][i1]-vx[i2-1][i1])-0.041666666666667f*(vx[i2+1][i1]-vx[i2-2][i1]);
		pz[i2][i1]=((1.-0.5*dt*d1z[i1])*pz[i2][i1]+dt*tmp*_dz*diff1)/(1.+0.5*dt*d1z[i1]);
		px[i2][i1]=((1.-0.5*dt*d2x[i2])*px[i2][i1]+dt*tmp*_dx*diff2)/(1.+0.5*dt*d2x[i2]);
		p[i2][i1]=px[i2][i1]+pz[i2][i1];
	}
}

int main(int argc, char* argv[])
{
	int jt, ft, it, i2, i1, sx, sz;
	float tmp, vmax;
	float *wlt, *d2x, *d1z;
	float **v0, **vv, **p, **pz, **px, **vz, **vx;
	sf_file Fv, Fw;

    	sf_init(argc,argv);
#ifdef _OPENMP
    	omp_init();
#endif

	Fv = sf_input("in");/* veloctiy model */
	Fw = sf_output("out");/* wavefield snaps */

    	if (!sf_histint(Fv,"n1",&nz)) sf_error("No n1= in input");/* veloctiy model: nz */
    	if (!sf_histint(Fv,"n2",&nx)) sf_error("No n2= in input");/* veloctiy model: nx */
    	if (!sf_histfloat(Fv,"d1",&dz)) sf_error("No d1= in input");/* veloctiy model: dz */
    	if (!sf_histfloat(Fv,"d2",&dx)) sf_error("No d2= in input");/* veloctiy model: dx */
    	if (!sf_getint("nb",&nb)) nb=30; /* thickness of PML ABC */
    	if (!sf_getint("nt",&nt)) sf_error("nt required");/* number of time steps */
    	if (!sf_getfloat("dt",&dt)) sf_error("dt required");/* time sampling interval */
    	if (!sf_getfloat("fm",&fm)) fm=20.0; /*dominant freq of Ricker wavelet */
   	if (!sf_getint("ft",&ft)) ft=0; /* first recorded time */
    	if (!sf_getint("jt",&jt)) jt=1;	/* time interval */

	sf_putint(Fw,"n1",nz);
	sf_putint(Fw,"n2",nx);
    	sf_putint(Fw,"n3",(nt-ft)/jt);
    	sf_putfloat(Fw,"d3",jt*dt);
    	sf_putfloat(Fw,"o3",ft*dt);

	_dx=1./dx;
	_dz=1./dz;
	nzpad=nz+2*nb;
	nxpad=nx+2*nb;
	sx=nxpad/2;
	sz=nzpad/2;

	wlt=sf_floatalloc(nt);
	v0=sf_floatalloc2(nz,nx); 	
	vv=sf_floatalloc2(nzpad, nxpad);
	p =sf_floatalloc2(nzpad, nxpad);
	pz=sf_floatalloc2(nzpad, nxpad);
	px=sf_floatalloc2(nzpad, nxpad);
	vz=sf_floatalloc2(nzpad, nxpad);
	vx=sf_floatalloc2(nzpad, nxpad);
	d1z=sf_floatalloc(nzpad);
	d2x=sf_floatalloc(nxpad);

	for(it=0;it<nt;it++){
		tmp=SF_PI*fm*(it*dt-1.0/fm);tmp*=tmp;
		wlt[it]=(1.0-2.0*tmp)*expf(-tmp);
	}
	sf_floatread(v0[0],nz*nx,Fv);
	expand2d(vv, v0);
	memset(p [0],0,nzpad*nxpad*sizeof(float));
	memset(px[0],0,nzpad*nxpad*sizeof(float));
	memset(pz[0],0,nzpad*nxpad*sizeof(float));
	memset(vx[0],0,nzpad*nxpad*sizeof(float));
	memset(vz[0],0,nzpad*nxpad*sizeof(float));
	vmax=v0[0][0];
	for(i2=0; i2<nx; i2++)
	for(i1=0; i1<nz; i1++)
		vmax=SF_MAX(v0[i2][i1],vmax);
	pmlcoeff_init(d1z, d2x, vmax);

	for(it=0; it<nt; it++)
	{
		if(it>=ft)
		{
			window2d(v0,p);
			sf_floatwrite(v0[0],nz*nx,Fw);
		}
		p[sx][sz]+=wlt[it];
		step_forward(p, pz, px, vz, vx, vv, d1z, d2x);
	}

	free(wlt);
	free(*v0); free(v0);
	free(*vv); free(vv);
	free(*p); free(p);
	free(*px); free(px);
	free(*pz); free(pz);
	free(*vx); free(vx);
	free(*vz); free(vz);
	free(d1z);
	free(d2x);

    	exit(0);
}

