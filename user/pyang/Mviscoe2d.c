/* 2D 4-th order visco-elastic wave propagation using sponge ABC
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

  Reference: JOA Robertsson, JO Blanch, WW Symes, Viscoelastic 
  finite-difference modeling, Geophysics 59 (9), 1444-1456
*/
#include <rsf.h>

#ifdef _OPENMP
#include <omp.h>
#endif

const static float c1=1.125, c2=-1./24.;
static int nb, nz, nx, nt, nzpad, nxpad;
static float dz, dx, _dz, _dx, dt, fm;

static void expand2d(float** b, float** a)
/*< expand domain of 'a' to 'b': source(a)-->target(b) >*/
{
    int iz,ix;

#ifdef _OPENMP
#pragma omp parallel for default(none)		\
    private(ix,iz)				\
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
/*< window 'b' to 'a': source(b)-->target(a) >*/
{
    int iz,ix;

#ifdef _OPENMP
#pragma omp parallel for default(none)		\
    private(ix,iz)				\
    shared(b,a,nb,nz,nx)
#endif
    for     (ix=0;ix<nx;ix++) {
	for (iz=0;iz<nz;iz++) {
	    a[ix][iz]=b[nb+ix][nb+iz] ;
	}
    }
}

void variable_transform(float **vp, float **vs, float **rho, 
			float **taup, float ** taus, float **tauo)
/*< variable transformation to reduce compuation >*/
{
    int ix, iz;

    for(ix=0; ix<nxpad; ix++)
	for(iz=0; iz<nzpad; iz++)
	{
	    vp[ix][iz]=rho[ix][iz]*vp[ix][iz]*vp[ix][iz];
	    vs[ix][iz]=rho[ix][iz]*vs[ix][iz]*vs[ix][iz];
	    rho[ix][iz]=1.0/rho[ix][iz];

	    taup[ix][iz]=taup[ix][iz]/tauo[ix][iz];
	    taus[ix][iz]=taus[ix][iz]/tauo[ix][iz];
	    tauo[ix][iz]=1.0/tauo[ix][iz];
	}
}

void average_variable(float **rho, float **taup, float **taus, float **aver)
/*< average the parameters, check guide sofi 2d >*/
{
    int i1, i2;

    memcpy(aver[0],taus[0],nzpad*nxpad*sizeof(float));
    for(i2=1; i2<nxpad-1; i2++)
	for(i1=0; i1<nzpad-1; i1++)
	    aver[i2][i1]=4./(1./taus[i2][i1]+1./taus[i2][i1+1]+1./taus[i2+1][i1]+1./taus[i2+1][i1+1]);
    memcpy(taus[0],aver[0],nxpad*nzpad*sizeof(float));

    memcpy(aver[0],rho[0],nzpad*nxpad*sizeof(float));
    for(i2=1; i2<nxpad-1; i2++)
	for(i1=0; i1<nzpad-1; i1++)
	    aver[i2][i1]=(rho[i2][i1]+rho[i2][i1+1]+rho[i2+1][i1]+rho[i2+1][i1+1])/4.;
    memcpy(rho[0],aver[0],nzpad*nxpad*sizeof(float));

    //......
}


/*****************************************************************************
z corresponds to y in Robertsson's paper. The state variables are capitalized.
	Txx, Tzz, Tzx: sigma_xx, sigma_zz, sigma_zx;
	Vx, Vz:	 v_x, v_y;
	Rxx, Rzz, Rzx： r_xx, r_zz, r_zx;
*******************************************************************************/
void forward_Txx_Tzz_Txz(float **Txx, float **Tzz, float **Txz, 
			 float **Rxx, float **Rzz, float **Rxz, 
			 float **Vx, float **Vz, 
			 float **taup, float **taus, float **tauo, 
			 float **vp, float **vs)
/*< forward step: update Txx, Tzz, Txz, Rxx, Rzz, Rxz >*/
{
    int i1, i2;
    float DiffVz_z,DiffVx_z,DiffVz_x,DiffVx_x, tmpxx, tmpxz, tmpzz,tmp;

#ifdef _OPENMP
#pragma omp parallel for default(none)					\
    private(i1,i2,DiffVz_z, DiffVx_z, DiffVz_x, DiffVx_x,tmpxx, tmpzz,tmpxz,tmp) \
    shared(Vx,Vz,Txx,Tzz,Txz,Rxx, Rzz, Rxz, taup,taus,tauo, vp,vs,nxpad,nzpad,dt,_dx,_dz)
#endif
    for(i2=1; i2<nxpad-2; i2++)
	for(i1=1; i1<nzpad-2; i1++)
	{
	    DiffVz_z = c1*(Vz[i2][i1+1]-Vz[i2][i1])+c2*(Vz[i2][i1+2]-Vz[i2][i1-1]);
	    DiffVx_z = c1*(Vx[i2][i1+1]-Vx[i2][i1])+c2*(Vx[i2][i1+2]-Vx[i2][i1-1]);
	    DiffVz_x = c1*(Vz[i2+1][i1]-Vz[i2][i1])+c2*(Vz[i2+2][i1]-Vz[i2-1][i1]);
	    DiffVx_x = c1*(Vx[i2+1][i1]-Vx[i2][i1])+c2*(Vx[i2+2][i1]-Vx[i2-1][i1]);

	    // Here, mu=rho*vs^2-->vs; pci=lambda+2mu=rho*vp^2-->vp; lambda=rho*(vp^2-vs^2)-->vp-2vs
	    tmpxz=(_dz*DiffVx_z+_dx*DiffVz_x);
	    tmpxx=(_dx*DiffVx_x+_dz*DiffVz_z);
	    tmpzz=(_dx*DiffVx_x+_dz*DiffVz_z);
	    Txz[i2][i1]+=dt*( vs[i2][i1]*taus[i2][i1]*tmpxz+ Rxz[i2][i1] );
	    Txx[i2][i1]+=dt*( vp[i2][i1]*taup[i2][i1]*tmpxx-2.*vs[i2][i1]*taus[i2][i1]*_dz*DiffVz_z+ Rxx[i2][i1] );
	    Tzz[i2][i1]+=dt*( vp[i2][i1]*taup[i2][i1]*tmpzz-2.*vs[i2][i1]*taus[i2][i1]*_dx*DiffVx_x+ Rzz[i2][i1] );

	    tmp=dt*tauo[i2][i1];
	    Rxz[i2][i1]=( (1.-0.5*tmp)*Rxz[i2][i1]-tmp*vs[i2][i1]*(taus[i2][i1]-1.)*tmpxz )/(1.+0.5*tmp);
	    Rxx[i2][i1]=( (1.-0.5*tmp)*Rxx[i2][i1]-tmp*( vp[i2][i1]*(taup[i2][i1]-1.)*tmpxx-2.*vs[i2][i1]*(taus[i2][i1]-1.)*_dz*DiffVz_z ) )/(1.+0.5*tmp);
	    Rzz[i2][i1]=( (1.-0.5*tmp)*Rzz[i2][i1]-tmp*( vp[i2][i1]*(taup[i2][i1]-1.)*tmpzz-2.*vs[i2][i1]*(taus[i2][i1]-1.)*_dx*DiffVx_x ) )/(1.+0.5*tmp);
	}
}

void forward_Vx_Vz(float **Vx, float **Vz, 
		   float **Txx, float **Tzz, float **Txz,
		   float **rho)
/*< forward step: update Vx, Vz >*/
{
    int i1, i2;
    float DiffTzz_z, DiffTzx_x, DiffTxx_x, DiffTzx_z;

#ifdef _OPENMP
#pragma omp parallel for default(none)				\
    private(i1,i2,DiffTzz_z, DiffTzx_x, DiffTxx_x, DiffTzx_z)	\
    shared(Vx,Vz,Txx,Tzz,Txz,rho, nxpad,nzpad,dt,_dx,_dz)
#endif
    for(i2=2; i2<nxpad-1; i2++)
	for(i1=2; i1<nzpad-1; i1++)
	{
	    DiffTzz_z = c1*(Tzz[i2][i1]-Tzz[i2][i1-1])+c2*(Tzz[i2][i1+1]-Tzz[i2][i1-2]);
	    DiffTzx_x = c1*(Txz[i2][i1]-Txz[i2-1][i1])+c2*(Txz[i2+1][i1]-Txz[i2-2][i1]);
	    DiffTxx_x = c1*(Txx[i2][i1]-Txx[i2-1][i1])+c2*(Txx[i2+1][i1]-Txx[i2-2][i1]);
	    DiffTzx_z = c1*(Txz[i2][i1]-Txz[i2][i1-1])+c2*(Txz[i2][i1+1]-Txz[i2][i1-2]);
	    Vz[i2][i1]+=dt*rho[i2][i1]*(_dz*DiffTzz_z+_dx*DiffTzx_x);//here, rho refers to 1/rho
	    Vx[i2][i1]+=dt*rho[i2][i1]*(_dz*DiffTzx_z+_dx*DiffTxx_x);//here, rho refers to 1/rho
	}
}



void apply_sponge(float **u, float *bndr)
/*< apply absorbing boundary condition >*/
{
    int ix,iz;

#ifdef _OPENMP
#pragma omp parallel for			\
    private(ix,iz)				\
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
#pragma omp parallel for			\
    private(ix,iz)				\
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
    int jt, ft, kt, it, ib, sx, sz;
    float a, *wlt, *bndr;
    float **tmp, **vp, **vs, **rho, **Vx, **Vz, **Txx, **Tzz, **Txz, **Rxx, **Rzz, **Rxz, **taup, **taus, **tauo, **aver;
    sf_file Fvp, Fvs, Frho, Ftaup, Ftaus, Ftauo, Fwavx, Fwavz;
    
    sf_init(argc,argv);
#ifdef _OPENMP
    omp_init();
#endif

    Fvp = sf_input("in");/* p-wave veloctiy */
    Fvs = sf_input("vs");/* s-wave veloctiy */
    Frho = sf_input("rho");/* density */
    Ftaup = sf_input("taup");/* p-wave delay */
    Ftaus = sf_input("taus");/* s-wave delay */
    Ftauo = sf_input("tauo");/* delay */
    Fwavz = sf_output("out");/* z-component of wavefield */
    Fwavx = sf_output("wavx");/* x-component of wavefield */

    if (!sf_getbool("verb",&verb)) verb=false;    /* verbosity */
    if (!sf_histint(Fvp,"n1",&nz)) sf_error("No n1= in input");/* veloctiy model: nz */
    if (!sf_histint(Fvp,"n2",&nx)) sf_error("No n2= in input");/* veloctiy model: nx */
    if (!sf_histfloat(Fvp,"d1",&dz)) sf_error("No d1= in input");/* veloctiy model: dz */
    if (!sf_histfloat(Fvp,"d2",&dx)) sf_error("No d2= in input");/* veloctiy model: dx */
    if (!sf_getint("nb",&nb)) nb=30; /* thickness of sponge ABC */
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
    tmp=sf_floatalloc2(nz,nx);
    vp=sf_floatalloc2(nzpad, nxpad);
    vs=sf_floatalloc2(nzpad, nxpad);
    rho=sf_floatalloc2(nzpad, nxpad);
    taup=sf_floatalloc2(nzpad, nxpad);
    taus=sf_floatalloc2(nzpad, nxpad);
    tauo=sf_floatalloc2(nzpad, nxpad);
    Vx=sf_floatalloc2(nzpad, nxpad);
    Vz=sf_floatalloc2(nzpad, nxpad);
    Txx=sf_floatalloc2(nzpad, nxpad);
    Tzz=sf_floatalloc2(nzpad, nxpad);
    Txz=sf_floatalloc2(nzpad, nxpad);
    Rxx=sf_floatalloc2(nzpad, nxpad);
    Rzz=sf_floatalloc2(nzpad, nxpad);
    Rxz=sf_floatalloc2(nzpad, nxpad);
    aver=sf_floatalloc2(nzpad, nxpad);

    /* initialization */
    for(it=0;it<nt;it++){
	a=SF_PI*fm*(it*dt-1.0/fm);a*=a;
	wlt[it]=(1.0-2.0*a)*expf(-a);
    }
    for(ib=0;ib<nb;ib++)	{
	a=0.015*(nb-ib);
	bndr[ib]=expf(-a*a);
    }
    sf_floatread(tmp[0], nz*nx, Fvp);
    expand2d(vp, tmp);
    sf_floatread(tmp[0], nz*nx, Fvs);
    expand2d(vs, tmp);
    sf_floatread(tmp[0], nz*nx, Frho);
    expand2d(rho, tmp);
    sf_floatread(tmp[0], nz*nx, Ftaup);
    expand2d(taup, tmp);
    sf_floatread(tmp[0], nz*nx, Ftaus);
    expand2d(taus, tmp);
    sf_floatread(tmp[0], nz*nx, Ftauo);
    expand2d(tauo, tmp);
    variable_transform(vp, vs, rho, taup, taus, tauo);
    memset(Vx[0],0,nzpad*nxpad*sizeof(float));
    memset(Vz[0],0,nzpad*nxpad*sizeof(float));
    memset(Txx[0],0,nzpad*nxpad*sizeof(float));
    memset(Tzz[0],0,nzpad*nxpad*sizeof(float));
    memset(Txz[0],0,nzpad*nxpad*sizeof(float));
    memset(Rxx[0],0,nzpad*nxpad*sizeof(float));
    memset(Rzz[0],0,nzpad*nxpad*sizeof(float));
    memset(Rxz[0],0,nzpad*nxpad*sizeof(float));


    for(it=0; it<nt; it++)
    {
	Txx[sx][sz]+=wlt[it];
	Tzz[sx][sz]+=wlt[it];

	forward_Txx_Tzz_Txz(Txx, Tzz, Txz, Rxx, Rzz, Rxz, Vx, Vz, taup, taus, tauo, vp, vs);
	forward_Vx_Vz(Vx, Vz, Txx, Tzz, Txz, rho);

	apply_sponge(Vz, bndr);
	apply_sponge(Vx, bndr);
	apply_sponge(Tzz, bndr);
	apply_sponge(Txx, bndr);
	apply_sponge(Txz, bndr);
	apply_sponge(Rzz, bndr);
	apply_sponge(Rxx, bndr);
	apply_sponge(Rxz, bndr);

	if (it==kt)
	{
	    window2d(tmp, Vx);
	    sf_floatwrite(tmp[0], nz*nx, Fwavx);
	    window2d(tmp, Vz);
	    sf_floatwrite(tmp[0], nz*nx, Fwavz);
	}
	if (verb) sf_warning("%d of %d;", it, nt);
    }

    free(wlt);
    free(bndr);
    free(*tmp); free(tmp);
    free(*vp); free(vp);
    free(*vs); free(vs);
    free(*rho); free(rho);
    free(*taup); free(taup);
    free(*taus); free(taus);
    free(*tauo); free(tauo);
    free(*Vx); free(Vx);
    free(*Vz); free(Vz);
    free(*Txx); free(Txx);
    free(*Tzz); free(Tzz);
    free(*Txz); free(Txz);
    free(*Rxx); free(Rxx);
    free(*Rzz); free(Rzz);
    free(*Rxz); free(Rxz);
    free(*aver); free(aver);

    exit(0);
}

