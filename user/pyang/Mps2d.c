/* 2-D attenuating wavefield simulation using Fourier Pseudo Spectral method 
   for computing fractional laplacian instead of fractional time derivative
*/
/*
  Copyright (C) 2016 Univ. Grenoble Alpes, Pengliang Yang

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

  Reference: 
  [1] Chen, W., and S. Holm, 2004, Fractional Laplacian time-space models for
  linear and nonlinear lossy media exhibiting arbitrary frequency power-law
  dependency: The Journal of the Acoustical Society of America, 115, no. 4
  [2] J.M. Carcione, 2010, A generalization of the Fourier pseudospectral method,
  GEOPHYSICS, VOL. 75, NO. 6
  [3] H. Igel, Pseudo Spectral tutorial
*/
#include <rsf.h>
#include <math.h>
#include <complex.h>

#ifdef SF_HAS_FFTW
#include <fftw3.h>

void expand2d(float** b, float** a, int nz, int nx, int nb)
/*< expand domain of 'a' to 'b': source(a)-->destination(b) >*/
{
    int nxpad,nzpad,iz,ix;

    nzpad=nz+2*nb;
    nxpad=nx+2*nb;

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


void window2d(float **a, float **b, int nz, int nx, int nb)
/*< window 'b' to 'a': source(b)-->destination(a) >*/
{
    int iz,ix;

    for     (ix=0;ix<nx;ix++) {
	for (iz=0;iz<nz;iz++) {
	    a[ix][iz]=b[nb+ix][nb+iz] ;
	}
    }
}

void apply_sponge(float **p0, float *bndr, int nzpad, int nxpad, int nb)
/* apply absorbing boundary condition */
{
    int ix,iz,ib,ibx,ibz;
    float w;

    for(ib=0; ib<nb; ib++) {
	w = bndr[ib];

	ibz = nzpad-ib-1;
	for(ix=0; ix<nxpad; ix++) {
	    p0[ix][ib ] *= w; /*    top sponge */
	    p0[ix][ibz] *= w; /* bottom sponge */
	}

	ibx = nxpad-ib-1;
	for(iz=0; iz<nzpad; iz++) {
	    p0[ib ][iz] *= w; /*   left sponge */
	    p0[ibx][iz] *= w; /*  right sponge */
	}
    }
}


void add_source(int *sxz, float **p, int ns, int nb, int nz, float *source, bool add)
/*< add source or subtract source >*/
{
    int is, sx, sz;
    if(add){
	for(is=0;is<ns; is++){
	    sx=sxz[is]/nz+nb;
	    sz=sxz[is]%nz+nb;
	    p[sx][sz]+=source[is];
	}
    }else{
	for(is=0;is<ns; is++){
	    sx=sxz[is]/nz+nb;
	    sz=sxz[is]%nz+nb;
	    p[sx][sz]-=source[is];
	}
    }
}

//PML-like Sponge ABC
void sponge_coeff(float *bndr, int nb,float dt,float dx, float dz, float vmax)
{
    int ib;
    float Rc=1.e-5;
    float x, L=nb* SF_MAX(dx,dz);
    float d0=-3.*vmax*logf(Rc)/(2.*L);

    for(ib=0;ib<nb;ib++){
	x=(nb-ib)*dx/L;    
	bndr[ib]=expf(-d0*x*x*dt);
    }
}


void step_forward(float **p0, float **p1, float **vv, float dz, float dx,int nzpad, int nxpad)
/*< forward modeling step >*/
{
    int ix,iz;
    float tmp;
    tmp = 1.0/(dz*dz);
    float c11 = 4.0*tmp/3.0;
    float c12= -tmp/12.0;
    tmp = 1.0/(dx*dx);
    float c21 = 4.0*tmp/3.0;
    float c22= -tmp/12.0;
    float c0=-2.0*(c11+c12+c21+c22);

    for (ix=2; ix < nxpad-2; ix++) 
	for (iz=2; iz < nzpad-2; iz++) 
	{
	    tmp =	c0*p1[ix][iz]+
		c11*(p1[ix][iz-1]+p1[ix][iz+1])+
		c12*(p1[ix][iz-2]+p1[ix][iz+2])+
		c21*(p1[ix-1][iz]+p1[ix+1][iz])+
		c22*(p1[ix-2][iz]+p1[ix+2][iz]);
	    p0[ix][iz]=2.*p1[ix][iz]-p0[ix][iz]+vv[ix][iz]*tmp;
	}
}

int main(int argc, char* argv[])
{
    int nt,nz,nx,nb,nzpad,nxpad;
    int ns,*sxz;
    int iz,ix,it, kt;
    float fm,dt,dz,dx,tmp,vmax,dkz,dkx;
    float *bndr,*wlt,**v0, **vv, **p0, **p1, **beta, **ptr;
    sf_file Fv,FQp, Fw; 
    float *kx, *kz;
    fftwf_plan plan_forw,plan_inv;
    fftwf_complex *fttmp;
	
    /* Madagascar initialization */
    sf_init(argc,argv);	

    /* setup I/O files */
    Fv=sf_input("in");	/* read the data to be interpolated */
    FQp=sf_input("Qp");	/* quality factor */
    Fw=sf_output("out"); 	/* output the reconstructed data */
 
    if (!sf_histint(Fv,"n1",&nz)) sf_error("No n1= in input");
    /* veloctiy model: nz */
    if (!sf_histint(Fv,"n2",&nx)) sf_error("No n2= in input");
    /* veloctiy model: nx */
    if (!sf_histfloat(Fv,"d1",&dz)) sf_error("No d1= in input");
    /* veloctiy model: dz */
    if (!sf_histfloat(Fv,"d2",&dx)) sf_error("No d2= in input");
    /* veloctiy model: dx */
    if (!sf_getint("nb",&nb)) nb=20;
    /* thickness of sponge ABC */
    if (!sf_getint("nt",&nt)) sf_error("nt required");
    /* number of time steps */
    if (!sf_getfloat("dt",&dt)) sf_error("dt required");
    /* time sampling interval */
    if (!sf_getfloat("fm",&fm)) fm=20.0;
    /*dominant freq of Ricker wavelet */
    if (!sf_getint("kt",&kt)) kt=nt;
    
    sf_putint(Fw,"n1",nz);	
    sf_putint(Fw,"n2",nx);


    nzpad=nz+2*nb;
    nxpad=nx+2*nb;
    ns=1;

    v0=sf_floatalloc2(nz, nx);
    beta=sf_floatalloc2(nzpad, nxpad);
    vv=sf_floatalloc2(nzpad, nxpad);
    p0=sf_floatalloc2(nzpad, nxpad);
    p1=sf_floatalloc2(nzpad, nxpad);
    bndr=sf_floatalloc(nb);
    wlt=sf_floatalloc(nt);
    kz=sf_floatalloc(nzpad);
    kx=sf_floatalloc(nxpad);
    fttmp=(fftwf_complex*)fftwf_malloc(sizeof(fftwf_complex)*nzpad*nxpad);
    plan_forw=fftwf_plan_dft_2d(nzpad,nxpad,fttmp,fttmp,FFTW_FORWARD,FFTW_MEASURE);	
    plan_inv=fftwf_plan_dft_2d(nzpad,nxpad,fttmp,fttmp,FFTW_BACKWARD,FFTW_MEASURE);

    sf_floatread(v0[0],nz*nx,Fv);
    expand2d(vv, v0, nz, nx, nb);
    sf_floatread(v0[0],nz*nx,FQp);
    expand2d(beta, v0, nz, nx, nb); //store Q in beta
    vmax=vv[0][0];
    for(ix=0;ix<nxpad;ix++){
	for(iz=0;iz<nzpad;iz++){
	    if(vmax<vv[ix][iz]) vmax=vv[ix][iz];
	    tmp=vv[ix][iz]*dt;
	    vv[ix][iz]=tmp*tmp;// vv=vv^2*dt^2

	    //J.M. Carcione, 2010, equations 9 and 12
	    tmp=atanf(1./beta[ix][iz])/SF_PI; //gamma
	    beta[ix][iz]=1./(1.-tmp);//gamma-->beta
	}
    }

    sponge_coeff(bndr, nb,dt, dx, dz, vmax);
    memset(p0[0],0,nzpad*nxpad*sizeof(float));
    memset(p1[0],0,nzpad*nxpad*sizeof(float));

    wlt=(float*)malloc(nt*sizeof(float));
    for(it=0;it<nt;it++){
	tmp=SF_PI*fm*(it*dt-1.0/fm);tmp=tmp*tmp;
	wlt[it]=(1.0-2.0*tmp)*expf(-tmp);
    }
    sxz=(int*)malloc(ns*sizeof(int));
    sxz[0]=nz/2+nz*(nx/2);

    //pre-compute the discrete wavenumbers: kx and kz
    dkx=2.*SF_PI/(dx*nxpad);
    kx[0]=0;
    for(ix=1; ix<(nxpad+1)/2; ix++) {
	kx[ix]=ix*dkx;
	kx[nxpad-ix]=-ix*dkx;
    }
    if(nxpad%2==0) kx[nxpad/2]=(nxpad/2)*dkx;//Nyquist freq
    dkz=2.*SF_PI/(dz*nzpad);
    kz[0]=0;
    for(iz=1; iz<(nzpad+1)/2; iz++) {
	kz[iz]=iz*dkz;
	kz[nzpad-iz]=-iz*dkz;
    }
    if(nzpad%2==0) kz[nzpad/2]=(nzpad/2)*dkz;//Nyquist freq


    for(it=0; it<nt; it++)    {
	sf_warning("it=%d;",it);

	add_source(sxz, p1, 1, nb,nz, &wlt[it], true);

	//step forward using Fourier pseudo spectral method
	for(ix=0; ix<nxpad; ix++){
	    for(iz=0; iz<nzpad; iz++){
		fttmp[iz+nzpad*ix]=p1[ix][iz];
	    }
	}
	fftwf_execute(plan_forw);// fft
	for(ix=0; ix<nxpad; ix++){
	    for(iz=0; iz<nzpad; iz++){
		//fractional laplacian
		tmp=-(kx[ix]*kx[ix]+kz[iz]*kz[iz]);
		fttmp[iz+nzpad*ix]*= cpowf(tmp,beta[ix][iz]);
	    }
	}
	fftwf_execute(plan_inv);//ifft
	for(ix=0; ix<nxpad; ix++){
	    for(iz=0; iz<nzpad; iz++){
		p0[ix][iz]=2.*p1[ix][iz]-p0[ix][iz]+vv[ix][iz]*crealf(fttmp[iz+nzpad*ix])/(nzpad*nxpad);
	    }
	}
    
	//step_forward(p0, p1, vv,dz, dx, nzpad, nxpad);
	ptr=p0; p0=p1; p1=ptr;

	apply_sponge(p0,bndr,nzpad,nxpad,nb);
	apply_sponge(p1,bndr,nzpad,nxpad,nb);

	if(it==kt){
	    window2d(v0,p1,nz,nx,nb);
	    sf_floatwrite(v0[0],nz*nx,Fw);
	}
    }
    sf_warning(".");
	
    fftwf_destroy_plan(plan_forw);
    fftwf_destroy_plan(plan_inv);
    fftwf_free(fttmp);


    free(sxz);
    free(wlt);
    free(*v0); free(v0);
    free(*beta); free(beta);
    free(*vv); free(vv);
    free(*p0); free(p0);
    free(*p1); free(p1);
    free(bndr);
    free(kx);
    free(kz);

    exit(0);
}
#endif
