/*
  Copyright (C) 2013  Xi'an Jiaotong University (Pengliang Yang)

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
#include <time.h>

#ifdef _OPENMP
#include <omp.h>
#endif


#ifndef MAX
#define	MAX(x,y) ((x) > (y) ? (x) : (y))
#endif
#ifndef MIN
#define	MIN(x,y) ((x) < (y) ? (x) : (y))
#endif


struct fdm2{
    int nb;
    int   nz,nzpad;
    int   nx,nxpad;
    float dz;
    float dx;
    bool verb;
    bool frsf;
};

typedef struct fdm2 *fdm2d;

void expand(float** a, float** b, fdm2d fdm)
/*< expand domain of 'a' to 'b' >*/
{
    int iz,ix;

    for     (ix=0;ix<fdm->nx;ix++) {
	for (iz=0;iz<fdm->nz;iz++) {
	    b[fdm->nb+ix][fdm->nb+iz] = a[ix][iz];
	}
    }

    for     (ix=0; ix<fdm->nxpad; ix++) {
	for (iz=0; iz<fdm->nb;    iz++) {
	    b[ix][           iz  ] = b[ix][           fdm->nb  ];
	    b[ix][fdm->nzpad-iz-1] = b[ix][fdm->nzpad-fdm->nb-1];
	}
    }

    for     (ix=0; ix<fdm->nb;    ix++) {
	for (iz=0; iz<fdm->nzpad; iz++) {
	    b[           ix  ][iz] = b[           fdm->nb  ][iz];
	    b[fdm->nxpad-ix-1][iz] = b[fdm->nxpad-fdm->nb-1][iz];
	}
    }
}


void window(float **a, float **b, fdm2d fdm)
/*< window 'b' to 'a'>*/
{
    int iz,ix;
    for     (ix=0;ix<fdm->nx;ix++) {
	for (iz=0;iz<fdm->nz;iz++) {
	    a[ix][iz]=b[fdm->nb+ix][fdm->nb+iz] ;
	}
    }
}

void init_abc(float **vv, float**bx1,float**bx2,float **bz1,float **bz2,float dt,fdm2d fdm)
{
    float Rc = 1.0e-5;
    float L = fdm->nb* MAX(fdm->dx,fdm->dz);
    float d0= -3.0*logf(Rc)/(2.0*L);

    int iz,ix;    
    for     (ix=0; ix<fdm->nxpad; ix++) {
	for (iz=0; iz<fdm->nb;    iz++) {
		float tmp1=fdm->nb-iz;
		float tmp2=tmp1-0.5;
		tmp1=tmp1/fdm->nb;
		tmp2=tmp2/fdm->nb;
		tmp1=tmp1*tmp1;
		tmp2=tmp2*tmp2;
	    	bz1[ix][           iz ] = expf(-dt*d0*vv[ix][fdm->nb		]*tmp1);
	    	bz1[ix][2*fdm->nb-iz-1] = expf(-dt*d0*vv[ix][fdm->nzpad-fdm->nb-1]*tmp2);
	    	bz2[ix][           iz ] = expf(-dt*d0*vv[ix][fdm->nb		]*tmp2);
	    	bz2[ix][2*fdm->nb-iz-1] = expf(-dt*d0*vv[ix][fdm->nzpad-fdm->nb-1]*tmp1);
	}
    }
    
    for     (ix=0; ix<fdm->nb;    ix++) {
	for (iz=0; iz<fdm->nzpad; iz++) {
		float tmp1=fdm->nb-ix;
		float tmp2=tmp1-0.5;
		tmp1=tmp1/fdm->nb;
		tmp2=tmp2/fdm->nb;
		tmp1=tmp1*tmp1;
		tmp2=tmp2*tmp2;
	    	bx1[           ix ][iz] = expf(-dt*d0*vv[fdm->nb		][iz]*tmp1);
	    	bx1[2*fdm->nb-ix-1][iz] = expf(-dt*d0*vv[2*fdm->nxpad-fdm->nb-1	][iz]*tmp2);
	    	bx2[           ix ][iz] = expf(-dt*d0*vv[fdm->nb		][iz]*tmp2);
	    	bx2[2*fdm->nb-ix-1][iz] = expf(-dt*d0*vv[2*fdm->nxpad-fdm->nb-1	][iz]*tmp1);
	}
    }

}
void  init_abc_coef(float **d1z, float **d1x, float **d2z, float **d2x, float **vv, fdm2d fdm)
{
    float Rc = 1.0e-5;
    float L = fdm->nb* MAX(fdm->dx,fdm->dz);
    float d0 = -3.0*log(Rc)/(2.0*L*L*L);

    int ix, iz;
#ifdef _OPENMP
#pragma omp parallel for	\
	private(ix,iz)		\
	shared(d1z,d2x,fdm)
#endif
    for     (ix=0;ix<fdm->nxpad;ix++) {
	for (iz=0;iz<fdm->nzpad;iz++) {
	    d1z[ix][iz] = 0.0;
	    d1x[ix][iz] = 0.0;
	    d2z[ix][iz] = 0.0;
	    d2x[ix][iz] = 0.0;
	}
    }

#ifdef _OPENMP
#pragma omp parallel for	\
	private(ix,iz)		\
	shared(d1z,d2x,vv,fdm)
#endif
    for     (ix=0; ix<fdm->nxpad; ix++) {
	for (iz=0; iz<fdm->nb;    iz++) {
	    float z1=(fdm->nb-iz)*fdm->dz;
	    float z2=(fdm->nb-iz-0.5)*fdm->dz;
	    d1z[ix][           iz  ] = d0*vv[ix][fdm->nb		]*z1*z1;
	    d1z[ix][fdm->nzpad-iz-1] = d0*vv[ix][fdm->nzpad-fdm->nb-1	]*z2*z2;
	    d2z[ix][           iz  ] = d0*vv[ix][fdm->nb		]*z2*z2;
	    d2z[ix][fdm->nzpad-iz-1] = d0*vv[ix][fdm->nzpad-fdm->nb-1	]*z1*z1;
	}
    }

#ifdef _OPENMP
#pragma omp parallel for	\
	private(ix,iz)		\
	shared(d1z,d2x,vv,fdm)
#endif
    for     (ix=0; ix<fdm->nb;    ix++) {
	for (iz=0; iz<fdm->nzpad; iz++) {
	    float x1=(fdm->nb-ix)*fdm->dx;
	    float x2=(fdm->nb-ix-0.5)*fdm->dx;
	    d1x[           ix  ][iz] = d0*vv[fdm->nb		][iz]*x1*x1;
	    d1x[fdm->nxpad-ix-1][iz] = d0*vv[fdm->nxpad-fdm->nb-1][iz]*x2*x2;
	    d2x[           ix  ][iz] = d0*vv[fdm->nb		][iz]*x2*x2;
	    d2x[fdm->nxpad-ix-1][iz] = d0*vv[fdm->nxpad-fdm->nb-1][iz]*x1*x1;
	}
    }
}


//step: forward modelling at time t=it*dt
void step_forward(float **u, float **ux,float **uz,float **ax,float **az,float **d1z,float **d1x,float **d2z,float **d2x,float **vv, fdm2d fdm,float dt)
{
    int ix,iz,frsf=0;
    if(fdm->frsf) frsf=1;

    // order:   NJ=8
	float c1 = 1.1962890625000;
	float c2 =-0.0797526041667;
	float c3 = 0.0095703125000;
	float c4 =-0.0006975446429;

    // updata ax and az, at time t=it*dt
#ifdef _OPENMP
#pragma omp parallel for	\
	private(iz,ix)		\
	shared(az,ax,ux,uz,u,d2x,d2z,fdm,c1,c2,c3,c4,dt)
#endif
    for(ix=3;ix<fdm->nxpad-4;ix++)    {
        for(iz=3+frsf*(fdm->nb-3);iz<fdm->nzpad-4;iz++){
            az[ix][iz] = (  (1.0-0.5*dt*d2z[ix][iz])*az[ix][iz]
                    +dt*( c1*(u[ix][iz+1]-u[ix][iz])
                        +c2*(u[ix][iz+2]-u[ix][iz-1])
                        +c3*(u[ix][iz+3]-u[ix][iz-2])
                        +c4*(u[ix][iz+4]-u[ix][iz-3])
                    )/fdm->dz  )/(1.0+0.5*dt*d1z[ix][iz]);
            ax[ix][iz] = ( (1.0-0.5*dt*d2x[ix][iz])*ax[ix][iz]
                     +dt*( c1*(u[ix+1][iz]-u[ix][iz])
                        +c2*(u[ix+2][iz]-u[ix-1][iz])
                        +c3*(u[ix+3][iz]-u[ix-2][iz])
                        +c4*(u[ix+4][iz]-u[ix-3][iz])
                    )/fdm->dx )/(1.0+0.5*dt*d2x[ix][iz]);
        }
    }

            // updata ux and uz, at time t=it*dt
#ifdef _OPENMP
#pragma omp parallel for	\
	private(iz,ix)		\
	shared(az,ax,ux,uz,u,vv,d1x,d1z,fdm,c1,c2,c3,c4,dt)
#endif
    for(ix=4;ix<fdm->nxpad-3;ix++)    {
        for(iz=4+frsf*(fdm->nb-4);iz<fdm->nzpad-3;iz++){
            float tmp = vv[ix][iz]*vv[ix][iz];
            uz[ix][iz] = (   (1.0-0.5*dt*d1z[ix][iz])*uz[ix][iz]+
                        tmp*dt*( c1*(az[ix][iz]-az[ix][iz-1])
                                +c2*(az[ix][iz+1]-az[ix][iz-2])
                                +c3*(az[ix][iz+2]-az[ix][iz-3])
                                +c4*(az[ix][iz+3]-az[ix][iz-4])  )/fdm->dz
                                )/(1.0+0.5*dt*d1z[ix][iz]);
            ux[ix][iz] = (   (1.0-0.5*dt*d1x[ix][iz])*ux[ix][iz]+
                        tmp*dt*( c1*(ax[ix][iz]-ax[ix-1][iz])
                                +c2*(ax[ix+1][iz]-ax[ix-2][iz])
                                +c3*(ax[ix+2][iz]-ax[ix-3][iz])
                                +c4*(ax[ix+3][iz]-ax[ix-4][iz]) )/fdm->dx
                                )/(1.0+0.5*dt*d1x[ix][iz]);
            u[ix][iz] = ux[ix][iz]+uz[ix][iz];
        }
    }
}


void wavefield_init(float **u,float**uz, float **ux, float **az, float **ax, fdm2d fdm)
{
	int iz,ix;
#ifdef _OPENMP
#pragma omp parallel for	\
	private(iz,ix)		\
	shared(az,ax,ux,uz,u,fdm)
#endif
	for(ix=0;ix<fdm->nxpad;ix++)	{
		for(iz=0;iz<fdm->nzpad;iz++){
		 	u[ix][iz]=0.0;
		 	uz[ix][iz]=0.0;
		 	ux[ix][iz]=0.0;
		 	az[ix][iz]=0.0;
		 	ax[ix][iz]=0.0;
		}
	}
}


int main(int argc, char* argv[])
{
    	fdm2d fdm;
    	fdm = (fdm2d) sf_alloc(1,sizeof(*fdm));
    	fdm->frsf=true;
    	fdm->verb=false;
    	fdm->nb=20;
    	fdm->nz=300;
    	fdm->nx=400;
    	fdm->dz=5;
    	fdm->dx=5;
    	fdm->nzpad=fdm->nz+2*fdm->nb;
    	fdm->nxpad=fdm->nx+2*fdm->nb;

	int nt=1800;
	float dt=0.001;
	float fm=25.0;
	bool csdgather=true;
	int ns=20;
	int ng=100;

	int jsx=15;
	int jsz=0;
	int jgx=1;
	int jgz=0;
	int sxbeg=50;//x-begin point of source, index starting from 0
	int szbeg=30;//z-begin point of source, index starting from 0
	int gxbeg=0;//x-begin point of geophone, index starting from 0
	int gzbeg=2;//z-begin point of geophone, index starting from 0

	int *sz,*sx,*gz,*gx;
	sz=sf_intalloc(ns);
	sx=sf_intalloc(ns);
	gz=sf_intalloc(ng);
	gx=sf_intalloc(ng);

	if (!(sxbeg>=0 && szbeg>=0 && sxbeg+(ns-1)*jsx<fdm->nx && szbeg+(ns-1)*jsz<fdm->nz))	
	{ printf("sources exceeds the computing zone!\n"); exit(1);}
	for(int is=0;is<ns;is++) {
		sz[is]=szbeg+is*jsz;
		sx[is]=sxbeg+is*jsx;
	}
	
	int distx=sxbeg-gxbeg;
	int distz=szbeg-gzbeg;
	if (csdgather)
	{
		//distance between source and geophone at the beginning
		if (!(gxbeg>=0 && gzbeg>=0 && gxbeg+(ng-1)*jgx<fdm->nx && gzbeg+(ng-1)*jgz<fdm->nz &&
		(sxbeg+(ns-1)*jsx)+(ng-1)*jgx-distx <fdm->nx  && (szbeg+(ns-1)*jsz)+(ng-1)*jgz-distz <fdm->nz))	
		{ printf("geophones exceeds the computing zone!\n"); exit(1);}		
		for(int is=0;is<ng;is++) {
			gz[is]=gzbeg+is*jgz;
			gx[is]=gxbeg+is*jgx;
		}
	}

	float *wlt;
	wlt=sf_floatalloc(nt);
	for(int it=0;it<nt;it++){
		float a=SF_PI*fm*(it*dt-1.0/fm);a*=a;
		wlt[it]=(1.0-2.0*a)*expf(-a);
	}


	float **v0, **Iss, **Isg, **Img, **seis;
	float **vv, **d1z, **d1x, **d2z, **d2x;
	float **su, **sux, **suz, **sax, **saz;
	float **gu, **gux, **guz, **gax, **gaz;

        v0 = sf_floatalloc2(fdm->nz,fdm->nx);
        Iss = sf_floatalloc2(fdm->nz,fdm->nx);
        Isg = sf_floatalloc2(fdm->nz,fdm->nx);
        Img = sf_floatalloc2(fdm->nz,fdm->nx);
        seis = sf_floatalloc2(ng,nt);
        vv = sf_floatalloc2(fdm->nzpad,fdm->nxpad);
        d1z = sf_floatalloc2(fdm->nzpad,fdm->nxpad);
        d1x = sf_floatalloc2(fdm->nzpad,fdm->nxpad);
        d2z = sf_floatalloc2(fdm->nzpad,fdm->nxpad);
        d2x = sf_floatalloc2(fdm->nzpad,fdm->nxpad);
        su  = sf_floatalloc2(fdm->nzpad,fdm->nxpad);
        sux = sf_floatalloc2(fdm->nzpad,fdm->nxpad);
        suz = sf_floatalloc2(fdm->nzpad,fdm->nxpad);
        sax = sf_floatalloc2(fdm->nzpad,fdm->nxpad);
        saz = sf_floatalloc2(fdm->nzpad,fdm->nxpad);
        gu  = sf_floatalloc2(fdm->nzpad,fdm->nxpad);
        gux = sf_floatalloc2(fdm->nzpad,fdm->nxpad);
        guz = sf_floatalloc2(fdm->nzpad,fdm->nxpad);
        gax = sf_floatalloc2(fdm->nzpad,fdm->nxpad);
        gaz = sf_floatalloc2(fdm->nzpad,fdm->nxpad);

	for(int ix=0;ix<fdm->nx;ix++){
		for(int iz=0;iz<fdm->nz;iz++){
			if (iz<150) v0[ix][iz]=1800;
			else	v0[ix][iz]=2000;	
			Iss[ix][iz]=0.0;
			Isg[ix][iz]=0.0;
			Img[ix][iz]=0.0;	
		}
	}
  	expand(v0, vv, fdm);
	init_abc_coef(d1z, d1x, d2z, d2x, vv, fdm);
	wavefield_init(su, suz, sux, saz, sax, fdm);
	wavefield_init(gu, guz, gux, gaz, gax, fdm);

	clock_t start_t, end_t;
	FILE *fp;
	for(int is=0; is<ns; is++)
	{
	   start_t = clock();   

		fp=fopen("wav.dat","wb");
		if (fp==NULL) { printf("cannot open the file"); exit(1);}

		if (csdgather){
			gxbeg=sx[is]-distx;		
			for(int ig=0;ig<ng;ig++) gx[ig]=gxbeg+ig*jgx;
		}
		wavefield_init(su, suz, sux, saz, sax, fdm);
		for(int it=0; it<nt; it++){
		    // add the wavelet at the source locations
		    su[sx[is]+fdm->nb][sz[is]+fdm->nb]+=wlt[it];
		    step_forward(su, suz, sux, saz, sax, d1z, d1x, d2z, d2x, vv, fdm, dt);

		    for(int ig=0; ig<ng; ig++) seis[it][ig]=su[gx[ig]+fdm->nb][gz[ig]+fdm->nb];
		    window(v0, su, fdm);
		    fwrite(v0[0], sizeof(float), fdm->nz*fdm->nx, fp);
		}
		fclose(fp);

		fp=fopen("wav.dat","rb");
		if (fp==NULL) { printf("cannot open the file"); exit(1);}
		wavefield_init(gu, guz, gux, gaz, gax, fdm);
		for(int it=nt-1; it>-1; it--)
		{
		    // add the wavelet at the source locations
		    for(int ig=0;ig<ng;ig++) gu[gx[ig]+fdm->nb][gz[ig]+fdm->nb]+=seis[it][ig];
		    step_forward(gu, guz, gux, gaz, gax, d1z, d1x, d2z, d2x, vv, fdm, dt);

		    fseek(fp,it*fdm->nz*fdm->nx*sizeof(float),SEEK_SET);
		    fread(v0[0], sizeof(float), fdm->nz*fdm->nx, fp);
		    for(int ix=0;ix<fdm->nx;ix++){
			for(int iz=0;iz<fdm->nz;iz++){
				Isg[ix][iz]+=v0[ix][iz]*gu[ix+fdm->nb][iz+fdm->nb];
				Iss[ix][iz]+=v0[ix][iz]*v0[ix][iz];
			}
		    }
		}
		fclose(fp);

		for(int ix=0;ix<fdm->nx;ix++){
			for(int iz=0;iz<fdm->nz;iz++){
					Img[ix][iz]+=Isg[ix][iz]/(Iss[ix][iz]+SF_EPS);
			}
		}

	   end_t = clock();
	   printf("CPU Time of this shot: %f s\n",  (double)(end_t - start_t) / CLOCKS_PER_SEC);
	}

	fp=fopen("img.dat","wb");
	if (fp==NULL) { printf("cannot open the file"); exit(1);}
	fwrite(Img[0], sizeof(float), fdm->nz*fdm->nx, fp);

    	exit(0);
}

