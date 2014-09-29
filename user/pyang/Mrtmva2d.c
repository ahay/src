/* RTM with checkpointing in 2D visco-acoustic media
The wavefield reconstruction method can not be utilized in visco-acoustic
 and visco-elastic wave equation due to the dissipation. The solution of
 computation without disk I/O is the use of checkpointing technique.
*/
/*
  Copyright (C) 2014 Xi'an Jiaotong University (Pengliang Yang)

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
  
  Reference: William Symes, Reverse time migration with optimal checkpointing,
	Geophysics, v. 72 no. 5 p. SM213-SM221 doi: 10.1190/1.2742686 
*/
#include <rsf.h>

#ifdef _OPENMP
#include <omp.h>
#endif

static int nb, nz, nx, nzpad, nxpad, ns, ng;
static float dz, dx, _dz, _dx, dt, vmute;

typedef struct checkpointing checkpoint;
/*^*/
struct  checkpointing{//store 2d arrays as 1d vectors
	float *p;
	float *r;
	float *vx;
	float *vz;
};

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



void step_forward(float **p, float **r, float **vz, float **vx, float **vv, float **rho, float **tau, float **tau0)
/*< forward modeling step >*/
{
	int i1, i2;
	float tmp, tmp2, diff1, diff2;

#ifdef _OPENMP
#pragma omp parallel for default(none)			\
	private(i1,i2,diff1,diff2)			\
	shared(nzpad, nxpad, rho, p, vz, vx, dt, _dz, _dx)
#endif
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
		vz[i2][i1]-=dt*_dz*diff1/rho[i2][i1];
		vx[i2][i1]-=dt*_dx*diff2/rho[i2][i1];
	}

#ifdef _OPENMP
#pragma omp parallel for default(none)			\
	private(i1,i2,diff1, diff2, tmp, tmp2)		\
	shared(nzpad, nxpad, rho, tau, tau0, vv, p, r, vz, vx, dt, _dz, _dx)
#endif
	for(i2=4; i2<nxpad-3; i2++)
	for(i1=4; i1<nzpad-3; i1++)
	{
		tmp=vv[i2][i1]; tmp=tmp*tmp;
		diff1=	 1.196289062500000*(vz[i2][i1]-vz[i2][i1-1])
			-0.079752604166667*(vz[i2][i1+1]-vz[i2][i1-2])
			+0.009570312500000*(vz[i2][i1+2]-vz[i2][i1-3])
			-0.000697544642857*(vz[i2][i1+3]-vz[i2][i1-4]);
		diff2=	 1.196289062500000*(vx[i2][i1]-vx[i2-1][i1])
			-0.079752604166667*(vx[i2+1][i1]-vx[i2-2][i1])
			+0.009570312500000*(vx[i2+2][i1]-vx[i2-3][i1])
			-0.000697544642857*(vx[i2+3][i1]-vx[i2-4][i1]);
		tmp=tmp*rho[i2][i1]*(_dz*diff1+_dx*diff2);
		tmp2=dt/tau0[i2][i1];
		r[i2][i1]=((1.-0.5*tmp2)*r[i2][i1]-tmp2*tau[i2][i1]*tmp)/(1.+0.5*tmp2);
		p[i2][i1]-=dt*((1.+tau[i2][i1])*tmp+r[i2][i1]);
	}
}

void sg_init(int *sxz, int szbeg, int sxbeg, int jsz, int jsx, int ns)
/*< shot/geophone position initialize
sxz/gxz; szbeg/gzbeg; sxbeg/gxbeg; jsz/jgz; jsx/jgx; ns/ng; >*/
{
	int is, sz, sx;

	for(is=0; is<ns; is++)
	{
		sz=szbeg+is*jsz;
		sx=sxbeg+is*jsx;
		sxz[is]=sz+nz*sx;
	}
}

void add_source(int *sxz, float **p, int ns, float *source, bool add)
/*< add seismic sources in grid >*/
{
	int is, sx, sz;

	if(add){/* add sources*/
#ifdef _OPENMP
#pragma omp parallel for default(none)	\
	private(is,sx,sz)		\
	shared(p,source,sxz,nb,ns,nz)
#endif
		for(is=0;is<ns; is++){
			sx=sxz[is]/nz+nb;
			sz=sxz[is]%nz+nb;
			p[sx][sz]+=source[is];
		}
	}else{ /* subtract sources */
#ifdef _OPENMP
#pragma omp parallel for default(none)	\
	private(is,sx,sz)		\
	shared(p,source,sxz,nb,ns,nz)
#endif
		for(is=0;is<ns; is++){
			sx=sxz[is]/nz+nb;
			sz=sxz[is]%nz+nb;
			p[sx][sz]-=source[is];
		}
	}
}

void cross_correlation(float **image, float **sp, float **gp)
/*< cross correlation >*/
{
	int i1, i2;

#ifdef _OPENMP
#pragma omp parallel for default(none)	\
	private(i1,i2)			\
	shared(nz, nx, nb, sp, gp, image)
#endif
	for(i2=0; i2<nx; i2++)
	for(i1=0; i1<nz; i1++)
	{
		image[i2][i1]+=sp[i2+nb][i1+nb]*gp[i2+nb][i1+nb];
	}
}

void record_seis(float *seis_it, int *gxz, float **p, int ng)
/*< record seismogram at time it into a vector length of ng >*/
{
	int ig, gx, gz;

#ifdef _OPENMP
#pragma omp parallel for default(none)	\
	private(ig,gx,gz)		\
	shared(seis_it,p,gxz,nb,ng,nz)
#endif
	for(ig=0;ig<ng; ig++)
	{
		gx=gxz[ig]/nz+nb;
		gz=gxz[ig]%nz+nb;
		seis_it[ig]=p[gx][gz];
	}
}

void muting(float *seis_kt, int gzbeg, int szbeg, int gxbeg, int sxc, int jgx, int it, int tdmute)
/*< muting the direct arrivals >*/
{
	int id, kt;
	float a,b,t0;

#ifdef _OPENMP
#pragma omp parallel for default(none)	\
	private(id,a,b,t0,kt)		\
	shared(seis_kt,vmute,ng,it,dt,dx,dz,tdmute,gzbeg,szbeg,gxbeg,sxc,jgx)
#endif
	for(id=0;id<ng;id++)
	{
		a=dx*abs(gxbeg+id*jgx-sxc);
		b=dz*(gzbeg-szbeg);
		t0=sqrtf(a*a+b*b)/vmute;
		kt=t0/dt+tdmute;// tdmute manually added to obtain the best muting effect.
    		if (it<kt) seis_kt[id]=0.;
	}
}


/*-----------------------------------------------------------------------------*/
void variable_inverse(float **rho, float **tau0)
/*< inverse of variables >*/
{
	int i1, i2;
	
	for(i2=0; i2<nxpad; i2++)
	for(i1=0; i1<nzpad; i1++)
	{
		rho[i2][i1]=1./rho[i2][i1];
		tau0[i2][i1]=1./tau0[i2][i1];
	}
}

float funI0(float a, float w)
/*< function I0: w (angular freq)>*/
{
	float tmp=w*a;
	return logf(1.+tmp*tmp)/(2.*a);
}

float funI1(float a, float w)
/*< function I1: w (angular freq)>*/
{
	float tmp1=w*a;
	float tmp2=tmp1/(1.+tmp1*tmp1);
	return (atanf(tmp1)-tmp2)/(2.*a);
}

float funI2lk(float a, float b, float w)
/*< function I2: w (angular freq)>*/
{
	float tmp1=atanf(w*a)/a-atanf(w*b)/b;
	float tmp2=a*b/(b*b-a*a);
	return tmp1*tmp2;
}

void compute_tau(float **tau, float **Q0, float ***tau_l, int L, float wa, float wb)
/*<compute tau according to Q0
Reference: Joakim O. Blanch, Johan O. A. Robertsson, William W. Symes, Modeling 
of a constant Q: Methodology and algorithm for an efficient and optimally 
inexpensive viscoelastic technique, GEOPHYSICS Jan 1995, Vol. 60, No. 1, pp. 176-184
 >*/
{
	int i1, i2, l, k;
	float s0, s1, s2;

	for(i2=0; i2<nxpad; i2++)
	for(i1=0; i1<nzpad; i1++)
	{
		s0=s1=s2=0.;
		for(l=0; l<L-1; l++)
		{
			s0+=funI0(tau_l[i2][i1][l], wb)-funI0(tau_l[i2][i1][l], wa);
		 	s1+=funI1(tau_l[i2][i1][l], wb)-funI1(tau_l[i2][i1][l], wa);
			for(k=l+1; k<L; k++)
			s2+=funI2lk(tau_l[i2][i1][l], tau_l[i2][i1][k], wb)
			   -funI2lk(tau_l[i2][i1][l], tau_l[i2][i1][k], wa);
		}

		tau[i2][i1]=s0/((s1+s2)*Q0[i2][i1]);
	}
}
/*-----------------------------------------------------------------------------*/

int main(int argc, char* argv[])
{
	bool verb, csdgather;
	int ib, is, it, ic, nt, nc, ntc, itc,jsx,jsz,jgx,jgz,sxbeg,szbeg,gxbeg,gzbeg, distx, distz, tdmute;
	int *sxz, *gxz;
	float tmp, fm;
	float *wlt, *bndr;
	float **dcal,**rho, **tau, **tau0, **v0, **vv, **sp, **sr, **svz, **svx, **gp, **gr, **gvz, **gvx, **image;
	float ***cp;
	sf_file Fv, Frho, Ftau, Ftau0, Fw;
	checkpoint *checkpoints;

    	sf_init(argc,argv);
#ifdef _OPENMP
    	omp_init();
#endif

	Fv = sf_input("in");/* veloctiy model */
	Frho=sf_input("rho");/* density */
	Ftau=sf_input("tau");/* tau, computed according to quality factor Q */
	Ftau0=sf_input("tau0");/* tau0, computed according to quality factor Q */
	Fw = sf_output("out");/* image */

    	if(!sf_getbool("verb",&verb)) verb=false;
	/* verbosity */
    	if (!sf_histint(Fv,"n1",&nz)) sf_error("No n1= in input");
	/* veloctiy model: nz */
    	if (!sf_histint(Fv,"n2",&nx)) sf_error("No n2= in input");
	/* veloctiy model: nx */
    	if (!sf_histfloat(Fv,"d1",&dz)) sf_error("No d1= in input");
	/* veloctiy model: dz */
    	if (!sf_histfloat(Fv,"d2",&dx)) sf_error("No d2= in input");
	/* veloctiy model: dx */
    	if (!sf_getint("nb",&nb)) nb=30;
	/* thickness of PML ABC */
    	if (!sf_getint("ns",&ns)) ns=1; 
	/* number of shots */
    	if (!sf_getint("ng",&ng)) ng=1; 
	/* number of receivers per shot */
    	if (!sf_getint("nt",&nt)) sf_error("nt required");
	/* number of time steps */
    	if (!sf_getfloat("dt",&dt)) sf_error("dt required");
	/* time sampling interval */
    	if (!sf_getfloat("fm",&fm)) fm=20.0; 
	/*dominant freq of Ricker wavelet */
	if (!sf_getint("ntc",&ntc)) ntc=20; 
	/* how many time steps every checkpoint */
	if (!sf_getint("jsx",&jsx))   sf_error("no jsx");
	/* source x-axis  jump interval  */
    	if (!sf_getint("jsz",&jsz))   jsz=0;
	/* source z-axis jump interval  */
    	if (!sf_getint("jgx",&jgx))   jgx=1;
	/* receiver x-axis jump interval */
    	if (!sf_getint("jgz",&jgz))   jgz=0;
	/* receiver z-axis jump interval */
    	if (!sf_getint("sxbeg",&sxbeg))   sf_error("no sxbeg");
	/* x-begining index of sources, starting from 0 */
    	if (!sf_getint("szbeg",&szbeg))   sf_error("no szbeg");
	/* z-begining index of sources, starting from 0 */
    	if (!sf_getint("gxbeg",&gxbeg))   sf_error("no gxbeg");
	/* x-begining index of receivers, starting from 0 */
    	if (!sf_getint("gzbeg",&gzbeg))   sf_error("no gzbeg");
	/* z-begining index of receivers, starting from 0 */
	if (!sf_getbool("csdgather",&csdgather)) csdgather=true;
	/* common shot gather or not */
	if (!sf_getfloat("vmute",&vmute))   vmute=1500;
	/* muting velocity to remove the low-freq noise, unit=m/s*/
	if (!sf_getint("tdmute",&tdmute))   tdmute=2./(fm*dt);
	/* number of deleyed time samples to mute */


	_dx=1./dx;
	_dz=1./dz;
	nzpad=nz+2*nb;
	nxpad=nx+2*nb;
	nc=(nt+ntc-1)/ntc;

	/* allocate variables */
	wlt=sf_floatalloc(nt);
	bndr=sf_floatalloc(nb);
	sxz=sf_intalloc(ns);
	gxz=sf_intalloc(ng);
	v0=sf_floatalloc2(nz,nx); 	
	rho=sf_floatalloc2(nzpad, nxpad);
	tau=sf_floatalloc2(nzpad, nxpad);
	tau0=sf_floatalloc2(nzpad, nxpad);
	vv=sf_floatalloc2(nzpad, nxpad);
	sp =sf_floatalloc2(nzpad, nxpad);
	sr =sf_floatalloc2(nzpad, nxpad);
	svz=sf_floatalloc2(nzpad, nxpad);
	svx=sf_floatalloc2(nzpad, nxpad);
	gp =sf_floatalloc2(nzpad, nxpad);
	gr =sf_floatalloc2(nzpad, nxpad);
	gvz=sf_floatalloc2(nzpad, nxpad);
	gvx=sf_floatalloc2(nzpad, nxpad);
	image=sf_floatalloc2(nz, nx);
	dcal=sf_floatalloc2(ng,nt);
	checkpoints=(checkpoint*)malloc(nc*sizeof(checkpoint));
	for(ic=0; ic<nc; ic++){
		checkpoints[ic].p =(float*)malloc(nzpad*nxpad*sizeof(float));
		checkpoints[ic].r =(float*)malloc(nzpad*nxpad*sizeof(float));
		checkpoints[ic].vz=(float*)malloc(nzpad*nxpad*sizeof(float));
		checkpoints[ic].vx=(float*)malloc(nzpad*nxpad*sizeof(float));
	}
	cp=sf_floatalloc3(nzpad, nxpad, ntc);/* ntc snapshots */

	/* initialization */
	for(it=0;it<nt;it++){
		tmp=SF_PI*fm*(it*dt-1.0/fm);tmp*=tmp;
		wlt[it]=(1.0-2.0*tmp)*expf(-tmp);
	}
	for(ib=0;ib<nb;ib++){
		tmp=0.015*(nb-ib);
		bndr[ib]=expf(-tmp*tmp);
	}
	if (!(sxbeg>=0 && szbeg>=0 && sxbeg+(ns-1)*jsx<nx && szbeg+(ns-1)*jsz<nz))	
	{ sf_error("sources exceeds the computing zone!"); exit(1);}
	sg_init(sxz, szbeg, sxbeg, jsz, jsx, ns);
	distx=sxbeg-gxbeg;
	distz=szbeg-gzbeg;
	if (!(gxbeg>=0 && gzbeg>=0 && gxbeg+(ng-1)*jgx<nx && gzbeg+(ng-1)*jgz<nz))	
	{ sf_error("geophones exceeds the computing zone!"); exit(1);}
	if (csdgather)	{
		if (!( (sxbeg+(ns-1)*jsx)+(ng-1)*jgx-distx <nx  && (szbeg+(ns-1)*jsz)+(ng-1)*jgz-distz <nz))	
		{ sf_error("geophones exceeds the computing zone!"); exit(1);}
	}
	sg_init(gxz, gzbeg, gxbeg, jgz, jgx, ng);
	sf_floatread(v0[0],nz*nx,Fv);
	expand2d(vv, v0);
	sf_floatread(v0[0],nz*nx,Frho);
	expand2d(rho, v0);
	sf_floatread(v0[0],nz*nx,Ftau);
	expand2d(tau, v0);
	sf_floatread(v0[0],nz*nx,Ftau0);
	expand2d(tau0, v0);
	memset(image[0],0,nz*nx*sizeof(float));

	for(is=0; is<ns; is++)
	{
		if (csdgather){
			gxbeg=sxbeg+is*jsx-distx;
			sg_init(gxz, gzbeg, gxbeg, jgz, jgx, ng);
		}
		memset(sp [0],0,nzpad*nxpad*sizeof(float));
		memset(sr [0],0,nzpad*nxpad*sizeof(float));
		memset(svx[0],0,nzpad*nxpad*sizeof(float));
		memset(svz[0],0,nzpad*nxpad*sizeof(float));
		ic=0;
		for(it=0; it<nt; it++)/* generate wavefields at checkpoints */
		{
			add_source(&sxz[is], sp, 1, &wlt[it], true);
			step_forward(sp, sr, svz, svx, vv, rho, tau, tau0);
			apply_sponge(sp, bndr);
			apply_sponge(sr, bndr);
			apply_sponge(svx, bndr);
			apply_sponge(svz, bndr);

			record_seis(dcal[it], gxz, sp, ng);
			muting(dcal[it], gzbeg, szbeg, gxbeg, sxbeg+is*jsx, jgx, it, tdmute);

			if(it==ic*ntc){/* record wavefields at checkpoints */
				memcpy(checkpoints[ic].p, sp[0], nzpad*nxpad*sizeof(float));
				memcpy(checkpoints[ic].r, sr[0], nzpad*nxpad*sizeof(float));
				memcpy(checkpoints[ic].vz, svz[0], nzpad*nxpad*sizeof(float));
				memcpy(checkpoints[ic].vx, svx[0], nzpad*nxpad*sizeof(float));
				ic++;
			}
		}

		memset(gp [0],0,nzpad*nxpad*sizeof(float));
		memset(gr [0],0,nzpad*nxpad*sizeof(float));
		memset(gvx[0],0,nzpad*nxpad*sizeof(float));
		memset(gvz[0],0,nzpad*nxpad*sizeof(float));
		ic=0;
		for(it=0; it<nt; it++)
		{
			if(it==ic*ntc){/* re-modeling ntc snapshots from the ic-th checkpoint */
				memcpy(sp[0], checkpoints[ic].p, nzpad*nxpad*sizeof(float));
				memcpy(sr[0], checkpoints[ic].r, nzpad*nxpad*sizeof(float));
				memcpy(svz[0], checkpoints[ic].vz, nzpad*nxpad*sizeof(float));
				memcpy(svx[0], checkpoints[ic].vx, nzpad*nxpad*sizeof(float));
				for(itc=0; itc<SF_MIN(ntc,nt-ic*ntc); itc++)
				{
					add_source(&sxz[is], sp, 1, &wlt[it], true);
					step_forward(sp, sr, svz, svx, vv, rho, tau, tau0);
					apply_sponge(sp, bndr);
					apply_sponge(sr, bndr);
					apply_sponge(svx, bndr);
					apply_sponge(svz, bndr);
					memcpy(cp[itc][0], sp[0], nzpad*nxpad*sizeof(float));
				}
				ic++;
			}

			add_source(gxz, gp, ng, dcal[nt-1-it], true);
			step_forward(gp, gr, gvz, gvx, vv, rho, tau, tau0);
			apply_sponge(gp, bndr);
			apply_sponge(gr, bndr);
			apply_sponge(gvx, bndr);
			apply_sponge(gvz, bndr);
		
			cross_correlation(image, gp, cp[it-ic*ntc]);
		}
	}
	sf_floatwrite(image[0], nz*nx,Fw);/* the image needs laplacian filtering to remove low-freq noise */

	/* free variables */
	free(wlt);
	free(bndr);
	free(sxz);
	free(gxz);
	free(*rho); free(rho);
	free(*tau); free(tau);
	free(*tau0); free(tau0);
	free(*v0); free(v0);
	free(*vv); free(vv);
	free(*sp); free(sp);
	free(*sr); free(sr);
	free(*svx); free(svx);
	free(*svz); free(svz);
	free(*gp); free(gp);
	free(*gr); free(gr);
	free(*gvx); free(gvx);
	free(*gvz); free(gvz);
	free(*image); free(image);
	free(*dcal); free(dcal);
	for(ic=0; ic<nc; ic++){
		free(checkpoints[ic].p);
		free(checkpoints[ic].r);
		free(checkpoints[ic].vz);
		free(checkpoints[ic].vx);
	}
	free(checkpoints);
	free(**cp); free(*cp); free(cp);

    	exit(0);
}

