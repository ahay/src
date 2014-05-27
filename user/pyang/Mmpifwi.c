/* Time domain full waveform inversion using MPI parallel programming 
Note: 	Here, the Enquist absorbing boundary condition is applied!
 */
/*
  Copyright (C) 2014  Xi'an Jiaotong University, UT Austin (Pengliang Yang)

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
#include <mpi.h>

#ifdef _OPENMP
#include <omp.h>
#endif

static bool csdgather;
static int nz,nx,nt,ns,ng;
static float dx, dz, fm, dt;
static int *sxz, *gxz;
static float *wlt, *bndr, *dobs, *dcal, *dres;
static float **vv, **sp0, **sp1, **sp2, **gp0, **gp1, **gp2, **g0, **g1, **cg, **lap, **vtmp;


void matrix_transpose(float *matrix, float *trans, int n1, int n2)
/*< matrix transpose: matrix tansposed to be trans >*/
{
	int i1, i2;

	for(i2=0; i2<n2; i2++)
	for(i1=0; i1<n1; i1++)
	    trans[i2+n2*i1]=matrix[i1+n1*i2];
}


void step_forward(float **p0, float **p1, float **p2, float **vv, float dtz, float dtx)
{
    int ix,iz;
    float v1,v2,diff1,diff2;
    

	for (ix=1; ix < nx-1; ix++) 
	for (iz=1; iz < nz-1; iz++) 
	{
	    	v1=vv[ix][iz]*dtz; v1=v1*v1;
	    	v2=vv[ix][iz]*dtx; v2=v2*v2;
		diff1 =	p1[ix][iz-1]-2.0*p1[ix][iz]+p1[ix][iz+1];
		diff2 =	p1[ix-1][iz]-2.0*p1[ix][iz]+p1[ix+1][iz];
	    	diff1*=v1;
	    	diff2*=v2;		
		p2[ix][iz]=2*p1[ix][iz]-p0[ix][iz]+diff1+diff2;
	}


/*
    for (ix=0; ix < nx; ix++) 
    for (iz=0; iz < nz; iz++) 
    {
	    v1=vv[ix][iz]*dtz; v1=v1*v1;
	    v2=vv[ix][iz]*dtx; v2=v2*v2;
	    diff1=diff2=-2.0*p1[ix][iz];
	    if (iz-1>=0) diff1+=p1[ix][iz-1];
	    if (iz+1<nz) diff1+=p1[ix][iz+1];
	    if (ix-1>=0) diff2+=p1[ix-1][iz];
	    if (ix+1<nx) diff2+=p1[ix+1][iz];
	    diff1*=v1;
	    diff2*=v2;
	    p2[ix][iz]=2.0*p1[ix][iz]-p0[ix][iz]+diff1+diff2;
    }
*/
    for (ix=1; ix < nx-1; ix++) { 
	//top boundary
/*
	if(iz==0){
	    	diff1=	(p1[ix][iz+1]-p1[ix][iz])-
			(p0[ix][iz+1]-p0[ix][iz]);
	    	diff2=	c21*(p1[ix-1][iz]+p1[ix+1][iz]) +
			c22*(p1[ix-2][iz]+p1[ix+2][iz]) +
			c20*p1[ix][iz];
		diff1*=sqrtf(vv[ix][iz])/dz;
 	    	diff2*=vv[ix][iz]/(2.0*dx*dx);
		p2[ix][iz]=2*p1[ix][iz]-p0[ix][iz]+diff1+diff2;
	} 
*/
	/* bottom boundary */
	if(iz==nz-1){
	    	v1=vv[ix][iz]*dtz; 
	    	v2=vv[ix][iz]*dtx;
	    	diff1=-(p1[ix][iz]-p1[ix][iz-1])+(p0[ix][iz]-p0[ix][iz-1]);
	    	diff2=p1[ix-1][iz]-2.0*p1[ix][iz]+p1[ix+1][iz];
		diff1*=v1;
 	    	diff2*=0.5*v2*v2;
		p2[ix][iz]=2.0*p1[ix][iz]-p0[ix][iz]+diff1+diff2;
	}
    }

    for (iz=1; iz <nz-1; iz++){ 
	/* left boundary */
    	if(ix==0){
	    	v1=vv[ix][iz]*dtz; 
	    	v2=vv[ix][iz]*dtx;
	    	diff1=p1[ix][iz-1]-2.0*p1[ix][iz]+p1[ix][iz+1];
	    	diff2=(p1[ix+1][iz]-p1[ix][iz])-(p0[ix+1][iz]-p0[ix][iz]);
 	    	diff1*=0.5*v1*v1;
		diff2*=v2;
		p2[ix][iz]=2.0*p1[ix][iz]-p0[ix][iz]+diff1+diff2;
	}
	/* right boundary */
    	if(ix==nx-1){
	    	v1=vv[ix][iz]*dtz; 
	    	v2=vv[ix][iz]*dtx;
	    	diff1=p1[ix][iz-1]-2.0*p1[ix][iz]+p1[ix][iz+1];
	    	diff2=-(p1[ix][iz]-p1[ix-1][iz])+(p0[ix][iz]-p0[ix-1][iz]);
 	    	diff1*=0.5*v1*v1;
		diff2*=v2;
		p2[ix][iz]=2.0*p1[ix][iz]-p0[ix][iz]+diff1+diff2;
	}
    }  
}

void step_backward(float **lap, float **p0, float **p1, float **p2, float **vv, float dtz, float dtx)
{
    int ix,iz;
    float v1,v2,diff1,diff2;
    

	for (ix=1; ix < nx-1; ix++) 
	for (iz=1; iz < nz-1; iz++) 
	{
	    	v1=vv[ix][iz]*dtz; v1=v1*v1;
	    	v2=vv[ix][iz]*dtx; v2=v2*v2;
		diff1 =	p1[ix][iz-1]-2.0*p1[ix][iz]+p1[ix][iz+1];
		diff2 =	p1[ix-1][iz]-2.0*p1[ix][iz]+p1[ix+1][iz];
	    	diff1*=v1;
	    	diff2*=v2;		
		p2[ix][iz]=2*p1[ix][iz]-p0[ix][iz]+diff1+diff2;
		lap[ix][iz]+=diff1+diff2;
	}
}

void add_source(float **p, float *source, int *sxz, int ns, bool add)
{
	int is, sx, sz;
	if(add){
		for(is=0;is<ns; is++){
			sx=sxz[is]/nz;
			sz=sxz[is]%nz;
			p[sx][sz]+=source[is];
		}
	}else{
		for(is=0;is<ns; is++){
			sx=sxz[is]/nz;
			sz=sxz[is]%nz;
			p[sx][sz]-=source[is];
		}
	}
}

void record_seis(float *seis_it, int *gxz, float **p, int ng)
/* record seismogram at time it into a vector length of ng */
{
	int ig, gx, gz;
	for(ig=0;ig<ng; ig++)
	{
		gx=gxz[ig]/nz;
		gz=gxz[ig]%nz;
		seis_it[ig]=p[gx][gz];
	}
}

void sg_init(int *sxz, int szbeg, int sxbeg, int jsz, int jsx, int ns)
/* shot/geophone position initialize */
{
	int is, sz, sx;
	for(is=0; is<ns; is++)
	{
		sz=szbeg+is*jsz;
		sx=sxbeg+is*jsx;
		sxz[is]=sz+nz*sx;
	}
}

void rw_bndr(float **p, float *bndr, bool read)
/*< if read==true, read boundaries into variables;
 else write/save boundaries (for 2nd order FD) >*/
{
	int i;
	if(read){
		for(i=0; i<nz; i++)
		{
			p[0][i]=bndr[i];
			p[nx-1][i]=bndr[i+nz];
		}
		for(i=0; i<nx; i++) p[i][nz-1]=bndr[i+2*nz];
	}else{
		for(i=0; i<nz; i++)
		{
			bndr[i]=p[0][i];
			bndr[i+nz]=p[nx-1][i];
		}
		for(i=0; i<nx; i++) bndr[i+2*nz]=p[i][nz-1];
	}
}

void cal_residual(float *dobs, float *dcal, float *dres, int ng)
/*< calculate residual >*/
{
  int ig;
  for(ig=0; ig<ng; ig++){
    dres[ig]=dcal[ig]-dobs[ig];
  }
}

void cal_gradient(float **grad, float **lap, float **gp)
{
  int ix, iz;
  for(ix=0; ix<nx; ix++){
    for(iz=0; iz<nz; iz++){
      grad[ix][iz]+=lap[ix][iz]*gp[ix][iz];
    }
  }
}

int main(int argc, char* argv[])
{
	bool verb, precon;
	int is, it, iter, niter, distx, distz, csd, rbell;
	int sxbeg,szbeg,gxbeg,gzbeg,jsx,jsz,jgx,jgz;
	float dtx, dtz,amp, tmp;
	float *trans, **ptr=NULL;
	sf_file vinit, shots, vupdates, grads;

	float tstart, tend, timer;
	int rank, size;
	MPI_Comm comm;
	comm=MPI_COMM_WORLD;

	MPI_Init(&argc,&argv);
	MPI_Comm_size(comm, &size);
	MPI_Comm_rank(comm, &rank);
	tstart=MPI_Wtime();
	if(rank==0) fprintf(stderr, "numprocs=%d\n",size);

    	/* initialize Madagascar */
    	sf_init(argc,argv);

    	/*< set up I/O files >*/
    	vinit=sf_input ("in");   /* initial velocity model, unit=m/s */
	shots=sf_input("shots"); /* recorded shots from exact velocity model */
    	vupdates=sf_output("out"); /* updated velocity in iterations */ 
    	grads=sf_output("grads");  /* gradient in iterations */ 

    	/* get parameters from velocity model and recorded shots */
	if (!sf_getbool("verb",&verb)) verb=true;
    	if (!sf_histint(vinit,"n1",&nz)) sf_error("no n1");
    	if (!sf_histint(vinit,"n2",&nx)) sf_error("no n2");
    	if (!sf_histfloat(vinit,"d1",&dz)) sf_error("no d1");
   	if (!sf_histfloat(vinit,"d2",&dx)) sf_error("no d2");

   	if (!sf_histint(shots,"n1",&nt)) sf_error("no nt");
	/* total modeling time steps */
   	if (!sf_histint(shots,"n2",&ng)) sf_error("no ng");
	/* total receivers in each shot */
   	if (!sf_histint(shots,"n3",&ns)) sf_error("no ns");
	/* number of shots */
   	if (!sf_histfloat(shots,"d1",&dt)) sf_error("no dt");
	/* time sampling interval */
   	if (!sf_histfloat(shots,"amp",&amp)) sf_error("no amp");
	/* maximum amplitude of ricker */
   	if (!sf_histfloat(shots,"fm",&fm)) sf_error("no fm");
	/* dominant freq of ricker */
   	if (!sf_histint(shots,"sxbeg",&sxbeg)) sf_error("no sxbeg");
	/* x-begining index of sources, starting from 0 */
   	if (!sf_histint(shots,"szbeg",&szbeg)) sf_error("no szbeg");
	/* x-begining index of sources, starting from 0 */
   	if (!sf_histint(shots,"gxbeg",&gxbeg)) sf_error("no gxbeg");
	/* x-begining index of receivers, starting from 0 */
   	if (!sf_histint(shots,"gzbeg",&gzbeg)) sf_error("no gzbeg");
	/* x-begining index of receivers, starting from 0 */
   	if (!sf_histint(shots,"jsx",&jsx)) sf_error("no jsx");
	/* source x-axis  jump interval  */
   	if (!sf_histint(shots,"jsz",&jsz)) sf_error("no jsz");
	/* source z-axis jump interval  */
   	if (!sf_histint(shots,"jgx",&jgx)) sf_error("no jgx");
	/* receiver x-axis jump interval  */
   	if (!sf_histint(shots,"jgz",&jgz)) sf_error("no jgz");
	/* receiver z-axis jump interval  */
   	if (!sf_histint(shots,"csdgather",&csd)) sf_error("csdgather or not required");
	/* default, common shot-gather; if n, record at every point*/
    	if (!sf_getint("niter",&niter))   niter=100;
	/* number of iterations */
	if (!sf_getint("rbell",&rbell))	  rbell=2;
	/* radius of bell smooth */
	if (!sf_getbool("precon",&precon)) precon=false;
	/* precondition or not */
	
	if(rank==0)       {
	  sf_putint(vupdates,"n1",nz);	
	  sf_putint(vupdates,"n2",nx);
	  sf_putfloat(vupdates,"d1",dz);
	  sf_putfloat(vupdates,"d2",dx);
	  sf_putstring(vupdates,"label1","Depth");
	  sf_putstring(vupdates,"label2","Distance");
	  sf_putstring(vupdates,"label3","Iteration");
	  sf_putint(vupdates,"n3",niter);
	  sf_putint(vupdates,"d3",1);
	  sf_putint(vupdates,"o3",1);
	  sf_putint(grads,"n1",nz);	
	  sf_putint(grads,"n2",nx);
	  sf_putint(grads,"n3",niter);
	  sf_putfloat(grads,"d1",dz);
	  sf_putfloat(grads,"d2",dx);
	  sf_putint(grads,"d3",1);
	  sf_putint(grads,"o3",1);
	  sf_putstring(grads,"label1","Depth");
	  sf_putstring(grads,"label2","Distance");
	  sf_putstring(grads,"label3","Iteration");
	}
	dtx=dt/dx; 
	dtz=dt/dz; 

	wlt=(float*)malloc(nt*sizeof(float));
	bndr=(float*)malloc(nt*(2*nz+nx)*sizeof(float));
	dcal=(float*)malloc(ng*sizeof(float));
	dobs=(float*)malloc(ng*nt*sizeof(float));
	dres=(float*)malloc(ns*ng*nt*sizeof(float));
	trans=(float*)malloc(ng*nt*sizeof(float));
	vv=sf_floatalloc2(nz, nx);
	sp0=sf_floatalloc2(nz, nx);
	sp1=sf_floatalloc2(nz, nx);
	sp2=sf_floatalloc2(nz, nx);
	gp0=sf_floatalloc2(nz, nx);
	gp1=sf_floatalloc2(nz, nx);
	gp2=sf_floatalloc2(nz, nx);
	lap=sf_floatalloc2(nz, nx);
	g0=sf_floatalloc2(nz, nx);
	g1=sf_floatalloc2(nz, nx);
	cg=sf_floatalloc2(nz, nx);
	vtmp=sf_floatalloc2(nz, nx);
	sxz=(int*)malloc(ns*sizeof(int));
	gxz=(int*)malloc(ng*sizeof(int));

	for(it=0; it<nt; it++){
		tmp=SF_PI*fm*(it*dt-1.0/fm);tmp=tmp*tmp;
		wlt[it]=amp*(1.0-2.0*tmp)*expf(-tmp);
	}
	memset(bndr,0,nt*(2*nz+nx)*sizeof(float));
	memset(dcal,0,ng*sizeof(float));
	memset(dobs,0,ng*nt*sizeof(float));
	memset(dres,0,ns*ng*nt*sizeof(float));
	memset(trans,0,ng*nt*sizeof(float));
	sf_floatread(vv[0],nz*nx,vinit);
	memset(sp0[0],0,nz*nx*sizeof(float));
	memset(sp1[0],0,nz*nx*sizeof(float));
	memset(sp2[0],0,nz*nx*sizeof(float));
	memset(gp0[0],0,nz*nx*sizeof(float));
	memset(gp1[0],0,nz*nx*sizeof(float));
	memset(gp2[0],0,nz*nx*sizeof(float));
	memset(lap[0],0,nz*nx*sizeof(float));
	memset(g0[0],0,nz*nx*sizeof(float));
	memset(g1[0],0,nz*nx*sizeof(float));
	memset(cg[0],0,nz*nx*sizeof(float));
	memset(vtmp[0],0,nz*nx*sizeof(float));
	if (!(sxbeg>=0 && szbeg>=0 && sxbeg+(ns-1)*jsx<nx && szbeg+(ns-1)*jsz<nz))	
	  { fprintf(stderr,"sources exceeds the computing zone!\n");
	    MPI_Finalize(); exit(1);}
 	sg_init(sxz, szbeg, sxbeg, jsz, jsx, ns);
	distx=sxbeg-gxbeg;
	distz=szbeg-gzbeg;
	if (csdgather)	{
	  if (!(gxbeg>=0 && gzbeg>=0 && gxbeg+(ng-1)*jgx<nx && gzbeg+(ng-1)*jgz<nz &&
		(sxbeg+(ns-1)*jsx)+(ng-1)*jgx-distx <nx  && (szbeg+(ns-1)*jsz)+(ng-1)*jgz-distz <nz))	
	    { fprintf(stderr,"geophones exceeds the computing zone!\n");
	      MPI_Finalize(); exit(1);}
	}else{
	  if (!(gxbeg>=0 && gzbeg>=0 && gxbeg+(ng-1)*jgx<nx && gzbeg+(ng-1)*jgz<nz))	
	    { fprintf(stderr,"geophones exceeds the computing zone!\n");
	      MPI_Finalize(); exit(1);}
	}
	sg_init(gxz, gzbeg, gxbeg, jgz, jgx, ng);

	for(iter=0; iter<niter; iter++){
	  memcpy(g0[0], g1[0], nz*nx*sizeof(float));
	  memset(g1[0], 0, nz*nx*sizeof(float));

	  for(is=rank; is<ns; is+=size)	    {
	      sf_seek(shots, is*ng*nt*sizeof(float), SEEK_SET);
	      sf_floatread(dobs, ng*nt, shots);
	      matrix_transpose(dobs, trans, nt, ng);

	      if (csdgather)	{
		gxbeg=sxbeg+is*jsx-distx;
		sg_init(gxz, gzbeg, gxbeg, jgz, jgx, ng);
	      }
	      memset(sp0[0],0,nz*nx*sizeof(float));
	      memset(sp1[0],0,nz*nx*sizeof(float));
	      memset(sp2[0],0,nz*nx*sizeof(float));
	      for(it=0; it<nt; it++)   {
		add_source(sp1, &wlt[it], &sxz[is], 1, true);
		step_forward(sp0, sp1, sp2, vv, dtz, dtx);
		ptr=sp0; sp0=sp1; sp1=sp2; sp2=ptr;
		rw_bndr(sp0, &bndr[it*(2*nz+nx)], false);		

		record_seis(dcal, gxz, sp0, ng);
		cal_residual(&dobs[it*ng], dcal, &dres[is*ng*nt+it*ng], ng);
	      }

	      ptr=sp0; sp0=sp1; sp1=ptr;
	      memset(lap[0],0,nz*nx*sizeof(float));
	      memset(gp0[0],0,nz*nx*sizeof(float));
	      memset(gp1[0],0,nz*nx*sizeof(float));
	      memset(gp2[0],0,nz*nx*sizeof(float));
	      for(it=nt-1; it>-1; it--)   {
		add_source(gp1, &dres[is*ng*nt+it*ng], &gxz[is], ng, true);
		step_forward(gp0, gp1, gp2, vv, dtz, dtx);
		ptr=gp0; gp0=gp1; gp1=gp2; gp2=ptr;

		rw_bndr(sp1, &bndr[it*(2*nz+nx)], true);
		step_backward(lap, sp0, sp1, sp2, vv, dtz, dtx);
		add_source(sp1, &wlt[it], &sxz[is], 1, false);
		ptr=sp0; sp0=sp1; sp1=sp2; sp2=ptr;
		
		cal_gradient(g1, lap, gp0);
	      }	      
	  }



	}


	free(sxz);
	free(gxz);
	free(bndr);
	free(dcal);
 	free(dobs);
	free(dres);
	free(trans);
	free(wlt);
	free(*vv); free(vv);
	free(*sp0); free(sp0);
	free(*sp1); free(sp1);
	free(*sp2); free(sp2);
	free(*gp0); free(gp0);
	free(*gp1); free(gp1);
	free(*gp2); free(gp2);
	free(*lap); free(lap);
	free(*g0); free(g0);
	free(*g1); free(g1);
	free(*cg); free(cg);
	free(*vtmp); free(vtmp);

	tend=MPI_Wtime();
	timer=tend-tstart;
	if(verb) fprintf(stderr,"rank=%d time=%g (s)\n",rank, timer);

	MPI_Finalize();

    	exit(0);
}

