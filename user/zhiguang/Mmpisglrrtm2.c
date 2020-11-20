/* Paralleled stagger-grid lowrank RTM modified based on sfsglfdrtm2 (serial program)
*/
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
#include <mpi.h>
#ifdef _OPENMP
#include <omp.h>
#endif
#include <math.h>
#include <time.h>

#define pml_fvxz {vxn1[ix][iz]=((1-dt*pmldx[ix]/2)*vxn0[ix][iz]-dt/fdenx[ix][iz]*ldx(txxn0, ix, iz))/(1+dt*pmldx[ix]/2);  vzn1[ix][iz]=((1-dt*pmldz[iz]/2)*vzn0[ix][iz]-dt/fdenz[ix][iz]*ldz(txxn0, ix, iz))/(1+dt*pmldz[iz]/2);}
#define pml_bvxz {vxn0[ix][iz]=((1-dt*pmldx[ix]/2)*vxn1[ix][iz]+dt/bdenx[ix][iz]*ldx(txxn0, ix, iz))/(1+dt*pmldx[ix]/2);  vzn0[ix][iz]=((1-dt*pmldz[iz]/2)*vzn1[ix][iz]+dt/bdenz[ix][iz]*ldz(txxn0, ix, iz))/(1+dt*pmldz[iz]/2);}
#define pml_ftxx {txxn1x[ix][iz]=((1-dt*pmldx[ix]/2)*txxn0x[ix][iz]-dt*fc11[ix][iz]*ldx(vxn1, ix-1, iz))/(1+dt*pmldx[ix]/2);  txxn1z[ix][iz]=((1-dt*pmldz[iz]/2)*txxn0z[ix][iz]-dt*fc11[ix][iz]*ldz(vzn1, ix, iz-1))/(1+dt*pmldz[iz]/2);  txxn1[ix][iz]=txxn1x[ix][iz]+txxn1z[ix][iz];}
#define pml_btxx {txxn0x[ix][iz]=((1-dt*pmldx[ix]/2)*txxn1x[ix][iz]+dt*bc11[ix][iz]*ldx(vxn1, ix-1, iz))/(1+dt*pmldx[ix]/2);  txxn0z[ix][iz]=((1-dt*pmldz[iz]/2)*txxn1z[ix][iz]+dt*bc11[ix][iz]*ldz(vzn1, ix, iz-1))/(1+dt*pmldz[iz]/2);  txxn0[ix][iz]=txxn0x[ix][iz]+txxn0z[ix][iz];}

static int nx, nz, nxz, nt, ng, wfnt, wfinv, ntau, nxb, nzb, nxzb;
static int spx, spz, ginv, gp;
static int lenx, lenz, nfd, pmlsize, srcrange;
static int it, is, itau, iturn;
static int *sxx, *sxz, *szx, *szz;

static float dt, wfdt, dtau, tau0;
static float pmld0, srctrunc, srcalpha, tmpfloat;
static float *wavelet;
static float **fden, **fc11, **bden, **bc11;
static float ***Gx, ***Gz;
static float **txxn1x, **txxn1z, **txxn0x, **txxn0z, **transp;
static float *pmldx, *pmldz;
static float **fdenx, **fdenz, **bdenx, **bdenz;
static bool srcdecay;

int sglfdfor2(float ***fwf, float **rcd, bool verb);
int sglfdback2(float ***mig1, float **mig2, float ***fwf, float **localrec, bool verb, bool wantwf, sf_file Ftmpbwf);
void init();
float ldx(float **data, int i, int j);
float ldz(float **data, int i, int j);
void explsource(float **data);
void pml_vxz(float **vxn1, float **vzn1, float **vxn0, float **vzn0, float **txxn0);
void pml_vxzb(float **vxn1, float **vzn1, float **vxn0, float **vzn0, float **txxn0);
void pml_txx(float **txxn1, float **vxn1, float **vzn1);
void pml_txxb(float **txxn0, float **vxn1, float **vzn1);
void zero3(float ***data, int n1, int n2, int n3);
void zero2(float **data, int n1, int n2);

int main(int argc, char *argv[])
{
	clock_t tstart, tend;
	double duration;

	int numprocs, rank;
	float *sendbuf, *recvbuf;
	MPI_Comm Comm=MPI_COMM_WORLD;

	bool verb, wantrecord, wantwf, onlyrecord;
	sf_file Ffvel, Ffden, Fbvel, Fbden;
	sf_file Fsrc, Frcd, Fimg1, Fimg2;
	sf_file FGx, FGz, Fsxx, Fsxz, Fszx, Fszz;
	sf_file Ftmpfwf=NULL, Ftmpbwf=NULL;

	sf_axis at, ax, az, atau;

	int shtbgn, shtinv, shtnmb, shtpad, shtnmb0;
	int snapturn, tmpint;
	int ix, iz;

	float **fvel, **bvel;
	float ***fwf, ***record=NULL, **localrec;
	float ***img1, **img2, ***mig1, **mig2;
	float *tmpsxx, *tmpsxz, *tmpszx, *tmpszz;

	sf_init(argc, argv);

	MPI_Init(&argc, &argv);
	MPI_Comm_size(Comm, &numprocs);
	MPI_Comm_rank(Comm, &rank);

	tstart=clock();
	if(rank==0) sf_warning("numprocs=%d", numprocs);

	if(!sf_getbool("verb", &verb)) verb=true;
	if(!sf_getbool("wantrecord", &wantrecord)) wantrecord=false;
	if(!sf_getbool("wantwf", &wantwf)) wantwf=false;
	if(!sf_getbool("onlyrecord", &onlyrecord)) onlyrecord=false;

	Fsrc=sf_input("-input");
	Fimg1=sf_output("-output");
	Fimg2=sf_output("img2");
	Ffvel=sf_input("fvel");
	Ffden=sf_input("fden");
	Fbvel=sf_input("bvel");
	Fbden=sf_input("bden");

	if(wantrecord)
		Frcd=sf_input("record");
	else
		Frcd=sf_output("record");

	if(wantwf){
		Ftmpfwf=sf_output("tmpfwf");
		Ftmpbwf=sf_output("tmpbwf");
	}

	FGx=sf_input("Gx");
	FGz=sf_input("Gz");
	Fsxx=sf_input("sxx");
	Fsxz=sf_input("sxz");
	Fszx=sf_input("szx");
	Fszz=sf_input("szz");
	
	at=sf_iaxa(Fsrc, 1); nt=sf_n(at); dt=sf_d(at);
	if(!sf_getbool("srcdecay", &srcdecay)) srcdecay=true;
	if(!sf_getint("srcrange", &srcrange)) srcrange=3;
	if(!sf_getfloat("srctrunc", &srctrunc)) srctrunc=0.2;
	if(!sf_getfloat("srcalpha", &srcalpha)) srcalpha=0.5;
	wavelet=sf_floatalloc(nt);
	sf_floatread(wavelet, nt, Fsrc);

	if(!sf_getint("pmlsize", &pmlsize)) pmlsize=30;
	if(!sf_getint("nfd", &nfd)) sf_error("Need half of the FD order!");
	if(!sf_getfloat("pmld0", &pmld0)) pmld0=200;

	if(!sf_getint("shtnmb", &shtnmb)) sf_error("Need shot number!");
	if(!sf_getint("shtinv", &shtinv)) sf_error("Need shot interval!");
	if(!sf_getint("shtbgn", &shtbgn)) shtbgn=0;
	shtpad=numprocs-shtnmb%numprocs;
	shtnmb0=shtnmb+shtpad;

	az=sf_iaxa(Ffvel, 1); nzb=sf_n(az);
	ax=sf_iaxa(Ffvel, 2); nxb=sf_n(ax);
	nxzb=nxb*nzb;
	nz=nzb-2*nfd-2*pmlsize;
	nx=nxb-2*nfd-2*pmlsize;

	if(!sf_getint("snapturn", &snapturn)) snapturn=1;
	if(!sf_getint("ginv", &ginv)) ginv=1;
	if(!sf_getint("wfinv", &wfinv)) wfinv=1;
	if(!sf_getint("spz", &spz)) spz=6;
	if(!sf_getint("gp", &gp)) gp=0;
	ng=(nx-1)/ginv+1;
	wfnt=(nt-1)/wfinv+1;
	wfdt=dt*wfinv;

	if(!sf_getint("ntau", &ntau)) ntau=1;
	if(!sf_getfloat("dtau", &dtau)) dtau=wfdt;
	if(!sf_getfloat("tau0", &tau0)) tau0=0;
	atau=sf_iaxa(Fsrc, 1);
	sf_setn(atau, ntau);
	sf_setd(atau, dtau);
	sf_seto(atau, tau0);

	if(!sf_histint(FGx, "n1", &nxz)) sf_error("No n1= in FGx!");
	if(nxz != nxzb) sf_error("Dimension error!");
	if(!sf_histint(FGx, "n2", &lenx)) sf_error("No n2= in FGx!");
	if(!sf_histint(FGz, "n2", &lenz)) sf_error("No n2= in FGz!");
	Gx=sf_floatalloc3(nzb, nxb, lenx);
	Gz=sf_floatalloc3(nzb, nxb, lenz);
	sxx=sf_intalloc(lenx);
	sxz=sf_intalloc(lenx);
	szx=sf_intalloc(lenz);
	szz=sf_intalloc(lenz);
	tmpsxx=sf_floatalloc(lenx);
	tmpsxz=sf_floatalloc(lenx);
	tmpszx=sf_floatalloc(lenz);
	tmpszz=sf_floatalloc(lenz);
	sf_floatread(Gx[0][0], nxzb*lenx, FGx);
	sf_floatread(Gz[0][0], nxzb*lenz, FGz);
	sf_floatread(tmpsxx, lenx, Fsxx);
	sf_floatread(tmpsxz, lenx, Fsxz);
	sf_floatread(tmpszx, lenz, Fszx);
	sf_floatread(tmpszz, lenz, Fszz);
	for (ix=0; ix<lenx; ix++){
		sxx[ix]=(int)tmpsxx[ix];
		sxz[ix]=(int)tmpsxz[ix];
	}
	for (iz=0; iz<lenz; iz++){
		szx[iz]=(int)tmpszx[iz];
		szz[iz]=(int)tmpszz[iz];
	}

	fvel=sf_floatalloc2(nzb, nxb);
	fden=sf_floatalloc2(nzb, nxb);
	fc11=sf_floatalloc2(nzb, nxb);
	bvel=sf_floatalloc2(nzb, nxb);
	bden=sf_floatalloc2(nzb, nxb);
	bc11=sf_floatalloc2(nzb, nxb);
	sf_floatread(fvel[0], nxzb, Ffvel);
	sf_floatread(fden[0], nxzb, Ffden);
	sf_floatread(bvel[0], nxzb, Fbvel);
	sf_floatread(bden[0], nxzb, Fbden);
	for (ix=0; ix<nxb; ix++){
		for (iz=0; iz<nzb; iz++){
			fc11[ix][iz]=fden[ix][iz]*fvel[ix][iz]*fvel[ix][iz];
			bc11[ix][iz]=bden[ix][iz]*bvel[ix][iz]*bvel[ix][iz];
		}
	}

	if(wantrecord){
		/* check record data */
		sf_histint(Frcd, "n1", &tmpint);
		if(tmpint != nt) sf_error("Not matched dimensions!");
		sf_histint(Frcd, "n2", &tmpint);
		if(tmpint != ng) sf_error("Not matched dimensions!");
		sf_histint(Frcd, "n3", &tmpint);
		if(tmpint != shtnmb) sf_error("Not matched dimensions!");
	}

	if(rank==0){
		record=sf_floatalloc3(nt, ng, shtnmb0);
		if(wantrecord){
			sf_floatread(record[0][0], nt*ng*shtnmb, Frcd);
			for(is=shtnmb; is<shtnmb0; is++)
				for(ix=0; ix<ng; ix++)
					for(it=0; it<nt; it++)
						record[is][ix][it]=0.0;
		}
	}

	img1=sf_floatalloc3(nz, nx, ntau);
	mig1=sf_floatalloc3(nz, nx, ntau);
	img2=sf_floatalloc2(nz, nx);
	mig2=sf_floatalloc2(nz, nx);
	zero3(img1, nz, nx, ntau);
	zero2(img2, nz, nx);

	sf_setn(az, nz);
	sf_setn(ax, ng);
	if(!wantrecord){
		sf_oaxa(Frcd, at, 1);
		sf_oaxa(Frcd, ax, 2);
		sf_putint(Frcd, "n3", shtnmb);
		sf_putint(Frcd, "d3", shtinv);
		sf_putint(Frcd, "o3", shtbgn);
	}

	sf_setn(ax, nx);
	if(wantwf){
		sf_setn(at, wfnt);
		sf_setd(at, wfdt);

		sf_oaxa(Ftmpfwf, az, 1);
		sf_oaxa(Ftmpfwf, ax, 2);
		sf_oaxa(Ftmpfwf, at, 3);

		sf_oaxa(Ftmpbwf, az, 1);
		sf_oaxa(Ftmpbwf, ax, 2);
		sf_oaxa(Ftmpbwf, at, 3);
	}

	sf_oaxa(Fimg1, az, 1);
	sf_oaxa(Fimg1, ax, 2);
	sf_oaxa(Fimg1, atau, 3);
	sf_oaxa(Fimg2, az, 1);
	sf_oaxa(Fimg2, ax, 2);

	fwf=sf_floatalloc3(nz, nx, wfnt);
	localrec=sf_floatalloc2(nt, ng);

	if(verb){
		sf_warning("==================================");
		sf_warning("nx=%d nz=%d nt=%d", nx, nz, nt);
		sf_warning("wfnt=%d wfdt=%f wfinv=%d dt=%f", wfnt, wfdt, wfinv, dt);
		sf_warning("nxb=%d nzb=%d pmlsize=%d nfd=%d", nxb, nzb, pmlsize, nfd);
		sf_warning("ntau=%d dtau=%f tau0=%f", ntau, dtau, tau0);
		sf_warning("shtnmb=%d shtbgn=%d shtinv=%d", shtnmb, shtbgn, shtinv);
		sf_warning("lenx=%d lenz=%d spz=%d gp=%d", lenx, lenz, spz, gp);
		sf_warning("==================================");
	}

	init();

	for(iturn=0; iturn*numprocs<shtnmb; iturn++){
		is=iturn*numprocs+rank;
		if(is<shtnmb){
			sf_warning("ishot/nshot: %d/%d", is+1, shtnmb);
			spx=is*shtinv+shtbgn;
			sglfdfor2(fwf, localrec, verb);
		}

		if(wantrecord){
			recvbuf=localrec[0];
			if(rank==0) sendbuf=record[iturn*numprocs][0];
			else sendbuf=NULL;
			MPI_Scatter(sendbuf, ng*nt, MPI_FLOAT, recvbuf, ng*nt, MPI_FLOAT, 0, Comm);
		}else{
			sendbuf=localrec[0];
			if(rank==0) recvbuf=record[iturn*numprocs][0];
			else recvbuf=NULL;
			MPI_Gather(sendbuf, ng*nt, MPI_FLOAT, recvbuf, ng*nt, MPI_FLOAT, 0, Comm);
		}

		if(wantwf && rank==0 && iturn==snapturn-1) wantwf=true;
		else wantwf=false;
		if(wantwf) sf_floatwrite(fwf[0][0], wfnt*nx*nz, Ftmpfwf);

		if(!onlyrecord && is<shtnmb){
			sglfdback2(mig1, mig2, fwf, localrec, verb, wantwf, Ftmpbwf);
			for(itau=0; itau<ntau; itau++){
				for(ix=0; ix<nx; ix++){
					for(iz=0; iz<nz; iz++){
						img1[itau][ix][iz]+=mig1[itau][ix][iz];
					}
				}
			}
			for(ix=0; ix<nx; ix++){
				for(iz=0; iz<nz; iz++){
					img2[ix][iz]+=mig2[ix][iz];
				}
			}
		}
		MPI_Barrier(Comm);
	} //end of iturn

	if(!onlyrecord){
	if(rank==0){
		sendbuf=(float *)MPI_IN_PLACE;
		recvbuf=img1[0][0];
	}else{
		sendbuf=img1[0][0];
		recvbuf=NULL;
	}
	MPI_Reduce(sendbuf, recvbuf, ntau*nx*nz, MPI_FLOAT, MPI_SUM, 0, Comm);

	if(rank==0){
		sendbuf=MPI_IN_PLACE;
		recvbuf=img2[0];
	}else{
		sendbuf=img2[0];
		recvbuf=NULL;
	}
	MPI_Reduce(sendbuf, recvbuf, nx*nz, MPI_FLOAT, MPI_SUM, 0, Comm);
	}

	if(rank==0){
		if(!wantrecord){
			sf_floatwrite(record[0][0], shtnmb*ng*nt, Frcd);
		}
		sf_floatwrite(img1[0][0], ntau*nx*nz, Fimg1);
		sf_floatwrite(img2[0], nx*nz, Fimg2);
	}

	tend=clock();
	duration=(double)(tend-tstart)/CLOCKS_PER_SEC;
	sf_warning(">>The CPU time of sfmpilfdrtm2 is: %f seconds<<", duration);
	MPI_Finalize();
	exit(0);
}

int sglfdfor2(float ***fwf, float **rcd, bool verb)
{
	float **txxn1, **txxn0, **vxn1, **vxn0, **vzn1, **vzn0;
	int wfit;
	int ix, iz;

	txxn1=sf_floatalloc2(nzb, nxb);
	txxn0=sf_floatalloc2(nzb, nxb);
	vxn1=sf_floatalloc2(nzb, nxb);
	vxn0=sf_floatalloc2(nzb, nxb);
	vzn1=sf_floatalloc2(nzb, nxb);
	vzn0=sf_floatalloc2(nzb, nxb);

	zero2(txxn1, nzb, nxb);
	zero2(txxn0, nzb, nxb);
	zero2(vxn1, nzb, nxb);
	zero2(vxn0, nzb, nxb);
	zero2(vzn1, nzb, nxb);
	zero2(vzn0, nzb, nxb);

	zero2(txxn1x, nzb, nxb);
	zero2(txxn1z, nzb, nxb);
	zero2(txxn0x, nzb, nxb);
	zero2(txxn0z, nzb, nxb);

	wfit=0;
	for(it=0; it<nt; it++){
//		sf_warning("test txxn1[801][30]=%d",txxn1[801][30])
		if(verb) sf_warning("Forward it=%d/%d;", it+1, nt);
#ifdef _OPENMP
#pragma omp parallel for private(ix,iz)
#endif
		for(ix=nfd+pmlsize; ix<nfd+pmlsize+nx; ix++){
			for(iz=nfd+pmlsize; iz<nfd+pmlsize+nz; iz++){
				vxn1[ix][iz]=vxn0[ix][iz]-dt/fdenx[ix][iz]*ldx(txxn0, ix, iz);
				vzn1[ix][iz]=vzn0[ix][iz]-dt/fdenz[ix][iz]*ldz(txxn0, ix, iz);
			}
		}

		pml_vxz(vxn1, vzn1, vxn0, vzn0, txxn0);
#ifdef _OPENMP
#pragma omp parallel for private(ix,iz)
#endif
		for(ix=nfd+pmlsize; ix<nfd+pmlsize+nx; ix++){
			for(iz=nfd+pmlsize; iz<nfd+pmlsize+nz; iz++){
				txxn1[ix][iz]=txxn0[ix][iz]-dt*fc11[ix][iz]*(ldx(vxn1, ix-1, iz) + ldz(vzn1, ix, iz-1));
			}
		}

		pml_txx(txxn1, vxn1, vzn1);

		if((it*dt)<srctrunc){
			explsource(txxn1);
		}

		if(it%wfinv==0){
#ifdef _OPENMP
#pragma omp parallel for private(ix,iz)
#endif
			for(ix=0; ix<nx; ix++){
				for(iz=0; iz<nz; iz++){
					fwf[wfit][ix][iz]=txxn0[ix+nfd+pmlsize][iz+nfd+pmlsize];
				}
			}
			wfit++;
		}
#ifdef _OPENMP
#pragma omp parallel for private(ix)
#endif
		for(ix=0; ix<ng; ix++){
			rcd[ix][it]=txxn0[ix*ginv+pmlsize+nfd][pmlsize+nfd+gp];
		}

		transp=txxn0; txxn0=txxn1; txxn1=transp;
		transp=vxn0;  vxn0=vxn1;   vxn1=transp;
		transp=vzn0;  vzn0=vzn1;   vzn1=transp;
	} // end of it
	if(verb) sf_warning(".");
	return 0;;
}

int sglfdback2(float ***mig1, float **mig2, float ***fwf, float **localrec, bool verb, bool wantwf, sf_file Ftmpbwf)
{
	float **txxn1, **txxn0, **vxn1, **vxn0, **vzn1, **vzn0;
	float **sill, **ccr, ***bwf;
	int wfit, htau;
	int ix, iz;
	float tau;

	sill=sf_floatalloc2(nz, nx);
	ccr=sf_floatalloc2(nz, nx);
	bwf=sf_floatalloc3(nz, nx, wfnt);
	zero2(sill, nz, nx);
	zero2(ccr, nz, nx);
	zero3(mig1, nz, nx, ntau);

	txxn1=sf_floatalloc2(nzb, nxb);
	txxn0=sf_floatalloc2(nzb, nxb);
	vxn1=sf_floatalloc2(nzb, nxb);
	vxn0=sf_floatalloc2(nzb, nxb);
	vzn1=sf_floatalloc2(nzb, nxb);
	vzn0=sf_floatalloc2(nzb, nxb);

	zero2(txxn1, nzb, nxb);
	zero2(txxn0, nzb, nxb);
	zero2(vxn1, nzb, nxb);
	zero2(vxn0, nzb, nxb);
	zero2(vzn1, nzb, nxb);
	zero2(vzn0, nzb, nxb);

	zero2(txxn1x, nzb, nxb);
	zero2(txxn1z, nzb, nxb);
	zero2(txxn0x, nzb, nxb);
	zero2(txxn0z, nzb, nxb);
	
	wfit=wfnt-1;
	for(it=nt-1; it>=0; it--){
		if(verb) sf_warning("Backward it=%d/%d;", it+1, nt);
#ifdef _OPENMP
#pragma omp parallel for private(ix,iz)
#endif
		for(ix=nfd+pmlsize; ix<nfd+pmlsize+nx; ix++){
			for(iz=nfd+pmlsize; iz<nfd+pmlsize+nz; iz++){
				txxn0[ix][iz]=txxn1[ix][iz]+dt*bc11[ix][iz]*(ldx(vxn1, ix-1, iz) +ldz(vzn1, ix, iz-1));
			}
		}

		pml_txxb(txxn0, vxn1, vzn1);
#ifdef _OPENMP
#pragma omp parallel for private(ix)
#endif
		for(ix=0; ix<ng; ix++){
			txxn0[ix*ginv+pmlsize+nfd][pmlsize+nfd+gp]+=localrec[ix][it];
		}
#ifdef _OPENMP
#pragma omp parallel for private(ix,iz)
#endif
		for(ix=nfd+pmlsize; ix<nfd+pmlsize+nx; ix++){
			for(iz=nfd+pmlsize; iz<nfd+pmlsize+nz; iz++){
				vxn0[ix][iz]=vxn1[ix][iz]+dt/bdenx[ix][iz]*ldx(txxn0, ix, iz);
				vzn0[ix][iz]=vzn1[ix][iz]+dt/bdenz[ix][iz]*ldz(txxn0, ix, iz);
			}
		}

		pml_vxzb(vxn1, vzn1, vxn0, vzn0, txxn0);

		transp=txxn1; txxn1=txxn0; txxn0=transp;
		transp=vxn1; vxn1=vxn0; vxn0=transp;
		transp=vzn1; vzn1=vzn0; vzn0=transp;

		if(it%wfinv==0){
			for(ix=0; ix<nx; ix++)
				for(iz=0; iz<nz; iz++){
					bwf[wfit][ix][iz]=txxn0[ix+pmlsize+nfd][iz+pmlsize+nfd];
					ccr[ix][iz]+=fwf[wfit][ix][iz]*bwf[wfit][ix][iz];
					sill[ix][iz]+=fwf[wfit][ix][iz]*fwf[wfit][ix][iz];
				}
			wfit--;
		}
	} //end of it
	if(verb) sf_warning(".");

	for(itau=0; itau<ntau; itau++){
		tau=itau*dtau+tau0;
		htau=tau/wfdt;
		for(it=abs(htau); it<wfnt-abs(htau); it++){
			for(ix=0; ix<nx; ix++){
				for(iz=0; iz<nz; iz++){
					mig1[itau][ix][iz]+=fwf[it+htau][ix][iz]*bwf[it-htau][ix][iz];
				}
			}
		}//end of it
	} // end of itau

	for(ix=0; ix<nx; ix++){
		for(iz=0; iz<nz; iz++){
			mig2[ix][iz]=ccr[ix][iz]/(sill[ix][iz]+SF_EPS);
		}
	}

	if(wantwf) sf_floatwrite(bwf[0][0], wfnt*nx*nz, Ftmpbwf);
	return 0;
}


void init()
{
	int ix, iz;
	fdenx = sf_floatalloc2(nzb, nxb);
	fdenz=sf_floatalloc2(nzb, nxb);

	for(ix=0; ix<nxb; ix++){
		for(iz=0; iz<nzb; iz++){
			fdenx[ix][iz]=fden[ix][iz];
			fdenz[ix][iz]=fden[ix][iz];
		}
	}

	for(ix=0; ix<nxb-1; ix++){
		for(iz=0; iz<nzb; iz++){
			fdenx[ix][iz]=(fden[ix][iz]+fden[ix+1][iz])*0.5;
		}
	}

	for(ix=0; ix<nxb; ix++){
		for (iz=0; iz<nzb-1; iz++){
			fdenz[ix][iz]=(fden[ix][iz]+fden[ix][iz+1])*0.5;
		}
	}

	bdenx=sf_floatalloc2(nzb, nxb);
	bdenz=sf_floatalloc2(nzb, nxb);

	for(ix=0; ix<nxb; ix++){
		for(iz=0; iz<nzb; iz++){
			bdenx[ix][iz]=bden[ix][iz];
			bdenz[ix][iz]=bden[ix][iz];
		}
	}

	for(ix=0; ix<nxb-1; ix++){
		for(iz=0; iz<nzb; iz++){
			bdenx[ix][iz]=(bden[ix][iz]+bden[ix+1][iz])*0.5;
		}
	}

	for(ix=0; ix<nxb; ix++){
		for (iz=0; iz<nzb-1; iz++){
			bdenz[ix][iz]=(bden[ix][iz]+bden[ix][iz+1])*0.5;
		}
	}

	pmldx=sf_floatalloc(nxb);
	pmldz=sf_floatalloc(nzb);

	txxn1x=sf_floatalloc2(nzb, nxb);
	txxn1z=sf_floatalloc2(nzb, nxb);
	txxn0x=sf_floatalloc2(nzb, nxb);
	txxn0z=sf_floatalloc2(nzb, nxb);

	for(ix=0; ix<nxb; ix++){
		pmldx[ix]=0.0;
	}
	for(iz=0; iz<nzb; iz++){
		pmldz[iz]=0.0;
	}

	for(ix=nfd; ix<nfd+pmlsize; ix++){
		tmpfloat=(float)(nfd+pmlsize-ix)/pmlsize;
		pmldx[ix]=pmld0*tmpfloat*tmpfloat;
	}
	for(ix=nfd+pmlsize+nx; ix<nfd+nx+2*pmlsize; ix++){
		tmpfloat=(float)(ix-nfd-pmlsize-nx+1)/pmlsize;
		pmldx[ix]=pmld0*tmpfloat*tmpfloat;
	}
	for(iz=nfd; iz<nfd+pmlsize; iz++){
		tmpfloat=(float)(nfd+pmlsize-iz)/pmlsize;
		pmldz[iz]=pmld0*tmpfloat*tmpfloat;
	}
	for(iz=nz+nfd+pmlsize; iz<nz+nfd+2*pmlsize; iz++){
		tmpfloat=(float)(iz-nz-nfd-pmlsize+1)/pmlsize;
		pmldz[iz]=pmld0*tmpfloat*tmpfloat;
	}
}

float ldx(float **data, int i, int j)
{
	float res=0.0;
	int il;
	for(il=0; il<lenx; il++){
		res += 0.5*(data[i+sxx[il]][j+sxz[il]] - data[i-sxx[il]+1][j-sxz[il]])*Gx[il][i][j];
	}
	return res;
}

float ldz(float **data, int i, int j)
{
	float res=0.0;
	int il;
	
	for(il=0; il<lenz; il++){
		res += 0.5*(data[i+szx[il]][j+szz[il]] - data[i-szx[il]][j-szz[il]+1])*Gz[il][i][j];
	}
	return res;
}

void explsource(float **data)
{
	int ix, iz;
	float phi;
	if(srcdecay){
		for(ix=-srcrange; ix<=srcrange; ix++){
			for(iz=-srcrange; iz<=srcrange; iz++){
				phi=exp(-srcalpha*srcalpha*(ix*ix+iz*iz));
				data[spx+ix+pmlsize+nfd][spz+iz+pmlsize+nfd]+=wavelet[it]*phi;
			}
		}
	}else{
		data[spx+pmlsize+nfd][spz+pmlsize+nfd]+=wavelet[it];
	}
}


void pml_vxz(float **vxn1, float **vzn1, float **vxn0, float **vzn0, float **txxn0)
{
	int ix, iz;
#ifdef _OPENMP
#pragma omp parallel for private(ix,iz)	
#endif
	for(ix=nfd; ix<nfd+2*pmlsize+nx; ix++){
		for(iz=nfd; iz<nfd+pmlsize; iz++){
			pml_fvxz
		}
	}
#ifdef _OPENMP
#pragma omp parallel for private(ix,iz)	
#endif
	for(ix=nfd; ix<nfd+pmlsize; ix++){
		for(iz=nfd+pmlsize; iz<nfd+pmlsize+nz; iz++){
			pml_fvxz
		}
	}
#ifdef _OPENMP
#pragma omp parallel for private(ix,iz)	
#endif
	for(ix=nfd+pmlsize+nx; ix<nfd+2*pmlsize+nx; ix++){
		for(iz=nfd+pmlsize; iz<nfd+pmlsize+nz; iz++){
			pml_fvxz
		}
	}
#ifdef _OPENMP
#pragma omp parallel for private(ix,iz)	
#endif
	for(ix=nfd; ix<nfd+2*pmlsize+nx; ix++){
		for(iz=nfd+pmlsize+nz; iz<nfd+2*pmlsize+nz; iz++){
			pml_fvxz
		}
	}
}

void pml_vxzb(float **vxn1, float **vzn1, float **vxn0, float **vzn0, float **txxn0)
{
	int ix, iz;
#ifdef _OPENMP
#pragma omp parallel for private(ix,iz)	
#endif
	for(ix=nfd; ix<nfd+2*pmlsize+nx; ix++){
		for(iz=nfd; iz<nfd+pmlsize; iz++){
			pml_bvxz
		}
	}
#ifdef _OPENMP
#pragma omp parallel for private(ix,iz)	
#endif
	for(ix=nfd; ix<nfd+pmlsize; ix++){
		for(iz=nfd+pmlsize; iz<nfd+pmlsize+nz; iz++){
			pml_bvxz
		}
	}
#ifdef _OPENMP
#pragma omp parallel for private(ix,iz)	
#endif
	for(ix=nfd+pmlsize+nx; ix<nfd+2*pmlsize+nx; ix++){
		for(iz=nfd+pmlsize; iz<nfd+pmlsize+nz; iz++){
			pml_bvxz
		}
	}
#ifdef _OPENMP
#pragma omp parallel for private(ix,iz)	
#endif
	for(ix=nfd; ix<nfd+2*pmlsize+nx; ix++){
		for(iz=nfd+pmlsize+nz; iz<nfd+2*pmlsize+nz; iz++){
			pml_bvxz
		}
	}
}

void pml_txx(float **txxn1, float **vxn1, float **vzn1)
{
	int ix, iz;
#ifdef _OPENMP
#pragma omp parallel for private(ix,iz)	
#endif
	for(ix=nfd; ix<nfd+2*pmlsize+nx; ix++){
		for(iz=nfd; iz<nfd+pmlsize; iz++){
			pml_ftxx
		}
	}

#pragma omp parallel for private(ix,iz)	
	for(ix=nfd; ix<nfd+pmlsize; ix++){
		for(iz=nfd+pmlsize; iz<nfd+pmlsize+nz; iz++){
			pml_ftxx
		}
	}
#ifdef _OPENMP
#pragma omp parallel for private(ix,iz)	
#endif
	for(ix=nfd+pmlsize+nx; ix<nfd+2*pmlsize+nx; ix++){
		for(iz=nfd+pmlsize; iz<nfd+pmlsize+nz; iz++){
			pml_ftxx
		}
	}
#ifdef _OPENMP
#pragma omp parallel for private(ix,iz)	
#endif
	for(ix=nfd; ix<nfd+2*pmlsize+nx; ix++){
		for(iz=nfd+pmlsize+nz; iz<nfd+2*pmlsize+nz; iz++){
			pml_ftxx
		}
	}

	transp=txxn0x; txxn0x=txxn1x; txxn1x=transp;
	transp=txxn0z; txxn0z=txxn1z; txxn1z=transp;
}

void pml_txxb(float **txxn0, float **vxn1, float **vzn1)
{
	int ix, iz;
#ifdef _OPENMP
#pragma omp parallel for private(ix,iz)	
#endif
	for(ix=nfd; ix<nfd+2*pmlsize+nx; ix++){
		for(iz=nfd; iz<nfd+pmlsize; iz++){
			pml_btxx
		}
	}
#ifdef _OPENMP
#pragma omp parallel for private(ix,iz)	
#endif
	for(ix=nfd; ix<nfd+pmlsize; ix++){
		for(iz=nfd+pmlsize; iz<nfd+pmlsize+nz; iz++){
			pml_btxx
		}
	}
#ifdef _OPENMP
#pragma omp parallel for private(ix,iz)	
#endif
	for(ix=nfd+pmlsize+nx; ix<nfd+2*pmlsize+nx; ix++){
		for(iz=nfd+pmlsize; iz<nfd+pmlsize+nz; iz++){
			pml_btxx
		}
	}
#ifdef _OPENMP
#pragma omp parallel for private(ix,iz)	
#endif
	for(ix=nfd; ix<nfd+2*pmlsize+nx; ix++){
		for(iz=nfd+pmlsize+nz; iz<nfd+2*pmlsize+nz; iz++){
			pml_btxx
		}
	}

	transp=txxn1x; txxn1x=txxn0x; txxn0x=transp;
	transp=txxn1z; txxn1z=txxn0z; txxn0z=transp;
}

void zero3(float ***data, int n1, int n2, int n3)
{
	int i,j,k;
	for(i=0; i<n3; i++){
		for(j=0; j<n2; j++){
			for(k=0; k<n1; k++){
				data[i][j][k]=0.0;
			}
		}
	}
}

void zero2(float **data, int n1, int n2)
{
	int i,j;
	for(i=0; i<n2; i++){
		for(j=0; j<n1; j++){
			data[i][j]=0.0;
		}
	}
}
