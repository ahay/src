/*  A k-space staggered-grid lowrank finite-difference for elastic and viscoelastic seismic-wave modeling
*/
/*
  Copyright (C)2014 Institue of Geology and Geophysics, Chinese Academy of Sciences (Jun Yan) 
					2009 University of Texas at Austin
  
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


/* low rank fd scheme */
static int *Sxx, *Sxz, *Szx, *Szz;

static float lrdx(float **data, int ix, int iz, int lenx, float ***gx);
/*<Low rank finite difference : d/dx>*/

static float lrdz(float **data, int ix, int iz, int lenz, float ***gz);
/*< Low rank finite difference : d/dz >*/


float **extmodel2d(float **init_model,int nz,int nx,int np);
/*< extend 2d model >*/

void explsourcet(float **txx/*@out@*/,
				 float fpeak,
		 		 int it, float dt,
                 int sx, int sz 
				);
/*< explosive source for stress txx >*/ 

int main(int argc, char **argv)
{
	int nx, nz, nt, nxb, nzb, npml, ix, iz, it, verb;
	int isrcx, isrcz,igz, snap, nsnap, ncourant, order, Nhalf, freesurface;
	float dx, dz, dt, ddx, ddxp, ddp, ddz, ddzp, dmin, dmax, o1, o2;
	float fpeak;
	float pmlthick, rcoef, dp0, dpleft, dpright, dptmp, dpcoe;

	float **vp_1, **vp_2, **vx_1, **vx_2, **vz_1, **vz_2;
	float **tpxx_1, **tpxx_2, **tpzz_1, **tpzz_2, **txx, **tzz, **txz_1, **txz_2;
	float *dpx1, *dpx2, *dpz1, *dpz2;

	/* input data */
	float **vp, **vs, **rho, **vppml, **vspml, **rhopml, **rhopmlx, **rhopmlz, **lambda, **mu, **lambda2mu, muxy, muxz, muyz,  vpmax, vpmin;
	float **datx, **datz, **datpx, **datpz, **datsx, **datsz, **data;

	/* low rank finite scheme */
	float ***Gpx, ***Gpz, ***Gsx, ***Gsz;
	float *Sxxtmp, *Sxztmp, *Szxtmp, *Szztmp;

	sf_file Gpxfp, Gpzfp, Gsxfp, Gszfp,sxxfp, sxzfp, szxfp, szzfp;
	/* output data and file */
	sf_file snapfpx, snapfpz, vpfp, vsfp, rhofp, datfp,snapfppx, snapfppz, snapfpsx, snapfpsz;
	/* check stability condition */
	double S2=1., S4=0.857142855, S6=0.805369132, S8=0.777417902, S10=0.76673391, S12=0.746790483;

    sf_init(argc,argv);

	/* vp */
	vpfp = sf_input("in");   
    vsfp = sf_input("vs");  
	rhofp = sf_input("rho"); 
	/* lowrank decompositon matrix */
	Gpxfp = sf_input("Gpx"); 
	Gpzfp = sf_input("Gpz");
	Gsxfp = sf_input("Gsx");
	Gszfp = sf_input("Gsz");
	/* FD stencil */
	sxxfp = sf_input("sxx"); 
	sxzfp = sf_input("sxz");
	szxfp = sf_input("szx");
	szzfp = sf_input("szz");
	
	/* output wave shotsnaps */
	snapfpx = sf_output("out");
	snapfpz = sf_output("snapz");
	snapfppx = sf_output("snappx");
	snapfppz = sf_output("snappz");
	snapfpsx = sf_output("snapsx");
	snapfpsz = sf_output("snapsz");
	
	if (SF_FLOAT != sf_gettype(vpfp)) sf_error("Need float input");
	if (!sf_histint(vpfp,"n1",&nz)) sf_error("No n1= in input");
	if (!sf_histfloat(vpfp,"d1",&dz)) sf_error("No d1= in input");
    if (!sf_histint(vpfp,"n2",&nx)) sf_error("No n2= in input");
    if (!sf_histfloat(vpfp,"d2",&dx)) sf_error("No d2= in input");
    if (!sf_histfloat(vpfp,"o1",&o1)) o1=0.0;
    if (!sf_histfloat(vpfp,"o2",&o2)) o2=0.0;

	/* fd scheme */   
	if (!sf_getint("order",&order)) order=12;
	/* source */
    if (!sf_getfloat("dt",&dt)) sf_error("Need dt input");
    if (!sf_getint("nt",&nt)) sf_error("Need nt input");
    if (!sf_getfloat("fpeak",&fpeak)) fpeak=20.0;
    if (!sf_getint("isrcx",&isrcx)) sf_error("Need isrcx ");
    if (!sf_getint("isrcz",&isrcz)) sf_error("Need isrcz ");

    if (!sf_getint("npml",&npml)) npml=30;
    if (!sf_getint("verb",&verb)) verb=1;

	/* freesurface */
    if (!sf_getint("freesurface",&freesurface)) freesurface=0;
	/* recerver */
    if (!sf_getint("igz",&igz)) igz=1;
	if (!sf_getint("snap",&snap)) snap=1;

	nsnap=0;
	for (it=0; it < nt; it++) {
			if (it%snap == 0) nsnap++;
	 }

	/* extended grid size */
	npml = 30;
	nxb = nx + 2*npml;
	nzb = nz + 2*npml;
	
	/* source */
	isrcx +=  npml;
	isrcz +=  npml;
	
	/* recerver */
	igz += npml;

	/* FD order */
	Nhalf = order/2;
    
	sf_putfloat(snapfpx,"d1",dz);
    sf_putfloat(snapfpx,"d2",dx);
    sf_putfloat(snapfpx,"d3",dt*snap);
    sf_putfloat(snapfpx,"o1",o1); 
    sf_putfloat(snapfpx,"o2",o2); 
    sf_putfloat(snapfpx,"o3",0.0); 

	sf_putint(snapfpx,"n1",nzb);
    sf_putint(snapfpx,"n2",nxb);
    sf_putint(snapfpx,"n3",nsnap);

    sf_putfloat(snapfpz,"d1",dz);
    sf_putfloat(snapfpz,"d2",dx);
    sf_putfloat(snapfpz,"d3",dt*snap);
    sf_putfloat(snapfpz,"o1",o1); 
    sf_putfloat(snapfpz,"o2",o2); 
    sf_putfloat(snapfpz,"o3",0.0); 

	sf_putint(snapfpz,"n1",nzb);
    sf_putint(snapfpz,"n2",nxb);
    sf_putint(snapfpz,"n3",nsnap);


    sf_putfloat(snapfppx,"d1",dz);
    sf_putfloat(snapfppx,"d2",dx);
    sf_putfloat(snapfppx,"d3",dt*snap);
    sf_putfloat(snapfppx,"o1",o1); 
    sf_putfloat(snapfppx,"o2",o2); 
    sf_putfloat(snapfppx,"o3",0.0); 

	sf_putint(snapfppx,"n1",nzb);
    sf_putint(snapfppx,"n2",nxb);
    sf_putint(snapfppx,"n3",nsnap);


    sf_putfloat(snapfppz,"d1",dz);
    sf_putfloat(snapfppz,"d2",dx);
    sf_putfloat(snapfppz,"d3",dt*snap);
    sf_putfloat(snapfppz,"o1",o1); 
    sf_putfloat(snapfppz,"o2",o2); 
    sf_putfloat(snapfppz,"o3",0.0); 

	sf_putint(snapfppz,"n1",nzb);
    sf_putint(snapfppz,"n2",nxb);
    sf_putint(snapfppz,"n3",nsnap);


    sf_putfloat(snapfpsx,"d1",dz);
    sf_putfloat(snapfpsx,"d2",dx);
    sf_putfloat(snapfpsx,"d3",dt*snap);
    sf_putfloat(snapfpsx,"o1",o1); 
    sf_putfloat(snapfpsx,"o2",o2); 
    sf_putfloat(snapfpsx,"o3",0.0); 

	sf_putint(snapfpsx,"n1",nzb);
    sf_putint(snapfpsx,"n2",nxb);
    sf_putint(snapfpsx,"n3",nsnap);


    sf_putfloat(snapfpsz,"d1",dz);
    sf_putfloat(snapfpsz,"d2",dx);
    sf_putfloat(snapfpsz,"d3",dt*snap);
    sf_putfloat(snapfpsz,"o1",o1); 
    sf_putfloat(snapfpsz,"o2",o2); 
    sf_putfloat(snapfpsz,"o3",0.0); 

	sf_putint(snapfpsz,"n1",nzb);
    sf_putint(snapfpsz,"n2",nxb);
    sf_putint(snapfpsz,"n3",nsnap);

	/* rho model */
	rho = sf_floatalloc2(nz, nx);
	sf_floatread(rho[0],nz*nx,rhofp);
	rhopml = extmodel2d(rho,nz,nx,npml);

	rhopmlx = sf_floatalloc2(nzb, nxb);
	rhopmlz = sf_floatalloc2(nzb, nxb);

	for (ix=0; ix < nxb-1; ix++) {
		for (iz=0; iz < nzb; iz++) {
			rhopmlx[ix][iz] = rhopml[ix][iz] + rhopml[ix+1][iz];
		}
	}
	for (iz=0; iz < nzb; iz++) {
		rhopmlx[nxb-1][iz] = rhopml[nxb-1][iz];
	}

	for (ix=0; ix < nxb; ix++) {
		for (iz=0; iz < nzb-1; iz++) {
			rhopmlz[ix][iz] = rhopml[ix][iz] + rhopml[ix][iz+1];
		}
	}
	for (ix=0; ix < nxb; ix++) {
		rhopmlx[ix][nzb-1] = rhopml[ix][nzb-1];
	}


	/* read velocity model */
	vp = sf_floatalloc2(nz,nx);
	vs = sf_floatalloc2(nz,nx);

    sf_floatread(vp[0],nz*nx,vpfp);
    sf_floatread(vs[0],nz*nx,vsfp);

	vppml = extmodel2d(vp,nz,nx,npml);
	vspml = extmodel2d(vs,nz,nx,npml);

	free(vp[0]);free(vp);
	free(vs[0]);free(vs);

	vpmin = vpmax = vppml[0][0];
	for (ix=0; ix < nxb; ix++) {
		for (iz=0; iz < nzb; iz++) {
			if (vpmax <= vppml[ix][iz])
				vpmax = vppml[ix][iz];
			if (vpmin >= vppml[ix][iz])
				vpmin = vppml[ix][iz];
		}
	} 
	
	if(dx <=dz){
		 dmin = dx;
		 dmax=dz;
	} else {
		 dmin = dz;
		 dmax = dx;
	}

	 if (verb) {
		 fprintf(stderr, "vpmax=%f, vpmin=%f\n", vpmax, vpmin);
		 fprintf(stderr, "dmax=%f, dmin=%f\n", dmax, dmin); }
	
	/* check stability condition */
	ncourant = (vpmax * dt) / dmin;

	switch (order) {
		case 12 : if (ncourant > S12/sqrt(2.)) sf_warning("The scheme is unstable");break; 
		case 10 : if (ncourant > S10/sqrt(2.)) sf_warning("The scheme is unstable");break; 
		case 8  : if (ncourant > S8/sqrt(2.)) sf_warning("The scheme is unstable");break; 
		case 6  : if (ncourant > S6/sqrt(2.)) sf_warning("The scheme is unstable");break; 
		case 4  : if (ncourant > S4/sqrt(2.)) sf_warning("The scheme is unstable");break; 
		case 2  : if (ncourant > S2/sqrt(2.)) sf_warning("The scheme is unstable");break; 
		default : sf_error("The order should be 2N, N=1, 2, 3, 4, 5, 6");
	}
	
	switch (order) {
		case 12 : if (vpmin/dmax/fpeak < 5) sf_warning("Frequency dispersion may happen");break; 
		case 10 : if (vpmin/dmax/fpeak < 6) sf_warning("Frequency dispersion may happen");break; 
		case 8  : if (vpmin/dmax/fpeak < 7) sf_warning("Frequency dispersion may happen");break; 
		case 6  : if (vpmin/dmax/fpeak < 8) sf_warning("Frequency dispersion may happen");break; 
		case 4  : if (vpmin/dmax/fpeak < 10) sf_warning("Frequency dispersion may happen");break; 
		case 2  : if (vpmin/dmax/fpeak < 20) sf_warning("Frequency dispersion may happen");break; 
		default : sf_error("The order should be 2N, N=1, 2, 3, 4, 5, 6");
	}

	/* effective media parameters */
	lambda = sf_floatalloc2(nzb,nxb);
	mu = sf_floatalloc2(nzb,nxb);
	lambda2mu = sf_floatalloc2(nzb,nxb);
	for (ix=0; ix < nxb; ix++) {
		for (iz=0; iz < nzb; iz++) {
			lambda2mu[ix][iz] = rhopml[ix][iz]*vppml[ix][iz]*vppml[ix][iz];
		}
	}
	for (ix=0; ix < nxb; ix++) {
		for (iz=0; iz < nzb; iz++) {
			mu[ix][iz] = rhopml[ix][iz]*vspml[ix][iz]*vspml[ix][iz];
		}
	}
	for (ix=0; ix < nxb; ix++) {
		for (iz=0; iz < nzb; iz++) {
			lambda[ix][iz] = lambda2mu[ix][iz]-2.*mu[ix][iz];
		}
	}
	
	free(vspml[0]);free(vspml);

	//PML x direction
	rcoef = 0.00001;
	pmlthick = npml * dx;
	dp0 = 3.* 1. * log(1./rcoef) / (2. * pmlthick);
	
	dpleft = pmlthick;
	dpright = (nxb - 1) * dx - pmlthick;

	dpx1 = sf_floatalloc(nxb);
	dpx2 = sf_floatalloc(nxb);

	for (ix=0; ix < nxb; ix++) {
		dptmp = ix * dx;
		if (dptmp < dpleft) {
			dpx1[ix] = dp0 * pow((dpleft - dptmp)/pmlthick, 2);
			dpx2[ix] = dp0 * pow((dpleft - dptmp-0.5*dx)/pmlthick, 2);
		} else if (dptmp >= 0.9999*dpright) {
			dpx1[ix] = dp0 * pow((dptmp - dpright)/pmlthick, 2);
			dpx2[ix] = dp0 * pow((dptmp - dpright+0.5*dx)/pmlthick, 2);
		} else {
			dpx1[ix] = 0.;
			dpx2[ix] = 0.;
		}
	}
	
	//PML z direction
	pmlthick = npml * dz;
	dp0 = 3.* 1. * log(1./rcoef) / (2. * pmlthick);
	
	dpleft = pmlthick;
	dpright = (nzb - 1) * dz - pmlthick;

	dpz1 = sf_floatalloc(nzb);
	dpz2 = sf_floatalloc(nzb);

	for (ix=0; ix < nzb; ix++) {
		dptmp = ix * dz;
		if (dptmp < dpleft) {
			dpz1[ix] = dp0 * pow((dpleft - dptmp)/pmlthick, 2);
			dpz2[ix] = dp0 * pow((dpleft - dptmp-0.5*dz)/pmlthick, 2);
		} else if (dptmp >= 0.9999*dpright) {
			dpz1[ix] = dp0 * pow((dptmp - dpright)/pmlthick, 2);
			dpz2[ix] = dp0 * pow((dptmp - dpright+0.5*dz)/pmlthick, 2);
		} else {
			dpz1[ix] = 0.;
			dpz2[ix] = 0.;
		}
	}


	/* read low rank fd scheme */
	Gpx = sf_floatalloc3(nzb, nxb, Nhalf);
	Gpz = sf_floatalloc3(nzb, nxb, Nhalf);
	Gsx = sf_floatalloc3(nzb, nxb, Nhalf);
	Gsz = sf_floatalloc3(nzb, nxb, Nhalf);

    sf_floatread(Gpx[0][0],nzb*nxb*Nhalf,Gpxfp);
    sf_floatread(Gpz[0][0],nzb*nxb*Nhalf,Gpzfp);
    sf_floatread(Gsx[0][0],nzb*nxb*Nhalf,Gsxfp);
    sf_floatread(Gsz[0][0],nzb*nxb*Nhalf,Gszfp);

	Sxxtmp = sf_floatalloc(Nhalf);
	Sxztmp = sf_floatalloc(Nhalf);
	Szxtmp = sf_floatalloc(Nhalf);
	Szztmp = sf_floatalloc(Nhalf);

	sf_floatread(Sxxtmp, Nhalf, sxxfp);
	sf_floatread(Sxztmp, Nhalf, sxzfp);
	sf_floatread(Szxtmp, Nhalf, szxfp);
	sf_floatread(Szztmp, Nhalf, szzfp);

	Sxx = sf_intalloc(Nhalf);
	Sxz = sf_intalloc(Nhalf);
	Szx = sf_intalloc(Nhalf);
	Szz = sf_intalloc(Nhalf);

	for (ix=0; ix < Nhalf; ix++) {
		Sxx[ix] = (int) Sxxtmp[ix];
		Sxz[ix] = (int) Sxztmp[ix];
		Szx[ix] = (int) Szxtmp[ix];
		Szz[ix] = (int) Szztmp[ix];
	}


	free(Sxxtmp);free(Sxztmp);
	free(Szxtmp);free(Szztmp);
	
	/* allocate workspace */
	vp_1 = sf_floatalloc2(nzb, nxb);
	vp_2 = sf_floatalloc2(nzb, nxb);
	vx_1 = sf_floatalloc2(nzb, nxb);
	vx_2 = sf_floatalloc2(nzb, nxb);
	vz_1 = sf_floatalloc2(nzb, nxb);
	vz_2 = sf_floatalloc2(nzb, nxb);

	tpxx_1 = sf_floatalloc2(nzb, nxb);
	tpxx_2 = sf_floatalloc2(nzb, nxb);
	tpzz_1 = sf_floatalloc2(nzb, nxb);
	tpzz_2 = sf_floatalloc2(nzb, nxb);
	txx = sf_floatalloc2(nzb, nxb);
	tzz = sf_floatalloc2(nzb, nxb);
	txz_1 = sf_floatalloc2(nzb, nxb);
	txz_2 = sf_floatalloc2(nzb, nxb);
	
	/* output p and s */
	datx = sf_floatalloc2(nzb, nxb);
	datz = sf_floatalloc2(nzb, nxb);
	datpx = sf_floatalloc2(nzb, nxb);
	datpz = sf_floatalloc2(nzb, nxb);
	datsx = sf_floatalloc2(nzb, nxb);
	datsz = sf_floatalloc2(nzb, nxb);

	/* initialize arrays */
#ifdef _OPENMP
#pragma omp parallel for schedule(dynamic,1)  private(ix,iz)
#endif
	for (ix=0; ix < nxb; ix++) {
		for (iz=0; iz < nzb; iz++) {
			vp_1[ix][iz] = 0.;
			vp_2[ix][iz] = 0.;
			vx_1[ix][iz] = 0.;
			vx_2[ix][iz] = 0.;
			vz_1[ix][iz] = 0.;
			vz_2[ix][iz] = 0.;

			tpxx_1[ix][iz] = 0.;
			tpxx_2[ix][iz] = 0.;
			tpzz_1[ix][iz] = 0.;
			tpzz_2[ix][iz] = 0.;
			txx[ix][iz] = 0.;
			tzz[ix][iz] = 0.;
			txz_1[ix][iz] = 0.;
			txz_2[ix][iz] = 0.;

			datx[ix][iz] = 0.;
			datz[ix][iz] = 0.;
			datpx[ix][iz] = 0.;
			datpz[ix][iz] = 0.;
			datsx[ix][iz] = 0.;
			datsz[ix][iz] = 0.;
		}
	}


	/* start time loop */
	for (it=0; it < nt; it++) {
		fprintf(stderr, "\b\b\b\b\b\b%d", it);


		//velocity

#ifdef _OPENMP
#pragma omp parallel for schedule(dynamic,1)  private(ix,iz,ddx,ddp,ddz,dpcoe) shared(vx_1,vx_2,vp_1,dpx1,dpz1,txx,txz_1,txz_2,tpxx_1,tpxx_2, vppml, rhopml)	
#endif
			for (ix=Nhalf; ix < nxb-Nhalf+1; ix++) {
				for (iz=Nhalf; iz < nzb-Nhalf+1; iz++) {

					ddx = lrdx(txx, ix, iz, Nhalf, Gsx);
					
					ddp = lrdx(tpxx_1, ix, iz, Nhalf, Gpx) + \
						  lrdx(tpxx_2, ix, iz, Nhalf, Gpx);

					ddz = lrdz(txz_1, ix, iz, Nhalf, Gsz) + \
						  lrdz(txz_2, ix, iz, Nhalf, Gsz);

					dpcoe = dpx1[ix];
					vx_1[ix][iz] = (vx_1[ix][iz]*(1./dt -dpcoe*vppml[ix][iz]/2.) + \
										ddx/rhopml[ix][iz])/(1./dt + dpcoe*vppml[ix][iz]/2.);
					
					vp_1[ix][iz] = (vp_1[ix][iz]*(1./dt -dpcoe*vppml[ix][iz]/2.) + \
										ddp/rhopml[ix][iz])/(1./dt + dpcoe*vppml[ix][iz]/2.);

					dpcoe = dpz1[iz];
					vx_2[ix][iz] = (vx_2[ix][iz]*(1./dt -dpcoe*vppml[ix][iz]/2.) + \
										ddz/rhopml[ix][iz])/(1./dt + dpcoe*vppml[ix][iz]/2.);
				}
			}


#ifdef _OPENMP
#pragma omp parallel for schedule(dynamic,1)  private(ix,iz,ddx,ddp,ddz,dpcoe) shared(vz_1,vz_2,vp_2,dpx2,dpz2,tpzz_1,tpzz_2,tzz,txz_1,txz_2, vppml, rhopml)	
#endif
			for (ix=Nhalf-1; ix < nxb-Nhalf; ix++) {
				for (iz=Nhalf-1; iz < nzb-Nhalf; iz++) {

					ddx = lrdx(txz_1, ix+1, iz, Nhalf, Gsx) + \
						  lrdx(txz_2, ix+1, iz, Nhalf, Gsx);

					ddp = lrdz(tpzz_1, ix, iz+1, Nhalf, Gpz) + \
						  lrdz(tpzz_2, ix, iz+1, Nhalf, Gpz);

					ddz = lrdz(tzz, ix, iz+1, Nhalf, Gsz);

					dpcoe = dpx2[ix];
					vz_1[ix][iz] = (vz_1[ix][iz]*(1./dt -dpcoe*vppml[ix][iz]/2.) + \
										ddx/rhopml[ix][iz])/(1./dt + dpcoe*vppml[ix][iz]/2.);
					dpcoe = dpz2[iz];
					vz_2[ix][iz] = (vz_2[ix][iz]*(1./dt -dpcoe*vppml[ix][iz]/2.) + \
										ddz/rhopml[ix][iz])/(1./dt + dpcoe*vppml[ix][iz]/2.);
					vp_2[ix][iz] = (vp_2[ix][iz]*(1./dt -dpcoe*vppml[ix][iz]/2.) + \
										ddp/rhopml[ix][iz])/(1./dt + dpcoe*vppml[ix][iz]/2.);

				}
			}

	//stress

#ifdef _OPENMP
#pragma omp parallel for schedule(dynamic,1)  private(ix,iz,ddx,ddxp,ddz,ddzp,dpcoe) shared(vx_1,vx_2,vz_1,vz_2,vp_1,vp_2,dpx2,dpz1,txx,tzz,tpxx_1,tpxx_2,tpzz_1,tpzz_2, vppml, lambda2mu, mu)	
#endif
			for (ix=Nhalf-1; ix < nxb-Nhalf; ix++) {//ix = Nhalf -1
				for (iz=Nhalf; iz < nzb-Nhalf+1; iz++) { // iz < nzb-Nhalf+1
					ddx = lrdx(vx_1, ix+1, iz, Nhalf, Gsx) + \
						  lrdx(vx_2, ix+1, iz, Nhalf, Gsx) + \
						  lrdx(vp_1, ix+1, iz, Nhalf, Gsx);
					ddxp = lrdx(vx_1, ix+1, iz, Nhalf, Gpx) + \
						  lrdx(vx_2, ix+1, iz, Nhalf, Gpx) + \
						  lrdx(vp_1, ix+1, iz, Nhalf, Gpx);

					ddz = lrdz(vz_1, ix, iz, Nhalf, Gsz) + \
						  lrdz(vz_2, ix, iz, Nhalf, Gsz) + \
						  lrdz(vp_2, ix, iz, Nhalf, Gsz);
					ddzp = lrdz(vz_1, ix, iz, Nhalf, Gpz) + \
						  lrdz(vz_2, ix, iz, Nhalf, Gpz) + \
						  lrdz(vp_2, ix, iz, Nhalf, Gpz);
				
					dpcoe = dpx2[ix];
					tzz[ix][iz] = (tzz[ix][iz]*(1./dt -dpcoe*vppml[ix][iz]/2.) - \
							 2.*mu[ix][iz]*ddx)/(1./dt + dpcoe*vppml[ix][iz]/2.);
					
					tpxx_1[ix][iz] = (tpxx_1[ix][iz]*(1./dt -dpcoe*vppml[ix][iz]/2.) + \
							 lambda2mu[ix][iz]*ddxp)/(1./dt + dpcoe*vppml[ix][iz]/2.);
					
					tpzz_1[ix][iz] = (tpzz_1[ix][iz]*(1./dt -dpcoe*vppml[ix][iz]/2.) + \
							 lambda2mu[ix][iz]*ddxp)/(1./dt + dpcoe*vppml[ix][iz]/2.);

					dpcoe = dpz1[iz];
					txx[ix][iz] = (txx[ix][iz]*(1./dt -dpcoe*vppml[ix][iz]/2.) - \
							 2.*mu[ix][iz]*ddz)/(1./dt + dpcoe*vppml[ix][iz]/2.);
					
					tpxx_2[ix][iz] = (tpxx_2[ix][iz]*(1./dt -dpcoe*vppml[ix][iz]/2.) + \
							 lambda2mu[ix][iz]*ddzp)/(1./dt + dpcoe*vppml[ix][iz]/2.);
					tpzz_2[ix][iz] = (tpzz_2[ix][iz]*(1./dt -dpcoe*vppml[ix][iz]/2.) + \
							 lambda2mu[ix][iz]*ddzp)/(1./dt + dpcoe*vppml[ix][iz]/2.);
				}
			}

#ifdef _OPENMP
#pragma omp parallel for schedule(dynamic,1)  private(ix,iz,ddx,ddz,dpcoe) shared(vx_1,vx_2,vz_1,vz_2,vp_1,vp_2,dpx1,dpz2,txz_1,txz_2, vppml, mu)	
#endif
			for (ix=Nhalf; ix < nxb-Nhalf+1; ix++) {
				for (iz=Nhalf-1; iz < nzb-Nhalf; iz++) { //Nhalf-1
					
					ddx = lrdx(vz_1, ix, iz, Nhalf, Gsx) + \
						  lrdx(vz_2, ix, iz, Nhalf, Gsx) + \
						  lrdx(vp_2, ix, iz, Nhalf, Gsx);

					ddz = lrdz(vx_1, ix, iz+1, Nhalf, Gsz) + \
						  lrdz(vx_2, ix, iz+1, Nhalf, Gsz) + \
						  lrdz(vp_1, ix, iz+1, Nhalf, Gsz);
					
					dpcoe = dpx1[ix];
					txz_1[ix][iz] = (txz_1[ix][iz]*(1./dt -dpcoe*vppml[ix][iz]/2.) + \
							mu[ix][iz]*ddx)/(1./dt + dpcoe*vppml[ix][iz]/2.);

					dpcoe = dpz2[iz];
					txz_2[ix][iz] = (txz_2[ix][iz]*(1./dt -dpcoe*vppml[ix][iz]/2.) + \
							mu[ix][iz]*ddz)/(1./dt + dpcoe*vppml[ix][iz]/2.);
				}
			}

		//add source

		explsourcet(vz_1, fpeak, it, dt, isrcx,  isrcz);
		explsourcet(vz_2, fpeak, it, dt, isrcx,  isrcz);
		explsourcet(vp_2, fpeak, it, dt, isrcx,  isrcz);


		//x
#ifdef _OPENMP
#pragma omp parallel for schedule(dynamic,1)  private(iz)
#endif
			for (iz=0; iz < nzb; iz++) {
				vp_1[0][iz] = 0.;
				vp_2[0][iz] = 0.;
				vx_1[0][iz] = 0.;
				vx_2[0][iz] = 0.;
				vz_1[0][iz] = 0.;
				vz_2[0][iz] = 0.;
			}


#ifdef _OPENMP
#pragma omp parallel for schedule(dynamic,1)  private(iz)
#endif
			for (iz=0; iz < nzb; iz++) {
				vp_1[nxb-1][iz] = 0.;
				vp_2[nxb-1][iz] = 0.;
				vx_1[nxb-1][iz] = 0.;
				vx_2[nxb-1][iz] = 0.;
				vz_1[nxb-1][iz] = 0.;
				vz_2[nxb-1][iz] = 0.;
			}
		
		//z
#ifdef _OPENMP
#pragma omp parallel for schedule(dynamic,1)  private(ix)
#endif
			for (ix=0; ix < nxb; ix++) {
				vp_1[ix][0] = 0.;
				vp_2[ix][0] = 0.;
				vx_1[ix][0] = 0.;
				vx_2[ix][0] = 0.;
				vz_1[ix][0] = 0.;
				vz_2[ix][0] = 0.;
			}

#ifdef _OPENMP
#pragma omp parallel for schedule(dynamic,1)  private(ix)
#endif
			for (ix=0; ix < nxb; ix++) {
				vp_1[ix][nzb-1] = 0.;
				vp_2[ix][nzb-1] = 0.;
				vx_1[ix][nzb-1] = 0.;
				vx_2[ix][nzb-1] = 0.;
				vz_1[ix][nzb-1] = 0.;
				vz_2[ix][nzb-1] = 0.;
			}

		/* output x, z components */
		if (it%snap == 0) {
			for (ix=0; ix < nx + 2*npml; ix++) {
				for (iz=0; iz < nz + 2*npml; iz++) {
					dptmp = vx_1[ix][iz] + vx_2[ix][iz] + vp_1[ix][iz];
            		sf_floatwrite(&dptmp,1,snapfpx);
					
					dptmp = vz_1[ix][iz] + vz_2[ix][iz] + vp_2[ix][iz];
            		sf_floatwrite(&dptmp,1,snapfpz);
				}
			}
		}

		/* output separated p and s wavefiled */
		if (it%snap == 0) {
			for (ix=0; ix < nx + 2*npml; ix++) {
				for (iz=0; iz < nz + 2*npml; iz++) {
					datpx[ix][iz] = vp_1[ix][iz];
					datpz[ix][iz] = vp_2[ix][iz];
					
					datsx[ix][iz] = vx_1[ix][iz] + vx_2[ix][iz];
					datsz[ix][iz] = vz_1[ix][iz] + vz_2[ix][iz];
					
				}
			}
      sf_floatwrite(datpx[0],nxb*nzb,snapfppx);
      sf_floatwrite(datpz[0],nxb*nzb,snapfppz);
      sf_floatwrite(datsx[0],nxb*nzb,snapfpsx);
      sf_floatwrite(datsz[0],nxb*nzb,snapfpsz);
		}

	} // end time loop

	free(vp_1[0]);free(vp_1);
	free(vp_2[0]);free(vp_2);
	free(vx_1[0]);free(vx_1);
	free(vx_2[0]);free(vx_2);
	free(vz_1[0]);free(vz_1);
	free(vz_2[0]);free(vz_2);

	free(tpxx_1[0]);free(tpxx_1);
	free(tpxx_2[0]);free(tpxx_2);
	free(tpzz_1[0]);free(tpzz_1);
	free(tpzz_2[0]);free(tpzz_2);
	free(txx[0]);free(txx);
	free(tzz[0]);free(tzz);
	free(txz_1[0]);free(txz_1);
	free(txz_2[0]);free(txz_2);
	
	free(dpx1);free(dpx2);
	free(dpz1);free(dpz2);

	
	free(vppml[0]);free(vppml);
	free(rhopmlx[0]);free(rhopmlx);
	free(rhopmlz[0]);free(rhopmlz);

	free(lambda[0]);free(lambda);
	free(mu[0]);free(mu);
	free(lambda2mu[0]);free(lambda2mu);

	
	free(datx[0]);free(datx);
	free(datz[0]);free(datz);
	
	free(datpx[0]);free(datpx);
	free(datpz[0]);free(datpz);
	free(datsx[0]);free(datsx);
	free(datsz[0]);free(datsz);

	free(Gpx[0][0]);free(Gpx[0]);free(Gpx);
	free(Gpz[0][0]);free(Gpz[0]);free(Gpz);
	free(Gsx[0][0]);free(Gsx[0]);free(Gsx);
	free(Gsz[0][0]);free(Gsz[0]);free(Gsz);

	free(Sxx);free(Sxz);
	free(Szx);free(Szz);
	
	fprintf(stderr, "\nDone\n");

	return 0;

}

static float lrdx(float **data, int ix, int iz, int lenx, float ***gx)
/*<Low rank finite difference : d/dx>*/
{
    float res = 0.0;
    int il;
    for (il = 0; il < lenx; il++) {
	res += 0.5*(data[ix+Sxx[il]-1][iz+Sxz[il]] - data[ix-Sxx[il]][iz-Sxz[il]])*gx[il][ix][iz];
    }
    return res;
}


static float lrdz(float **data, int ix, int iz, int lenz, float ***gz)
/*<Low rank finite difference : d/dz>*/
{
    float res = 0.0;
    int il;
    for (il = 0; il < lenz; il++) {
	res += 0.5*(data[ix+Szx[il]][iz+Szz[il]-1] - data[ix-Szx[il]][iz-Szz[il]])*gz[il][ix][iz];
    }
    return res;
}


float **extmodel2d(float **init_model,int nz,int nx,int np)
{
	float **p;
	int ix,iz;
	int nx2=nx+2*np; 
	int nz2=nz+2*np;

	p=sf_floatalloc2(nz2,nx2);	
	
	for (ix=0; ix < np; ix++) {
		for (iz=np; iz < nz+np; iz++) {
			p[ix][iz] = init_model[0][iz-np];
		}
	}
	
	for (ix=nx+np; ix < nx2; ix++) {
		for (iz=np; iz < nz+np; iz++) {
			p[ix][iz] = init_model[nx-1][iz-np];
		}
	}
	
	for (ix=np; ix < nx+np; ix++) {
		for (iz=nz+np; iz < nz2; iz++) {
			p[ix][iz] = init_model[ix-np][nz-1];
		}
	}

	for (ix=np; ix < np+nx; ix++) {
		for (iz=0; iz < np; iz++) {
			p[ix][iz] = init_model[ix-np][0];
		}
	}

	for (ix=0; ix < np; ix++) {
		for (iz=0; iz < np; iz++) {
			p[ix][iz] = init_model[0][0];
		}
	}
			
	for (ix=nx+np; ix < nx2; ix++) {
		for (iz=0; iz < np; iz++) {
			p[ix][iz] = init_model[nx-1][0];
		}
	}						

	for (ix=0; ix < np; ix++) {
		for (iz=nz+np; iz < nz2; iz++) {
			p[ix][iz] = init_model[0][nz-1];
		}
	}
	
	for (ix=nx+np; ix < nx2; ix++) {
		for (iz=nz+np; iz < nz2; iz++) {
			p[ix][iz] = init_model[nx-1][nz-1];
		}
	}
	
	for (ix=np; ix < nx+np; ix++) {
		for (iz=np; iz < nz+np; iz++) {
			p[ix][iz] = init_model[ix-np][iz-np];
		}
	}
			
	return p;
}	

void explsourcet(float **txx/*@out@*/,
				 float fpeak,
		 		 int it, float dt,
                 int sx, int sz 
				)
/*<explosive source for stress txx>*/ 
{
	float wavelet, tmp;
	tmp = SF_PI*fpeak*(it*dt-1.2/fpeak);
	wavelet = (1.-2.*tmp*tmp)*exp(-tmp*tmp);
  
	txx[sx][sz] += wavelet * dt;
     
}

