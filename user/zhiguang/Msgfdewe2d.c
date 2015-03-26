/* 2-D staggered-grid elastic time-domain FD modeling 

10th order in space, 2th order in time. 
*/
/*
  Copyright (C) 2013 China University of Petroleum
  
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
#include <math.h>
/*^*/

#ifdef _OPENMP
#include <omp.h>
#endif

#define TS (2.0/favg)
#define NS ((int)(TS/dt)+1)

static double coef[5]={1.211243f, -0.08972168f, 0.01384277f, -0.00176566f, 0.0001186795f};
static float *source;
static float *absx, *absz;
static int npd, rangex, rangez, sx, sz;
static float dtx, dtz;

// function declaration**************
void get_source(int ns, float dt, float ts, float favg);
void get_absorb_par(int nx, int nz, int nabs);
void zero2d(int n1, int n2, float **array);
void get_current_vel_rho(int is, int ds, int nx, int nz, int nvx, float **vp0, 
              float **vs0, float **rho0, float **vp, float **vs, float **rho);
void pad_vel_rho(int nx, int nz, float **in, float **out);
void update_vel(int it, int hsz, float **u, float **w, float **txx, float **tzz, 
                         float **txz, float **rho, float **datax, float **dataz);
void update_stress(float **u, float **w, float **txx, float **tzz, float **txz, 
                                                        float **vp, float **vs);
void ptsrc(int it, int nx, int nz, float **txx, float **tzz);

// main function************
int main(int argc, char* argv[])
{  
    int ishot, nshot, dshot;
    int it, ix, iz, nt, nvx, nx, nz, nabs;
    int sx_ini, sz_ini, hsz;
    float dt, tmax, dx, dz, favg;
    float **vp0, **vs0, **rho0;
    float **vp, **vs, **rho;
    float **u, **w, **txx, **tzz, **txz;
    float **datax, **dataz;
    sf_file dataxf, datazf, vpf, vsf, rhof;
    
    sf_init(argc, argv);
    
#ifdef _OPENMP
    omp_init();
#endif

    dataxf= sf_output("datax");
    /* horizontal component of records */
    datazf= sf_output("dataz");
    /* vertical component of records */
    vpf   = sf_input("vp");
    /* p-wave velocity file */
    vsf   = sf_input("vs");
    /* s-wave velocity file */
    rhof  = sf_input("rho");
    /* density file */
    
    if(!sf_histint  (vpf, "n1", &nz))  sf_error("No n1= in vpf");
    if(!sf_histint  (vpf, "n2", &nvx)) sf_error("No n2= in vpf");
    if(!sf_histfloat(vpf, "d1", &dz))  sf_error("No d1= in vpf");
    if(!sf_histfloat(vpf, "d2", &dx))  sf_error("No d2= in vpf");
    
    if(!sf_getint("nshot", &nshot)) nshot=1;
    /* number of shots */
    if(!sf_getint("dshot", &dshot)) dshot=1;
    /* shot interval, multiple of receiver intervals */
    if(!sf_getint("nx",    &nx))    nx=nvx;
    /* coverage area for each shot */
    if(!sf_getint("nabs",  &nabs))  nabs=50;
    /* width of padded boundary */
    
    if(!sf_getint("sx_ini", &sx_ini)) 
    sf_error("Need input the coordinates of the shot source");
    /* horizontal position of shot point */
    if(!sf_getint("sz_ini", &sz_ini)) sz_ini=5;
    /* vertical position of shot point */
    if(!sf_getint("hsz", &hsz)) hsz=5;
    /* vertical position of recerivers */
    
    if(!sf_getfloat("tmax", &tmax)) sf_error("Need input record time");
    /* record length */
    if(!sf_getfloat("dt",   &dt  )) sf_error("Need input time interval");
    /* time interval */
    if(!sf_getfloat("peak", &favg)) favg=20;
    /* peak frequency for Ricker wavelet (in Hz) */
    
    nt=(int)(tmax/dt)+1;
    sf_warning("tmax=%f s, dt=%f s, nt=%d", tmax, dt, nt);
    
    npd=nabs+5;
    rangex=nx+2*npd;
    rangez=nz+2*npd;
    sx=sx_ini+npd;
    sz=sz_ini+npd;
    
    dtx=dt/dx;
    dtz=dt/dz;
    
    sf_putint  (dataxf, "n1", nt);
    sf_putfloat(dataxf, "d1", dt);
    sf_putfloat(dataxf, "o1", 0.);
    sf_putstring(dataxf,"unit1","s");
    sf_putint  (dataxf, "n2", nx);
    sf_putfloat(dataxf, "d2", dx);
    sf_putstring(dataxf, "unit2", "m");
    sf_putfloat(dataxf, "o2", -dx*(sx_ini-1));
    sf_putint  (dataxf, "n3", nshot);
    sf_putint  (dataxf, "d3", 1);
    sf_putint  (dataxf, "o3", 1);
    
    sf_putint  (datazf, "n1", nt);
    sf_putfloat(datazf, "d1", dt);
    sf_putfloat(datazf, "o1", 0.);
    sf_putstring(datazf,"unit1","s");
    sf_putint  (datazf, "n2", nx);
    sf_putfloat(datazf, "d2", dx);
    sf_putfloat(datazf, "o2", -dx*(sx_ini-1));
    sf_putstring(datazf,"unit2","m");
    sf_putint  (datazf, "n3", nshot);
    sf_putint  (datazf, "d3", 1);
    sf_putint  (datazf, "o3", 1);
     
    vp0 = sf_floatalloc2(nz, nvx);
    vs0 = sf_floatalloc2(nz, nvx);
    rho0= sf_floatalloc2(nz, nvx);
    
    sf_floatread(vp0[0], nz*nvx, vpf);
    sf_floatread(vs0[0], nz*nvx, vsf);
    sf_floatread(rho0[0], nz*nvx, rhof);
    
    for(ix=0; ix<nvx; ix++){
        for(iz=0; iz<nz; iz++){
            vp0[ix][iz]=rho0[ix][iz]*vp0[ix][iz]*vp0[ix][iz];
            vs0[ix][iz]=rho0[ix][iz]*vs0[ix][iz]*vs0[ix][iz];
            rho0[ix][iz]=1.0/rho0[ix][iz];
        }
    }
    
    vp=sf_floatalloc2(rangez, rangex);
    vs=sf_floatalloc2(rangez, rangex);
    rho=sf_floatalloc2(rangez, rangex);
    
    u=sf_floatalloc2(rangez, rangex);
    w=sf_floatalloc2(rangez, rangex);
    txx=sf_floatalloc2(rangez, rangex);
    tzz=sf_floatalloc2(rangez, rangex);
    txz=sf_floatalloc2(rangez, rangex);
    
    datax=sf_floatalloc2(nt, nx);
    dataz=sf_floatalloc2(nt, nx);
    
    source=sf_floatalloc(NS);
    get_source(NS, dt, TS, favg);
    
    absx=sf_floatalloc(rangex);
    absz=sf_floatalloc(rangez);
    get_absorb_par(nx, nz, nabs);
    
    for(ishot=1; ishot<=nshot; ishot++){
        sf_warning("---%dth shot of %d is beginning!---", ishot, nshot);
        
        zero2d(rangez, rangex, u);
        zero2d(rangez, rangex, w);
        zero2d(rangez, rangex, txx);
        zero2d(rangez, rangex, tzz);
        zero2d(rangez, rangex, txz);
        
        zero2d(nt, nx, datax);
        zero2d(nt, nx, dataz);
        
        get_current_vel_rho(ishot, dshot, nx, nz, nvx, vp0, vs0, rho0, vp, vs, rho);
        
        for(it=0; it<nt; it++){
            sf_warning("@@@@:is=%d, it=%d;",ishot, it+1);
            update_vel(it, hsz, u, w, txx, tzz, txz, rho, datax, dataz);
            update_stress(u, w, txx, tzz, txz, vp, vs);
            if(it<NS) ptsrc(it, nx, nz, txx, tzz);
        }
        
        sf_floatwrite(datax[0], nt*nx, dataxf);
        sf_floatwrite(dataz[0], nt*nx, datazf);
    }
    
    free(*vp0); free(vp0);
    free(*vs0); free(vs0);
    free(*rho0); free(rho0);
    free(*vp); free(vp);
    free(*vs); free(vs);
    free(*rho); free(rho);
    free(*u); free(*w); free(*txx); free(*tzz); free(*txz);
    free(u); free(w); free(txx); free(tzz); free(txz);
    free(*datax); free(datax); free(*dataz); free(dataz);
    free(source); free(absx); free(absz);
    
    exit(0);
}

// subroutines*************
void get_source(int ns, float dt, float ts, float favg)
{
    int it; 
    float t, tmp;
    t=-ts/2.;
    for(it=0; it<ns; it++){
      tmp=favg*t*SF_PI;
      tmp*=tmp;
      source[it]=(1-2*tmp)*expf(-tmp);
      t+=dt;
    }
    return;
}

void get_absorb_par(int nx, int nz, int nabs)
{
    int ix, iz;
    
    for(ix=0; ix<rangex; ix++)
        absx[ix]=1.0;
    for(iz=0; iz<rangez; iz++)
        absz[iz]=1.0;
    
    for(ix=0; ix<nabs; ix++){
        absx[nabs+4-ix]=1.0-0.15*(1.0*(ix+1)/nabs)*(1.0*(ix+1)/nabs);
        absx[nabs+5+nx+ix]=1.0-0.15*(1.0*(ix+1)/nabs)*(1.0*(ix+1)/nabs);
    }
    for(iz=0; iz<nabs; iz++){
        absz[nabs+4-iz]=1.0-0.15*(1.0*(iz+1)/nabs)*(1.0*(iz+1)/nabs);
        absz[nabs+5+nz+iz]=1.0-0.15*(1.0*(iz+1)/nabs)*(1.0*(iz+1)/nabs);
    }
    return;
}

void zero2d(int n1, int n2, float **array)
{
    int i1, i2;
    for(i2=0; i2<n2; i2++){
        for(i1=0; i1<n1; i1++){
            array[i2][i1]=0.0;
        }
    }
    return;
}

void get_current_vel_rho(int is, int ds, int nx, int nz, int nvx, float **vp0, 
               float **vs0, float **rho0, float **vp, float **vs, float **rho)
{
    int ix, iz;
    int vstart, vend;    
    float **vptemp, **vstemp, **rhotemp;
    
    vstart=(is-1)*ds;
    vend=vstart+nx-1;
    
    if(vstart<0) sf_error("Vstart is less than zero!");
    if(vend>=nvx) sf_error("Vend is bigger than nvx!");
    
    vptemp=sf_floatalloc2(nz, nx);
    vstemp=sf_floatalloc2(nz, nx);
    rhotemp=sf_floatalloc2(nz, nx);
    
    for(ix=vstart; ix<=vend; ix++){
        for(iz=0; iz<nz; iz++){
            vptemp[ix-vstart][iz]=vp0[ix][iz];
            vstemp[ix-vstart][iz]=vs0[ix][iz];
            rhotemp[ix-vstart][iz]=rho0[ix][iz];
        }
    }
    
    pad_vel_rho(nx, nz, vptemp, vp);
    pad_vel_rho(nx, nz, vstemp, vs);
    pad_vel_rho(nx, nz, rhotemp, rho);
    
    free(*vptemp); free(vptemp);
    free(*vstemp); free(vstemp);
    free(*rhotemp); free(rhotemp);
    return;
}

void pad_vel_rho(int nx, int nz, float **in, float **out)
{
    int ix, iz;
    for(ix=npd; ix<nx+npd; ix++)
    for(iz=npd; iz<nz+npd; iz++)
        out[ix][iz]=in[ix-npd][iz-npd];
    
    for(ix=0; ix<npd; ix++)
    for(iz=npd; iz<nz+npd; iz++)
        out[ix][iz]=out[npd][iz];
        
    for(ix=nx+npd; ix<rangex; ix++)
    for(iz=npd; iz<nz+npd; iz++)
        out[ix][iz]=out[nx+npd-1][iz];
        
    for(iz=0; iz<npd; iz++)
    for(ix=0; ix<rangex; ix++)
        out[ix][iz]=out[ix][npd];
        
    for(iz=nz+npd; iz<rangez; iz++)
    for(ix=0; ix<rangex; ix++)
        out[ix][iz]=out[ix][nz+npd-1];
        
    return;
}
        
void update_vel(int it, int hsz, float **u, float **w, float **txx, float **tzz, 
                         float **txz, float **rho, float **datax, float **dataz)
{
    int ix, jz, ic;
    float dtxx, dtzz, dtxz, dtzx;

#ifdef _OPENMP
#pragma omp parallel for                 \
    private(ix,jz,dtxx,dtzz,dtxz,dtzx,ic) \
    shared(it,rangex,rangez,coef,txx,tzz,txz,u,w,rho,dtx,dtz,absx,absz,npd,hsz,datax,dataz)
#endif
    for(ix=5; ix<rangex-5; ix++){
	for(jz=5; jz<rangez-5; jz++){
	    dtxx=0.0; dtzz=0.0; dtxz=0.0; dtzx=0.0;
			
	    for(ic=0; ic<5; ic++){
		dtxx+=coef[ic]*(txx[ix+ic][jz]-txx[ix-ic-1][jz]);
		dtzz+=coef[ic]*(tzz[ix][jz+ic+1]-tzz[ix][jz-ic]);
		dtxz+=coef[ic]*(txz[ix][jz+ic]-txz[ix][jz-ic-1]);
		dtzx+=coef[ic]*(txz[ix+ic+1][jz]-txz[ix-ic][jz]);
	    }
			
	    u[ix][jz]=(u[ix][jz]+rho[ix][jz]*(dtx*dtxx+dtz*dtxz))*absx[ix]*absz[jz];
	    w[ix][jz]=(w[ix][jz]+rho[ix][jz]*(dtx*dtzx+dtz*dtzz))*absx[ix]*absz[jz];

	    if(ix>=npd && ix<rangex-npd && jz==hsz+npd){
		datax[ix-npd][it]=u[ix][jz];
		dataz[ix-npd][it]=w[ix][jz];
	    }
	}
    }
    return;
}

void update_stress(float **u, float **w, float **txx, float **tzz, 
                              float **txz, float **vp, float **vs)
{
    int ix, jz, ic;
    float dux, duz, dwx, dwz;

#ifdef _OPENMP
#pragma omp parallel for                  \
    private(ix,jz,dux,duz,dwx,dwz,ic)     \
    shared(rangex,rangez,coef,u,w,txx,tzz,txz,dtx,dtz,vp,vs,absx,absz)
#endif
    for(ix=5; ix<rangex-5; ix++){
	for(jz=5; jz<rangez-5; jz++){
	    dux=0.0; duz=0.0; dwx=0.0; dwz=0.0;

	    for(ic=0; ic<5; ic++){
		dux+=coef[ic]*(u[ix+ic+1][jz]-u[ix-ic][jz]);
		duz+=coef[ic]*(u[ix][jz+ic+1]-u[ix][jz-ic]);
		dwx+=coef[ic]*(w[ix+ic][jz]-w[ix-ic-1][jz]);
		dwz+=coef[ic]*(w[ix][jz+ic]-w[ix][jz-ic-1]);
	    }

	txx[ix][jz]=(txx[ix][jz]+dtx*vp[ix][jz]*dux+dtz*(vp[ix][jz]-2*vs[ix][jz])*dwz)*absx[ix]*absz[jz];
	tzz[ix][jz]=(tzz[ix][jz]+dtz*vp[ix][jz]*dwz+dtx*(vp[ix][jz]-2*vs[ix][jz])*dux)*absx[ix]*absz[jz];
	txz[ix][jz]=(txz[ix][jz]+vs[ix][jz]*(dtz*duz+dtx*dwx))*absx[ix]*absz[jz];
	}
    }
    return;
}        
            
void ptsrc(int it, int nx, int nz, float **txx, float **tzz)
{
    int ix, jz;
    int ifx, ilx, jfz, jlz;
    float dex, dez;

    ifx=(npd>sx-3)?npd:sx-3;
    ilx=(nx+npd-1<sx+3)?nx+npd-1:sx+3;
    jfz=(npd>sz-3)?npd:sz-3;
    jlz=(nz+npd-1<sz+3)?nz+npd-1:sz+3;

    for(ix=ifx; ix<=ilx; ix++){
	for(jz=jfz; jz<=jlz; jz++){
	    dex=ix-sx;
	    dez=jz-sz;
	    txx[ix][jz]+=1000.0*source[it]*expf(-dex*dex-dez*dez);
	    tzz[ix][jz]+=1000.0*source[it]*expf(-dex*dex-dez*dez);
	}
    }
    return;
}
