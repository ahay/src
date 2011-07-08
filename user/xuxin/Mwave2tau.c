/* 2-D time-extrapolation in pseudo-depth domain */

/*
  Copyright (C) 2011 KAUST
  
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

/* Pseudo-depth wave equation: P_{tt} = v^2 *(P_{xx} + 2*s*P_{xz} + s^2*P_{zz}) + v0^2*P_{zz} + f */

#include <rsf.h>
#ifdef _OPENMP
#include <omp.h>
#endif
#include "wave.h"
#include "interp.h"

/* FD stencils */
#define C0 -2.500000000000000000000 /* -30/12 */
#define CA  1.333333333333333333333 /*  16/12 */
#define CB -0.083333333333333333333 /*  -1/12 */
#define DXDX(a,ix,iz,idx2) (C0*a[iz][ix]+\
                            CA*(a[iz][ix-1]+a[iz][ix+1])+\
                            CB*(a[iz][ix-2]+a[iz][ix+2]))*idx2
#define DZDZ(a,ix,iz,idz2) (C0*a[iz][ix]+\
                            CA*(a[iz-1][ix]+a[iz+1][ix])+\
                            CB*(a[iz-2][ix]+a[iz+2][ix]))*idz2
#define CX  0.250000000000000000000
#define DXDZ(a,ix,iz,idxz) (a[iz-1][ix-1]-a[iz-1][ix+1]-a[iz+1][ix-1]+a[iz+1][ix+1])*idxz*CX
#define C1  0.666666666666666666667 /*  2/3 */
#define C2 -0.083333333333333333333 /* -1/12 */
#define DX(a,ix,iz,idx) (C1*(a[iz][ix+1]-a[iz][ix-1])+\
                         C2*(a[iz][ix+2]-a[iz][ix-2]))*idx
#define DZ(a,ix,iz,idz) (C1*(a[iz+1][ix]-a[iz-1][ix])+\
                         C2*(a[iz+2][ix]-a[iz-2][ix]))*idz

int main(int argc, char* argv[])
{
    bool fsrf,wrtdat,dens;
    sf_file Fi=NULL,Fo=NULL,Fv=NULL,Fro=NULL,Fss=NULL,Frr=NULL,Fd=NULL,Fs=NULL;
    float **v=NULL,**ro=NULL,**rox=NULL,**roz=NULL,**um=NULL,**uo=NULL,**up=NULL,**ua=NULL,**vt2=NULL,**tt=NULL,**s2=NULL,**sv=NULL,**s=NULL;
    float *ww=NULL,*dd=NULL,*xs=NULL,*zs=NULL,*xr=NULL,*zr=NULL;
    int nx,nz,ix,iz,nbx,nbz,npx,npz,x,z,nt,it,ns,nr,n;
    int jsnap,ompnth,ompath,ompchunk,srf;
    float ox,oz,dx,dz,idx,idz,idxz,idx2,idz2,ot,dt,dt2,v0,v2;    
    sf_axis ax,az,at,aj;
    wave2dp wave=NULL;
    abcone2dp abc=NULL;
    sponge2dp spo=NULL;
    lint2dp lints=NULL,lintr=NULL;

    sf_init(argc,argv);

    /* files */
    Fi  = sf_input("in");    /* source data */
    Fv  = sf_input("vel");   /* velocity */
    Fs  = sf_input("sig");   /* sigma */
    Fss = sf_input("sou");   /* source coord */
    Fo  = sf_output("out");  /* wavefield */
    if (NULL != sf_getstring("den")) {
	dens = true;
	Fro= sf_input("den");   /* density */
    } else {
	dens = false;
    }
    if (NULL != sf_getstring("rec") && NULL != sf_getstring("dat")) {
	wrtdat = true;
	Frr = sf_input("rec");  /* receiver coord */
	Fd  = sf_output("dat"); /* seismogram */ 
    } else {
	wrtdat = false;
    }

    /* parameters */
    if (!sf_getbool("free",&fsrf)) fsrf=false;
    /* if y, free surface */
    if (!sf_getint("jsnap",&jsnap)) jsnap=1;
    /* jump in wavefield snapshots */
    if (!sf_getint("nbx",&nbx)) nbx=20;
    /* num of boundary grid */
    if (!sf_getint("nbz",&nbz)) nbz=20;
    /* num of boundary grid */
    if (!sf_getint("ompnth",&ompnth)) ompnth=1;
    /* OpenMP num of threads */
    if (!sf_getint("ompchunk",&ompchunk)) ompchunk=1;
    /* OpenMP chunksize */
    if (!sf_getfloat("v0",&v0)) v0=1.0f;
    /* scaling velocity */
#ifdef _OPENMP
#pragma omp parallel
    ompath = omp_get_num_threads();
    omp_set_num_threads(ompnth);
    sf_warning("using %d threads of a total of %d",ompnth,ompath);
#endif

    /* axis */
    at = sf_iaxa(Fi,2);
    ax = sf_iaxa(Fv,1);
    az = sf_iaxa(Fv,2);
    nt = sf_n(at); ot = sf_o(at); dt = sf_d(at);
    nx = sf_n(ax); ox = sf_o(ax); dx = sf_d(ax); 
    nz = sf_n(az); oz = sf_o(az); dz = sf_d(az);
    if (dens) {
	if (!sf_histint(Fro,"n1",&n) || n != nx) sf_error("Need n1=%d in den",nx);
	if (!sf_histint(Fro,"n2",&n) || n != nz) sf_error("Need n1=%d in den",nz);
    }
    if (!sf_histint(Fs,"n1",&n) || n != nx) sf_error("Need n1=%d in sig",nx);
    if (!sf_histint(Fs,"n2",&n) || n != nz) sf_error("Need n1=%d in sig",nz);
    for (n=0, it=0; it<nt; it++)
	if (!(it % jsnap)) n++;
    aj = sf_maxa(n,ot,dt*jsnap);
    sf_setlabel(at,"t"); sf_setunit(at,"");
    sf_setlabel(ax,"x"); sf_setunit(at,"");
    sf_setlabel(az,"z"); sf_setunit(at,"");
    sf_setlabel(aj,"t"); sf_setunit(aj,"");
    sf_oaxa(Fo,ax,1);
    sf_oaxa(Fo,az,2);
    sf_oaxa(Fo,aj,3);
    idx = 1.0f/dx; idx2 = idx*idx;
    idz = 1.0f/dz; idz2 = idz*idz;
    idxz = idx*idz;
    dt2 = dt*dt;

    /* grids */
    wave = wave2d_init(ax,az,at,nbx,nbz,fsrf);
    npx = sf_n(wave->px);
    npz = sf_n(wave->pz);
    srf = fsrf ? 0 : nbz;
    sf_raxa(ax); sf_raxa(wave->px);    
    sf_raxa(az); sf_raxa(wave->pz);
    sf_raxa(at); sf_raxa(aj);

    /* source coord */
    if (!sf_histint(Fi,"n1",&ns)) sf_error("No n1= in input");
    if (!sf_histint(Fss,"n1",&n) || n != ns)  sf_error("Need n1=%d in sou",ns);
    xs = sf_floatalloc(ns);
    zs = sf_floatalloc(ns);
    sf_floatread(xs,ns,Fss);
    sf_floatread(zs,ns,Fss);
    lints = lint2d_init(xs,zs,ns,wave->px,wave->pz);
    sf_warning("ns=%d",ns);
   
    /* receiver coord */
    if (wrtdat) {
	if (!sf_histint(Frr,"n1",&nr)) sf_error("No n1= in rec");
	xr = sf_floatalloc(nr);
	zr = sf_floatalloc(nr);
	sf_floatread(xr,nr,Frr);
	sf_floatread(zr,nr,Frr);
	lintr = lint2d_init(xr,zr,nr,wave->px,wave->pz);
	sf_warning("nr=%d",nr);

	sf_putint(Fd,"n1",nr);
	sf_putstring(Fd,"label1","r");
	sf_putstring(Fd,"unit1","");
	sf_oaxa(Fd,at,2);
    }    

    /* read velocity */
    tt = sf_floatalloc2(nx ,nz );
    v  = sf_floatalloc2(npx,npz);
    vt2= sf_floatalloc2(npx,npz);
    sf_floatread(tt[0],nx*nz,Fv);
    expand2(tt,v,wave);
    for     (iz=0; iz<npz; iz++) {
	for (ix=0; ix<npx; ix++) {
	    vt2[iz][ix] = v[iz][ix]*v[iz][ix] *dt*dt;
	}
    }

    /* read density */
    if (dens) {
	ro =sf_floatalloc2(npx,npz);
	rox=sf_floatalloc2(npx,npz);
	roz=sf_floatalloc2(npx,npz);
	sf_floatread(tt[0],nx*nz,Fro);
	expand2(tt,ro,wave);
	for     (iz=0; iz<nz; iz++) {
	    for (ix=0; ix<nx; ix++) {
		x = ix + nbx;
		z = iz + srf;
		rox[z][x] = DX(ro,x,z,idx)/ro[z][x];
		if (fsrf && iz<2) 
		    roz[z][x]=0.0f;
		else
		    roz[z][x] = DZ(ro,x,z,idz)/ro[z][x];
	    }
	}
    }

    /* read sigma */
    s  =sf_floatalloc2(npx,npz);
    s2 =sf_floatalloc2(npx,npz); /* 2*s */
    sv =sf_floatalloc2(npx,npz); /* s^2+(v0/v)^2 */
    sf_floatread(tt[0],nx*nz,Fs);
    expand2(tt,s,wave);
    for     (iz=0; iz<npz; iz++) {
	for (ix=0; ix<npx; ix++) {
	    v2= v[iz][ix]*v[iz][ix];
	    s2[iz][ix] = s[iz][ix]*2.0f;
	    sv[iz][ix] = s[iz][ix]*s[iz][ix] + v0/v[iz][ix] * v0/v[iz][ix];
	}
    }

    um =sf_floatalloc2(npx,npz);
    uo =sf_floatalloc2(npx,npz);
    up =sf_floatalloc2(npx,npz);
    ua =sf_floatalloc2(npx,npz);
    ww =sf_floatalloc(ns);
    if (wrtdat)
	dd =sf_floatalloc(nr);

    for (ix=0; ix<npx*npz; ix++) {
	um[0][ix] = 0.0f;
	uo[0][ix] = 0.0f;
	up[0][ix] = 0.0f;
	ua[0][ix] = 0.0f;
    }

    abc=abcone2d_tau_init(v,v0,wave);
    spo=sponge2d_init(0.00005,wave);

    /* MAIN */
    for (it=0; it<nt; it++) {
	sf_warning("it=%d;",it);

	/* write wavefield */
	if (!(it % jsnap)) {
	    for(iz=srf; iz<npz-nbz; iz++)
		sf_floatwrite(&uo[iz][nbx],nx,Fo);
	}

	/* write seismogram */
	if (wrtdat) {
	    lint2d_extract(uo,dd,lintr);
	    sf_floatwrite(dd,nr,Fd);
	}

	/* inject source */
	sf_floatread(ww,ns,Fi);
	lint2d_inject(uo,ww,lints);

	/* laplacian */
#ifdef _OPENMP
#pragma omp parallel for schedule(dynamic,ompchunk) \
    private(ix,iz,x,z)				    \
    shared(ua,uo,s2,sv,nx,nz,nbx,srf,idx2,idz2,idxz,fsrf)
#endif
	for     (iz=0; iz < nz; iz++) {
	    for (ix=0; ix < nx; ix++) {
		x = ix + nbx;
		z = iz + srf;
		if (fsrf && iz<2)
		    ua[z][x]=0.0f;
		else
		    ua[z][x] = DXDX(uo,x,z,idx2) 
                	     + DXDZ(uo,x,z,idxz)*s2[z][x]
			     + DZDZ(uo,x,z,idz2)*sv[z][x];
		if (dens) {
		    /* density */
		}
	    }
	}

	/* time stepping */
#ifdef _OPENMP
#pragma omp parallel for schedule(dynamic,ompchunk) \
    private(ix,iz,x,z)				    \
    shared(ua,uo,um,up,vt2,nx,nz,nbx,srf)
#endif
	for     (iz=0; iz<nz; iz++) {
	    for (ix=0; ix<nx; ix++) {
		x = ix + nbx;
		z = iz + srf;
		if (fsrf && iz<2)
		    up[z][x]=0.0f;
		else 
		    up[z][x]=uo[z][x]*2.0f-um[z][x]+ua[z][x]*vt2[z][x];
	    }
	}
	/* absorbing BC */
	abcone2d_absorb(up,uo,wave,abc);
	/* sponge */
	sponge2d_sponge(up,wave,spo);

	/* circulate pointer */
	um=uo;
	uo=up;
	up=um;
    } /* it */

    sf_warning(".");
    exit(0);
}
