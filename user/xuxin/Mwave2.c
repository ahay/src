/* 2-D acoustic time-extrapolation */

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
#define DZDZ(a,ix,iz,idz2) (C0* a[ix][iz]+\
                            CA*(a[ix][iz-1]+a[ix][iz+1])+\
                            CB*(a[ix][iz-2]+a[ix][iz+2]))*idz2
#define DXDX(a,ix,iz,idx2) (C0* a[ix][iz]+\
                            CA*(a[ix-1][iz]+a[ix+1][iz])+\
                            CB*(a[ix-2][iz]+a[ix+2][iz]))*idx2
#define C1  0.666666666666666666667 /*  2/3 */
#define C2 -0.083333333333333333333 /* -1/12 */
#define DZ(a,ix,iz,idz) (C1*(a[ix][iz+1]-a[ix][iz-1])+\
                         C2*(a[ix][iz+2]-a[ix][iz-2]))*idz
#define DX(a,ix,iz,idx) (C1*(a[ix+1][iz]-a[ix-1][iz])+\
                         C2*(a[ix+2][iz]-a[ix-2][iz]))*idx

int main(int argc, char* argv[])
{
    bool fsrf,wrtdat,dens;
    sf_file Fin=NULL,Fout=NULL,Fvel=NULL,Fden=NULL,Fsou=NULL,Frec=NULL,Fdat=NULL;
    float **vel=NULL,**den=NULL,**rox=NULL,**roz=NULL,**um=NULL,**uo=NULL,**up=NULL,**ua=NULL,**vt2=NULL,**tt=NULL;
    float *ww=NULL,*dd=NULL,*xs=NULL,*zs=NULL,*xr=NULL,*zr=NULL;
    int nx,nz,ix,iz,nbx,nbz,npx,npz,x,z,nt,it,ns,nr,n;
    int jsnap,ompnth,ompath,ompchunk,srf;
    float ox,oz,dx,dz,idx,idz,idx2,idz2,ot,dt;    
    sf_axis ax,az,at,aj;
    wave2dp wave=NULL;
    abcone2dp abc=NULL;
    sponge2dp spo=NULL;
    lint2dp lints=NULL,lintr=NULL;

    sf_init(argc,argv);
    /* files */
    Fin  = sf_input("in");    /* source data */
    Fvel = sf_input("vel");   /* velocity */
    Fsou = sf_input("sou");   /* source coord */
    Fout = sf_output("out");  /* wavefield */
    if (NULL != sf_getstring("den")) {
	dens = true;
	Fden=sf_input("den");   /* density */
    } else {
	dens = false;
    }
    if (NULL != sf_getstring("rec") && NULL != sf_getstring("dat")) {
	wrtdat = true;
	Frec = sf_input("rec");  /* receiver coord */
	Fdat = sf_output("dat"); /* seismogram */ 
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
#ifdef _OPENMP
#pragma omp parallel
    ompath=omp_get_num_threads();
    omp_set_num_threads(ompnth);
    sf_warning("using %d threads of a total of %d",ompnth,ompath);
#endif

    /* axis */
    at = sf_iaxa(Fin,2);
    az = sf_iaxa(Fvel,1);
    ax = sf_iaxa(Fvel,2);
    nt = sf_n(at); ot = sf_o(at); dt = sf_d(at);
    nx = sf_n(ax); ox = sf_o(ax); dx = sf_d(ax); 
    nz = sf_n(az); oz = sf_o(az); dz = sf_d(az);
    if (dens) {
	if (!sf_histint(Fden,"n1",&n) || n != nz) sf_error("Need n1=%d in den",nz);
	if (!sf_histint(Fden,"n2",&n) || n != nx) sf_error("Need n1=%d in den",nx);
    }
    for (n=0, it=0; it<nt; it++)
	if (!(it % jsnap)) n++;
    aj = sf_maxa(n,ot,dt*jsnap);
    sf_setlabel(at,"t"); sf_setunit(at,"");
    sf_setlabel(ax,"x"); sf_setunit(at,"");
    sf_setlabel(az,"z"); sf_setunit(at,"");
    sf_setlabel(aj,"t"); sf_setunit(aj,"");
    sf_oaxa(Fout,az,1);
    sf_oaxa(Fout,ax,2);
    sf_oaxa(Fout,aj,3);
    idx = 1.0f/dx; idx2 = idx*idx;
    idz = 1.0f/dz; idz2 = idz*idz;

    /* grids */
    wave = wave2d_init(ax,az,at,nbx,nbz,fsrf);
    npx = sf_n(wave->px);
    npz = sf_n(wave->pz);
    /* srf = wave->fsrf ? 0 : wave->nbz; */
    srf = fsrf ? 0 : nbz;
    sf_raxa(ax); sf_raxa(wave->px);    
    sf_raxa(az); sf_raxa(wave->pz);
    sf_raxa(at); sf_raxa(aj);

    /* source coord */
    if (!sf_histint(Fin,"n1",&ns))    sf_error("No n1= in input");
    if (!sf_histint(Fsou,"n1",&n) || n != ns)  sf_error("Need n1=%d in sou",ns);
    xs = sf_floatalloc(ns);
    zs = sf_floatalloc(ns);
    sf_floatread(xs,ns,Fsou);
    sf_floatread(zs,ns,Fsou);
    lints = lint2d_init(zs,xs,ns,wave->pz,wave->px);
    sf_warning("ns=%d",ns);
   
    /* receiver coord */
    if (wrtdat) {
	if (!sf_histint(Frec,"n1",&nr))   sf_error("No n1= in rec");
	xr = sf_floatalloc(nr);
	zr = sf_floatalloc(nr);
	sf_floatread(xr,nr,Frec);
	sf_floatread(zr,nr,Frec);
	lintr = lint2d_init(zr,xr,nr,wave->pz,wave->px);
	sf_warning("nr=%d",nr);

	sf_putint(Fdat,"n1",nr);
	sf_putstring(Fdat,"label1","r");
	sf_putstring(Fdat,"unit1","");
	sf_oaxa(Fdat,at,2);
    }    

    /* read velocity */
    tt  = sf_floatalloc2(nz ,nx );
    vel = sf_floatalloc2(npz,npx);
    vt2 = sf_floatalloc2(npz,npx);
    sf_floatread(tt[0],nx*nz,Fvel);
    expand2(tt,vel,wave);
    for     (iz=0; iz<npz; iz++) {
	for (ix=0; ix<npx; ix++) {
	    vt2[ix][iz] = vel[ix][iz]*vel[ix][iz] *dt*dt;
	}
    }

    /* read density */
    if (dens) {
	den=sf_floatalloc2(npz,npx);
	rox=sf_floatalloc2(npz,npx);
	roz=sf_floatalloc2(npz,npx);
	sf_floatread(tt[0],nx*nz,Fden);
	expand2(tt,den,wave);
	for     (iz=0; iz<nz; iz++) {
	    for (ix=0; ix<nx; ix++) {
		x = ix + nbx;
		z = iz + srf;
		rox[x][z] = DX(den,x,z,idx)/den[x][z];
		if (fsrf && iz<2) 
		    roz[x][z]=0.0f;
		else
		    roz[x][z] = DZ(den,x,z,idz)/den[x][z];
	    }
	}
    }

    um =sf_floatalloc2(npz,npx);
    uo =sf_floatalloc2(npz,npx);
    up =sf_floatalloc2(npz,npx);
    ua =sf_floatalloc2(npx,npx);
    ww =sf_floatalloc(ns);
    if (wrtdat)
	dd =sf_floatalloc(nr);

    for (ix=0; ix<npx*npz; ix++) {
	um[0][ix] = 0.0f;
	uo[0][ix] = 0.0f;
	up[0][ix] = 0.0f;
	ua[0][ix] = 0.0f;
    }

    abc=abcone2d_init(vel,wave);
    spo=sponge2d_init(0.00005,wave);

    /* MAIN */
    for (it=0; it<nt; it++) {
	sf_warning("it=%d;",it);

	/* write wavefield */
	if (!(it % jsnap)) {
	    /* for(iz=srf; iz<npz-nbz; iz++) */
	    /* 	sf_floatwrite(&uo[iz][nbx],nx,Fout); */
	    for (ix=nbx; ix < npx-nbx; ix++)
		sf_floatwrite(&uo[ix][srf],nz,Fout);
	}

	/* write seismogram */
	if (wrtdat) {
	    lint2d_extract(uo,dd,lintr);
	    sf_floatwrite(dd,nr,Fdat);
	}

	/* inject source */
	sf_floatread(ww,ns,Fin);
	lint2d_inject(uo,ww,lints);

	/* laplacian */
#ifdef _OPENMP
#pragma omp parallel for schedule(dynamic,ompchunk) \
    private(ix,iz,x,z)				    \
    shared(uo,ua,nx,nz,nbx,srf,idx2,idz2,fsrf,rox,roz,dens)
#endif	
	for     (iz=0; iz<nz; iz++) {
	    for (ix=0; ix<nx; ix++) {
		x = ix + nbx;
		z = iz + srf;
		if (fsrf && iz<2)
		    ua[x][z]=0.0f;
		else
		    ua[x][z] = DXDX(uo,x,z,idx2) + DZDZ(uo,x,z,idz2);
		if (dens) {
		    ua[x][z] -= rox[x][z]*DX(uo,x,z,idx);
                    ua[x][z] -= roz[x][z]*DZ(uo,x,z,idz);
		} 
	    }
	}
	/* time stepping */
#ifdef _OPENMP
#pragma omp parallel for schedule(dynamic,ompchunk) \
    private(ix,iz,x,z)				    \
    shared(ua,uo,um,up,vt2,nx,nz,nbx,srf,fsrf)
#endif
	for     (iz=0; iz<nz; iz++) {
	    for (ix=0; ix<nx; ix++) {
		x = ix + nbx; 
		z = iz + srf;
		if (fsrf && iz<2)
		    up[x][z]=0.0f;
		else
		    up[x][z]=2.0f*uo[x][z]-um[x][z]+vt2[x][z]*ua[x][z];
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
