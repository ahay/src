/* Time-domain acoustic FD modeling */
/*
  Copyright (C) 2006 Colorado School of Mines
  
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

#define C1  0.66666666666666666666 /*  2/3  */	
#define C2 -0.08333333333333333333 /* -1/12 */
#define D1(a,i2,i1,s) (C2*(a[i2  ][i1+2] - a[i2  ][i1-2]) +  \
                       C1*(a[i2  ][i1+1] - a[i2  ][i1-1])  )*s
#define D2(a,i2,i1,s) (C2*(a[i2+2][i1  ] - a[i2-2][i1  ]) +  \
                       C1*(a[i2+1][i1  ] - a[i2-1][i1  ])  )*s

int main(int argc, char* argv[])
{
    bool verb; /* verbosity flag */
    bool abc;  /* absorbing boundary conditions flag */
    bool free; /* free surface flag*/
    bool snap; /* wavefield snapshots flag */
    bool dens;
    int  jsnap;/* save wavefield every *jsnap* time steps */

    /* I/O files */
    sf_file Fw,Fs,Fr;
    sf_file Fd,Fu;
    sf_file Fv=NULL; /* velocity */
    sf_file Fe=NULL; /* density */

    /* cube axes */
    sf_axis at,az,ax,as,ar;
    int it,iz,ix,is,ir, iop;
    int nt,nz,nx,ns,nr,nz2,nx2;
    float z0,dz,x0,dx,idx,idz,dt,dt2;

    /* arrays */
    pt2d   *ss, *rr; /* source/receiver locations */
    float  *ww=NULL; /* wavelet */
    float  *dd=NULL; /* data */
    float **vv=NULL; /* velocity */
    float **ee=NULL; /* density  */

    float *fzs,*fxs,    *fzr,*fxr;
    int   *jzs,*jxs,    *jzr,*jxr;

    float *ws00,*ws01,*ws10,*ws11;
    float *wr00,*wr01,*wr10,*wr11;

    float **um,**uo,**up,**ud,**vp,**ro,**tt,**ut;
    float  *bzl,*bzh,*bxl,*bxh;  /* boundary */

    int   nop=2;       /* Laplacian operator size */
    float c0, c1, c2;  /* Laplacian operator coefficients */
    float co,c1x,c2x,c1z,c2z;

    int  nbz,nbx; /* boundary size */
    float tz, tx; /* sponge boundary decay coefficients */
    float dp;
    float ws;     /* injected data */

    int ompchunk;  /* OpenMP data chunk size */

/*------------------------------------------------------------*/

    /* init RSF */
    sf_init(argc,argv);

    if(! sf_getint("ompchunk",&ompchunk)) ompchunk=1;

    if(! sf_getbool("verb",&verb)) verb=false;
    if(! sf_getbool( "abc",&abc ))  abc=false;
    if(! sf_getbool("snap",&snap)) snap=false;
    if(! sf_getbool("free",&free)) free=false;
    if(! sf_getbool("dens",&dens)) dens=false;

    Fw = sf_input ("in" ); /* wavelet */
    Fv = sf_input ("vel"); /* velocity */
    Fs = sf_input ("sou"); /* sources */
    Fr = sf_input ("rec"); /* receivers */
    Fu = sf_output("wfl"); /* wavefield */
    Fd = sf_output("out"); /* data */

    if(dens) Fe = sf_input("den"); /* density */

    /* read axes*/
    at=sf_iaxa(Fw,1); sf_setlabel(at,"t"); if(verb) sf_raxa(at); /* time */
    az=sf_iaxa(Fv,1); sf_setlabel(az,"z"); if(verb) sf_raxa(az); /* depth */
    ax=sf_iaxa(Fv,2); sf_setlabel(ax,"x"); if(verb) sf_raxa(ax); /* space */
    as=sf_iaxa(Fs,2); sf_setlabel(as,"s"); if(verb) sf_raxa(as); /* source */
    ar=sf_iaxa(Fr,2); sf_setlabel(ar,"r"); if(verb) sf_raxa(ar); /* receiver */

    nt=sf_n(at); dt=sf_d(at);
    nz=sf_n(az);
    nx=sf_n(ax);
    ns=sf_n(as);
    nr=sf_n(ar);

    /* configure wavefield snapshots */
    if(snap) {
	if(! sf_getint("jsnap",&jsnap)) jsnap=nt;
    }

/*------------------------------------------------------------*/

    /* expand domain for absorbing boundary conditions */
    if(abc) {
	if(! sf_getint("nbz",&nbz)) nbz=nop; if(nbz<nop) nbz=nop;
	if(! sf_getint("nbx",&nbx)) nbx=nop; if(nbx<nop) nbx=nop;
	
	if(! sf_getfloat("tz",&tz)) tz=0.025;
	if(! sf_getfloat("tx",&tx)) tx=0.025;
    } else {
	nbz=nop;
	nbx=nop;
    }
    /* expanded domain ( az+2 nz, ax+2 nx ) */
    nz2=nz+2*nbz; dz=sf_d(az); z0=sf_o(az)-nbz*dz; 
    nx2=nx+2*nbx; dx=sf_d(ax); x0=sf_o(ax)-nbx*dx; 

    sf_setn(az,nz2); sf_seto(az,z0); if(verb) sf_raxa(az);
    sf_setn(ax,nx2); sf_seto(ax,x0); if(verb) sf_raxa(ax);
    
/*------------------------------------------------------------*/

    /* setup output data header */
    sf_oaxa(Fd,ar,1);
    sf_oaxa(Fd,at,2);

    /* setup output wavefield header */
    if(snap) {
	sf_setn(at,nt/jsnap);
	sf_setd(at,dt*jsnap);

	sf_oaxa(Fu,az,1);
	sf_oaxa(Fu,ax,2);
	sf_oaxa(Fu,at,3);
    }

    /* Laplacian coefficients */
    c0=-30./12.; 
    c1=+16./12.;
    c2=- 1./12.;

    dt2 = dt*dt;
    idz = 1/dz;
    idx = 1/dx;

    co = c0 * (idx*idx+idz*idz);
    c1x= c1 *  idx*idx;
    c2x= c2 *  idx*idx;
    c1z= c1 *          idz*idz;
    c2z= c2 *          idz*idz;

/*------------------------------------------------------------*/
     
    /* allocate arrays */
    ww=sf_floatalloc (nt);    sf_floatread(ww   ,nt   ,Fw);
    vv=sf_floatalloc2(nz,nx); sf_floatread(vv[0],nz*nx,Fv);

    ee=sf_floatalloc2(nz,nx); 
    if(dens) {
	sf_floatread(ee[0],nz*nx,Fe); 
    } else {
	for (iz=0; iz<nz; iz++) {
	    for (ix=0; ix<nx; ix++) {
		ee[ix][iz]=1;
	    }
	}
    }
    
    /* allocate source/receiver point arrays */
    ss = (pt2d*) sf_alloc(ns,sizeof(*ss)); 
    rr = (pt2d*) sf_alloc(nr,sizeof(*rr)); 

    pt2dread1(Fs,ss,ns,3); /* read 3 elements (x,z,v) */
    pt2dread1(Fr,rr,nr,2); /* read 2 elements (x,z)   */

    dd=sf_floatalloc(nr);
    for(ir=0;ir<nr;ir++) {
	dd[ir]=0;
    }
    
    jzs=sf_intalloc(ns); fzs=sf_floatalloc(ns); 
    jzr=sf_intalloc(nr); fzr=sf_floatalloc(nr);
    jxs=sf_intalloc(ns); fxs=sf_floatalloc(ns);
    jxr=sf_intalloc(nr); fxr=sf_floatalloc(nr);

    ws00 = sf_floatalloc(ns); wr00 = sf_floatalloc(nr); 
    ws01 = sf_floatalloc(ns); wr01 = sf_floatalloc(nr);
    ws10 = sf_floatalloc(ns); wr10 = sf_floatalloc(nr);
    ws11 = sf_floatalloc(ns); wr11 = sf_floatalloc(nr);
/*------------------------------------------------------------*/

    for (is=0;is<ns;is++) {

	if(ss[is].z >= z0 && 
	   ss[is].z <  z0 + (nz2-1)*dz &&
	   ss[is].x >= x0 && 
	   ss[is].x <  x0 + (nx2-1)*dx) {
	    
	    jzs[is] = (int)( (ss[is].z-z0)/dz);
	    fzs[is] =        (ss[is].z-z0)/dz - jzs[is];	    
	    jxs[is] = (int)( (ss[is].x-x0)/dx);
	    fxs[is] =        (ss[is].x-x0)/dx - jxs[is];
	} else {
	    jzs[is] = 0; jxs[is] = 0;
	    fzs[is] = 1; fxs[is] = 0;
	    ss[is].v= 0;
	}

	ws00[is] = (1-fzs[is])*(1-fxs[is]);
	ws01[is] = (  fzs[is])*(1-fxs[is]);
	ws10[is] = (1-fzs[is])*(  fxs[is]);
	ws11[is] = (  fzs[is])*(  fxs[is]);

    }

    for (ir=0;ir<nr;ir++) {

	if(rr[ir].z >= z0 && 
	   rr[ir].z < z0 + (nz2-1)*dz &&
	   rr[ir].x >= x0 && 
	   rr[ir].x < x0 + (nx2-1)*dx) {
	    
	    jzr[ir] = (int)( (rr[ir].z-z0)/dz);
	    fzr[ir] =        (rr[ir].z-z0)/dz - jzr[ir];
	    jxr[ir] = (int)( (rr[ir].x-x0)/dx);
	    fxr[ir] =        (rr[ir].x-x0)/dx - jxr[ir];

	    rr[ir].v=1;
	} else {
	    jzr[ir] = 0;
	    fzr[ir] = 1;
	    rr[ir].v= 0;
	}

	wr00[ir] = (1-fzr[ir])*(1-fxr[ir]);
	wr01[ir] = (  fzr[ir])*(1-fxr[ir]);
	wr10[ir] = (1-fzr[ir])*(  fxr[ir]);
	wr11[ir] = (  fzr[ir])*(  fxr[ir]);
    }
    
/*------------------------------------------------------------*/
    
    /* allocate temporary arrays */
    um=sf_floatalloc2(nz2,nx2);
    uo=sf_floatalloc2(nz2,nx2);
    up=sf_floatalloc2(nz2,nx2);
    ud=sf_floatalloc2(nz2,nx2);
    tt=sf_floatalloc2(nz2,nx2);

#ifdef _OPENMP
#pragma omp parallel for schedule(dynamic,ompchunk) private(iz,ix) shared(nx2,nz2,um,uo,up,ud,tt)
#endif
    for (iz=0; iz<nz2; iz++) {
	for (ix=0; ix<nx2; ix++) {
	    um[ix][iz]=0;
	    uo[ix][iz]=0;
	    up[ix][iz]=0;
	    ud[ix][iz]=0;
	    tt[ix][iz]=1;
	}
    }

/*------------------------------------------------------------*/

    /* velocity in the expanded domain (vp=vv^2)*/
    vp=sf_floatalloc2(nz2,nx2);
    ro=sf_floatalloc2(nz2,nx2);
    
#ifdef _OPENMP
#pragma omp parallel for schedule(dynamic,ompchunk) private(iz,ix) shared(nz,nx,vp,ro,vv,ee)
#endif
    for (iz=0; iz<nz; iz++) {
	for (ix=0; ix<nx; ix++) {
	    vp[nbx+ix][nbz+iz] = vv[ix][iz] * vv[ix][iz];
	    ro[nbx+ix][nbz+iz] = ee[ix][iz];
	}
    }
    /* fill boundaries */
    for (iz=0; iz<nbz; iz++) {
	for (ix=0; ix<nx2; ix++) {
	    vp[ix][    iz  ] = vp[ix][    nbz  ];
	    vp[ix][nz2-iz-1] = vp[ix][nz2-nbz-1];
	    
	    ro[ix][    iz  ] = ro[ix][    nbz  ];
	    ro[ix][nz2-iz-1] = ro[ix][nz2-nbz-1];
	}
    }
    for (iz=0; iz<nz2; iz++) {
	for (ix=0; ix<nbx; ix++) {
	    vp[    ix  ][iz] = vp[    nbx  ][iz];
	    vp[nx2-ix-1][iz] = vp[nx2-nbx-1][iz];

	    ro[    ix  ][iz] = ro[    nbx  ][iz];
	    ro[nx2-ix-1][iz] = ro[nx2-nbx-1][iz];
	}
    }

/*------------------------------------------------------------*/

    /* free surface */
    if(abc && free) {
	for (iz=0; iz<nbz; iz++) {
	    for (ix=0; ix<nx2; ix++) {
		vp[ix][iz]=0;
	    }
	}
    }

/*------------------------------------------------------------*/

    /* sponge ABC setup */
    if(abc) {
	for (iz=0; iz<nbz; iz++) {
	    for (ix=0; ix<nx2; ix++) {
		tt[ix][    iz  ] = exp( - (tz*(nbz-iz))*(tz*(nbz-iz)) );
		tt[ix][nz2-iz-1] = tt[ix][iz];
	    }
	}
	for (iz=0; iz<nz2; iz++) {
	    for (ix=0; ix<nbx; ix++) {
		tt[    ix  ][iz] = exp( - (tx*(nbx-ix))*(tx*(nbx-ix)) );
		tt[nx2-ix-1][iz] = tt[ix][iz];
	    }
	}
    }

    /* one-way ABC setup */
    bzl=sf_floatalloc(nx2);
    bzh=sf_floatalloc(nx2);
    bxl=sf_floatalloc(nz2);
    bxh=sf_floatalloc(nz2);
    
    for (ix=0;ix<nx2;ix++) {
	dp = vp[ix][    nop  ] *dt/dz; bzl[ix] = (1-dp)/(1+dp);
	dp = vp[ix][nz2-nop-1] *dt/dz; bzh[ix] = (1-dp)/(1+dp);
    }
    for (iz=0;iz<nz2;iz++) {
	dp = vp[    nop  ][iz] *dt/dx; bxl[iz] = (1-dp)/(1+dp);
	dp = vp[nx2-nop-1][iz] *dt/dx; bxh[iz] = (1-dp)/(1+dp);
    }
/*------------------------------------------------------------*/
    /* 
     *  MAIN LOOP
     */
    if(verb) fprintf(stderr,"\n");
    for (it=0; it<nt; it++) {
	if(verb) fprintf(stderr,"\b\b\b\b\b%d",it);
	
	if(dens) { 	/* variable density */
#ifdef _OPENMP
#pragma omp parallel for schedule(dynamic,ompchunk) private(iz,ix) shared(nop,nx2,nz2,ud,uo,ro,co,c1x,c1z,c2x,c2z,idx,idz)
#endif
	    for(    ix=nop; ix<nx2-nop; ix++) {
		for(iz=nop; iz<nz2-nop; iz++) {

		    /* 4th order Laplacian operator */
		    ud[ix][iz] = 
			co * uo[ix  ][iz  ] + 
			c1x*(uo[ix-1][iz  ] + uo[ix+1][iz  ]) +
			c2x*(uo[ix-2][iz  ] + uo[ix+2][iz  ]) +
			c1z*(uo[ix  ][iz-1] + uo[ix  ][iz+1]) +
			c2z*(uo[ix  ][iz-2] + uo[ix  ][iz+2]);	  

		    /* density terms */
		    ud[ix][iz] -= (
			D1(uo,ix,iz,idz) * D1(ro,ix,iz,idz) +
			D2(uo,ix,iz,idx) * D2(ro,ix,iz,idx) ) / ro[ix][iz];
		}
	    }   

	} else {	/* constant density */

#ifdef _OPENMP
#pragma omp parallel for schedule(dynamic,ompchunk) private(iz,ix) shared(nop,nx2,nz2,ud,uo,co,c1x,c1z,c2x,c2z)
#endif
	    for(    ix=nop; ix<nx2-nop; ix++) {
		for(iz=nop; iz<nz2-nop; iz++) {

		    /* 4th order Laplacian operator */
		    ud[ix][iz] = 
			co * uo[ix  ][iz  ] + 
			c1x*(uo[ix-1][iz  ] + uo[ix+1][iz  ]) +
			c2x*(uo[ix-2][iz  ] + uo[ix+2][iz  ]) +
			c1z*(uo[ix  ][iz-1] + uo[ix  ][iz+1]) +
			c2z*(uo[ix  ][iz-2] + uo[ix  ][iz+2]);	  
		}
	    }
	}
	
	/* inject wavelet */
	for (is=0;is<ns;is++) {
	    ws = ww[it] * ss[is].v;
	    ud[ jxs[is]  ][ jzs[is]  ] -= ws * ws00[is];
	    ud[ jxs[is]  ][ jzs[is]+1] -= ws * ws01[is];
	    ud[ jxs[is]+1][ jzs[is]  ] -= ws * ws10[is];
	    ud[ jxs[is]+1][ jzs[is]+1] -= ws * ws11[is];
	}

#ifdef _OPENMP
#pragma omp parallel for schedule(dynamic,ompchunk) private(ix,iz) shared(nx2,nz2,ud,uo,um,up,vp,dt2)
#endif
	for(    ix=0; ix<nx2; ix++) {
	    for(iz=0; iz<nz2; iz++) {

		/* time step and velocity scale*/
		up[ix][iz] = 2*uo[ix][iz] - 
		               um[ix][iz] + 
		               ud[ix][iz] * vp[ix][iz] * dt2; 
	    }
	}
	/* circulate arrays */
	ut=um;
	um=uo;
	uo=up;
	up=ut;


	/* one-way ABC apply */	
	if(abc) {
#ifdef _OPENMP
#pragma omp parallel for schedule(dynamic,ompchunk) private(ix,iz,iop) shared(nx2,nz2,nop,uo,um,bzl,bzh)
#endif
	    for(ix=0;ix<nx2;ix++) {
		for(iop=0;iop<nop;iop++) {
		    iz = nop-iop;
		    uo      [ix][iz  ] 
			= um[ix][iz+1] 
			+(um[ix][iz  ]
			- uo[ix][iz+1]) * bzl[ix];
		    
		    iz = nz2-nop+iop-1;
		    uo      [ix][iz  ] 
			= um[ix][iz-1]
			+(um[ix][iz  ]
			- uo[ix][iz-1]) * bzh[ix];
		}
	    }

#ifdef _OPENMP
#pragma omp parallel for schedule(dynamic,1) private(ix,iz,iop) shared(nx2,nz2,nop,uo,um,bzl,bzh)
#endif
	    for(iop=0;iop<nop;iop++) {
		for(iz=0;iz<nz2;iz++) {
		    ix = nop-iop;
		    uo      [ix  ][iz] 
			= um[ix+1][iz] 
			+(um[ix  ][iz]
			- uo[ix+1][iz]) * bxl[iz];
		    
		    ix = nx2-nop+iop-1;
		    uo      [ix  ][iz] 
			= um[ix-1][iz]
			+(um[ix  ][iz]
			- uo[ix-1][iz]) * bxh[iz];
		}
	    }
	}
	
	/* sponge ABC apply */
	if(abc) {
#ifdef _OPENMP
#pragma omp parallel for schedule(dynamic,ompchunk) private(ix,iz) shared(nx2,nz2,uo,um,ud,tt)
#endif
	    for(    ix=0; ix<nx2; ix++) {
		for(iz=0; iz<nz2; iz++) {
		    uo[ix][iz] *= tt[ix][iz];
		    um[ix][iz] *= tt[ix][iz];
		    ud[ix][iz] *= tt[ix][iz];
		}
	    }
	}
	
	/* write wavefield */
	if(snap && it%jsnap==0) {
	    sf_floatwrite(uo[0],nz2*nx2,Fu);
	}

	/* collect data */
#ifdef _OPENMP
#pragma omp parallel for schedule(dynamic,1) private(ir) shared(dd,rr,uo,jzr,wr00,wr01,wr10,wr11)
#endif
	for (ir=0;ir<nr;ir++) {
	    dd[ir] =
		uo[ jxr[ir]  ][ jzr[ir]  ] * wr00[ir] +
		uo[ jxr[ir]  ][ jzr[ir]+1] * wr01[ir] +
		uo[ jxr[ir]+1][ jzr[ir]  ] * wr10[ir] +
		uo[ jxr[ir]+1][ jzr[ir]+1] * wr11[ir];
	    dd[ir] *= rr[ir].v;
	}
	/* write data */
	sf_floatwrite(dd,nr,Fd);
    }
    if(verb) fprintf(stderr,"\n");

    exit (0);
}
