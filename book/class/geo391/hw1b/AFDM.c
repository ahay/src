/* Time-domain acoustic FD modeling */
/*
  Copyright (C) 2004 University of Texas at Austin
  
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

#define C1  0.66666666666666666666 //  2/3
#define C2 -0.08333333333333333333 // -1/12
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
    bool dens; /* variable density flag */

    int  jsnap;/* save wavefield every *jsnap* time steps */

    /* I/O files */
    sf_file Fw=NULL; /* wavelet   */
    sf_file Fs=NULL; /* sources   */
    sf_file Fr=NULL; /* receivers */
    sf_file Fd=NULL; /* data      */
    sf_file Fv=NULL; /* velocity  */
    sf_file Fe=NULL; /* density   */
    sf_file Fu=NULL; /* wavefield */

    /* I/O arrays */
    float  *ww=NULL; /* wavelet   */
    pt2d   *ss=NULL; /* sources   */
    pt2d   *rr=NULL; /* receivers */
    float  *dd=NULL; /* data      */

    float **vv=NULL; /* velocity  */
    float **vp=NULL; /* velocity in expanded domain */

    float **ee=NULL; /* density   */
    float **ro=NULL; /* density in expanded domain */

    float **um,**uo,**up,**ud; /* wavefield */
    /* 
       um = U @ t-1
       uo = U @ t
       up = U @ t+1
    */

    float **tt=NULL; /* temporary */

    /* cube axes */
    axa at,az,ax,as,ar;
    axa bt,bz,bx;
    int it,iz,ix,is,ir, iop;
    float idx,idz,dt2; /* 1/dx, 1/dz, 1/dt2*/

    /* linear interpolation weights/indices */
    float *fzs,*fxs,    *fzr,*fxr;
    int   *jzs,*jxs,    *jzr,*jxr;
    float *ws00,*ws01,*ws10,*ws11;
    float *wr00,*wr01,*wr10,*wr11;

    /* Laplacian */
    int   nop=2;       /* Laplacian operator size */
    float c0, c1, c2;  /* Laplacian operator coefficients */
    float co,c1x,c2x,c1z,c2z;

    /* boundary arrays */
    float  *bzl,*bzh,*bxl,*bxh;

    int  nbz,nbx; /* boundary size */
    float tz, tx; /* sponge boundary decay coefficients */
    float dp;
    float ws;     /* injected data */

/*------------------------------------------------------------*/

    /* init RSF */
    sf_init(argc,argv);

    if(! sf_getbool("verb",&verb)) verb=false;
    if(! sf_getbool( "abc",&abc ))  abc=false;
    if(! sf_getbool("snap",&snap)) snap=false;
    if(! sf_getbool("free",&free)) free=false;
    if(! sf_getbool("dens",&dens)) dens=false;

    ;        Fw = sf_input ("in" ); /* wavelet   */
    ;        Fv = sf_input ("vel"); /* velocity  */
    ;        Fs = sf_input ("sou"); /* sources   */
    ;        Fr = sf_input ("rec"); /* receivers */
    ;        Fu = sf_output("wfl"); /* wavefield */
    ;        Fd = sf_output("out"); /* data      */
    if(dens) Fe = sf_input ("den"); /* density   */

    /* read axes*/
    iaxa(Fw,&at,1); at.l="t"; if(verb) raxa(at); /* time */
    iaxa(Fv,&az,1); az.l="z"; if(verb) raxa(az); /* depth */
    iaxa(Fv,&ax,2); ax.l="x"; if(verb) raxa(ax); /* space */
    iaxa(Fs,&as,2); as.l="s"; if(verb) raxa(as); /* sources */
    iaxa(Fr,&ar,2); ar.l="r"; if(verb) raxa(ar); /* receivers */

    /* configure wavefield snapshots */
    if(snap) {
	if(! sf_getint("jsnap",&jsnap)) jsnap=at.n;
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
    /* expanded domain ( az+2 nbz, ax+2 nbx ) */
    bz.n=az.n+2*nbz; bz.o=az.o-nbz*az.d; bz.d=az.d; bz.l="z";
    bx.n=ax.n+2*nbx; bx.o=ax.o-nbx*ax.d; bx.d=ax.d; bx.l="x";

    if(verb) raxa(bz);
    if(verb) raxa(bx);
/*------------------------------------------------------------*/

    /* setup output wavefield header */
    if(snap) {
	bt.n = at.n / jsnap;
	bt.o = at.o;
	bt.d = at.d * jsnap;
	bt.l = "t";

	oaxa(Fu,&bz,1);
	oaxa(Fu,&bx,2);
	oaxa(Fu,&bt,3);
    }

    /* setup output data header */
    oaxa(Fd,&ar,1);
    oaxa(Fd,&at,2);

/*------------------------------------------------------------*/

    /* Laplacian coefficients */
    c0=-30./12.; 
    c1=+16./12.;
    c2=- 1./12.;

    dt2 =    at.d*at.d;
    idz = 1/az.d;
    idx = 1/ax.d;

    co = c0 * (idx*idx+idz*idz);
    c1x= c1 *  idx*idx;
    c2x= c2 *  idx*idx;
    c1z= c1 *          idz*idz;
    c2z= c2 *          idz*idz;

/*------------------------------------------------------------*/
     
    /* allocate arrays */
    ww=sf_floatalloc (at.n);      sf_floatread(ww   ,at.n     ,Fw);
    vv=sf_floatalloc2(az.n,ax.n); sf_floatread(vv[0],az.n*ax.n,Fv);

    ee=sf_floatalloc2(az.n,ax.n); /* input density */
    if(dens) {
	sf_floatread(ee[0],az.n*ax.n,Fe); 
    } else {
	for (iz=0; iz<az.n; iz++) {
	    for (ix=0; ix<ax.n; ix++) {
		ee[ix][iz]=1;
	    }
	}
    }

    /* allocate source/receiver point arrays */
    ss = (pt2d*) sf_alloc(as.n,sizeof(*ss)); 
    rr = (pt2d*) sf_alloc(ar.n,sizeof(*rr)); 

    pt2dread1(Fs,ss,as.n,3); /* read 3 elements (x,z,v) */
    pt2dread1(Fr,rr,ar.n,2); /* read 2 elements (x,z)   */

/*------------------------------------------------------------*/
    /* linear interpolation */

    dd=sf_floatalloc(ar.n);
    jzs=sf_intalloc(as.n); fzs=sf_floatalloc(as.n); 
    jzr=sf_intalloc(ar.n); fzr=sf_floatalloc(ar.n);
    jxs=sf_intalloc(as.n); fxs=sf_floatalloc(as.n);
    jxr=sf_intalloc(ar.n); fxr=sf_floatalloc(ar.n);

    ws00 = sf_floatalloc(as.n); wr00 = sf_floatalloc(ar.n); 
    ws01 = sf_floatalloc(as.n); wr01 = sf_floatalloc(ar.n);
    ws10 = sf_floatalloc(as.n); wr10 = sf_floatalloc(ar.n);
    ws11 = sf_floatalloc(as.n); wr11 = sf_floatalloc(ar.n);

    /* precompute sources coefficients */
    for (is=0;is<as.n;is++) {

	if(ss[is].z >= bz.o && 
	   ss[is].z <  bz.o + (bz.n-1)*bz.d &&
	   ss[is].x >= bx.o && 
	   ss[is].x <  bx.o + (bx.n-1)*bx.d) {
	    
	    jzs[is] = (int)( (ss[is].z-bz.o)/bz.d);
	    fzs[is] =        (ss[is].z-bz.o)/bz.d - jzs[is];	    
	    jxs[is] = (int)( (ss[is].x-bx.o)/bx.d);
	    fxs[is] =        (ss[is].x-bx.o)/bx.d - jxs[is];
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

    /* precompute receivers coefficients */
    for (ir=0;ir<ar.n;ir++) {

	if(rr[ir].z >= bz.o && 
	   rr[ir].z < bz.o + (bz.n-1)*bz.d &&
	   rr[ir].x >= bx.o && 
	   rr[ir].x < bx.o + (bx.n-1)*bx.d) {
	    
	    jzr[ir] = (int)( (rr[ir].z-bz.o)/bz.d);
	    fzr[ir] =        (rr[ir].z-bz.o)/bz.d - jzr[ir];
	    jxr[ir] = (int)( (rr[ir].x-bx.o)/bx.d);
	    fxr[ir] =        (rr[ir].x-bx.o)/bx.d - jxr[ir];

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
    
    /* allocate wavefield/temporary arrays */
    um=sf_floatalloc2(bz.n,bx.n);
    uo=sf_floatalloc2(bz.n,bx.n);
    up=sf_floatalloc2(bz.n,bx.n);
    ud=sf_floatalloc2(bz.n,bx.n);
    tt=sf_floatalloc2(bz.n,bx.n);

    for (iz=0; iz<bz.n; iz++) {
	for (ix=0; ix<bx.n; ix++) {
	    um[ix][iz]=0;
	    uo[ix][iz]=0;
	    up[ix][iz]=0;
	    ud[ix][iz]=0;
	    tt[ix][iz]=1;
	}
    }

/*------------------------------------------------------------*/

    /* velocity/density in expanded domain */
    vp=sf_floatalloc2(bz.n,bx.n);
    ro=sf_floatalloc2(bz.n,bx.n);
    
    for (iz=0; iz<az.n; iz++) {
	for (ix=0; ix<ax.n; ix++) {
	    vp[nbx+ix][nbz+iz] = vv[ix][iz] * vv[ix][iz]; /* vp = vv^2 */
	    ro[nbx+ix][nbz+iz] = ee[ix][iz];              /* density   */
	}
    }
    /* fill boundaries of expanded domain */
    for (iz=0; iz<nbz; iz++) {
	for (ix=0; ix<bx.n; ix++) {
	    vp[ix][     iz  ] = vp[ix][     nbz  ];
	    vp[ix][bz.n-iz-1] = vp[ix][bz.n-nbz-1];

	    ro[ix][     iz  ] = ro[ix][     nbz  ];
	    ro[ix][bz.n-iz-1] = ro[ix][bz.n-nbz-1];
	}
    }
    for (iz=0; iz<bz.n; iz++) {
	for (ix=0; ix<nbx; ix++) {
	    vp[     ix  ][iz] = vp[     nbx  ][iz];
	    vp[bx.n-ix-1][iz] = vp[bx.n-nbx-1][iz];

	    ro[     ix  ][iz] = ro[     nbx  ][iz];
	    ro[bx.n-ix-1][iz] = ro[bx.n-nbx-1][iz];
	}
    }

/*------------------------------------------------------------*/

    /* force free surface if free=y */
    if(abc && free) {
	for (iz=0; iz<nbz; iz++) {
	    for (ix=0; ix<bx.n; ix++) {
		vp[ix][iz]=0;
	    }
	}
    }

/*------------------------------------------------------------*/

    /* sponge ABC setup */
    if(abc) {
	for (iz=0; iz<nbz; iz++) {
	    for (ix=0; ix<bx.n; ix++) {
		tt[ix][     iz  ] = exp( - (tz*(nbz-iz))*(tz*(nbz-iz)) );
		tt[ix][bz.n-iz-1] = tt[ix][iz];
	    }
	}
	for (iz=0; iz<bz.n; iz++) {
	    for (ix=0; ix<nbx; ix++) {
		tt[     ix  ][iz] = exp( - (tx*(nbx-ix))*(tx*(nbx-ix)) );
		tt[bx.n-ix-1][iz] = tt[ix][iz];
	    }
	}
    }

    /* one-way ABC setup */
    bzl=sf_floatalloc(bx.n);
    bzh=sf_floatalloc(bx.n);
    bxl=sf_floatalloc(bz.n);
    bxh=sf_floatalloc(bz.n);
    
    for (ix=0;ix<bx.n;ix++) {
	dp = vp[ix][     nop  ] *at.d/bz.d; bzl[ix] = (1-dp)/(1+dp);
	dp = vp[ix][bz.n-nop-1] *at.d/bz.d; bzh[ix] = (1-dp)/(1+dp);
    }
    for (iz=0;iz<bz.n;iz++) {
	dp = vp[     nop  ][iz] *at.d/bx.d; bxl[iz] = (1-dp)/(1+dp);
	dp = vp[bx.n-nop-1][iz] *at.d/bx.d; bxh[iz] = (1-dp)/(1+dp);
    }
/*------------------------------------------------------------*/

    /* 
     *  MAIN LOOP
     */
    if(verb) fprintf(stderr,"\n");
    for (it=0; it<at.n; it++) {
	if(verb) fprintf(stderr,"\b\b\b\b\b%d",it);
	
	/* 4th order Laplacian operator */
	for (ix=nop; ix<bx.n-nop; ix++) {
	    for (iz=nop; iz<bz.n-nop; iz++) {
		ud[ix][iz] = 
		    co * uo[ix  ][iz  ] + 
		    c1x*(uo[ix-1][iz  ] + uo[ix+1][iz  ]) +
		    c2x*(uo[ix-2][iz  ] + uo[ix+2][iz  ]) +
		    c1z*(uo[ix  ][iz-1] + uo[ix  ][iz+1]) +
		    c2z*(uo[ix  ][iz-2] + uo[ix  ][iz+2]);	  
	    }
	}

	/* variable density */
	if(dens) {
	    for (ix=nop; ix<bx.n-nop; ix++) {
		for (iz=nop; iz<bz.n-nop; iz++) {
		    ud[ix][iz] -= 
			D1(uo,ix,iz,idz) * 
			D1(ro,ix,iz,idz)/ro[ix][iz];

		    ud[ix][iz] -= 
			D2(uo,ix,iz,idx) * 
			D2(ro,ix,iz,idx)/ro[ix][iz];
		}
	    }   
	}

	/* inject wavelet */
	for (is=0;is<as.n;is++) {
	    ws = ww[it] * ss[is].v;
	    ud[ jxs[is]  ][ jzs[is]  ] -= ws * ws00[is];
	    ud[ jxs[is]  ][ jzs[is]+1] -= ws * ws01[is];
	    ud[ jxs[is]+1][ jzs[is]  ] -= ws * ws10[is];
	    ud[ jxs[is]+1][ jzs[is]+1] -= ws * ws11[is];
	}

	/* velocity scale */
	for (ix=0; ix<bx.n; ix++) {
	    for (iz=0; iz<bz.n; iz++) {
		ud[ix][iz] *= vp[ix][iz];
	    }
	}
	
	/* time step */
	for (ix=0; ix<bx.n; ix++) {
	    for (iz=0; iz<bz.n; iz++) {
		up[ix][iz] = 2*uo[ix][iz] - um[ix][iz] + ud[ix][iz] * dt2; 
		um[ix][iz] =   uo[ix][iz];
		uo[ix][iz] =   up[ix][iz];
	    }
	}
	
	/* one-way ABC apply */
	if(abc) {
	    for(ix=0;ix<bx.n;ix++) {
		for(iop=0;iop<nop;iop++) {
		    iz = nop-iop;
		    uo      [ix][iz  ] 
			= um[ix][iz+1] 
			+(um[ix][iz  ]
			- uo[ix][iz+1]) * bzl[ix];
		    
		    iz = bz.n-nop+iop-1;
		    uo      [ix][iz  ] 
			= um[ix][iz-1]
			+(um[ix][iz  ]
			- uo[ix][iz-1]) * bzh[ix];
		}
	    }

	    for(iop=0;iop<nop;iop++) {
		for(iz=0;iz<bz.n;iz++) {
		    ix = nop-iop;
		    uo      [ix  ][iz] 
			= um[ix+1][iz] 
			+(um[ix  ][iz]
			- uo[ix+1][iz]) * bxl[iz];
		    
		    ix = bx.n-nop+iop-1;
		    uo      [ix  ][iz] 
			= um[ix-1][iz]
			+(um[ix  ][iz]
			- uo[ix-1][iz]) * bxh[iz];
		}
	    }
	}
	
	/* sponge ABC apply */
	if(abc) {
	    for (ix=0; ix<bx.n; ix++) {
		for (iz=0; iz<bz.n; iz++) {
		    uo[ix][iz] *= tt[ix][iz];
		    um[ix][iz] *= tt[ix][iz];
		    ud[ix][iz] *= tt[ix][iz];
		}
	    }
	}
	
	/* write wavefield */
	if(snap && it%jsnap==0) {
	    sf_floatwrite(uo[0],bz.n*bx.n,Fu);
	}

	/* write data */
	for (ir=0;ir<ar.n;ir++) {
	    dd[ir] =
		uo[ jxr[ir]  ][ jzr[ir]  ] * wr00[ir] +
		uo[ jxr[ir]  ][ jzr[ir]+1] * wr01[ir] +
		uo[ jxr[ir]+1][ jzr[ir]  ] * wr10[ir] +
		uo[ jxr[ir]+1][ jzr[ir]+1] * wr11[ir];
	    dd[ir] *= rr[ir].v;
	}
	sf_floatwrite(dd,ar.n,Fd);

    }
    if(verb) fprintf(stderr,"\n");    

    exit (0);
}
