/* 
 * exploding reflector time-domain acoustic FD modeling
 */
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

int main(int argc, char* argv[])
{
    bool verb; /* verbosity flag */
    bool abc;  /* absorbing boundary conditions flag */
    bool free; /* free surface flag*/
    bool snap; /* wavefield snapshots flag */
    int  jsnap;/* save wavefield every *jsnap* time steps */

    /* I/O files */
    sf_file Fw,Fv,Fs,Fr;
    sf_file Fd,Fo;

    /* cube axes */
    axa at,az,ax,as,ar;
    axa bt,bz,bx;
    int it,iz,ix,is,ir, iop;
    float idx,idz,dt2;

    /* arrays */
    float  *ww,**vv; /* wavelet, velocity */
    pt2d   *ss, *rr; /* source/receiver locations */
    float  *dd;      /* data */

    float *fzs,*fxs,    *fzr,*fxr;
    int   *jzs,*jxs,    *jzr,*jxr;

    float *ws00,*ws01,*ws10,*ws11;
    float *wr00,*wr01,*wr10,*wr11;

    float **um,**uo,**up,**ud,**vp,**tt;
    float  *bzl,*bzh,*bxl,*bxh;  /* boundary */

    int   nop=2;       /* Laplacian operator size */
    float c0, c1, c2;  /* Laplacian operator coefficients */
    float co,c1x,c2x,c1z,c2z;

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

    Fw = sf_input ("in" ); /* wavelet */
    Fv = sf_input ("vel"); /* velocity */
    Fs = sf_input ("sou"); /* sources */
    Fr = sf_input ("rec"); /* receivers */
    Fo = sf_output("wfl"); /* wavefield */
    Fd = sf_output("out"); /* data */

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
    /* expanded domain ( az+2 nz, ax+2 nx ) */
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

	oaxa(Fo,&bz,1);
	oaxa(Fo,&bx,2);
	oaxa(Fo,&bt,3);
    }

    /* setup output data header */
    oaxa(Fd,&ar,1);
    oaxa(Fd,&at,2);

    dt2 =    at.d*at.d;
    idz = 1/(az.d*az.d);
    idx = 1/(ax.d*ax.d);

    /* Laplacian coefficients */
    c0=-30./12.; 
    c1=+16./12.;
    c2=- 1./12.;

    co = c0 * (idx+idz);
    c1x= c1 *  idx;
    c2x= c2 *  idx;
    c1z= c1 *      idz;
    c2z= c2 *      idz;

/*------------------------------------------------------------*/
     
    /* allocate arrays */
    ww=sf_floatalloc (at.n);      sf_floatread(ww   ,at.n     ,Fw);
    vv=sf_floatalloc2(az.n,ax.n); sf_floatread(vv[0],az.n*ax.n,Fv);

    /* allocate source/receiver point arrays */
    ss = (pt2d*) sf_alloc(as.n,sizeof(*ss)); 
    rr = (pt2d*) sf_alloc(ar.n,sizeof(*rr)); 

    readpt2d(Fs,ss,as.n,3); /* read 3 elements (x,z,v) */
    readpt2d(Fr,rr,ar.n,2); /* read 2 elements (x,z)   */

    dd=sf_floatalloc(ar.n);
    jzs=sf_intalloc(as.n); fzs=sf_floatalloc(as.n); 
    jzr=sf_intalloc(ar.n); fzr=sf_floatalloc(ar.n);
    jxs=sf_intalloc(as.n); fxs=sf_floatalloc(as.n);
    jxr=sf_intalloc(ar.n); fxr=sf_floatalloc(ar.n);

    ws00 = sf_floatalloc(as.n); wr00 = sf_floatalloc(ar.n); 
    ws01 = sf_floatalloc(as.n); wr01 = sf_floatalloc(ar.n);
    ws10 = sf_floatalloc(as.n); wr10 = sf_floatalloc(ar.n);
    ws11 = sf_floatalloc(as.n); wr11 = sf_floatalloc(ar.n);
/*------------------------------------------------------------*/

    for (is=0;is<as.n;is++) {

	if(ss[is].z >= bz.o && 
	   ss[is].z <  bz.o + (bz.n-1)*bz.d &&
	   ss[is].x >= bx.o && 
	   ss[is].x < bx.o + (bx.n-1)*bx.d) {
	    
	    jzs[is] = (int)( (ss[is].z-bz.o)/bz.d);
	    fzs[is] =        (ss[is].z-bz.o)/bz.d - jzs[is];	    
	    jxs[is] = (int)( (ss[is].x-bx.o)/bx.d);
	    fxs[is] =        (ss[is].x-bx.o)/bx.d - jxs[is];
	} else {
	    jzs[is] = 0;
	    fzs[is] = 1;
	    ss[is].v= 0;
	}

	ws00[is] = (1-fzs[is])*(1-fxs[is]);
	ws01[is] = (  fzs[is])*(1-fxs[is]);
	ws10[is] = (1-fzs[is])*(  fxs[is]);
	ws11[is] = (  fzs[is])*(  fxs[is]);

    }

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
    
    /* allocate temporary arrays */
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

    /* velocity in the expanded domain (vp=vv^2)*/
    vp=sf_floatalloc2(bz.n,bx.n);
    for (iz=0; iz<az.n; iz++) {
	for (ix=0; ix<ax.n; ix++) {
	    vp[nbx+ix][nbz+iz] = vv[ix][iz] * vv[ix][iz];
	}
    }
    /* fill boundaries */
    for (iz=0; iz<nbz; iz++) {
	for (ix=0; ix<bx.n; ix++) {
	    vp[ix][     iz  ] = vp[ix][     nbz  ];
	    vp[ix][bz.n-iz-1] = vp[ix][bz.n-nbz-1];
	}
    }
    for (iz=0; iz<bz.n; iz++) {
	for (ix=0; ix<nbx; ix++) {
	    vp[     ix  ][iz] = vp[     nbx  ][iz];
	    vp[bx.n-ix-1][iz] = vp[bx.n-nbx-1][iz];
	}
    }

/*------------------------------------------------------------*/

    /* free surface */
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

	/* velocity scale */
	for (ix=0; ix<bx.n; ix++) {
	    for (iz=0; iz<bz.n; iz+=2) {
		ud[ix][iz  ] *= vp[ix][iz  ];
		ud[ix][iz+1] *= vp[ix][iz+1];
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

	/* time step */
	for (ix=0; ix<bx.n; ix++) {
	    for (iz=0; iz<bz.n; iz+=2) {
		up[ix][iz  ] = 2*uo[ix][iz  ] - um[ix][iz  ] + ud[ix][iz  ] * dt2; 
		um[ix][iz  ] =   uo[ix][iz  ];
		uo[ix][iz  ] =   up[ix][iz  ];
		
		up[ix][iz+1] = 2*uo[ix][iz+1] - um[ix][iz+1] + ud[ix][iz+1] * dt2;
		um[ix][iz+1] =   uo[ix][iz+1];
		uo[ix][iz+1] =   up[ix][iz+1];
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
		for (iz=0; iz<bz.n; iz+=2) {
		    uo[ix][iz  ] *= tt[ix][iz  ];
		    um[ix][iz  ] *= tt[ix][iz  ];
		    ud[ix][iz  ] *= tt[ix][iz  ];

		    uo[ix][iz+1] *= tt[ix][iz+1];
		    um[ix][iz+1] *= tt[ix][iz+1];
		    ud[ix][iz+1] *= tt[ix][iz+1];
		}
	    }
	}
	
	/* write wavefield */
	if(snap && it%jsnap==0) {
	    sf_floatwrite(uo[0],bz.n*bx.n,Fo);
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
