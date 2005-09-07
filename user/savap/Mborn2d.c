/* 
 * Born modeling
 * pcs 2005
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

    /* cube axes */
    axa at,az,ax,as,ar;
    axa bt,bz,bx;
    int it,iz,ix,is,ir, iop;
    float idx,idz,dt2;

    /* Laplacian */
    int   nop=2;       /* Laplacian operator size */
    float c0, c1, c2;  /* Laplacian operator coefficients */
    float co,c1x,c2x,c1z,c2z;

    int  nbz,nbx; /* boundary size */
    float tz, tx; /* sponge boundary decay coefficients */
    float dp;
    float ws;     /* injected data */

    /* linear interpolation */
    float *fzs,*fxs,    *fzr,*fxr;
    int   *jzs,*jxs,    *jzr,*jxr;

    float *ws00,*ws01,*ws10,*ws11;
    float *wr00,*wr01,*wr10,*wr11;

    /* boundary */
    float  *bzl,*bzh,*bxl,*bxh;

    /* I/O files */
    sf_file Fw,Fs,Fr;
    float  *ww;      /* wavelet */
    pt2d   *ss, *rr; /* source/receiver locations */

    float **tt; /* taper */

    /* background */
    sf_file Bv,Bd,Bu; /* velocity, data, wavefield */
    float **bvv,**bvo;               /* velocity   */
    float  *bdd;                     /* data       */
    float **bum,**buo,**bup,**bud;   /* wavefields */

    /* perturbation */
    sf_file Pv,Pd,Pu;
    float **pvv,**pvo;
    float  *pdd;
    float **pum,**puo,**pup,**pud;

/*------------------------------------------------------------*/

    /* init RSF */
    sf_init(argc,argv);

    if(! sf_getbool("verb",&verb)) verb=false;
    if(! sf_getbool( "abc",&abc ))  abc=false;
    if(! sf_getbool("snap",&snap)) snap=false;
    if(! sf_getbool("free",&free)) free=false;

    Fw = sf_input ("in" ); /* wavelet */
    Fs = sf_input ("sou"); /* sources */
    Fr = sf_input ("rec"); /* receivers */

    Bv = sf_input ("vel"); /* velocity */
    Bu = sf_output("wfl"); /* wavefield */
    Bd = sf_output("out"); /* data */

    Pv = sf_input ("ref"); /* velocity */
    Pu = sf_output("liw"); /* linearized wavefield */
    Pd = sf_output("lid"); /* linearized data */

    /* read axes*/
    iaxa(Fw,&at,1); at.l="t"; if(verb) raxa(at); /* time */
    iaxa(Fs,&as,2); as.l="s"; if(verb) raxa(as); /* sources */
    iaxa(Fr,&ar,2); ar.l="r"; if(verb) raxa(ar); /* receivers */

    iaxa(Bv,&az,1); az.l="z"; if(verb) raxa(az); /* depth */
    iaxa(Bv,&ax,2); ax.l="x"; if(verb) raxa(ax); /* space */

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

	oaxa(Bu,&bz,1);
	oaxa(Bu,&bx,2);
	oaxa(Bu,&bt,3);

	oaxa(Pu,&bz,1);
	oaxa(Pu,&bx,2);
	oaxa(Pu,&bt,3);
    }

    /* setup output data header */
    oaxa(Bd,&ar,1);
    oaxa(Bd,&at,2);

    oaxa(Pd,&ar,1);
    oaxa(Pd,&at,2);

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

    bvv=sf_floatalloc2(az.n,ax.n); sf_floatread(bvv[0],az.n*ax.n,Bv);
    pvv=sf_floatalloc2(az.n,ax.n); sf_floatread(pvv[0],az.n*ax.n,Pv);

    /* allocate source/receiver point arrays */
    ss = (pt2d*) sf_alloc(as.n,sizeof(*ss)); 
    rr = (pt2d*) sf_alloc(ar.n,sizeof(*rr)); 

    pt2dread1(Fs,ss,as.n,3); /* read 3 elements (x,z,v) */
    pt2dread1(Fr,rr,ar.n,2); /* read 2 elements (x,z)   */

    bdd=sf_floatalloc(ar.n);
    pdd=sf_floatalloc(ar.n);

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
    bum=sf_floatalloc2(bz.n,bx.n);
    buo=sf_floatalloc2(bz.n,bx.n);
    bup=sf_floatalloc2(bz.n,bx.n);
    bud=sf_floatalloc2(bz.n,bx.n);

    pum=sf_floatalloc2(bz.n,bx.n);
    puo=sf_floatalloc2(bz.n,bx.n);
    pup=sf_floatalloc2(bz.n,bx.n);
    pud=sf_floatalloc2(bz.n,bx.n);

    tt=sf_floatalloc2(bz.n,bx.n);

    for (iz=0; iz<bz.n; iz++) {
	for (ix=0; ix<bx.n; ix++) {
	    bum[ix][iz]=pum[ix][iz]=0;
	    buo[ix][iz]=puo[ix][iz]=0;
	    bup[ix][iz]=pup[ix][iz]=0;
	    bud[ix][iz]=pud[ix][iz]=0;
	    tt[ix][iz]=1;
	}
    }

/*------------------------------------------------------------*/

    /* velocity in the expanded domain (vo=vv^2)*/
    bvo=sf_floatalloc2(bz.n,bx.n);
    pvo=sf_floatalloc2(bz.n,bx.n);

    for (iz=0; iz<az.n; iz++) {
	for (ix=0; ix<ax.n; ix++) {
	    bvo[nbx+ix][nbz+iz] = bvv[ix][iz] * bvv[ix][iz];
	    pvo[nbx+ix][nbz+iz] = pvv[ix][iz];
	}
    }
    /* fill boundaries */
    for (iz=0; iz<nbz; iz++) {
	for (ix=0; ix<bx.n; ix++) {
	    bvo[ix][     iz  ] = bvo[ix][     nbz  ];
	    bvo[ix][bz.n-iz-1] = bvo[ix][bz.n-nbz-1];

	    pvo[ix][     iz  ] = pvo[ix][     nbz  ];
	    pvo[ix][bz.n-iz-1] = pvo[ix][bz.n-nbz-1];
	}
    }
    for (iz=0; iz<bz.n; iz++) {
	for (ix=0; ix<nbx; ix++) {
	    bvo[     ix  ][iz] = bvo[     nbx  ][iz];
	    bvo[bx.n-ix-1][iz] = bvo[bx.n-nbx-1][iz];

	    pvo[     ix  ][iz] = pvo[     nbx  ][iz];
	    pvo[bx.n-ix-1][iz] = pvo[bx.n-nbx-1][iz];
	}
    }
 
/*------------------------------------------------------------*/

    /* free surface */
    if(abc && free) {
	for (iz=0; iz<nbz; iz++) {
	    for (ix=0; ix<bx.n; ix++) {
		bvo[ix][iz]=0;
		pvo[ix][iz]=0;
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
	dp = bvo[ix][     nop  ] *at.d/bz.d; bzl[ix] = (1-dp)/(1+dp);
	dp = bvo[ix][bz.n-nop-1] *at.d/bz.d; bzh[ix] = (1-dp)/(1+dp);
    }
    for (iz=0;iz<bz.n;iz++) {
	dp = bvo[     nop  ][iz] *at.d/bx.d; bxl[iz] = (1-dp)/(1+dp);
	dp = bvo[bx.n-nop-1][iz] *at.d/bx.d; bxh[iz] = (1-dp)/(1+dp);
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
		bud[ix][iz] = 
		    co * buo[ix  ][iz  ] + 
		    c1x*(buo[ix-1][iz  ] + buo[ix+1][iz  ]) +
		    c2x*(buo[ix-2][iz  ] + buo[ix+2][iz  ]) +
		    c1z*(buo[ix  ][iz-1] + buo[ix  ][iz+1]) +
		    c2z*(buo[ix  ][iz-2] + buo[ix  ][iz+2]);	  

		pud[ix][iz] = 
		    co * puo[ix  ][iz  ] + 
		    c1x*(puo[ix-1][iz  ] + puo[ix+1][iz  ]) +
		    c2x*(puo[ix-2][iz  ] + puo[ix+2][iz  ]) +
		    c1z*(puo[ix  ][iz-1] + puo[ix  ][iz+1]) +
		    c2z*(puo[ix  ][iz-2] + puo[ix  ][iz+2]);	 
	    }
	}

	/* inject source */
	for (is=0;is<as.n;is++) {
	    ws = ww[it] * ss[is].v;
	    bud[ jxs[is]  ][ jzs[is]  ] -= ws * ws00[is];
	    bud[ jxs[is]  ][ jzs[is]+1] -= ws * ws01[is];
	    bud[ jxs[is]+1][ jzs[is]  ] -= ws * ws10[is];
	    bud[ jxs[is]+1][ jzs[is]+1] -= ws * ws11[is];
	}

	for (ix=0; ix<bx.n; ix++) {
	    for (iz=0; iz<bz.n; iz++) {
		pud[ix][iz] -= bud[ix][iz] * 2*pvo[ix][iz];
	    }
	}

	/* velocity scale */
	for (ix=0; ix<bx.n; ix++) {
	    for (iz=0; iz<bz.n; iz++) {
		bud[ix][iz] *= bvo[ix][iz];
		pud[ix][iz] *= bvo[ix][iz];
	    }
	}
	
	/* time step */
	for (ix=0; ix<bx.n; ix++) {
	    for (iz=0; iz<bz.n; iz++) {
		bup[ix][iz] = 2*buo[ix][iz] - bum[ix][iz] + bud[ix][iz] * dt2; 
		bum[ix][iz] =   buo[ix][iz];
		buo[ix][iz] =   bup[ix][iz];

		pup[ix][iz] = 2*puo[ix][iz] - pum[ix][iz] + pud[ix][iz] * dt2; 
		pum[ix][iz] =   puo[ix][iz];
		puo[ix][iz] =   pup[ix][iz];
	    }
	}
	
	/* one-way ABC apply */
	if(abc) {
	    for(ix=0;ix<bx.n;ix++) {
		for(iop=0;iop<nop;iop++) {
		    iz = nop-iop;
		    buo      [ix][iz  ] 
			= bum[ix][iz+1] 
			+(bum[ix][iz  ]
			- buo[ix][iz+1]) * bzl[ix];
		    puo      [ix][iz  ] 
			= pum[ix][iz+1] 
			+(pum[ix][iz  ]
			- puo[ix][iz+1]) * bzl[ix];
		    
		    iz = bz.n-nop+iop-1;
		    buo      [ix][iz  ] 
			= bum[ix][iz-1]
			+(bum[ix][iz  ]
			- buo[ix][iz-1]) * bzh[ix];
		    puo      [ix][iz  ] 
			= pum[ix][iz-1]
			+(pum[ix][iz  ]
			- puo[ix][iz-1]) * bzh[ix];
		}
	    }

	    for(iop=0;iop<nop;iop++) {
		for(iz=0;iz<bz.n;iz++) {
		    ix = nop-iop;
		    buo      [ix  ][iz] 
			= bum[ix+1][iz] 
			+(bum[ix  ][iz]
			- buo[ix+1][iz]) * bxl[iz];
		    puo      [ix  ][iz] 
			= pum[ix+1][iz] 
			+(pum[ix  ][iz]
			- puo[ix+1][iz]) * bxl[iz];
		    
		    ix = bx.n-nop+iop-1;
		    buo      [ix  ][iz] 
			= bum[ix-1][iz]
			+(bum[ix  ][iz]
			- buo[ix-1][iz]) * bxh[iz];
		    puo      [ix  ][iz] 
			= pum[ix-1][iz]
			+(pum[ix  ][iz]
			- puo[ix-1][iz]) * bxh[iz];
		}
	    }
	}
	
	/* sponge ABC apply */
	if(abc) {
	    for (ix=0; ix<bx.n; ix++) {
		for (iz=0; iz<bz.n; iz++) {
		    buo[ix][iz] *= tt[ix][iz];
		    bum[ix][iz] *= tt[ix][iz];
		    bud[ix][iz] *= tt[ix][iz];

		    puo[ix][iz] *= tt[ix][iz];
		    pum[ix][iz] *= tt[ix][iz];
		    bud[ix][iz] *= tt[ix][iz];
		}
	    }
	}
	
	/* write wavefield */
	if(snap && it%jsnap==0) {
	    sf_floatwrite(buo[0],bz.n*bx.n,Bu);
	    sf_floatwrite(puo[0],bz.n*bx.n,Pu);
	}

	/* write data */
	for (ir=0;ir<ar.n;ir++) {
	    bdd[ir] =
		buo[ jxr[ir]  ][ jzr[ir]  ] * wr00[ir] +
		buo[ jxr[ir]  ][ jzr[ir]+1] * wr01[ir] +
		buo[ jxr[ir]+1][ jzr[ir]  ] * wr10[ir] +
		buo[ jxr[ir]+1][ jzr[ir]+1] * wr11[ir];
	    bdd[ir] *= rr[ir].v;

	    pdd[ir] =
		puo[ jxr[ir]  ][ jzr[ir]  ] * wr00[ir] +
		puo[ jxr[ir]  ][ jzr[ir]+1] * wr01[ir] +
		puo[ jxr[ir]+1][ jzr[ir]  ] * wr10[ir] +
		puo[ jxr[ir]+1][ jzr[ir]+1] * wr11[ir];
	    pdd[ir] *= rr[ir].v;
	}
	sf_floatwrite(bdd,ar.n,Bd);
	sf_floatwrite(pdd,ar.n,Pd);

    }
    if(verb) fprintf(stderr,"\n");    

    exit (0);
}
