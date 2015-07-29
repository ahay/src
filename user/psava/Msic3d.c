/* Local slant stacks I.C. */
/*
  Copyright (C) 2007 Colorado School of Mines
  
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

int main(int argc, char* argv[])
{
    bool verb;
    bool stack;
    float sig;
    
    sf_file Fs=NULL; /* input   */
    sf_file Fr=NULL; /* input   */
    sf_file Fi=NULL; /* output  */

    float ***us=NULL;
    float ***ur=NULL;
    float  **ii=NULL;

    float   *gg=NULL; /*  taper */
    float ***ts=NULL;
    float ***tr=NULL;

    /* cube axes */
    sf_axis at,az,ax,aa,bb,ll,ak; /* cube axes */
    int     nt,nz,nx,na,nb,nl;
    int     it,iz,ix,ia,ib,il;
    int     jt,jz,jx;
    float oa,da,a;
    float ob,db,b;
    float ol,dl,l,lt,lx,lz;
/*    float ot,dt;
      float ox,dx; */

    int   ***kt=NULL,*pkt;
    int   ***kx=NULL,*pkx;
    int   ***kz=NULL,*pkz;

    int   ***ht=NULL,*pht;
    int   ***hx=NULL,*phx;
    int   ***hz=NULL,*phz;

    int   ***gt=NULL,*pgt;
    int   ***gx=NULL,*pgx;
    int   ***gz=NULL,*pgz;

    int     lkt, lkx, lkz;
    int     lgt, lgx, lgz;
    int     lht, lhx, lhz;

    float   lgg;

    int ompchunk; 

    /*------------------------------------------------------------*/
    /* init RSF */
    sf_init(argc,argv);

    if(! sf_getint("ompchunk",&ompchunk)) ompchunk=1;  /* OpenMP data chunk size */
    if(! sf_getbool("verb",   &verb))     verb=false;  /* verbosity flag */
    if(! sf_getbool("stack",  &stack))   stack=false;
    if(! sf_getfloat("sig",   &sig))     sig=1.0;

    Fs = sf_input ("in" );
    Fr = sf_input ("ur" );
    Fi = sf_output("out");

    /* angle axis (in degrees) */
    if(! sf_getint  ("nanga",&na)) na=1;
    if(! sf_getfloat("oanga",&oa)) oa=0.0;
    if(! sf_getfloat("danga",&da)) da=1.0;
    aa=sf_maxa(na,oa,da);
    sf_setlabel(aa,"a"); 
    sf_setunit (aa,"");

    /* angle axis (in degrees) */
    if(! sf_getint  ("nangb",&nb)) nb=1;
    if(! sf_getfloat("oangb",&ob)) ob=0.0;
    if(! sf_getfloat("dangb",&db)) db=1.0;
    bb=sf_maxa(nb,ob,db);
    sf_setlabel(bb,"b"); 
    sf_setunit (bb,"");

    /* length axis (in samples!!!) */
    if(! sf_getint  ("nl",&nl)) nl=1;
    if(! sf_getfloat("dl",&dl)) dl=1.;
    if(! sf_getfloat("ol",&ol)) ol=0.;
    ll=sf_maxa(nl,ol,dl);
    sf_setlabel(ll,"l"); 
    sf_setunit (ll,"");

    /* input axes */
    at = sf_iaxa(Fs,1); sf_setlabel(at,"t"); sf_setunit(at,""); 
    ax = sf_iaxa(Fs,2); sf_setlabel(ax,"x"); sf_setunit(ax,""); 
    az = sf_iaxa(Fs,3); sf_setlabel(az,"z"); sf_setunit(az,"");

    nt = sf_n(at); /* ot=sf_o(at); dt=sf_d(at); */
    nx = sf_n(ax); /* ox=sf_o(ax); dx=sf_d(ax); */
    nz = sf_n(az); 

    if(verb) {
	sf_raxa(az);
	sf_raxa(ax);
	sf_raxa(at);
	sf_raxa(aa);
	sf_raxa(bb);
	sf_raxa(ll);
    }    

    /* set output axes */
    if(stack)
	ak=sf_maxa(1,0,1); 
    else
	ak=sf_maxa(nb,0,1);
    
    /* setup output header */
    sf_oaxa(Fi,ax,1);
    sf_oaxa(Fi,az,2);
    sf_oaxa(Fi,ak,3);

    nz = sf_n(az);
    nx = sf_n(ax);
    nt = sf_n(at);
    
    /*------------------------------------------------------------*/
    /* allocate arrays */
    us=sf_floatalloc3(nt,nx,nz);  /*   source wavefield */
    ur=sf_floatalloc3(nt,nx,nz);  /* receiver wavefield */
    ii=sf_floatalloc2(   nx,nz);  /* image */

    ts=sf_floatalloc3(nt,nx,nz);
    tr=sf_floatalloc3(nt,nx,nz);
    
    gg=sf_floatalloc(2*nl+1);

    kt=sf_intalloc3 (2*nl+1,na,nb);
    kx=sf_intalloc3 (2*nl+1,na,nb);
    kz=sf_intalloc3 (2*nl+1,na,nb);

    ht=sf_intalloc3 (2*nl+1,na,nb);
    hx=sf_intalloc3 (2*nl+1,na,nb);
    hz=sf_intalloc3 (2*nl+1,na,nb);

    gt=sf_intalloc3 (2*nl+1,na,nb);
    gx=sf_intalloc3 (2*nl+1,na,nb);
    gz=sf_intalloc3 (2*nl+1,na,nb);
    
    /*------------------------------------------------------------*/
    /* taper */
    for(il=0;il<2*nl+1;il++) {
	l = ol + (il-nl)*dl;
	l /= (nl/2);
	l /= sig;
	gg[il] = exp(-l*l);
	gg[il]=1.;
    }

    /*------------------------------------------------------------*/
    /* nearest neighbor interpolation */
    for(ib=0;ib<nb;ib++) {
	b  = ob + ib * db;
	b *= SF_PI/180.;
	

	for(ia=0;ia<na;ia++) {
	    a  = oa + ia * da;
	    a *= SF_PI/180.;
	    
	    for(il=0;il<2*nl+1;il++){
		l = ol + (il-nl)*dl;   
		
/*		lt = l*cos(b)*sin(a);*/
/*		lx = l*cos(b)*cos(a);*/
/*		lz = l*sin(b);*/

		lz = l*cos(a) * sin(b);
		lx = l*cos(a) * cos(b);
		lt = l*sin(a);

		kt[ib][ia][il] = (floor)(lt);	    
		kx[ib][ia][il] = (floor)(lx);
		kz[ib][ia][il] = (floor)(lz);
		
		ht[ib][ia][il] = SF_ABS(kt[ib][ia][il]); 
		hx[ib][ia][il] = SF_ABS(kx[ib][ia][il]); 
		hz[ib][ia][il] = SF_ABS(kz[ib][ia][il]); 
		
		gt[ib][ia][il] = nt - ht[ib][ia][il];
		gx[ib][ia][il] = nx - hx[ib][ia][il];
		gz[ib][ia][il] = nz - hz[ib][ia][il];
	    }
	}
    }

    /*------------------------------------------------------------*/
    for(    iz=0; iz<nz; iz++) {
	for(ix=0; ix<nx; ix++) {         /* init image */
	    ii[iz][ix] = 0;
	}
    }
    if(stack)
	sf_floatwrite(ii[0],nx*nz,Fi);      
    else
	for (ib=0; ib<nb; ib++) {
	    sf_floatwrite(ii[0],nx*nz,Fi);      
	}
    sf_seek(Fi,0,SEEK_SET);

    sf_floatread(us[0][0],nt*nx*nz,Fs);  /* read   source wavefield */
    sf_floatread(ur[0][0],nt*nx*nz,Fr);  /* read receiver wavefield */

    /*------------------------------------------------------------*/    
    if(verb) fprintf(stderr,"  b   a\n");
    if(verb) fprintf(stderr,"%3d %3d\n",nb-1,na-1);    
    for(ib=0;ib<nb;ib++) {    

	for(    iz=0; iz<nz; iz++) {
	    for(ix=0; ix<nx; ix++) {         /* init image */
		ii[iz][ix] = 0;
	    }
	}

	for(ia=0;ia<na;ia++) {
	    if(verb) fprintf(stderr,"%3d %3d",ib,ia);

#ifdef _OPENMP
#pragma omp parallel for schedule(dynamic,ompchunk) private(iz,ix,it) shared(nz,nx,nt,ts,tr)
#endif	
	    for(        iz=0; iz<nz; iz++) {
		for(    ix=0; ix<nx; ix++) {
		    for(it=0; it<nt; it++) {	
			ts[iz][ix][it] = 0;
			tr[iz][ix][it] = 0;
		    }     /* t */
		}         /* x */
	    }             /* z */

	    pkt = kt[ib][ia];
	    pkx = kx[ib][ia];
	    pkz = kz[ib][ia];

	    pht = ht[ib][ia];	    
	    phx = hx[ib][ia];
	    phz = hz[ib][ia];

	    pgt = gt[ib][ia];
	    pgx = gx[ib][ia];
	    pgz = gz[ib][ia];

#ifdef _OPENMP
#pragma omp parallel for schedule(dynamic,ompchunk) private(il,iz,ix,it,lkz,lkx,lkt,lhz,lhx,lht,lgz,lgx,lgt,jz,jt,jx,lgg) shared(ia,ib,pkz,pkx,pkt,phz,phx,pht,pgz,pgx,pgt,nl,nz,nx,nt,ts,tr,gg,us,ur)
#endif	
	    for(il=0;il<2*nl+1;il++) {
		lgg=gg[il];
		
		lkt=pkt[il];
		lkx=pkx[il];
		lkz=pkz[il];
		
		lht=pht[il];
		lhx=phx[il];
		lhz=phz[il];

		lgt=pgt[il];
		lgx=pgx[il];
		lgz=pgz[il];
		
		for(        iz=lhz; iz<lgz; iz++) { jz = iz+lkz;			    
		    for(    ix=lhx; ix<lgx; ix++) { jx = ix+lkx;
			for(it=lht; it<lgt; it++) { jt = it+lkt;
			    ts[ iz ][ ix ][ it ] += lgg * us[ jz ][ jx ][ jt ];
			    tr[ iz ][ ix ][ it ] += lgg * ur[ jz ][ jx ][ jt ];
			} /* ht loop */
		    }     /* hx loop */		
		}         /* hz loop */

	    }	          /*  l loop */
	
#ifdef _OPENMP
#pragma omp parallel for schedule(dynamic,ompchunk) private(iz,ix,it) shared(nz,nx,nt,ts,tr,ii)
#endif	
	    for(        iz=0; iz<nz; iz++) {
		for(    ix=0; ix<nx; ix++) {
		    for(it=0; it<nt; it++) {
			ii[iz][ix] += ts[iz][ix][it] * tr[iz][ix][it];
		    }     /* t */
		}         /* x */
	    }             /* z */
	    
	    if(verb) fprintf(stderr,"\b\b\b\b\b\b\b\b\b\b\b\b");
	}                 /* a loop */

	sf_floatwrite(ii[0],nx*nz,Fi);      
	if(stack) sf_seek(Fi,0,SEEK_SET);
    }                     /* b loop */
    
    if(verb) fprintf(stderr,"\n");

    exit (0);
}
