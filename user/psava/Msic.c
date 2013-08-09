/* Local slant stacks I.C. */
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

int main(int argc, char* argv[])
{
    bool verb;
    float sig;
    
    sf_file Fs=NULL; /* input   */
    sf_file Fr=NULL; /* input   */
    sf_file Fi=NULL; /* output  */

    float ***us=NULL;
    float ***ur=NULL;
    float   *ii=NULL;

    float   *gg=NULL; /*  taper */
    float ***ww=NULL; /* weight */
    int   ***kt=NULL; /* slant-stack index on axis 1*/
    int   ***kx=NULL; /* slant-stack index on axis 2*/

    float **ts=NULL;
    float **tr=NULL;

    /* cube axes */
    sf_axis at,az,ax,aa,ll,ak; /* cube axes */
    int     nt,nz,nx,na,nl,nb;
    int        it,ix,ia,il,ib;
    int        jt,jx;
    int        mt,mx;
    int        ht,hx;
    int        gt,gx;
    float oa,da,a;
    float ol, dl,l,lt,lx;
    float ft;
    float fx;
    float wo;

    int ic;
    int ompchunk; 

/*------------------------------------------------------------*/

    /* init RSF */
    sf_init(argc,argv);

    if(! sf_getint("ompchunk",&ompchunk)) ompchunk=1;  /* OpenMP data chunk size */
    if(! sf_getbool("verb",   &verb))     verb=false;  /* verbosity flag */
    if(! sf_getint("nbuf",    &nb))       nb=1;        /* buffer size */
    if(! sf_getfloat("sig",   &sig))     sig=1.0;

    Fs = sf_input ("in" );
    Fr = sf_input ("ur" );
    Fi = sf_output("out");

    /* angle axis (in degrees) */
    if(! sf_getint  ("na",&na)) na=1;
    if(! sf_getfloat("oa",&oa)) oa=0.0;
    if(! sf_getfloat("da",&da)) da=1.0;
    aa=sf_maxa(na,oa,da);
    sf_setlabel(aa,""); 
    sf_setunit (aa,"");

    /* length axis (in samples) */
    if(! sf_getint  ("nl",&nl)) nl=1;
    if(! sf_getfloat("dl",&dl)) dl=1.;
    if(! sf_getfloat("ol",&ol)) ol=0.;
    ll=sf_maxa(nl,ol,dl);
    sf_setlabel(ll,""); 
    sf_setunit (ll,"");

    /* input axes */
    at = sf_iaxa(Fs,1); sf_setlabel(at,"t"); sf_setunit(at,""); 
    ax = sf_iaxa(Fs,2); sf_setlabel(ax,"x"); sf_setunit(ax,""); 
    az = sf_iaxa(Fs,3); sf_setlabel(az,"z"); sf_setunit(az,"");

    nt = sf_n(at);  
    nx = sf_n(ax);  
    nz = sf_n(az); 

    if(verb) {
	sf_raxa(az);
	sf_raxa(ax);
	sf_raxa(at);
	sf_raxa(aa);
	sf_raxa(ll);
    }    

    /* set output axes */
    ak=sf_maxa(1,0,1); 

    /* setup output header */
    sf_oaxa(Fi,ax,1);
    sf_oaxa(Fi,az,2);
    sf_oaxa(Fi,ak,3);

    nz = sf_n(az); nb = SF_MIN(nb,nz);
    nx = sf_n(ax);
    nt = sf_n(at);
    
/*------------------------------------------------------------*/

    /* allocate arrays */
    us=sf_floatalloc3(nt,nx,nb);  /*   source wavefield */
    ur=sf_floatalloc3(nt,nx,nb);  /* receiver wavefield */
    ii=sf_floatalloc (   nx);     /* image */

    ts=sf_floatalloc2(nt,nx);
    tr=sf_floatalloc2(nt,nx);
    
    gg=sf_floatalloc (  2*nl+1);
    ww=sf_floatalloc3(4,2*nl+1,na);
    kt=sf_intalloc3  (4,2*nl+1,na);
    kx=sf_intalloc3  (4,2*nl+1,na);
    
/*------------------------------------------------------------*/
    /* taper */
    for(il=0;il<2*nl+1;il++) {
	l = ol + (il-nl)*dl;
	l /= (nl/2);
	l /= sig;
	gg[il] = exp(-l*l);
    }

/*------------------------------------------------------------*/
    /* compute bilinear interpolation indices and weights */
    for(ia=0;ia<na;ia++) {
	a  = oa + ia * da;
	a *= SF_PI/180.;
	
	for(il=0;il<2*nl+1;il++){
	    l = ol + (il-nl)*dl;   
	    
	    lt = l*sin(a);
	    lx = l*cos(a);	    
	    
	    kt[ia][il][0] = (floor)(lt);
	    kt[ia][il][1] = kt[ia][il][0] +1;
	    kt[ia][il][2] = kt[ia][il][0];
	    kt[ia][il][3] = kt[ia][il][0] +1;
	    
	    kx[ia][il][0] = (floor)(lx);
	    kx[ia][il][1] = kx[ia][il][0];
	    kx[ia][il][2] = kx[ia][il][0] +1;
	    kx[ia][il][3] = kx[ia][il][0] +1;
	    
	    ft = lt-kt[ia][il][0];
	    fx = lx-kx[ia][il][0];
	    
	    ww[ia][il][0] = (1-ft)*(1-fx);
	    ww[ia][il][1] =    ft *(1-fx);
	    ww[ia][il][2] = (1-ft)*   fx ;
	    ww[ia][il][3] =    ft *   fx ;

	    for(ic=0;ic<4;ic++) {
		ww[ia][il][ic] *= gg[il];
	    }
	}
    }
    
/*------------------------------------------------------------*/

    for (; nz > 0; nz -= nb) {
	if (nb > nz) nb=nz;

	sf_floatread(us[0][0],nt*nx*nb,Fs);  /* read   source wavefield */
	sf_floatread(ur[0][0],nt*nx*nb,Fr);  /* read receiver wavefield */

	if(verb) fprintf(stderr,"  z   a\n");
	if(verb) fprintf(stderr,"%3d %3d\n",nb-1,na-1);

	for(        ib=0; ib<nb; ib++) {	    
	    for(    ix=0; ix<nx; ix++) {
		ii[ix] = 0;
	    }

	    for(        ia=0;ia<na;    ia++) {	
		if(verb) fprintf(stderr,"%3d %3d",ib,ia);
		
		for(    ix=0; ix<nx; ix++) {
		    for(it=0; it<nt; it++) {	
			ts[ix][it] = 0;
			tr[ix][it] = 0;
		    }
		}

#ifdef _OPENMP
#pragma omp parallel for schedule(dynamic,ompchunk) private(il,ic,ix,it,mx,mt,ht,hx,gt,gx,jt,jx,wo) shared(ia,kt,kx,nl,nx,nt,ts,tr,ww,us,ur)
#endif	
		for(    il=0;il<2*nl+1;il++) {
		    for(ic=0;ic<4;     ic++) {
			mt=kt[ia][il][ic]; ht = SF_ABS(mt); gt = nt - ht;
			mx=kx[ia][il][ic]; hx = SF_ABS(mx); gx = nx - hx;
			wo=ww[ia][il][ic];

			for(    ix=hx; ix<gx; ix++) { jx = ix+mx;
			    for(it=ht; it<gt; it++) { jt = it+mt;
				ts[ ix ][ it ] += wo * us[ib][ jx ][ jt ];
				tr[ ix ][ it ] += wo * ur[ib][ jx ][ jt ];
			    } /* x loop */
			}     /* t loop */

		    }         /* c loop */
		}             /* l loop */

		for(    ix=0; ix<nx; ix++) {
		    for(it=0; it<nt; it++) {
			ii[ix] += ts[ix][it] * tr[ix][it];
		    }
		}
		if(verb) fprintf(stderr,"\b\b\b\b\b\b\b\b\b\b\b\b");
	    }                 /* a loop */

	    sf_floatwrite(ii,nx,Fi);

	}                     /* z loop (in block) */
    }                         /* z loop (blocks) */

    exit (0);
}
