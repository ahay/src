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

    sf_file Fs=NULL; /* input   */
    sf_file Fr=NULL; /* input   */
    sf_file Fi=NULL; /* output  */

    float ***us=NULL;
    float ***ur=NULL;
    float  **ii=NULL;

    float   *gg=NULL; /*  taper */
    float ***ww=NULL; /* weight */
    int   ***kz=NULL; /* slant-stack index on axis 1*/
    int   ***kx=NULL; /* slant-stack index on axis 2*/

    /* cube axes */
    sf_axis at,az,ax,aa,ll,ak; /* cube axes */
    int     nt,nz,nx,na,nl,nb;
    int        iz,ix,ia,il,ib;
    int        jz,jx;
    float oa,da,a;
    float ol, dl,l,lz,lx;
    float oz,dz,fz;
    float ox,dx,fx;

    int ic;
    int ompchunk; 

    float ts,tr;

/*------------------------------------------------------------*/

    /* init RSF */
    sf_init(argc,argv);

    if(! sf_getint("ompchunk",&ompchunk)) ompchunk=1;  /* OpenMP data chunk size */
    if(! sf_getbool("verb",   &verb))     verb=false;  /* verbosity flag */
    if(! sf_getint("nbuf",    &nb))       nb=1;        /* buffer size */

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
    az = sf_iaxa(Fs,1); sf_setlabel(az,"z"); sf_setunit(az,"");
    ax = sf_iaxa(Fs,2); sf_setlabel(ax,"x"); sf_setunit(ax,""); 
    at = sf_iaxa(Fs,3); sf_setlabel(at,"t"); sf_setunit(at,""); 
    nz = sf_n(az); oz=sf_o(az); dz=sf_d(az);
    nx = sf_n(ax); ox=sf_o(ax); dx=sf_d(ax);

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
    sf_oaxa(Fi,az,1);
    sf_oaxa(Fi,ax,2);
    sf_oaxa(Fi,ak,3);

    nz = sf_n(az);
    nx = sf_n(ax);
    nt = sf_n(at);
    
    nb = SF_MIN(nb,nt);

/*------------------------------------------------------------*/

    /* allocate arrays */
    us=sf_floatalloc3(nz,nx,nb);  //   source wavefield
    ur=sf_floatalloc3(nz,nx,nb);  // receiver wavefield
    ii=sf_floatalloc2(nz,nx);     // image

    /* init output */
    for (ix=0; ix<nx; ix++) {
	for (iz=0; iz<nz; iz++) {
	    ii[ix][iz]=0;
	}
    }

    gg=sf_floatalloc (  2*nl);
    ww=sf_floatalloc3(4,2*nl,na);
    kz=sf_intalloc3  (4,2*nl,na);
    kx=sf_intalloc3  (4,2*nl,na);
    
/*------------------------------------------------------------*/
    /* taper */
    for(il=0;il<2*nl;il++) {
	l = ol + (il-nl)*dl;
	l /= (nl/2);
	gg[il] = exp(-l*l);
    }

/*------------------------------------------------------------*/
    /* compute indices and weight */
    for(ia=0;ia<na;ia++) {
	a  = oa + ia * da;
	a *= SF_PI/180.;
	
	for(il=0;il<2*nl;il++){
	    l = ol + (il-nl)*dl;   
	    
	    lz = l*sin(a);
	    lx = l*cos(a);	    
	    
	    kz[ia][il][0] = (floor)(lz);
	    kz[ia][il][1] = kz[ia][il][0] +1;
	    kz[ia][il][2] = kz[ia][il][0];
	    kz[ia][il][3] = kz[ia][il][0] +1;
	    
	    kx[ia][il][0] = (floor)(lx);
	    kx[ia][il][1] = kx[ia][il][0];
	    kx[ia][il][2] = kx[ia][il][0] +1;
	    kx[ia][il][3] = kx[ia][il][0] +1;
	    
	    fz = lz-kz[ia][il][0];
	    fx = lx-kx[ia][il][0];
	    
	    ww[ia][il][0] = (1-fz)*(1-fx);
	    ww[ia][il][1] =    fz *(1-fx);
	    ww[ia][il][2] = (1-fz)*   fx ;
	    ww[ia][il][3] =    fz *   fx ;

	    for(ic=0;ic<4;ic++) {
		ww[ia][il][ic] *= gg[il];
	    }
	}
    }
    
/*------------------------------------------------------------*/

    for (; nt > 0; nt -= nb) {
	if (nb > nt) nb=nt;
	if(verb) sf_warning("nsiz=%ld nbuf=%ld",nt,nb);

	sf_floatread(us[0][0],nz*nx*nb,Fs);
	sf_floatread(ur[0][0],nz*nx*nb,Fr);

	for(            ib=0;  ib<nb;    ib++) {
	    if(verb) fprintf(stderr,"\b\b\b\b\b%d",ib);
#ifdef _OPENMP
#pragma omp parallel for schedule(dynamic,ompchunk) private(ia,ix,iz,il,ic,ts,tr) shared(na,nz,nx,nl,ii,us,ur,ww,kx,kz,)
#endif
	    for(        ix=nl; ix<nx-nl; ix++) {
		for(    iz=nl; iz<nz-nl; iz++) {

		    for(ia=0;  ia<na;    ia++) {			
			ts=0;
			tr=0;
			for(il=0;il<2*nl;il++) {
			    for(ic=0;ic<4;ic++) {
				jz = iz+kz[ia][il][ic];
				jx = ix+kx[ia][il][ic];

				ts += ww[ia][il][ic] * us[ib][ jx ][ jz ];
				tr += ww[ia][il][ic] * ur[ib][ jx ][ jz ];
			    }
			}
			ii[ix][iz] += ts*tr;
		    } // ia

		} // iz
	    } // ix
	} // ib
	if(verb) fprintf(stderr,"\b\b\b\b\b");

    }
    sf_floatwrite(ii[0],nz*nx,Fi);
    
    exit (0);
}
