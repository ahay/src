/* Imaging condition */
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
    int version;
    float eps;

    sf_file Fi,Fs,Fr;    /* I/O files */

    sf_axis at,az,ax,aa; /* cube axes */
    int     nt,nz,nx,nb;
    int        iz,ix,ib;

    float  **ii=NULL;

    float  **t1=NULL;
    float  **t2=NULL;
    float  **t3=NULL;

    float ***us=NULL;
    float ***ur=NULL;

    float ts,tr;

    int ompchunk; 

/*------------------------------------------------------------*/
    /* init RSF */
    sf_init(argc,argv);

    if(! sf_getint("ompchunk",&ompchunk)) ompchunk=1;  /* OpenMP data chunk size */
    if(! sf_getbool("verb",   &verb))     verb=false;  /* verbosity flag */
    if(! sf_getint("nbuf",    &nb))       nb=1;        /* buffer size */
    if(! sf_getint("version", &version))  version=0;   /* I.C. version (see paper) */
    if(! sf_getfloat("eps",   &eps))      eps=1e-6;    /* epsilon */

    Fs = sf_input ("in" ); /*   source wavefield */
    Fr = sf_input ("ur" ); /* receiver wavefield */
    Fi = sf_output("out"); /* image */

    /* read axes */
    az=sf_iaxa(Fs,1); sf_setlabel(az,"z"); if(verb) sf_raxa(az); /* depth */
    ax=sf_iaxa(Fs,2); sf_setlabel(ax,"x"); if(verb) sf_raxa(ax); /* position */
    at=sf_iaxa(Fs,3); sf_setlabel(at,"t"); if(verb) sf_raxa(at); /* time */

    /* set output axes */
    aa=sf_maxa(1,0,1); 
    sf_oaxa(Fi,aa,3);

    nz = sf_n(az);
    nx = sf_n(ax);
    nt = sf_n(at);

    nb = SF_MIN(nb,nt);

/*------------------------------------------------------------*/
    /* allocate arrays */
    us=sf_floatalloc3(nz,nx,nb);  /*  source wavefield  */
    ur=sf_floatalloc3(nz,nx,nb);  /* receiver wavefield */
    ii=sf_floatalloc2(nz,nx);     /* image */
    if(version>1) {
	t3=sf_floatalloc2(nz,nx);
	t2=sf_floatalloc2(nz,nx);
	t1=sf_floatalloc2(nz,nx);
    }

    /* init output */
    for (ix=0; ix<nx; ix++) {
	for (iz=0; iz<nz; iz++) {
	    ii[ix][iz]=0;
	}
    }
    if(version>1) {
	for (ix=0; ix<nx; ix++) {
	    for (iz=0; iz<nz; iz++) {
		t3[ix][iz]=0;
		t2[ix][iz]=0;
		t1[ix][iz]=0;
	    }
	}
    }

    switch (version){
	case 3: /* SUM(Us Ur) / sqrt( SUM(Us Us) SUM(Ur Ur) ) */
	    for (; nt > 0; nt -= nb) {
		if (nb > nt) nb=nt;
		if(verb) sf_warning("nsiz=%ld nbuf=%ld",nt,nb);
		
		sf_floatread(us[0][0],nz*nx*nb,Fs);
		sf_floatread(ur[0][0],nz*nx*nb,Fr);	

		for(ib=0; ib<nb; ib++) {
		    if(verb) fprintf(stderr,"\b\b\b\b\b%d",ib);
#ifdef _OPENMP
#pragma omp parallel for schedule(dynamic,ompchunk) private(iz,ix) shared(ib,nz,nx,ii,us,ur)
#endif
		    for(    ix=0; ix<nx; ix++) {
			for(iz=0; iz<nz; iz++) {
			    ts = us[ib][ix][iz];
			    tr = ur[ib][ix][iz];			    
			    t1[ix][iz] += (ts * tr);
			    t2[ix][iz] += (ts * ts);
			    t3[ix][iz] += (tr * tr);
			} /* nz */
		    } /* nx */ 
		} /* nb */
		if(verb) fprintf(stderr,"\b\b\b\b\b");
	    } /* nt */
	    for (ix=0; ix<nx; ix++) {
		for (iz=0; iz<nz; iz++) {
		    ii[ix][iz] = t1[ix][iz] / ( sqrt( t2[ix][iz] * t3[ix][iz] ) + (nz*nx*nb)*eps);
		}
	    }
	    sf_floatwrite(ii[0],nz*nx,Fi);	
	    
	    break;
	    
	    /*------------------------------------------------------------*/
	case 2: /* SUM( Us Ur) / SUM (Us Us) */
	    for (; nt > 0; nt -= nb) {
		if (nb > nt) nb=nt;
		if(verb) sf_warning("nsiz=%ld nbuf=%ld",nt,nb);
		
		sf_floatread(us[0][0],nz*nx*nb,Fs);
		sf_floatread(ur[0][0],nz*nx*nb,Fr);	

		for(ib=0; ib<nb; ib++) {
		    if(verb) fprintf(stderr,"\b\b\b\b\b%d",ib);
#ifdef _OPENMP
#pragma omp parallel for schedule(dynamic,ompchunk) private(iz,ix) shared(ib,nz,nx,ii,us,ur)
#endif
		    for(    ix=0; ix<nx; ix++) {
			for(iz=0; iz<nz; iz++) {
			    ts = us[ib][ix][iz];
			    tr = ur[ib][ix][iz];			    
			    t1[ix][iz] += (ts * tr);
			    t2[ix][iz] += (ts * ts);
			} /* nz */
		    } /* nx */ 
		} /* nb */
		if(verb) fprintf(stderr,"\b\b\b\b\b");
	    } /* nt */
	    for (ix=0; ix<nx; ix++) {
		for (iz=0; iz<nz; iz++) {
		    ii[ix][iz] = t1[ix][iz] / ( t2[ix][iz] + (nz*nx*nb)*eps*(nz*nx*nb)*eps);
		}
	    }
	    sf_floatwrite(ii[0],nz*nx,Fi);	    
	    break;
	    
	    /*------------------------------------------------------------*/
	case 1: /* SUM( (Us Ur)/(Us Us) ) */
	    for (; nt > 0; nt -= nb) {
		if (nb > nt) nb=nt;
		if(verb) sf_warning("nsiz=%ld nbuf=%ld",nt,nb);
		
		sf_floatread(us[0][0],nz*nx*nb,Fs);
		sf_floatread(ur[0][0],nz*nx*nb,Fr);	

		for(ib=0; ib<nb; ib++) {
		    if(verb) fprintf(stderr,"\b\b\b\b\b%d",ib);
#ifdef _OPENMP
#pragma omp parallel for schedule(dynamic,ompchunk) private(iz,ix) shared(ib,nz,nx,ii,us,ur)
#endif
		    for(    ix=0; ix<nx; ix++) {
			for(iz=0; iz<nz; iz++) {
			    ts = us[ib][ix][iz];
			    tr = ur[ib][ix][iz];			    
			    ii[ix][iz] += (ts * tr) / (ts * ts + eps*eps);
			} /* nz */
		    } /* nx */ 
		} /* nb */
		if(verb) fprintf(stderr,"\b\b\b\b\b");
	    } /* nt */
	    sf_floatwrite(ii[0],nz*nx,Fi);
	    break;
	    
	    /*------------------------------------------------------------*/
	case 0: /* Us * Ur */
	    for (; nt > 0; nt -= nb) {
		if (nb > nt) nb=nt;
		if(verb) sf_warning("nsiz=%ld nbuf=%ld",nt,nb);
		
		sf_floatread(us[0][0],nz*nx*nb,Fs);
		sf_floatread(ur[0][0],nz*nx*nb,Fr);

		for(ib=0; ib<nb; ib++) {
		    if(verb) fprintf(stderr,"\b\b\b\b\b%d",ib);
#ifdef _OPENMP
#pragma omp parallel for schedule(dynamic,ompchunk) private(iz,ix) shared(ib,nz,nx,ii,us,ur)
#endif
		    for(    ix=0; ix<nx; ix++) {
			for(iz=0; iz<nz; iz++) {    
			    ii[ix][iz] += us[ib][ix][iz]
				*         ur[ib][ix][iz];
			} /* nz */
		    } /* nx */ 
		} /* nb */
		if(verb) fprintf(stderr,"\b\b\b\b\b");
	    } /* nt */
	    sf_floatwrite(ii[0],nz*nx,Fi);	    
	    break;
    } /* switch */

    exit (0);
}
