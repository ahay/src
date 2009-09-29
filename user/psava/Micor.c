/* Interferometric cross-correlation of time series (zero-lag output) */
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

    sf_file Fi,Fs,Fr;    /* I/O files */
    sf_axis at,az,ax,aa; /* cube axes */

    int     nt,nz,nx, nhz,nhx,nht;
    int     it,iz,ix, ihz,ihx,iht;
    int     jt,jz,jx;
    int     kt,kz,kx;

    float **ii=NULL, ***us=NULL,***ur=NULL; /* arrays */

    int ompchunk; 

    int lox,hix;
    int loz,hiz;
    int lot,hit;

/*------------------------------------------------------------*/

    /* init RSF */
    sf_init(argc,argv);

    if(! sf_getint("ompchunk",&ompchunk)) ompchunk=1;  /* OpenMP data chunk size */
    if(! sf_getbool("verb",&verb)) verb=false;         /* verbosity flag */

    if(! sf_getint("nhz",&nhz)) nhz=0;
    if(! sf_getint("nhx",&nhx)) nhx=0;
    if(! sf_getint("nht",&nht)) nht=1;
    sf_warning("nht=%d nhx=%d nhz=%d",2*nht+1,2*nhx+1,2*nhz+1);

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

    /* allocate work arrays */
    us=sf_floatalloc3(nz,nx,nt);
    ur=sf_floatalloc3(nz,nx,nt);
    ii=sf_floatalloc2(nz,nx);

    /* init output */
    for (iz=0; iz<nz; iz++) {
	for (ix=0; ix<nx; ix++) {
	    ii[ix][iz]=0;
	}
    }

    sf_floatread(us[0][0],nz*nx*nt,Fs);
    sf_floatread(ur[0][0],nz*nx*nt,Fr);

    if(verb) fprintf(stderr," ht  hx\n");
    if(verb) fprintf(stderr,"%3d %3d\n",2*nht,2*nhx);
    for(        iht=-nht; iht<nht+1; iht++) { lot=SF_ABS(iht); hit=nt-SF_ABS(iht);
	for(    ihx=-nhx; ihx<nhx+1; ihx++) { lox=SF_ABS(ihx); hix=nx-SF_ABS(ihx);
	    if(verb) fprintf(stderr,"%3d %3d",nht+iht,nhx+ihx);
#ifdef _OPENMP
#pragma omp parallel for schedule(dynamic,ompchunk) private(ihz,it,iz,ix,loz,hiz,jt,jx,jz,kt,kx,kz) shared(iht,ihx,lot,hit,lox,hix,nhz,ii,us,ur)
#endif		
	    for(ihz=-nhz; ihz<nhz+1; ihz++) { loz=SF_ABS(ihz); hiz=nz-SF_ABS(ihz);
		for(        it=lot; it<hit; it++) { jt=it-iht; kt=it+iht;
		    for(    ix=lox; ix<hix; ix++) { jx=ix-ihx; kx=ix+ihx;
			for(iz=loz; iz<hiz; iz++) { jz=iz-ihz; kz=iz+ihz;
			    ii[ix][iz] += us[jt][jx][jz] 
				*         ur[kt][kx][kz];
			} /* nz */
		    } /* nx */
		} /* nt */
	    } /* nhz */
	    if(verb) fprintf(stderr,"\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b");
	} /* nhx */
    } /* nht */
    if(verb) fprintf(stderr,"\n");

    sf_floatwrite(ii[0],nz*nx,Fi);    

    exit (0);
}
