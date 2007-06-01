/* Assymptotic Wigner distribution in space-time */
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

    sf_file Fu,Fw; /* I/O files */
    sf_axis at,az,ax; /* cube axes */

    float ***uu=NULL, ***ww=NULL;
    
    int nhz,nhx,nht;
    int ihz,ihx,iht;
    int  nz, nx, nt;
    int  iz, ix, it;
    int  jz, jx, jt;
    int  kz, kx, kt;

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
    if(! sf_getint("nht",&nht)) nht=0;
    sf_warning("nht=%d nhx=%d nhz=%d",2*nht+1,2*nhx+1,2*nhz+1);
    
    Fu = sf_input ("in" ); /*   input field */
    Fw = sf_output("out"); /* wigner distribution */

    /* read axes */
    az=sf_iaxa(Fu,1); sf_setlabel(az,"z"); if(verb) sf_raxa(az); /* depth */
    ax=sf_iaxa(Fu,2); sf_setlabel(ax,"x"); if(verb) sf_raxa(ax); /* position */
    at=sf_iaxa(Fu,3); sf_setlabel(at,"t"); if(verb) sf_raxa(at); /* time */

    nz = sf_n(az);
    nx = sf_n(ax);
    nt = sf_n(at);

    /* allocate work arrays */
    uu=sf_floatalloc3(nz,nx,nt);
    ww=sf_floatalloc3(nz,nx,nt);

    /* read input */
    sf_floatread(uu[0][0],nz*nx*nt,Fu);

    if(verb) fprintf(stderr," ht  hx  hz\n");
    if(verb) fprintf(stderr,"%3d %3d %3d\n",2*nht,2*nhx,2*nhz);
    for(        iht=-nht; iht<nht+1; iht++) { lot=SF_ABS(iht); hit=nt-SF_ABS(iht);
	for(    ihx=-nhx; ihx<nhx+1; ihx++) { lox=SF_ABS(ihx); hix=nx-SF_ABS(ihx);
	    for(ihz=-nhz; ihz<nhz+1; ihz++) { loz=SF_ABS(ihz); hiz=nz-SF_ABS(ihz);

		if(verb) fprintf(stderr,"%3d %3d %3d",nht+iht,nhx+ihx,nhz+ihz);

#ifdef _OPENMP
#pragma omp parallel for schedule(dynamic,ompchunk) private(it,iz,ix,jt,jx,jz,kt,kx,kz) shared(iht,ihx,ihz,lot,hit,lox,hix,loz,hiz,uu,ww)
#endif	
		for(        it=lot; it<hit; it++) { jt=it-iht; kt=it+iht;
		    for(    ix=lox; ix<hix; ix++) { jx=ix-ihx; kx=ix+ihx;
			for(iz=loz; iz<hiz; iz++) { jz=iz-ihz; kz=iz+ihz;
			    ww[it][ix][iz] += uu[jt][jx][jz] 
				*             uu[kt][kx][kz];
			} /* nz */
		    }     /* nx */
		}         /* nt */

		if(verb) fprintf(stderr,"\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b");

	    } /* nhz */
	}     /* nhx */
    }         /* nht */
    
    /* write output */
    sf_floatwrite(ww[0][0],nz*nx*nt,Fw);  

    exit (0);
}
