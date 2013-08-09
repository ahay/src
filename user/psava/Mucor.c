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

    if(verb) fprintf(stderr,"   t    x\n");
    if(verb) fprintf(stderr,"%4d %4d\n",nt,nx);
    for(        it=nht; it<nt-nht; it++) { lot=-nht; hit=nht+1;
	for(    ix=nhx; ix<nx-nhx; ix++) { lox=-nhx; hix=nhx+1;
	    if(verb) fprintf(stderr,"%4d %4d",it,ix);
#ifdef _OPENMP
#pragma omp parallel for schedule(dynamic,ompchunk) private(iz,iht,ihx,ihz,loz,hiz,jt,jx) shared(ur,nz,nhz,lot,hit,lox,hix)
#endif		
	    for(iz=nhz; iz<nz-nhz; iz++) { loz=-nhz; hiz=nhz+1;
		for(        iht=lot; iht<hit; iht++) { jt=it-iht; 
		    for(    ihx=lox; ihx<hix; ihx++) { jx=ix-ihx; 
			for(ihz=loz; ihz<hiz; ihz++) { jz=iz-ihz;
			    ii[ix][iz] += ur[jt][jx][jz]
				*         us[jt][jx][jz];
			} /* nhz */
		    } /* nhx */
		} /* nht */
	    } /* nz */
	    if(verb) fprintf(stderr,"\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b");
	} /* nx */
    } /* nt */
    if(verb) fprintf(stderr,"\n");

    sf_floatwrite(ii[0],nz*nx,Fi);    

    exit (0);
}
