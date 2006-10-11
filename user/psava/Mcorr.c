/* Cross-correlation of time series */
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

    /* I/O files */
    sf_file Fi,Fs,Fr;

    /* cube axes */
    sf_axis at,az,ax,aa;
    int     nt,nz,nx, nhz,nhx;
    int     it,iz,ix, ihz,ihx;

    /* arrays */
    float **ii=NULL, **us=NULL,**ur=NULL;

    int ompchunk; 

/*------------------------------------------------------------*/

    /* init RSF */
    sf_init(argc,argv);

    if(! sf_getint("ompchunk",&ompchunk)) ompchunk=1;  /* OpenMP data chunk size */
    if(! sf_getbool("verb",&verb)) verb=false;         /* verbosity flag */

    if(! sf_getint("nhz",&nhz)) nhz=0;
    if(! sf_getint("nhx",&nhx)) nhx=0;

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
    us=sf_floatalloc2(nz,nx);
    ur=sf_floatalloc2(nz,nx);
    ii=sf_floatalloc2(nz,nx);

    /* init output */
    for (iz=0; iz<nz; iz++) {
	for (ix=0; ix<nx; ix++) {
	    ii[ix][iz]=0;
	}
    }

    if(verb) fprintf(stderr,"\n");
    for(it=0; it<nt; it++) {
	if(verb) fprintf(stderr,"\b\b\b\b\b%d",it);

	sf_floatread(us[0],nz*nx,Fs);
	sf_floatread(ur[0],nz*nx,Fr);

#ifdef _OPENMP
#pragma omp parallel for schedule(dynamic,ompchunk) private(iz,ix,ihz,ihx) shared(nz,nx,ii,us,ur)
#endif
	for(    ihz=-nhz; ihz<nhz+1; ihz++) {
	    for(ihx=-nhx; ihx<nhx+1; ihx++) {
		for(    iz=0+SF_ABS(ihz); iz<nz-SF_ABS(ihz); iz++) {
		    for(ix=0+SF_ABS(ihx); ix<nx-SF_ABS(ihx); ix++) {
			ii[ix][iz] += us[ix-ihx][iz-ihz] * ur[ix+ihx][iz+ihz];
		    } // nx
		} // nz
	    } // nhx
	} // nhz

    }
    if(verb) fprintf(stderr,"\n");

    sf_floatwrite(ii[0],nz*nx,Fi);

    exit (0);
}
