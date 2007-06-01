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

    float st,stp,stm;
    float sx,sxp,sxm;
    float sz,szp,szm;

    int ompchunk;

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

    st = 0;
    for        (iht=-nht; iht<=nht; iht++) {
	for    (ihx=-nhx; ihx<=nhx; ihx++) {
	    for(ihz=-nhz; ihz<=nhz; ihz++) {
		st += 
		    uu[nht-iht][nhx-ihx][nhz-ihz] *
		    uu[nht+iht][nhx+ihx][nhz+ihz];
	    }
	}
    }

    if(verb) fprintf(stderr,"\n");
    for(it=nht; it<nt-nht-1; it++) {
	sx = st;

	for(ix=nhx; ix<nx-nhx-1; ix++) {
	    sz = sx;

	    for(iz=nhz; iz<nz-nhz-1; iz++) {
		if(verb) fprintf(stderr,"%3d %3d %3d",it,ix,iz);

		szp=0;
		szm=0;
		for    (iht=-nht; iht<=nht; iht++) {
		    for(ihx=-nhx; ihx<=nhx; ihx++) {
			szp += 
			    uu[it-iht][ix-ihx][iz+nhz+1] *
			    uu[it+iht][ix+ihx][iz+nhz+1];
			szm += 
			    uu[it-iht][ix-ihx][iz-nhz  ] *
			    uu[it+iht][ix+ihx][iz-nhz  ];
		    }
		}
		sz += szp;
		sz -= szm;

		ww[it][ix][iz] = sz;
		if(verb) fprintf(stderr,"\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b");
	    }

	    sxp=0;
	    sxm=0;
	    for    (iht=-nht; iht<=nht; iht++) {
		for(ihz=-nhz; ihz<=nhz; ihz++) {
		    sxp += 
			uu[it-iht][ix+nhx+1][iz-ihz] *
			uu[it+iht][ix+nhx+1][iz+ihz];
		    sxm += 
			uu[it-iht][ix-nhx  ][iz-ihz] *
			uu[it+iht][ix-nhx  ][iz+ihz];
		}
	    }
	    sx += sxp;
	    sx -= sxm;
	}

	stp=0;
	stm=0;
	for    (ihz=-nhz; ihz<=nhz; ihz++) {
	    for(ihx=-nhx; ihx<=nhx; ihx++) {
		stp += 
		    uu[it-nht  ][ix-ihx][iz-ihz] *
		    uu[it-nht  ][ix+ihx][iz+ihz];
		stm += 
		    uu[it+nht+1][ix-ihx][iz-ihz] *
		    uu[it+nht+1][ix+ihx][iz+ihz];
	    }
	}
	st += stp;
	st -= stm;
    }
    if(verb) fprintf(stderr,"\n");
    
    /* write output */
    sf_floatwrite(ww[0][0],nz*nx*nt,Fw);  

    exit (0);
}
