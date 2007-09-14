/* Compute lagged-products */

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

/* 
 * inputs are wavefields organized as z-x-t
 * output are lagged products organized as z-x-hz-hx-ht
 */

#ifdef _OPENMP
#include <omp.h>
#endif

int main(int argc, char* argv[])
{
    bool verb;

    sf_file Fs,Fr,Fi;    /* I/O files */
    sf_axis az,ax,at,aa; /* cube axes */

    int     nz,nx,nt, nhz,nhx,nht;
    int     iz,ix,it, ihz,ihx,iht;
    int     izs,ixs,    jhz,jhx;
    int     izr,ixr;

    float ****ii=NULL,**us=NULL,**ur=NULL; /* arrays */

    int ompchunk; 

    int loz,hiz;
    int lox,hix;
    int lot,hit;

    sf_axis acz,acx;
    int     jcz,jcx;

    /*------------------------------------------------------------*/
    /* init RSF */
    sf_init(argc,argv);

    if(! sf_getint("ompchunk",&ompchunk)) ompchunk=1;  /* OpenMP data chunk size */
    if(! sf_getbool("verb",&verb)) verb=false;         /* verbosity flag */
    
    Fs = sf_input ("in" ); /*   source wavefield */
    Fr = sf_input ("ur" ); /* receiver wavefield */
    Fi = sf_output("out"); /* image */

    /* read axes */
    az=sf_iaxa(Fs,1); sf_setlabel(az,"z"); if(verb) sf_raxa(az);
    ax=sf_iaxa(Fs,2); sf_setlabel(ax,"x"); if(verb) sf_raxa(ax);
    at=sf_iaxa(Fs,3); sf_setlabel(at,"t"); if(verb) sf_raxa(at);

    if(!sf_getint ("jcx",&jcx) || nx==1) jcx=1;
    if(!sf_getint ("jcz",&jcz) || nz==1) jcz=1;

    acx = sf_maxa(SF_MAX(1,sf_n(ax)/jcx),sf_o(ax),sf_d(ax)*jcx); sf_setlabel(acx,"cx");
    acz = sf_maxa(SF_MAX(1,sf_n(az)/jcz),sf_o(az),sf_d(az)*jcz); sf_setlabel(acz,"cz");

    sf_oaxa(Fi,az,1);
    sf_oaxa(Fi,ax,2);

    nz = sf_n(az);
    nx = sf_n(ax);
    nt = sf_n(at);

    if(! sf_getint("nhz",&nhz)) nhz=0; /* number of lags on the z axis */
    if(! sf_getint("nhx",&nhx)) nhx=0; /* number of lags on the x axis */
    if(! sf_getint("nht",&nht)) nht=0; /* number of lags on the t axis */
    sf_warning("nhz=%d nhx=%d nht=%d",2*nhz+1,2*nhx+1,2*nht+1);

    /* set output axes */
    aa=sf_maxa(2*nhz+1,-nhz*sf_d(az),sf_d(az));
    sf_setlabel(aa,"hz");
    if(verb) sf_raxa(aa);
    sf_oaxa(Fi,aa,3);
    
    aa=sf_maxa(2*nhx+1,-nhx*sf_d(ax),sf_d(ax)); 
    sf_setlabel(aa,"hx");
    if(verb) sf_raxa(aa);
    sf_oaxa(Fi,aa,4);

    aa=sf_maxa(2*nht+1,-nht*sf_d(at),sf_d(at)); 
    sf_setlabel(aa,"ht");
    if(verb) sf_raxa(aa);
    sf_oaxa(Fi,aa,5);

    /* allocate work arrays */
    us=sf_floatalloc2(nz,nx);
    ur=sf_floatalloc2(nz,nx);
    ii=sf_floatalloc4(nz,nx,2*nhz+1,2*nhx+1);

    if(verb) fprintf(stderr," nt   ht\n");
    if(verb) fprintf(stderr,"%4d %3d \n",nt-1,2*nht);

    /* t-lag loop */
    for(iht=-nht; iht<nht+1; iht++) { lot=SF_ABS(iht); hit=nt-SF_ABS(iht);

	/* seek in input */
	sf_seek(Fs,(lot+iht)*nz*nx*sizeof(float),SEEK_SET);
	sf_seek(Fr,(lot-iht)*nz*nx*sizeof(float),SEEK_SET);
	
	/* zero output */
	for(        ihx=-nhx; ihx<nhx+1; ihx++) { jhx=nhx+ihx;
	    for(    ihz=-nhz; ihz<nhz+1; ihz++) { jhz=nhz+ihz;
		for(    ix=0; ix<nx; ix++) { 
		    for(iz=0; iz<nz; iz++) { 
			ii[jhx][jhz][ix][iz] = 0;
		    } /* nz */
		} /* nx */
	    } /* nhz */
	} /* nhx */

	/* t loop */
	for(it=lot; it<hit; it++) {
	    if(verb) fprintf(stderr,"%4d %3d",it,nht+iht);
	    
	    /* read input */
	    sf_floatread(us[0],nz*nx,Fs);
	    sf_floatread(ur[0],nz*nx,Fr);
	    
#ifdef _OPENMP
#pragma omp parallel for schedule(dynamic,ompchunk) \
    private(ihx,lox,hix,jhx, ihz,loz,hiz,jhz, ix,iz,ixs,izs,ixr,izr) \
    shared(nhx,nhz,ii,us,ur)
#endif		
	    for(          ihx=-nhx; ihx<nhx+1; ihx++) { lox=SF_ABS(ihx); hix=nx-lox; jhx=nhx+ihx;
		for(      ihz=-nhz; ihz<nhz+1; ihz++) { loz=SF_ABS(ihz); hiz=nz-loz; jhz=nhz+ihz;  

		    for(    ix=lox;  ix<hix;    ix++) { ixs=ix-ihx; ixr=ix+ihx;
			for(iz=loz;  iz<hiz;    iz++) { izs=iz-ihz; izr=iz+ihz;

			    ii[jhx][jhz][ix][iz] += us[ixs][izs] 
				*                   ur[ixr][izr];

			} /* nz */
		    } /* nx */  

		} /* nhz */
	    } /* nhx */
	    
	    if(verb) fprintf(stderr,"\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b");
	} /* nt */
	
	/* write output */
	sf_floatwrite(ii[0][0][0],nz*nx*(2*nhz+1)*(2*nhx+1),Fi);    

    } /* nht */
    if(verb) fprintf(stderr,"\n");
    
    exit (0);
}
