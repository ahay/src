/* Wilson-Burg factorization  */
/*
 Copyright (C) 2008 Colorado School of Mines
 
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
#include "wilson.h"
#include "wilson2.h"

#ifdef _OPENMP
#include <omp.h>
#endif

int main(int argc, char* argv[])
{
    bool verb,stable;
    int  ompchunk; 
	
    sf_file Fa=NULL; /* autocorrelation     */
    sf_file Ff=NULL; /* filter coefficients */
    sf_file Fl=NULL; /* filter lags         */
	
    /* cube axes */
    sf_axis a1,a2,af,aj;
    int     i1,i2;
    int     k1,k2;
    int     j1,j2;
	
    /* H.F. */
    sf_filter aa,ff;
	
    float **ain; /* autocorrelation array */
    int     nn,kk;
    int     na,ja;
    int     nf,jf;
    int    maxlag;
    int     niter;
    char    *file;
	
    /*------------------------------------------------------------*/
    /* init RSF */
    sf_init(argc,argv);
	
    if(! sf_getint ("ompchunk",&ompchunk)) ompchunk=1;     /* OMP chunk size */
    if(! sf_getbool("verb",    &verb))         verb=false; /* verbosity flag */
    if (!sf_getbool("stable",  &stable))     stable=false; /* stability flag */
	
    Fa = sf_input ("in" );
    Ff = sf_output("out");
    Fl = sf_output("lag");
    
    /* input axes */
    a1 = sf_iaxa(Fa,1);
    a2 = sf_iaxa(Fa,2);
    if(verb) { 
		sf_raxa(a1); 
		sf_raxa(a2); 
    }
	
    if(! sf_getint("niter",&niter)) niter=20; /* Wilson iterations */
    if(! sf_getint("nn",   &nn))    nn=1000; /* Helix diameter */
    if(! sf_getint("nf",   &nf))    nf=32;   /* factor coefficients */
	
    /*------------------------------------------------------------*/
    /* allocate autocorrelation */
    ain=sf_floatalloc2(sf_n(a1),sf_n(a2));
	
    /* read autocorrelation */
    sf_floatread(ain[0],sf_n(a1)*sf_n(a2),Fa);
	
    /* display autocorrelation */
    if(verb) {
		sf_warning("=======================================");
		for    (i2=0; i2<sf_n(a2); i2++) {
			for(i1=0; i1<sf_n(a1); i1++) {
				sf_warning("%2d %2d: %3.6g",i2,i1,ain[i2][i1]);
			}
		}
		sf_warning("=======================================");
    }
	
    /*------------------------------------------------------------*/
    /* allocate H.F. */
    na = (sf_n(a1)*sf_n(a2) - 1) / 2;
    aa = sf_allocatehelix(na);
    
    ff = sf_allocatehelix(nf);
    for(jf=0;jf<nf/2;jf++) {
		ff->lag[jf]=jf+1;
		ff->lag[nf/2+jf]=nn-nf/2+jf+1;
    }
    ff->h0=0.;
    for(jf=0;jf<nf/2;jf++) {
		ff->flt[jf]=0.;
    }
	
    /*------------------------------------------------------------*/
    /* transfer autocorrelation into H.F. */
    k1=(sf_n(a1)+1)/2-1;
    k2=(sf_n(a2)+1)/2-1;
	
    aa->h0=ain[k2][k1]/2.;
    kk=0;
    for    (i2=0; i2<sf_n(a2); i2++) {
		j2=(i2-k2)*nn;
		
		for(i1=0; i1<sf_n(a1); i1++) {
			j1=i1-k1;
			
			if(j2+j1>0) {
				aa->lag[kk]=j2+j1;
				aa->flt[kk]=ain[i2][i1];
				kk++;
			}
		}
    }
	
    /* display autocorrelation */
    if(verb) sf_displayhelix(aa);
	
    maxlag = 0;
    for( ja=0; ja<na; ja++) {
		if (aa->lag[ja] > maxlag) maxlag = aa->lag[ja];
    }
	
    /*------------------------------------------------------------*/
    /* factorize by Wilson-Burg */
    if (stable) {
		wilson2_init( maxlag*10);
		ff->h0 = wilson2_factor(niter, 2.*aa->h0, aa, 1, ff, verb, 1.e-6);
    } else {
		wilson_init(  maxlag*10);
		ff->h0 = wilson_factor (niter, 2.*aa->h0, aa,    ff, verb, 1.e-6);
    }
	
    /*------------------------------------------------------------*/
    /* scale factor */
    for( jf=0; jf<nf; jf++) {
		ff->flt[jf] *= ff->h0;
    }
    
    /* display factor */
    if(verb) sf_displayhelix(ff);
	
    /*------------------------------------------------------------*/
    /* write factor */
	
    /* filter */
    sf_putfloat(Ff,"a0",ff->h0);
    af=sf_maxa(nf,0,1); sf_oaxa(Ff,af,1);
    aj=sf_maxa( 1,0,1); sf_oaxa(Ff,aj,2);
	
    file = sf_getstring("lag"); 
    sf_putstring(Ff,"lag",file);
	
    sf_floatwrite(ff->flt,ff->nh,Ff);
	
    /* lags */
    sf_putint(Fl,"n",nn);
    sf_settype(Fl,SF_INT);
    af=sf_maxa(nf,0,1); sf_oaxa(Fl,af,1);
    aj=sf_maxa( 1,0,1); sf_oaxa(Fl,aj,2);
    sf_intwrite(ff->lag,ff->nh,Fl);
	
    /*------------------------------------------------------------*/
    /* deallocate H.F. */
    sf_deallocatehelix(ff);
    sf_deallocatehelix(aa);
	
	
    exit (0);
}
