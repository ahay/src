/* OpenMP lagged-products in the time- or freq-domain */

/*
  Copyright (C) 2009 Colorado School of Mines
  
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

/* 
 * inputs:      wavefields organized as z-x-t(w)
 * output: lagged products organized as hz-hx-ht @ irregular points
 */

#include <rsf.h>

#ifdef _OPENMP
#include <omp.h>
#include "omputil.h"
#endif

#define rCOR(a,b) (a*b)

#ifdef SF_HAS_COMPLEX_H
#define cWGH(a,b,c) (-1.*(conj(a)*b)*c)
#else
#define cWGH(a,b,c) (-1.*sf_cmul((sf_cmul(conj(a),b)),c))
#endif

int main(int argc, char* argv[])
{
    bool verb,rflg;

    sf_file Fs,Fr,Fi,Fc;           /* I/O files */
    sf_axis az,ax,at,ac,aw,aa;     /* cube axes */
    int     nz,nx,nt,nw,nhz,nhx,nht,nc;
    int           it,iw,ihz,ihx,iht,ic;
    
    float oht,dht,ht,w;

    float       ***ii=NULL;
    sf_complex  ***jj=NULL;
    float      **r_us=NULL,**r_ur=NULL; 
    sf_complex **c_us=NULL,**c_ur=NULL; 
    sf_complex   **tt=NULL;
    sf_complex   wt;

    pt2d *cc=NULL;
    int  *ccix=NULL, *cciz=NULL;
    bool *ccin=NULL;
    float czmin,cxmin;
    float czmax,cxmax;
    int   icz,  icx;
    int   mcz,  mcx;
    int   pcz,  pcx;
    int   ocz,  ocx;
    
#ifdef _OPENMP
    int ompnth=0;
#endif

    /*------------------------------------------------------------*/
    /* init RSF */
    sf_init(argc,argv);

    /* OMP parameters */
#ifdef _OPENMP
    ompnth=omp_init();
#endif

    if(! sf_getbool("verb",&verb)) verb=false;         /* verbosity flag */
    
    Fs = sf_input ("in" ); /*   source wavefield */
    Fr = sf_input ("ur" ); /* receiver wavefield */
    Fc = sf_input ("cc" ); /* CIP coordinates    */
    Fi = sf_output("out"); /* image */

    rflg = (sf_gettype(Fs) == SF_COMPLEX)?false:true;
    if(rflg) {
	sf_warning("real input");
    } else {
	sf_warning("complex input");
	sf_settype(Fi,SF_FLOAT);
    }

    /* read axes */
    az=sf_iaxa(Fs,1); sf_setlabel(az,"z"); nz = sf_n(az);
    ax=sf_iaxa(Fs,2); sf_setlabel(ax,"x"); nx = sf_n(ax);

    if(rflg){
	at=sf_iaxa(Fs,3); sf_setlabel(at,"t"); nt = sf_n(at);
    } else {
	aw=sf_iaxa(Fs,3); sf_setlabel(aw,"w"); nw = sf_n(aw);
    }

    /* CIP coordinates */
    ac = sf_iaxa(Fc,2); sf_setlabel(ac,"cc"); sf_setunit(ac,"");
    if(verb) sf_raxa(ac); 
    nc = sf_n(ac);

    if(! sf_getint("nhz",&nhz)) nhz=0; /* number of lags on the z axis */
    if(! sf_getint("nhx",&nhx)) nhx=0; /* number of lags on the x axis */
    if(! sf_getint("nht",&nht)) nht=0; /* number of lags on the t axis */

    if(verb) {
	sf_warning("nhz=%3d nhx=%3d nht=%3d",2*nhz+1,2*nhx+1,2*nht+1);

	sf_raxa(az);
	sf_raxa(ax);
	if(rflg){
	    sf_raxa(at);
	} else {
	    sf_raxa(aw);
	}
	sf_raxa(ac);
    }

    /* set output axes */
    sf_oaxa(Fi,ac,1);

    aa=sf_maxa(2*nhz+1,-nhz*sf_d(az),sf_d(az));
    sf_setlabel(aa,"hz"); sf_setunit(aa,"");
    if(verb) sf_raxa(aa);
    sf_oaxa(Fi,aa,2);
    
    aa=sf_maxa(2*nhx+1,-nhx*sf_d(ax),sf_d(ax)); 
    sf_setlabel(aa,"hx"); sf_setunit(aa,"");
    if(verb) sf_raxa(aa);
    sf_oaxa(Fi,aa,3);

    if(rflg) {
	aa=sf_maxa(2*nht+1,-nht*sf_d(at),sf_d(at));
	sf_setlabel(aa,"ht"); sf_setunit(aa,"");
    } else {
	if(! sf_getfloat("dht",&dht)) sf_error("need dht");
	oht=-nht*dht;
	aa=sf_maxa(2*nht+1,oht,dht);
	sf_setlabel(aa,"ht"); sf_setunit(aa,"s");    
    }
    if(verb) sf_raxa(aa);
    sf_oaxa(Fi,aa,4);

    /* work arrays */
    if(rflg) {
	r_us=sf_floatalloc2(nz,nx);
	r_ur=sf_floatalloc2(nz,nx);
    } else {
	c_us=sf_complexalloc2(nz,nx);
	c_ur=sf_complexalloc2(nz,nx);
	tt  =sf_complexalloc2(2*nht+1,nw);
	jj  =sf_complexalloc3(nc,2*nhz+1,2*nhx+1);

	/* precompute phase shift */
	for(iht=0;iht<2*nht+1;iht++) {
	    ht = oht+iht*dht;
	    
	    for(iw=0;iw<nw;iw++) {
		w = 2*SF_PI*( sf_o(aw)+iw*sf_d(aw) );
		tt[iht][iw] = sf_cmplx(cosf(2*w*ht),sinf(2*w*ht)); 
/*		sf_warning("%f %f %f+i%f",w,ht,cosf(-w*ht),sinf(-w*ht));*/
	    }
	}
    }

    ii=sf_floatalloc3(nc,2*nhz+1,2*nhx+1);

    /*------------------------------------------------------------*/
    /* CIP coordinates */
    cc= (pt2d*) sf_alloc(nc,sizeof(*cc));
    pt2dread1(Fc,cc,nc,2); /* read coordinates */

    ccix=sf_intalloc(nc);
    cciz=sf_intalloc(nc);
    ccin=sf_boolalloc(nc);

    cxmin = sf_o(ax) +             nhx *sf_d(ax);
    cxmax = sf_o(ax) + (sf_n(ax)-1-nhx)*sf_d(ax);
    czmin = sf_o(az) +             nhz *sf_d(az);
    czmax = sf_o(az) + (sf_n(az)-1-nhz)*sf_d(az);

    for(ic=0; ic<nc; ic++) {
	ccix[ic]= (cc[ic].x-sf_o(ax))/sf_d(ax);
	cciz[ic]= (cc[ic].z-sf_o(az))/sf_d(az);

	ccin[ic]=(cc[ic].x>=cxmin && cc[ic].x<=cxmax &&
		  cc[ic].z>=czmin && cc[ic].z<=czmax)?true:false;
    }

    /*------------------------------------------------------------*/
    if(rflg) {

	if(verb) fprintf(stderr,"nht   nt\n");
	if(verb) fprintf(stderr,"%03d %04d\n",2*nht,nt-2*nht);
	/* loop over time lags */
	for(iht=-nht; iht<nht+1; iht++) {
	    
	    /* zero output */
	    for    (ihx=0; ihx<2*nhx+1; ihx++) {
		for(ihz=0; ihz<2*nhz+1; ihz++) {
		    for    (ic=0; ic<nc; ic++) {
			ii[ihx][ihz][ic] = 0;
		    }
		}
	    }
	    
	    sf_seek(Fs,(nht-iht)*nz*nx*sizeof(float),SEEK_SET);
	    sf_seek(Fr,(nht+iht)*nz*nx*sizeof(float),SEEK_SET);
	    
	    /* loop over time */
	    for(it=0; it<nt-2*nht; it++) {
		if(verb) fprintf(stderr,"\b\b\b\b\b\b\b\b\b\b%03d %04d",nht+iht,it);
		
		/* read wavefield */
		sf_floatread  (r_us[0],nz*nx,Fs);
		sf_floatread  (r_ur[0],nz*nx,Fr);
		
#ifdef _OPENMP
#pragma omp parallel for schedule(dynamic)		\
    private(ihx,ihz,mcx,mcz,pcx,pcz,ocx,ocz,ic,icx,icz)	\
    shared(nhx,nhz,ii,r_us,r_ur,ccix,cciz,ccin)
#endif
		
		/* loop over CIPs */		    
		for(ic=0; ic<nc; ic++) {
		    
		    if(ccin[ic]) {
			icx=ccix[ic];
			icz=cciz[ic];		     
			
			/* loop over space lags */
			for    (ihx=-nhx; ihx<nhx+1; ihx++) { mcx=icx-ihx; pcx=icx+ihx; ocx=nhx+ihx;
			    for(ihz=-nhz; ihz<nhz+1; ihz++) { mcz=icz-ihz; pcz=icz+ihz; ocz=nhz+ihz;
				
				ii[ocx][ocz][ic] += rCOR(r_us[mcx][mcz],r_ur[pcx][pcz]);
			    } /* ihz */
			}     /* ihx */
			
		    } /* end if */
		    
		} /* ic */	
	    } /* it */
	    
	    /* write image */
	    sf_floatwrite(ii[0][0],nc*(2*nhz+1)*(2*nhx+1),Fi);
	    
	} /* iht */
	if(verb) fprintf(stderr,"\n");
    
    } else {

	if(verb) fprintf(stderr,"nht   nw\n");
	if(verb) fprintf(stderr,"%03d %04d\n",2*nht,nw-1);
	/* loop over time lags */
	for(iht=0; iht<2*nht+1; iht++) {

	    /* zero output */
	    for    (ihx=0; ihx<2*nhx+1; ihx++) {
		for(ihz=0; ihz<2*nhz+1; ihz++) {
		    for    (ic=0; ic<nc; ic++) {
			jj[ihx][ihz][ic] = 0;
		    }
		}
	    }
	    
	    sf_seek(Fs,0,SEEK_SET);
	    sf_seek(Fr,0,SEEK_SET);

	    /* loop over frequency */
	    for(iw=0; iw<nw; iw++) {
		if(verb) fprintf(stderr,"\b\b\b\b\b\b\b\b\b\b%03d %04d",iht,iw);

		wt = tt[iht][iw];

		/* read wavefield */
		sf_complexread  (c_us[0],nz*nx,Fs);
		sf_complexread  (c_ur[0],nz*nx,Fr);

#ifdef _OPENMP
#pragma omp parallel for schedule(dynamic)		\
    private(ihx,ihz,mcx,mcz,pcx,pcz,ocx,ocz,ic,icx,icz)	\
    shared(nhx,nhz,ii,c_us,c_ur,ccix,cciz,ccin,wt)
#endif

		/* loop over CIPs */		    
		for(ic=0; ic<nc; ic++) {
		    
		    if(ccin[ic]) {
			icx=ccix[ic];
			icz=cciz[ic];		     
			
			/* loop over space lags */
			for    (ihx=-nhx; ihx<nhx+1; ihx++) { mcx=icx-ihx; pcx=icx+ihx; ocx=nhx+ihx;
			    for(ihz=-nhz; ihz<nhz+1; ihz++) { mcz=icz-ihz; pcz=icz+ihz; ocz=nhz+ihz;
				
				jj[ocx][ocz][ic] += cWGH(c_us[mcx][mcz],c_ur[pcx][pcz],wt);

			    } /* ihz */
			}     /* ihx */
			
		    } /* end if */
		    
		} /* ic */
		
	    } /* iw */
	    
	    /* write image */
	    for    (ihx=0; ihx<2*nhx+1; ihx++) {
		for(ihz=0; ihz<2*nhz+1; ihz++) {
		    for    (ic=0; ic<nc; ic++) {
			ii[ihx][ihz][ic] = crealf(jj[ihx][ihz][ic]);
		    }
		}
	    }
	    sf_floatwrite(ii[0][0],nc*(2*nhz+1)*(2*nhx+1),Fi);
	    
	} /* iht */
	if(verb) fprintf(stderr,"\n");

    }

    /*------------------------------------------------------------*/
    free(**ii); free(*ii); free(ii);

     if(rflg) {
	free(*r_us); free(r_us);
	free(*r_ur); free(r_ur);
    } else {
	free(*c_us); free(c_us);
	free(*c_ur); free(c_ur);
	free(*tt);   free(tt);
	free(**jj);  free(*jj); free(jj);
    }

     free(cc);
     free(ccix);
     free(cciz);
     free(ccin);

    /*------------------------------------------------------------*/
    exit (0);
}
