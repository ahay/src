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
#define cWGH(a,b,c) (+1.*(conjf(a)*b)*c)
#define cCOR(a,b) (conjf(a)*b) 
#define cMUL(a,b) (a*b) 
#else
#define cWGH(a,b,c) (+1.*sf_cmul((sf_cmul(conjf(a),b)),c))
#define cCOR(a,b) (sf_cmul(conjf(a),b))
#define cMUL(a,b) (sf_cmul(a,b))
#endif

int main(int argc, char* argv[])
{
    bool verb,rflg,buf;

    sf_file Fs,Fr,Fi,Fc;           /* I/O files */
    sf_axis az,ax,at,ac,aw,aa;     /* cube axes */
    int     nz,nx,nt,nw,nhz,nhx,nht,nc;
    int           it,iw,ihz,ihx,iht,ic;
    
    float oht,dht,ht,w;
    int   nhx2,nhz2,nht2;

    float      ****ii=NULL;
    sf_complex ****jf=NULL;
    float      ***r_us=NULL,***r_ur=NULL; 
    sf_complex **c_us=NULL,**c_ur=NULL; 
    sf_complex   **tt=NULL;
    sf_complex   uu;
    sf_complex   wt;

    pt2d *cc=NULL;
    bool *ccin=NULL;
    float czmin,cxmin;
    float czmax,cxmax;
    int   icz,  icx;
    int   mcz,  mcx, mct;
    int   pcz,  pcx, pct;

    int **mcxall, **pcxall;
    int **mczall, **pczall;
    int  *mctall,  *pctall;
    int lht;

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
    if(! sf_getbool("buf", &buf)) buf=false;
    
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
	aw=sf_maxa(0,0,0); nw=0;
    } else {
	at=sf_maxa(0,0,0); nt=0;
	aw=sf_iaxa(Fs,3); sf_setlabel(aw,"w"); nw = sf_n(aw);
    }

    /* CIP coordinates */
    ac = sf_iaxa(Fc,2); sf_setlabel(ac,"cc"); sf_setunit(ac,"");
    nc = sf_n(ac);

    if(! sf_getint("nhz",&nhz)) nhz=0; /* number of lags on the z axis */
    if(! sf_getint("nhx",&nhx)) nhx=0; /* number of lags on the x axis */
    if(! sf_getint("nht",&nht)) nht=0; /* number of lags on the t axis */

    nhx2=2*nhx+1;
    nhz2=2*nhz+1;
    nht2=2*nht+1;

    if(verb) {
	sf_warning("nhz=%3d nhx=%3d nht=%3d",nhz2,nhx2,nht2);

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
    aa=sf_maxa(nhz2,-nhz*sf_d(az),sf_d(az));
    sf_setlabel(aa,"hz"); sf_setunit(aa,"");
    if(verb) sf_raxa(aa);
    sf_oaxa(Fi,aa,1);
    
    aa=sf_maxa(nhx2,-nhx*sf_d(ax),sf_d(ax)); 
    sf_setlabel(aa,"hx"); sf_setunit(aa,"");
    if(verb) sf_raxa(aa);
    sf_oaxa(Fi,aa,2);

    if(rflg) {
	aa=sf_maxa(nht2,-nht*sf_d(at),sf_d(at));
	sf_setlabel(aa,"ht"); sf_setunit(aa,"");
	oht=-nht*sf_d(at);
    } else {
	if(! sf_getfloat("dht",&dht)) sf_error("need dht");
	oht=-nht*dht;
	aa=sf_maxa(nht2,oht,dht);
	sf_setlabel(aa,"ht"); sf_setunit(aa,"s");    
    }
    if(verb) sf_raxa(aa);
    sf_oaxa(Fi,aa,3);

    sf_oaxa(Fi,ac,4);

    /* work arrays */
    if(rflg) {

	r_us=sf_floatalloc3(nz,nx,nht2);
	r_ur=sf_floatalloc3(nz,nx,nht2);

    } else {

	jf  =sf_complexalloc4(nht2,nhz2,nhx2,nc);

	c_us=sf_complexalloc2(nz,nx);
	c_ur=sf_complexalloc2(nz,nx);

	/* precompute phase for time delay */
	tt  =sf_complexalloc2(nht2,nw);
	for(iw=0;iw<nw;iw++) {
	    w = -2*SF_PI*( sf_o(aw)+iw*sf_d(aw) );
	    
	    for(iht=0;iht<nht2;iht++) {
		ht = oht+iht*dht;
		
		tt[iw][iht] = sf_cmplx(cosf(2*w*ht),sinf(2*w*ht)); 
	    }
	}
    }

    ii=sf_floatalloc4(nhz2,nhx2,nht2,nc);

    /*------------------------------------------------------------*/
    /* CIP coordinates */
    cc= (pt2d*) sf_alloc(nc,sizeof(*cc));
    pt2dread1(Fc,cc,nc,2); /* read coordinates */

    mcxall=sf_intalloc2(nhx2,nc);
    pcxall=sf_intalloc2(nhx2,nc);
    mczall=sf_intalloc2(nhz2,nc);
    pczall=sf_intalloc2(nhz2,nc);

    ccin=sf_boolalloc(nc);

    cxmin = sf_o(ax) +             nhx *sf_d(ax);
    cxmax = sf_o(ax) + (sf_n(ax)-1-nhx)*sf_d(ax);
    czmin = sf_o(az) +             nhz *sf_d(az);
    czmax = sf_o(az) + (sf_n(az)-1-nhz)*sf_d(az);

    for(ic=0; ic<nc; ic++) {
	ccin[ic]=(cc[ic].x>=cxmin && cc[ic].x<=cxmax &&
		  cc[ic].z>=czmin && cc[ic].z<=czmax)?true:false;

	if(ccin[ic]) {

	    icx = (cc[ic].x-sf_o(ax))/sf_d(ax);
	    for(ihx=-nhx; ihx<nhx+1; ihx++) {
		mcxall[ic][nhx+ihx] = icx-ihx;
		pcxall[ic][nhx+ihx] = icx+ihx;
	    }

	    icz = (cc[ic].z-sf_o(az))/sf_d(az);
	    for(ihz=-nhz; ihz<nhz+1; ihz++) {
		mczall[ic][nhz+ihz] = icz-ihz;
		pczall[ic][nhz+ihz] = icz+ihz;
	    }

	}
    }
   
    if(rflg) {
	lht=2*nht;

	mctall=sf_intalloc(nht2);
	pctall=sf_intalloc(nht2);
	for (iht=0; iht<nht2; iht++) { 
	    mctall[iht]=iht;
	    pctall[iht]=2*nht-iht;
	}
    }


    if(rflg) {
	/*------------------------------------------------------------*/
	/* t-domain */
	/*------------------------------------------------------------*/

#ifdef _OPENMP
#pragma omp parallel for schedule(dynamic)	\
    private(ic,iht, ihx, ihz)			\
    shared (nc,nht2,nhx2,nhz2,ii)
#endif	
	/* zero output */
	for(ic=0; ic<nc; ic++) {
	    for        (iht=0; iht<nht2; iht++) {
		for    (ihx=0; ihx<nhx2; ihx++) {
		    for(ihz=0; ihz<nhz2; ihz++) {
			ii[ic][iht][ihx][ihz] = 0;
		    }
		}
	    }
        }

	/* read wavefield @ [0...2nht-1]*/
	sf_floatread  (r_us[0][0],nz*nx*(2*nht),Fs);
	sf_floatread  (r_ur[0][0],nz*nx*(2*nht),Fr);

	if(verb) fprintf(stderr,"nt\n");
	for(it=nht;it<nt-nht;it++) { /* loop over time */
	    if(verb) fprintf(stderr,"\b\b\b\b\b\b\b\b\b\b%04d",nt-nht-1-it);
	    
		/* read wavefield @ [2nht]*/
		sf_floatread(r_us[ lht ][0],nz*nx,Fs);
		sf_floatread(r_ur[ lht ][0],nz*nx,Fr);
		
#ifdef _OPENMP
#pragma omp parallel for schedule(dynamic)				\
    private(ic,iht, ihx, ihz,              mcx,   pcx,   mcz,   pcz,   mct,   pct   ) \
    shared (nc,nht2,nhx2,nhz2,ii,r_us,r_ur,mcxall,pcxall,mczall,pczall,mctall,pctall,ccin)
#endif
		for(ic=0; ic<nc; ic++) { /* loop over CIPs */
		    if(ccin[ic]) {
			
			for        (iht=0; iht<nht2; iht++) { mct=mctall[iht];     pct=pctall[iht];
			    for    (ihx=0; ihx<nhx2; ihx++) { mcx=mcxall[ic][ihx]; pcx=pcxall[ic][ihx];
				for(ihz=0; ihz<nhz2; ihz++) { mcz=mczall[ic][ihz]; pcz=pczall[ic][ihz];
				    
				    ii[ic][iht][ihx][ihz] += rCOR(r_us[mct][mcx][mcz],
								  r_ur[pct][pcx][pcz]);
				    
				} /* ihz */
			    }     /* ihx */
			}         /* iht */
			
		    } /* end if */
		} /* ic */
		
		/* update wavefield index (cycle) */
		for(iht=0;iht<nht2;iht++) {
		    mctall[iht] = (mctall[iht]+1) % nht2;
		    pctall[iht] = (pctall[iht]+1) % nht2;
		}
		lht = (lht+1) % nht2;
		
	} /* it */
	if(verb) fprintf(stderr,"\n");
	
        /* write image */
        sf_floatwrite(ii[0][0][0],nc*(nht2)*(nhz2)*(nhx2),Fi);

    } else {
	/*------------------------------------------------------------*/
	/* f-domain */
	/*------------------------------------------------------------*/

#ifdef _OPENMP
#pragma omp parallel for schedule(dynamic)	\
    private(ic,iht, ihx, ihz)			\
    shared (nc,nht2,nhx2,nhz2,jf)
#endif			
	/* zero output */
	for(ic=0; ic<nc; ic++) {
	    for        (ihx=0; ihx<nhx2; ihx++) {
		for    (ihz=0; ihz<nhz2; ihz++) {
		    for(iht=0; iht<nht2; iht++) {
			jf[ic][ihx][ihz][iht] = sf_cmplx(0.,0.);
		    }
		}
	    }
	}

	if(verb) fprintf(stderr,"nw\n");
	for(iw=0; iw<nw; iw++) { /* loop over frequency */
	    if(verb) fprintf(stderr,"\b\b\b\b\b\b\b\b\b\b%04d",nw-1-iw);
	    
	    /* read wavefield */
	    sf_complexread  (c_us[0],nz*nx,Fs);
	    sf_complexread  (c_ur[0],nz*nx,Fr);
	    
#ifdef _OPENMP
#pragma omp parallel for schedule(dynamic)				\
    private(ic,iht, ihx, ihz,              mcx,   pcx,   mcz,   pcz   ,wt,uu) \
    shared (nc,nht2,nhx2,nhz2,jf,c_us,c_ur,mcxall,pcxall,mczall,pczall,ccin)
#endif			
	    for(ic=0; ic<nc; ic++) { /* loop over CIPs */		    
		if(ccin[ic]) {

		    for        (ihx=0; ihx<nhx2; ihx++) { mcx=mcxall[ic][ihx]; pcx=pcxall[ic][ihx];
			for    (ihz=0; ihz<nhz2; ihz++) { mcz=mczall[ic][ihz]; pcz=pczall[ic][ihz];

			    uu = cCOR( c_us[mcx][mcz], c_ur[pcx][pcz]);

			    for(iht=0; iht<nht2; iht++) { wt = tt[iw][iht];
#ifdef SF_HAS_COMPLEX_H
				jf[ic][ihx][ihz][iht] += cMUL(uu,wt);
#else
				jf[ic][ihx][ihz][iht] = sf_cadd(jf[ic][ihx][ihz][iht],cMUL(uu,wt));
#endif
			    }
			}
		    }

		} /* end if */
	    } /* ic */
	    
	} /* iw */
	if(verb) fprintf(stderr,"\n");
	
#ifdef _OPENMP
#pragma omp parallel for schedule(dynamic)	\
    private(ic,iht, ihx, ihz)			\
    shared (nc,nht2,nhx2,nhz2,jf,ii)
#endif			
	/* write image */
	for(ic=0; ic<nc; ic++) {
	    if(ccin[ic]) {
		for        (ihx=0; ihx<nhx2; ihx++) {
		    for    (ihz=0; ihz<nhz2; ihz++) {
			for(iht=0; iht<nht2; iht++) {
			    ii[ic][iht][ihx][ihz] = -crealf(jf[ic][ihx][ihz][iht]);
			}
		    }
		}
	    }
	}
	sf_floatwrite(ii[0][0][0],nc*(nht2)*(nhz2)*(nhx2),Fi);
    }

    /*------------------------------------------------------------*/
    free(***ii); free(**ii); free(*ii);free(ii);

	
     if(rflg) {
	free(*r_us); free(r_us);
	free(*r_ur); free(r_ur);
    } else {
	free(*c_us); free(c_us);
	free(*c_ur); free(c_ur);
	free(*tt);   free(tt);
	free(***jf); free(**jf); free(*jf);free(jf);
    }

     free(cc);
     free(ccin);

     free(*mcxall); free(mcxall);
     free(*pcxall); free(pcxall);
     free(*mczall); free(mczall);
     free(*pczall); free(pczall);

     /*------------------------------------------------------------*/

 
     exit (0);
}
