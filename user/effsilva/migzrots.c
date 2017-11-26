#include <math.h>
#include <stdlib.h>
#include <stdio.h>
#include <rsf.h>

#ifndef  _aarots_h
#include "aarots.h"
#endif

#include "migzrots.h"

/*--------------------------------------------------------------------------*/
void spreadSR( int nf,      float fmax,  float df,    
               float tgmax, float tgtap,
               int nt,      float ot,    float dt,
               int nx,      float ox,    float dx,
               int nz,      float oz,    float dz,
               float oztt,  float oxtt,
               float **px,  float **pz,
               float **tr,  float **ts, 
               float **datf, float **image )
/*< spread >*/
{
    int ix,iz,lx,lz,indf;
    float t,tts,ttr,tmin,tmax,w,tg,fatg,zero=0.0001,tg90=20,fc,p2;
    
    tmin = ot;
    tmax = ot+(nt-1)*dt;
    lx   = floor((ox-oxtt)/dx);
    lz   = floor((oz-oztt)/dz);
    fatg = 1/(tgmax-tgtap);
    
    for(ix=0;ix<nx;ix++) {
	for(iz=0;iz<nz;iz++) {
	    
	    tts = ts[ix+lx][iz+lz];
	    ttr = tr[ix+lx][iz+lz];
	    t   = tts+ttr;
	    if(t<tmax) {
		if(t>tmin) {
		    
		    if(fabs(pz[ix+lx][iz+lz])>zero) {
			tg   = fabs(px[ix+lx][iz+lz]/pz[ix+lx][iz+lz]);
		    }else{
			tg   = tg90;
		    }
		    
		    if(tg<tgmax) {
			
			p2   = px[ix+lx][iz+lz]*px[ix+lx][iz+lz] +
			    pz[ix+lx][iz+lz]*pz[ix+lx][iz+lz];
			
			w    = p2*sqrt(tts+ttr)*sqrt(tts)/sqrt(ttr*ttr*ttr);           
			
			if(tg>tgtap) w*=(tgmax-tg)*fatg;
			fc   = .25/(fabs(px[ix+lx][iz+lz])*dx);
			
			indf = floor((fmax-fc)/df); if(indf<0) indf=0;
			
			if(indf>(nf-1)) indf=nf-1;
			image[ix][iz] += w*ampf(datf[indf],dt,ot,t);
			
		    }
		}
	    }
	}
    }
}

/*...........................................................................*/
void derive_1( int n1, int n2, float d1, 
               float **y, float **dyd1 )
/*< derive 1 >*/
{
    int n=6,i1,i2;
    
    sf_deriv_init(n1, n, 0.);
    
    for (i2=0;i2<n2;i2++) {
	sf_deriv(y[i2],dyd1[i2]);
	for(i1=0;i1<n1;i1++) dyd1[i2][i1] /= d1;
    }
    
    sf_deriv_close();
    
}

/*...........................................................................*/
void derive_2( int n1, int n2, float d2,
               float **y, float **dyd2 )
/*< derive 2 >*/
{
    int n=6,i1,i2;
    float *yy=NULL,*der2=NULL;
/*
  hilbert_init(n2, n, 0.);
  der2 = sf_floatalloc(n2);
  
  for (i1=0;i1<n1;i1++) {
  for(i2=0;i2<n2;i2++) taux[i2] = tsg[i2*n1+i1];
  deriv(taux,der2);
  for(i2=0;i2<n2;i2++) dtdx[i2*n1+i1] = der2[i2]/d2;
  }
*/
    
    sf_deriv_init(n2, n, 0.);
    yy   = sf_floatalloc(n2);
    der2 = sf_floatalloc(n2);
    
    for (i1=0;i1<n1;i1++) {
	for(i2=0;i2<n2;i2++) yy[i2] = y[i2][i1];
	sf_deriv(yy,der2);
	for(i2=0;i2<n2;i2++) dyd2[i2][i1] = der2[i2]/d2;
    }
    
    free(yy); free(der2);
    sf_deriv_close();
    
}

