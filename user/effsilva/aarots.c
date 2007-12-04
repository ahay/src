#include <math.h>
#include <stdlib.h>
#include <stdio.h>
#include <assert.h>
#include <rsf.h>

#include "aarots.h"

/****************************************************************************/
#ifndef _aarots_h

typedef struct aaliasstruct aalias;
/*^*/

struct aaliasstruct{ 
    int nf;
    float fmax;
    float df;
    float *fcut;
    float **ftap;
    float dt;
    int nt;
    int na;
    int hsign; 
};
/*^*/

#endif

#ifndef SWAP
#define SWAP(a,b)     tempr=(a);(a)=(b);(b)=tempr
#endif


/*----------------------------------------------------------------------------*/
double ampd( double *tr, double dt, double ot, double t)
/*< ampd >*/
{
    int it;
    double a;
    
    it = floor((t-ot)/dt);
    a  = tr[it] + (tr[it+1]-tr[it])*(t-it*dt-ot)/dt;
    return a;    
}

/*----------------------------------------------------------------------------*/
double ampf( float *tr, float dt, float ot, float t )
/*< ampf >*/
{
    int it;
    float a;
    
    it = floor((t-ot)/dt);
    a  = tr[it] + (tr[it+1]-tr[it])*(t-it*dt-ot)/dt;
    return a;   
}

/*----------------------------------------------------------------------------*/
int initAalias( int hsign, 
		bool verb, 
		float fmax, 
		float df, 
                int nt,    
		float dt,  
		aalias *a )
/*< init >*/
{
    int na, na2, nf, i;
    float ILOG2=1.442695040888963, **ftap=NULL;
    float *fcut;
    
    a->hsign=hsign;
    nf = floor(fmax/df);  if(nf<1) nf=1;
    a->nf=nf;
    
    a->fcut = fcut = (float *) malloc(nf*sizeof(float));
    assert(fcut!=NULL);
    
    for(i=0;i<nf;i++) {
	fcut[i] = fmax-i*df;
	if(verb) fprintf(stderr,"i=%d freq=%f\n",i,fcut[i]);
    }
    
    na = (int) pow(2.0, floor(log(nt-0.01)*ILOG2) + 1.0);
    na2 = na/2;
    a->na = na; 
    a->nt = nt;
    a->dt = dt;
    
    if(verb) fprintf(stderr,"na=%d na2=%d\n",na,na2);
    a->ftap = ftap = (float **) malloc(nf*sizeof(float)); assert(ftap!=NULL);
    for(i=0;i<nf;i++) {
	ftap[i] = (float *) malloc(na2*sizeof(float));
     assert(ftap[i]!=NULL); }
    
    loadftap(na,dt,nf,fcut,ftap);
    
    return nf;
}

/*---------------------------------------------------------------------------*/
void loadftap( int na, float dt, int nf, float *cutf, float **ftap )
/*< load >*/
{
    int jf, jw, nn2;
    float cf1, cf2, w1, dw, fnyq, PI=3.141592653589793238;
    
    nn2  = na/2;
    fnyq = .5/dt;
    dw   = 2.0*fnyq/(na-1);
    
    for(jf=0;jf<nf;jf++){
	cf1 = 0.75*cutf[jf];
	cf2 = cutf[jf];
	for(jw=0;jw<nn2;jw++){
	    w1  = jw*dw;
	    if(w1<cf1){
		ftap[jf][jw] = 1;
	    }else{
		if(w1>cf2) {
		    ftap[jf][jw] = 0;
		}else{
		    ftap[jf][jw] =  0.5*(1+cos(PI*(1-4.0*(cf2-w1)/cf2)));
		}
	    }
	}
    }
}

/*----------------------------------------------------------------------------*/
void loadBank( aalias aa, float *dat, float **bank )
/*< load >*/
{
    bankHalfd(aa.nt,dat,aa.dt,aa.nf,aa.fcut,aa.na,aa.ftap,bank,aa.hsign);
}

/*---------------------------------------------------------------------------*/
void bankHalfd( int nt, float *data,  float dt, int nfreq, float *cutf,
                int na, float **ftap, float **bank, int ihalfsign )
/*< bank >*/
{
    int it, jf, jw, ifftsign, i, nn2, nn;
    float dw, fnyq, a, tempr;
    float *y, *xx;
    float isq2=0.7071067811865475244;
    
    nn2 = na/2; nn=na;
    
    y  = (float *) calloc(2*na,sizeof(float) );
    xx = (float *) calloc(2*na,sizeof(float) );
    
    for(it=0;it<nt;it++) y[2*it]=data[it];
    ifftsign = -1;
    four1(y-1,(unsigned long) na, ifftsign);
    
    fnyq = .5/dt;
    dw = 2.0 * fnyq / (na-1);
    
    y[0] = 0.0;   y[1] = 0.0;
    
    /******* Positive frequencies  */
    for(i=1;i<nn2+1;i++){
	a        = sqrt(i*dw) * isq2;
	tempr    = y[2*i];
	y[2*i]   = a*(y[2*i] - ihalfsign*y[2*i+1]);
	y[2*i+1] = a*(y[2*i+1] + ihalfsign*tempr);
    }
    
    /******* Negative frequencies  */
    for(i=nn2+1;i<na;i++){
	a        = sqrt((na-i)*dw) * isq2;
	tempr    = y[2*i];
	y[2*i]   = a*(y[2*i] + ihalfsign*y[2*i+1]);
	y[2*i+1] = a*(y[2*i+1] - ihalfsign*tempr);
    }
    
    for(jf=0;jf<nfreq;jf++) {
	for(jw=0;jw<nn2;jw++) {
	    xx[2*jw]      = ftap[jf][jw] * y[2*jw];   /* real part */
	    xx[2*jw+1]    = ftap[jf][jw] * y[2*jw+1]; /* imaginary part */
	    xx[2*(nn-jw)-2] = ftap[jf][jw] * y[2*(nn-jw)-2]; /* real part*/
	    xx[2*(nn-jw)-1] = ftap[jf][jw] * y[2*(nn-jw)-1]; /* imag part*/
	}
	ifftsign = 1;
	four1(xx-1, (unsigned long) nn, ifftsign);
	for(it=0;it<nt;it++) bank[jf][it] = (float) xx[2*it]/nn;
    }
    free(y);
    free(xx);   
}
/*************************************************************************
**                                                                      **
**  AUTHORS                                                             **
**                                                                      **
**    Rodrigo de Souza Portugal ( C version )                           **
**    Lucio Tunes dos Santos ( Matlab version )                         **
**                                                                      **
**    Mailing address: DMA - IMECC - UNICAMP                            **
**                     Caixa Postal 6065                                **
**                     13081-970 Campinas, SP, Brasil                   **
**                                                                      **
**    email: portugal@ime.unicamp.br                                    **
**           lucio@ime.unicamp.br                                       **
**                                                                      **
**                                                                      **
**  SUPPORT                                                             **
**                                                                      **
**    Research Foundation of the State of Sao Paulo (FAPESP-Brazil),    **
**    National Council for Scientific and Technological Development     **
**    (CNPq - Brazil)                                                   **
**    Ministry of Science and Technology (PRONEX - Brazil).             **
**                                                                      **
**                                                                      **
**  LAST REVISION: 08/10/2001                                           **
**                                                                      **
*************************************************************************/

void four1(float data[], unsigned long nn, int isign)
/*< four >*/
{
    unsigned long n,mmax,m,j,istep,i;
    double wtemp,wr,wpr,wpi,wi,theta;
    float tempr,tempi;
    
    n=nn << 1;
    j=1;
    for (i = 1; i < n; i+=2) {
	if (j > i) {
	    SWAP(data[j],data[i]);
	    SWAP(data[j+1],data[i+1]);
	}
	m=n >> 1;
	while (m >= 2 && j > m) {
	    j -= m;
	    m >>= 1;
	}
	j += m;
    }
    mmax=2;
    while (n > mmax) {
	istep = mmax << 1;
	theta = isign * (6.28318530717959/mmax);
	wtemp = sin(0.5 * theta);
	wpr   = -2.0 * wtemp * wtemp;
	wpi   = sin(theta);
	wr    = 1.0;
	wi    = 0.0;
	for (m=1;m<mmax;m+=2) {
	    for (i=m;i<=n;i+=istep) {
		j         =  i+mmax;
		tempr     =  wr * data[j]   - wi * data[j+1];
		tempi     =  wr * data[j+1] + wi * data[j];
		data[j]   =  data[i]   - tempr;
		data[j+1] =  data[i+1] - tempi;
		data[i]   += tempr;
		data[i+1] += tempi;
	    }
	    wr = (wtemp=wr) * wpr - wi    * wpi + wr;
	    wi = wi         * wpr + wtemp * wpi + wi;
	}
	mmax=istep;
    }
}
