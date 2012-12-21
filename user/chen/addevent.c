/* Add a analog event with dispersion and attenuation */

/*
  Copyright (C) 2012 Zhonghuan Chen, UT Austin, Tsinghua University
  
  This program is free software; you can redistribute it and/or modify
  it under the terms of the GNU General Public License as published by
  the Free Software Foundation; either version 2 of the License, or
  (at your option) any later version.
  
  This program is distributed in the hope that it will be useful,
  but WITHOUT ANY WA:RRANTY; without even the implied warranty of
  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
  GNU General Public License for more details.
  
  You should have received a copy of the GNU General Public License
  along with this program; if not, write to the Free Software
  Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA
*/


#include <rsf.h>
#include "traveltime.h"
#include "wavelet.h"


typedef struct tag_addevent
{
    float	df;	/* time sample rate */
    int 	of;	/* time samples */
    int		nf;	/* frequency numbers (half) */
    sf_complex 	*wvlt;	/* frequency domain wavelet */
    int 	type;	/* event type */
    float 	t0;	/* event traveltime at zero offset (x=0) */
    float 	a0;	/* event amplitude at reference frequency at zero offset (x=0,f=f0) */
    float 	*vw; 	/* velocity over frequencies */
    float 	qa;	/* event amplitude attenuation */
    float	f0;	/* center frequency */
    int 	a0ref;  /* reference point for a0 */
}addevent;


void* sf_addevent_init(int nfft, float of, float df, 	/* time domain sampling parameters */
		       int wvtype,float *wvp,						/* wavelet parameters */
		       int event, float t0, float v0, float a0,	/* event initial parameters */
		       float qv, float qa, float f0, int a0ref )
/*< initialize >*/
{
    addevent * pp;
    int ifreq;
    float freq;

    pp = (addevent*) malloc(sizeof(addevent));

    pp->nf=nfft/2+1;
    pp->of=of;
    pp->df=df;
    pp->t0=t0;
    pp->type=event;
    pp->f0=f0;	
    pp->a0=M_2_SQRTPI*a0;
    pp->qa=qa;
    pp->a0ref=a0ref;

    pp->vw=sf_floatalloc( pp->nf );
    pp->wvlt=sf_complexalloc( pp->nf );

    for(ifreq=0; ifreq<pp->nf; ifreq++)
    {
	freq=fabs(ifreq*df+of)/f0+0.001;
	if(qv<=0) pp->vw[ifreq] = v0;		/* dispersion velocity */
	else {
	    pp->vw[ifreq]=v0*pow(freq,1.0/(SF_PI*qv));
	}
    }

    switch(wvtype)
    {
	case 0:
	    for(ifreq=0;ifreq<pp->nf;ifreq++)
	    {
		pp->wvlt[ifreq]=sf_cmplx(ifreq*df+of,0.f);
	    }
	    sf_wvlt_frck(pp->nf, pp->wvlt, wvp);
	    break;
	case 1:
	    for(ifreq=0;ifreq<pp->nf;ifreq++)
	    {
		if(fabs(ifreq*df+of) < wvp[0]) pp->wvlt[ifreq] = sf_cmplx(1.0f,0.0f);
		else pp->wvlt[ifreq] = sf_cmplx(0.0f,0.0f);
	    }
	    break;
	default:
	    break;
    }
    return ( pp);
}

void sf_addevent(void* p, float x, sf_complex *ftr)
/*< modeling for one trace >*/
{
    int ifreq;
    float txw,freq,fa,deltat;
    addevent *pp;
    sf_complex fc1;

    pp = (addevent*) p;
    for(ifreq=0; ifreq<pp->nf; ifreq++)
    {
	freq = ifreq*pp->df+pp->of;
	txw  = sf_traveltime(pp->type, pp->t0, x, pp->vw[ifreq]);
	if(pp->qa > 0)
	{
	    if(pp->a0ref == 0) deltat = txw;
	    else deltat = txw-pp->t0;
	    fa = pp->a0*exp(-SF_PI/fabs(freq*deltat));
	}else{
	    fa = pp->a0;
	}
	fc1 = sf_cmplx(fa*cos(2.0*SF_PI*freq*txw),-fa*sin(2.0*SF_PI*freq*txw));
#ifdef SF_HAS_COMPLEX_H
	ftr[ifreq] += fc1*pp->wvlt[ifreq];
#else
	ftr[ifreq] = sf_cadd(ftr[ifreq],sf_cmul(fc1,pp->wvlt[ifreq]));
#endif
    }
}

void sf_addevent_close(void* p)
/*< release memory >*/
{
    addevent *pp;
    pp = (addevent*) p;
    if(pp->vw != NULL) free(pp->vw);
    if(pp->wvlt != NULL) free(pp->wvlt);
    free(p);
}


