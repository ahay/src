/* Add a analog event with dispersion and attenuation */

#include <rsf.h>
#include "traveltime.h"
#include "wavelet.h"


typedef struct tag_addevent
{
	float	df;		// time sample rate
	int 	of;		// time samples
	int		nf;		// frequency numbers (half)
	float complex 	*wvlt;	// frequency domain wavelet
	int 	type;	// event type
	float 	t0;		// event traveltime at zero offset (x=0)
	float 	a0;		// event amplitude at reference frequency at zero offset (x=0,f=f0)
	float 	*vw; 	// velocity over frequencies
	float 	qa;		// event amplitude attenuation
	float	f0;		// center frequency
	int 	a0ref;	// reference point for a0 
}addevent;


void* sf_addevent_init(int nfft, float of, float df, 	// time domain sampling parameters
	int wvtype,float *wvp,						// wavelet parameters
	int event, float t0, float v0, float a0,	// event initial parameters
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
		if(qv<=0) pp->vw[ifreq] = v0;		// dispersion velocity
		else {
			pp->vw[ifreq]=v0*pow(freq,1.0/(M_PI*qv));
		}
	}

	switch(wvtype)
	{
	case 0:
		for(ifreq=0;ifreq<pp->nf;ifreq++)
		{
			pp->wvlt[ifreq]=ifreq*df+of;
		}
		sf_wvlt_frck(pp->nf, pp->wvlt, wvp);
		break;
	default:
		break;
	}
	return ( pp);
}

void sf_addevent(void* p, float x, float complex *ftr)
/*< modeling for one trace >*/
{
	int ifreq;
	float txw,freq,fa,deltat;
	addevent *pp;
	float complex fc1;

	pp = (addevent*) p;
	for(ifreq=0; ifreq<pp->nf; ifreq++)
	{
		freq = ifreq*pp->df+pp->of;
		txw  = sf_traveltime(pp->type, pp->t0, x, pp->vw[ifreq]);
		if(pp->qa > 0)
		{
			if(pp->a0ref == 0) deltat = txw;
			else deltat = txw-pp->t0;
			fa = pp->a0*exp(-M_PI/fabs(freq*deltat));
		}else{
			fa = pp->a0;
		}
		fc1 = fa*cos(2.0*M_PI*freq*txw) - I*fa*sin(2.0*M_PI*freq*txw);
		ftr[ifreq] += fc1*pp->wvlt[ifreq];
	}
}

void sf_addevent_release(void* p)
/*< release memory >*/
{
	addevent *pp;
	pp = (addevent*) p;
	if(pp->vw != NULL) free(pp->vw);
	if(pp->wvlt != NULL) free(pp->wvlt);
	free(p);
}


