#include <float.h>

#include <rsf.h>

#include "shotfill.h"
#include "cbanded.h"

static int nh;
static float complex **a, *offd[2], **cc;
static float **c, *diag, s, h1, b0[3];
static void filt2matrix(const float complex *filt,
			float* d, float complex *dd);

void shotfill_init (int nh_in, float h0, float dh, float ds)
{
    nh = nh_in;
    offd[0] = sf_complexalloc(nh-1);
    offd[1] = sf_complexalloc(nh-2);
    diag = sf_floatalloc(nh);

    a = sf_complexalloc2(3,nh);
    c = sf_floatalloc2(3,nh);
    cc = sf_complexalloc2(3,nh);

    cbanded_init(nh, 2);
    
    s = -0.5*ds/dh;
    h1 = h0/dh;

    b0[0] = (2.-s)*(1.-s)/12.;
    b0[1] = (2.-s)*(2.+s)/6.;
    b0[2] = (2.+s)*(1.+s)/12.;
}

void shotfill_close (void)
{
    free(*a);
    free(a);
    free(offd[0]);
    free(offd[1]);
    free(diag);
    free(*c);
    free(c);
    free(*cc);
    free(cc);

    cbanded_close();
}

static void filt2matrix(const float complex *filt,
			float* d, float complex *dd)
{
    d[0] = crealf (filt[0] * conjf (filt[0]));
    d[1] = crealf (filt[1] * conjf (filt[1]));
    d[2] = crealf (filt[2] * conjf (filt[2]));
    dd[0] = filt[1] * conjf (filt[0]);
    dd[1] = filt[2] * conjf (filt[1]);
    dd[2] = filt[2] * conjf (filt[0]);
}

void shotfill_define (float w, const float complex* pef)
{
    float den, h, d[3];
    float complex *b, dd[3];
    int ih;

    filt2matrix(pef,d,dd);
  
    for (ih=1; ih < nh-1; ih++) {
	h = h1 + ih;

	den = 12.*h*h + (s-1.)*(s+1.)*w*w;
	
	b = a[ih];

	b[0] = b0[0] - 0.5*h*(I*(s-2.)/w+(2.*h+I*(s-1.)*w)/den);
	b[1] = b0[1] + h*(I*s/w+(2.*h+I*s*w)/den);
	b[2] = b0[2] - 0.5*h*(I*(s+2.)/w+(2.*h+I*(s+1.)*w)/den);

	filt2matrix(b,c[ih],cc[ih]);
    }

    diag[0] = c[1][0] + c[1][2] + d[2];
    diag[1] = 2.*c[1][1] + c[2][0] + c[2][2] + d[1] + d[2];
    for (ih=2; ih < nh-2; ih++) {
	diag[ih] = c[ih-1][0] + c[ih-1][2] + 2.*c[ih][1]
	    +      c[ih+1][0] + c[ih+1][2] + d[0] + d[1] + d[2];
    }
    diag[nh-2] = c[nh-3][0] + c[nh-3][2] + 2.*c[nh-2][1] + d[0] + d[1];
    diag[nh-1] = c[nh-2][0] + c[nh-2][2] + d[0];

    offd[0][0] = cc[1][0] + cc[1][1] + dd[1];
    offd[1][0] = 2.*cc[1][2] + dd[2];
    for (ih=1; ih < nh-2; ih++) {
	offd[0][ih] = 
	    cc[ih][0] + cc[ih+1][0] +  
	    cc[ih][1] + cc[ih+1][1] + dd[0] + dd[1];
	offd[1][ih] = 2.*cc[ih+1][2] + dd[2];
    }
    offd[0][nh-2] = cc[nh-2][0] + cc[nh-2][1] + dd[0];
    
    cbanded_define(diag,offd);
}

void shotfill_apply (const float complex *s1, const float complex *s2, 
		     float complex *s)
{
    float complex *b, sum1, sum2; 
    int ih;

    for (ih=0; ih < nh; ih++) {
	s[ih] = 0.;
    }

    for (ih=1; ih < nh-1; ih++) {
	b = a[ih];
	sum1 = s1[ih-1]*b[2] + s1[ih]*b[1] + s1[ih+1]*b[0];
	sum2 = 
	    s2[ih-1]*conjf(b[0]) + 
	    s2[ih  ]*conjf(b[1]) +
	    s2[ih+1]*conjf(b[2]);
	s[ih-1] += sum1*b[0] + sum2*conjf(b[2]);
	s[ih  ] += sum1*b[1] + sum2*conjf(b[1]);
	s[ih+1] += sum1*b[2] + sum2*conjf(b[0]);
    }

    cbanded_solve(s);
}

