#include <float.h>

#include <rsf.h>

#include "shotfill.h"
#include "cbanded.h"

static int nh;
static float complex **a, *offd[2], **cc;
static float **c, *diag, s, h1, b0[3];

void shotfill_init (int nh_in, float h0, float dh, float ds)
{
    nh = nh_in;
    a = sf_complexalloc2(nh,3);
    offd[0] = sf_complexalloc(nh-1);
    offd[1] = sf_complexalloc(nh-2);
    diag = sf_floatalloc(nh);
    c = sf_floatalloc2(nh,3);
    cc = sf_complexalloc2(nh,3);

    cbanded_init(nh, 2);
    
    s = ds/dh;
    h1 = h0/dh;

    b0[0] = (s-2.)*(s-1.)/12.;
    b0[1] = (s-2.)*(s+2.)/6.;
    b0[2] = (s+2.)*(s+1.)/12.;
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

void shotfill_define (float w, const float complex *pef)
{
    float d[3], den, h;
    float complex dd[3], b[3];
    int ih;

    d[0] = crealf (pef[0] * conjf (pef[0]));
    d[1] = crealf (pef[1] * conjf (pef[1]));
    d[2] = crealf (pef[2] * conjf (pef[2]));
    dd[0] = pef[1] * conjf (pef[0]);
    dd[1] = pef[2] * conjf (pef[1]);
    dd[2] = pef[2] * conjf (pef[0]);
  
    for (ih=1; ih < nh-1; ih++) {
	h = h1 + ih;

	den = 12.*h*h + (s-1.)*(s+1.)*w*w;
	
	b[0] = b0[0] - 0.5*h*(I*(s-2.)/w+(2.*h+I*(s-1.)*w)/den);
	b[1] = b0[1] + h*(I*s/w+(2.*h+I*s*w)/den);
	b[2] = b0[2] - 0.5*h*(I*(s+2.)/w+(2.*h+I*(s+1.)*w)/den);

	a[0][ih] = b[0];
	a[1][ih] = b[1];
	a[2][ih] = b[2];
	c[0][ih] = crealf(b[0] * conjf(b[0]));
	c[1][ih] = crealf(b[1] * conjf(b[1]));
	c[2][ih] = crealf(b[2] * conjf(b[2]));
	cc[0][ih] = b[0] * conjf(b[1]);
	cc[1][ih] = b[1] * conjf(b[2]);
	cc[2][ih] = b[0] * conjf(b[2]);
    }

    diag[0] = c[0][1] + c[2][1] + d[2];
    diag[1] = 2.*c[1][1] + c[0][2] + c[2][2] + d[1] + d[2];
    for (ih=2; ih < nh-2; ih++) {
	diag[ih] = c[0][ih-1] + c[2][ih-1] + 2.*c[1][ih]
	    +      c[0][ih+1] + c[2][ih+1] + d[0] + d[1] + d[2];
    }
    diag[nh-2] = c[0][nh-3] + c[2][nh-3] + 2.*c[1][nh-2] + d[0] + d[1];
    diag[nh-1] = c[0][nh-2] + c[2][nh-2] + d[0];

    offd[0][0] = cc[0][1] + cc[1][1] + dd[1];
    offd[1][0] = 2.*cc[2][1] + dd[2];
    for (ih=1; ih < nh-2; ih++) {
	offd[0][ih] = cc[0][ih] + cc[0][ih+1] + 
	    cc[1][ih] + cc[1][ih+1] + dd[0] + dd[1];
	offd[1][ih] = 2.*cc[2][ih+1] + dd[2];
    }
    offd[0][nh-2] = cc[0][nh-2] + cc[1][nh-2] + dd[0];
    
    cbanded_define(diag,offd);
}

void shotfill_apply (const float complex *s1, const float complex *s2, 
		     float complex *s)
{
    float complex b[3], sum1, sum2; 
    int ih;

    for (ih=0; ih < nh; ih++) {
	s[ih] = 0.;
    }

    for (ih=1; ih < nh-1; ih++) {
	b[0] = a[0][ih];
	b[1] = a[1][ih];
	b[2] = a[2][ih];
	sum1 = s2[ih-1]*b[0] + s2[ih]*b[1] + s2[ih+1]*b[2];
	sum2 = 
	    s1[ih-1]*conjf(b[2]) + 
	    s1[ih  ]*conjf(b[1]) +
	    s1[ih+1]*conjf(b[0]);
	s[ih-1] += sum1*b[2] + sum2*conjf(b[0]);
	s[ih  ] += sum1*b[1] + sum2*conjf(b[1]);
	s[ih+1] += sum1*b[0] + sum2*conjf(b[2]);
    }

    cbanded_solve(s);
}

