#include "spline.h"
#include "banded.h"
#include "tridiagonal.h"

static const float s3 = 0.75, s4 = 2./3., s6 = 11./20., s8 = 151./315.;
static const float m4 = 0.625, m6 = 77./144., m8 = 151./320.;  
static const float flt3[] = {0.125}; 
static const float flt4[] = {1./6.};
static const float mom4[] = {0.1875};
static const float mom6[] = {2./9., 1./96.};
static const float flt6[] = {26./120., 1./120.};
static const float flt8[] = {1191./5040., 120./5040., 1./5040.};
static const float mom8[] = {2741./11520., 298./11520., 3./11520.};

bands spline_init (int nw, int nd)
{
    bands slv;
    int   na;
    
    na = (nw>0)?(nw-1)/2:(-nw-1)/2;
    slv = banded_init(nd,na);
    
    switch (nw) {
	case -8:
	    banded_const_define (slv, m8, mom8);
	    break;
	case -6:
	    banded_const_define (slv, m6, mom6);
	    break;
	case -4:
	    banded_const_define (slv, m4, mom4);
	    break;
	case 3:
	    banded_const_define (slv, s3, flt3);
	    break;
	case 4:
	    banded_const_define (slv, s4, flt4);
	    break;
	case 6:
	    banded_const_define (slv, s6, flt6);
	    break;
	case 8:
	    banded_const_define (slv, s8, flt8);
	    break;
	default:
	    break;
    }
    
    return slv;
}

tris spline4_init (int nd)
{
    tris slv;
    
    slv = tridiagonal_init(nd);
    tridiagonal_const_define (slv, s4, flt4[0]);
    
    return slv;
}

void spline2 (bands slv1, bands slv2, int n1, int n2, float** dat, float* tmp)
{
    int i1, i2;
    for (i2 = 0; i2 < n2; i2++) {
	banded_solve (slv1, dat[i2]);
    }
    for (i1 = 0; i1 < n1; i1++) {
	for (i2 = 0; i2 < n2; i2++) {
	    tmp[i2] = dat[i2][i1];
	}
	banded_solve (slv2, tmp);
	for (i2 = 0; i2 < n2; i2++) {
	    dat[i2][i1] = tmp[i2];
	}
    }
}

/* 	$Id: spline.c,v 1.2 2003/10/01 22:45:56 fomels Exp $	 */
