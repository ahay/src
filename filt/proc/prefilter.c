/* 
   File: prefilter.c
   -----------------
   Convert data to B-spline coefficients
   by fast B-spline transform
*/

#include <rsf.h>

#include "prefilter.h"
#include "recfilt.h"

static float flt2[1] = {0.0};
static float flt3[1] = {0.1715729};
static float flt4[1] = {0.2679492};
static float flt6[2] = {0.4736716, 0.01855619};
static float flt8[3] = {0.6669837, 0.07161941, 0.0006001561};

static float a2 = 1.; 
static float a3 = 1.372583;
static float a4 = 1.607695;
static float a6 = 2.226744;
static float a8 = 3.024829;

static float *tmp1, *tmp2 /* temporary storage */;
static float a0 /* normalization */;
static int nt, pad;

/*
  Function: prefilter_init
  ------------------------
  Initialize internal storage
  nw  - spline order (interpolator length)
  nt  - length for temporary storage
  pad - length for padding
*/
void prefilter_init (int nw, int nt_in, int pad_in)
{
    pad = pad_in;
    nt = nt_in + 2*pad;

    tmp1 = sf_floatalloc (nt);
    tmp2 = sf_floatalloc (nt);

    switch (nw) {    
	case 8: 
	    a0 = a8;
	    recfilt_init (nt, 3, flt8);
	    break;
	case 6:
	    a0 = a6;
	    recfilt_init (nt, 2, flt6);
	    break;
	case 4:
	    a0 = a4;
	    recfilt_init (nt, 1, flt4);
	    break;
	case 3:
	    a0 = a3;
	    recfilt_init (nt, 1, flt3);
	    break;
	case 2:
	    a0 = a2;
	    recfilt_init (nt, 1, flt2);
	    break;
	default:
	    sf_error ("%s: spline length %d not implemented",__FILE__,nw); 
	    break;
    }
}

/*
  Function: prefilter_apply
  -------------------------
  Convert 1-D data to spline coefficients
  nd      - data length
  dat[nd] - data/coefficients (input/output)
*/
void prefilter_apply (int nd, float* dat)
{
    int i;

    /* copy data and pad */
    for (i = 0; i < nt; i++) {
	if (i-pad >= 0 && i -pad < nd)
	    tmp1[i] = dat[i-pad];
	else
	    tmp1[i] = 0.0;
    }

    /* recursive filtering in two directions */
    recfilt_lop (false,false,nt,nt,tmp1,tmp2);
    recfilt_lop (true, false,nt,nt,tmp1,tmp2);

    /* extract the output */
    for (i = pad; i < nd+pad; i++) {
	dat[i-pad] = tmp1[i]*a0;
    }
}

/*
  Function: prefilter
  -------------------
  Convert N-D data to spline coefficients
  dim    - dimensionality
  n[dim] - data dimensions
  dat    - data/coefficients as 1-D array (input/output)
*/
void prefilter (int dim, int* n, float* dat)
{
    int i, i2, j, kk, k, n12;
    int m[SF_MAX_DIM], m1[SF_MAX_DIM], n1[SF_MAX_DIM];

    if (dim == 1) {
	prefilter_apply (n[0], dat);
	return;
    }
  
    n12 = 1;
    for (k=0; k < dim; k++) {
	n12 *= n[k];
    
    }

    for (j=0; j < dim; j++) {
	for (k=0, kk=0; k < dim; k++) {
	    if (k != j) n1[kk++]=n[k];
	}

	for (i2 =0; i2 < n12/n[j]; i2++) {
	    sf_line2cart(dim-1, n1, i2, m1);
	    for (k=0, kk=0; k < dim; k++) {
		if (k != j) m[k] = m1[kk++];
	    }

	    for (i = 0; i < nt; i++) {
		if (i-pad >= 0 && i -pad < n[j]) {
		    m[j] = i-pad;
		    tmp1[i] = dat[sf_cart2line(dim, n, m)];
		} else {
		    tmp1[i] = 0.0;
		}
	    }

	    recfilt_lop (false,false,nt,nt,tmp1,tmp2);
	    recfilt_lop (true, false,nt,nt,tmp1,tmp2);
    
	    for (i = 0; i < nt; i++) {
		if (i-pad >= 0 && i -pad < n[j]) {
		    m[j] = i-pad;
		    dat[sf_cart2line(dim, n, m)] = tmp1[i]*a0;
		}
	    }
	}
    }
}

/*
  Function: prefilter_close
  -------------------------
  Free temporary storage
*/
void prefilter_close( void)
{
    free (tmp1);
    free (tmp2);
    recfilt_close();
}
      
/* 	$Id: prefilter.c,v 1.2 2003/10/01 22:45:56 fomels Exp $	 */

