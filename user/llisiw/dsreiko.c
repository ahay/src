/* Double square-root eikonal solver interface */
/*
  Copyright (C) 2011 University of Texas at Austin
  
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

#include <rsf.h>

#include "dsreiko.h"

struct Upd {
    double stencil, value, delta;
    int label;
};

static const double tol = 1.e-6;

static float *o, *v, *d;
static int *n, *in, s[3];
static float **x, **xn, **x1;
static int *offsets;
static float *t;

void pqueue_insert(float* v1);
float* pqueue_extract(void);
void pqueue_update(int index);
int neighbors_default();
int neighbours(float* time, int i);
int update(float value, float* time, int i);
float qsolve(float* time, int i);
bool updaten(int i, int m, float* res, struct Upd *vv[], struct Upd* xj[], double vr, double vs);
double newton(double a,double b,double c,double d,double e, double guess);
double ferrari(double a,double b,double c,double d,double e);

void dsreiko_init(int *n_in   /* length */,
		  float *o_in /* origin */,
		  float *d_in /* increment */)
/*< initialize >*/
{
    int maxband;

    n = n_in;
    o = o_in;
    d = d_in;

    s[0] = 1; s[1] = n[0]; s[2] = n[0]*n[1];

    /* maximum length of heap */
    maxband = 0;

    if (n[0] > 1) maxband += 2*n[1]*n[2];
    if (n[1] > 1) maxband += 2*n[0]*n[2];
    if (n[2] > 1) maxband += 2*n[0]*n[1];

    x = (float **) sf_alloc ((10*maxband+1),sizeof (float *));
    in = sf_intalloc(n[0]*n[1]*n[2]);

    offsets = (int *) sf_alloc (n[0]*n[1]*n[2],sizeof (int));
}

void dsreiko_fastmarch(float *time /* time */,
		       float *v_in /* slowness squared */)
/*< fast-marching solver >*/
{
    float *p;
    int npoints, i;

    t = time;
    v = v_in;

    xn = x;
    x1 = x+1;

    /* initialize from zero-offset plane */
    for (npoints =  neighbors_default();
	 npoints > 0;
	 npoints -= neighbours(t,i)) {
	/* Pick smallest value in the NarrowBand
	   mark as good, decrease points_left */

	p = pqueue_extract();

	if (p == NULL) {
	    sf_warning("%s: heap exausted!",__FILE__);
	    break;
	}
	
	i = p-t;

	in[i] = SF_IN;
    }
}

void dsreiko_mirror(float *time /*time*/)
/*< source-receiver reciprocity >*/
{
    int i, j, k;

    for (k=1; k < n[1]; k++) {
	for (j=0; j < k; j++) {
	    for (i=0; i < n[0]; i++) {
		time[k*s[2]+j*s[1]+i] = time[j*s[2]+k*s[1]+i];
	    }
	}
    }
}

void dsreiko_close()
/*< finalize >*/
{
    free(x);
    free(in);
}

void pqueue_insert(float* v1)
/* insert an element (smallest first) */
{
    int newOffset, tIndex;
    float **xi, **xq;
    unsigned int q;
    
    xi = ++xn;

    q = (unsigned int) (xn-x);
    for (q >>= 1; q > 0; q >>= 1) {
	xq = x + q;
	if (*v1 > **xq) break;
	*xi = *xq;

	/* moved down, update its offset */
	newOffset = xi-x;
	tIndex = (*xi)-t;
	offsets[tIndex] = newOffset;

	xi = xq;
    }
    *xi = v1; 
    
    /* far enough up, record the offset */
    newOffset = xi-x;
    tIndex = v1-t;
    offsets[tIndex] = newOffset;
}

float* pqueue_extract(void)
/* extract the smallest element */
{
    int newOffset, tIndex;
    unsigned int c;
    int nn;
    float *vv, *formerlyLast;
    float **xi, **xc;
    
    /* check if queue is empty */
    nn = (int) (xn-x);
    if (nn == 0) return NULL;

    vv = *x1;
    tIndex = vv-t;
    offsets[tIndex] = -1;

    *(xi = x1) = formerlyLast = *(xn--);
    nn--;
    for (c = 2; c <= (unsigned int) nn; c <<= 1) {
	xc = x + c;
	if (c < (unsigned int) nn && **xc > **(xc+1)) {
	    c++; xc++;
	}
	if (*formerlyLast <= **xc) break;
	*xi = *xc; 

	/* moved up, update its offset */
	newOffset = xi-x;
	tIndex = (*xi)-t;
	offsets[tIndex] = newOffset;

	xi = xc;
    }
    *xi = formerlyLast;

    /* far enough down, record the offset */
    newOffset = xi-x;
    tIndex = *xi-t;
    offsets[tIndex] = newOffset;

    return vv;
}

void pqueue_update(int index)
/* restore the heap */
{
    int newOffset, tIndex;
    unsigned int c;
    float **xc, **xi;
    
    c = offsets[index];
    xi = x+c;

    for (c >>= 1; c > 0; c >>= 1) {
	xc = x + c;
	if (t[index] > **xc) break;
	*xi = *xc;
	
	/* moved down, update its offset */
	newOffset = xi-x;
	tIndex = (*xi)-t;
	offsets[tIndex] = newOffset;
	
	xi = xc; 
    }
    *xi = t+index; 

    /* far enough up, record the offset */
    newOffset = xi-x;
    tIndex = *xi-t;
    offsets[tIndex] = newOffset;
}

int neighbors_default()
/* initialize source */
{
    int i, j, k, nxy;
    double vr, vs;

    /* total number of points */
    nxy = n[0]*n[1]*n[2];

    /* set all points to be out */
    for (i=0; i < nxy; i++) {
	in[i] = SF_OUT;
	t[i] = SF_HUGE;
	offsets[i] = -1;
    }
    
    /* zero-offset and h = s-r */
    for (j=0; j < n[1]; j++) {
	for (k=0; k <= j; k++) {
	    for (i=0; i < n[0]; i++) {
		in[j*s[2]+k*s[1]+i] = SF_IN;
		t[j*s[2]+k*s[1]+i] = 0.;
	    }
	}
    }

    /* h = r-s */
    for (j=0; j < n[1]-1; j++) {
	for (i=0; i < n[0]; i++) {
	    vr = (double)v[(j+1)*s[1]+i];
	    vs = (double)v[j*s[1]+i];

	    in[j*s[2]+(j+1)*s[1]+i] = SF_FRONT;
	    t[j*s[2]+(j+1)*s[1]+i] = (sqrt(vr)*d[1]+sqrt(vs)*d[2])/2.;

	    pqueue_insert(t+j*s[2]+(j+1)*s[1]+i);
	}
    }
    
    return (nxy-n[0]*(n[1]*n[1]+3*n[1]-2)/2);
}

int neighbours(float* time, int i)
/* update neighbors of gridpoint i, return number of updated points */
{
    int j, k, ix, np;
    
    np = 0;
    for (j=0; j < 3; j++) {
	ix = (i/s[j])%n[j];

	/* try both directions */
	if (ix+1 <= n[j]-1) {
	    k = i+s[j]; 
	    if (in[k] != SF_IN) np += update(qsolve(time,k),time,k);
	}
	if (ix-1 >= 0 ) {
	    k = i-s[j];
	    if (in[k] != SF_IN) np += update(qsolve(time,k),time,k);
	}
    }
    return np;
}

int update(float value, float* time, int i)
/* update gridpoint i with new value and modify wave front */
{
    /* only update when smaller than current value */
    if (value < time[i]) {
	time[i] = value;
	if (in[i] == SF_OUT) { 
	    in[i] = SF_FRONT;      
	    pqueue_insert(time+i);
	    return 1;
	} else {
	    pqueue_update(i);
	}
    }
    
    return 0;
}

float qsolve(float* time, int i)
/* find new traveltime at gridpoint i */
{
    int j, k, ix[3];
    float a, b, res;
    double vr, vs;
    struct Upd *vv[3], xx[3], *xj[3];

    for (j=0; j<3; j++) 
	ix[j] = (i/s[j])%n[j];

    vr = v[ix[1]*s[1]+ix[0]];
    vs = v[ix[2]*s[1]+ix[0]];

    for (j=0; j<3; j++) {
	if (ix[j] > 0) { 
	    k = i-s[j];
	    a = time[k];
	} else {
	    a = SF_HUGE;
	}
	
	if (ix[j] < n[j]-1) {
	    k = i+s[j];
	    b = time[k];
	} else {
	    b = SF_HUGE;
	}
	
	xj[j] = xx+j;
	xj[j]->label = j;
	xj[j]->delta = 1./(d[j]*d[j]);

	if (a < b) {
	    xj[j]->value = a;
	} else {
	    xj[j]->value = b;
	}
    }

    xj[0]->stencil = vr+vs;
    xj[1]->stencil = vr-vs;
    xj[2]->stencil = vs-vr;

    /* sort from smallest to largest */
    if (xx[0].value <= xx[1].value) {
	if (xx[1].value <= xx[2].value) {
	    vv[0] = xx; vv[1] = xx+1; vv[2] = xx+2;
	} else if (xx[2].value <= xx[0].value) {
	    vv[0] = xx+2; vv[1] = xx; vv[2] = xx+1;
	} else {
	    vv[0] = xx; vv[1] = xx+2; vv[2] = xx+1;
	}
    } else {
	if (xx[0].value <= xx[2].value) {
	    vv[0] = xx+1; vv[1] = xx; vv[2] = xx+2;
	} else if (xx[2].value <= xx[1].value) {
	    vv[0] = xx+2; vv[1] = xx+1; vv[2] = xx;
	} else {
	    vv[0] = xx+1; vv[1] = xx+2; vv[2] = xx;
	}
    }

    if(vv[2]->value < SF_HUGE) {   /* update from all three directions */
	if (updaten(i,3,&res,vv,xj,vr,vs) || 
	    updaten(i,2,&res,vv,xj,vr,vs) || 
	    updaten(i,1,&res,vv,xj,vr,vs)) return res;
    } else if(vv[1]->value < SF_HUGE) { /* update from two directions */
	if (updaten(i,2,&res,vv,xj,vr,vs) || 
	    updaten(i,1,&res,vv,xj,vr,vs)) return res;
    } else if(vv[0]->value < SF_HUGE) { /* update from only one direction */
	if (updaten(i,1,&res,vv,xj,vr,vs)) return res;
    }
	
    return SF_HUGE;
}

bool updaten(int i, int m, float* res, struct Upd *vv[], struct Upd *xj[], double vr, double vs)
/* calculate new traveltime */
{
    double a, b, c, d, e, t, discr;
    int j;    

    /* a*t^4 + b*t^3 + c*t^2 + d*t + e = 0. */
    a = b = c = d = e = 0.;

    /* one-sided */
    if (m == 1 && vv[0]->label == 0) {
	t = (sqrt(vs)+sqrt(vr))/sqrt(xj[0]->delta)+xj[0]->value;
	
	*res = t;
	return true;
    }
    
    if (m == 1 && vv[0]->label == 1) {
	t = sqrt(vr/xj[1]->delta)+xj[1]->value;

	*res = t;
	return true;
    }

    if (m == 1 && vv[0]->label == 2) {
	t = sqrt(vs/xj[2]->delta)+xj[2]->value;

	*res = t;
	return true;
    }

    /* two-sided */
    if (m == 2 && vv[2]->label == 0) {
	t = ((sqrt(vr/xj[1]->delta)+xj[1]->value)+(sqrt(vs/xj[2]->delta)+xj[2]->value))/2.;
	
	*res = t;
	return true;
    }

    if (m == 2 && vv[2]->label == 1) {
	a = xj[0]->delta+xj[2]->delta;
	b = xj[0]->value*xj[0]->delta+xj[2]->value*xj[2]->delta+sqrt(vr*xj[0]->delta);
	c = pow(xj[0]->value,2.)*xj[0]->delta+pow(xj[2]->value,2.)*xj[2]->delta+2.*sqrt(vr*xj[0]->delta)*xj[0]->value+(vr-vs);
	b /= a;

	discr=b*b-c/a;

	if (discr < 0.) return false;

	t = b + sqrt(discr);
	if (t <= vv[m-1]->value) return false;

	*res = t;
	return true;
    }

    if (m == 2 && vv[2]->label == 2) {
	a = xj[0]->delta+xj[1]->delta;
	b = xj[0]->value*xj[0]->delta+xj[1]->value*xj[1]->delta+sqrt(vs*xj[0]->delta);
	c = pow(xj[0]->value,2.)*xj[0]->delta+pow(xj[1]->value,2.)*xj[1]->delta+2.*sqrt(vs*xj[0]->delta)*xj[0]->value+(vs-vr);
	b /= a;

	discr=b*b-c/a;

	if (discr < 0.) return false;

	t = b + sqrt(discr);
	if (t <= vv[m-1]->value) return false;

	*res = t;
	return true;
    }

    /* three-sided */
    for (j=0; j<m; j++) {
	a += pow(vv[j]->delta,2.);
	b += -4.*vv[j]->value*pow(vv[j]->delta,2.);
	c += 6.*pow(vv[j]->value,2.)*pow(vv[j]->delta,2.)-2.*vv[j]->stencil*vv[j]->delta;
	d += -4.*pow(vv[j]->value,3.)*pow(vv[j]->delta,2.)+4.*vv[j]->stencil*vv[j]->value*vv[j]->delta;
	e += pow(vv[j]->value,4.)*pow(vv[j]->delta,2.)-2.*vv[j]->stencil*pow(vv[j]->value,2.)*vv[j]->delta;
    }    
    e += pow(vs,2.)+pow(vr,2.)-2.*vs*vr;


    a += -2.*xj[1]->delta*xj[2]->delta
	 +2.*xj[0]->delta*xj[2]->delta
	 +2.*xj[1]->delta*xj[0]->delta;
    
    b +=  4.*xj[1]->delta*xj[2]->delta*(xj[1]->value+xj[2]->value)
	 -4.*xj[0]->delta*xj[2]->delta*(xj[0]->value+xj[2]->value)
	 -4.*xj[1]->delta*xj[0]->delta*(xj[1]->value+xj[0]->value);
    
    c += -2.*xj[1]->delta*xj[2]->delta*(pow(xj[1]->value,2.)+4.*xj[1]->value*xj[2]->value+pow(xj[2]->value,2.))
	 +2.*xj[0]->delta*xj[2]->delta*(pow(xj[0]->value,2.)+4.*xj[0]->value*xj[2]->value+pow(xj[2]->value,2.))
	 +2.*xj[1]->delta*xj[0]->delta*(pow(xj[1]->value,2.)+4.*xj[1]->value*xj[0]->value+pow(xj[0]->value,2.));
    
    d +=  4.*xj[1]->delta*xj[2]->delta*(xj[1]->value*pow(xj[2]->value,2.)+xj[2]->value*pow(xj[1]->value,2.))
	 -4.*xj[0]->delta*xj[2]->delta*(xj[0]->value*pow(xj[2]->value,2.)+xj[2]->value*pow(xj[0]->value,2.))
	 -4.*xj[1]->delta*xj[0]->delta*(xj[1]->value*pow(xj[0]->value,2.)+xj[0]->value*pow(xj[1]->value,2.));
    
    e += -2.*xj[1]->delta*xj[2]->delta*pow(xj[1]->value,2.)*pow(xj[2]->value,2.)
	 +2.*xj[0]->delta*xj[2]->delta*pow(xj[0]->value,2.)*pow(xj[2]->value,2.)
	 +2.*xj[1]->delta*xj[0]->delta*pow(xj[1]->value,2.)*pow(xj[0]->value,2.);
    
    /*
    t = newton(a,b,c,d,e,vv[m-1]->value);
    */
    t = ferrari(a,b,c,d,e);

    if (t <= vv[m-1]->value) return false;

    *res = t;
    return true;
}

double newton(double a,double b,double c,double d,double e /* coefficients */,
	      double guess /* initial guess */)
/* quartic solve (Newton's method) */
{
    double val, der, root;

    root = guess;

    /* evaluate functional */
    val = a*pow(root,4.)+b*pow(root,3.)+c*pow(root,2.)+d*root+e;
    while (fabs(val) > tol) {

	/* evaluate derivative */
	der = 4.*a*pow(root,3.)+3.*b*pow(root,2.)+2.*c*root+d*root+d;

	if (fabs(der) <= tol) {
	    sf_warning("FAILURE: Newton's method meets local minimum.");
	    return -1.;
	}

	root -= val/der;
	val = a*pow(root,4.)+b*pow(root,3.)+c*pow(root,2.)+d*root+e;
    }

    return root;
}

double ferrari(double a,double b,double c,double d,double e /* coefficients */)
/* quartic solve (Ferrari's method) */
{
    double alpha, beta, gama, P, Q, y, W;
    double delta, root, temp;
    double complex R, U;

    alpha = -3./8.*pow(b,2.)/pow(a,2.)+c/a;
    beta  = 1./8.*pow(b,3.)/pow(a,3.)-1./2.*b*c/pow(a,2.)+d/a;
    gama  = -3./256.*pow(b,4.)/pow(a,4.)+1./16.*c*pow(b,2.)/pow(a,3.)-1./4.*b*d/pow(a,2.)+e/a;

    P = -1./12.*pow(alpha,2.)-gama;
    Q = -1./108.*pow(alpha,3.)+1./3.*alpha*gama-1./8.*pow(beta,2.);

    delta = 1./4.*pow(Q,2.)+1./27.*pow(P,3.);

    /* R could be complex, any complex root will work */
    if (delta >= 0.)
	R = -1./2.*Q+sqrt(delta)+I*0.;
    else
	R = -1./2.*Q+I*sqrt(-delta);
    
    U = cpow(R,1./3.);

    /* rotate U by exp(i*2*pi/3) until W is real*/
    if (2.*creal(U) <= alpha/3.)
	U *= -0.5+I*(sqrt(3.)/2.);
    if (2.*creal(U) <= alpha/3.)
	U *= -0.5+I*(sqrt(3.)/2.);

    /* y must be real since a,b,c,d,e are all real */
    if (cabs(U) <= 1.e-15)
	if (Q >= 0.)
	    y = -5./6.*alpha-pow(Q,1./3.);
	else
	    y = -5./6.*alpha+pow(-Q,1./3.);
    else
	y = -5./6.*alpha+2.*creal(U);

    W = sqrt(alpha+2.*y);

    root = -1.;

    /* return the largest real root */
    delta = -3.*alpha-2.*y-2.*beta/W;    
    if (delta >= 0.) {
	temp = -1./4.*b/a+1./2.*(W+sqrt(delta));
	if (temp > root) root = temp;
	
	temp = -1./4.*b/a+1./2.*(W-sqrt(delta));
	if (temp > root) root = temp;
    }

    delta = -3.*alpha-2.*y+2.*beta/W;    
    if (delta >= 0.) {
	temp = -1./4.*b/a+1./2.*(-W+sqrt(delta));
	if (temp > root) root = temp;
	temp = -1./4.*b/a+1./2.*(-W-sqrt(delta));
	if (temp > root) root = temp;
    }

    return root;
}
