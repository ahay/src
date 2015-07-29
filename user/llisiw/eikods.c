/* Source perturbation eikonal interface */
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
/*^*/

#include "eikods.h"

static int ndim, nt, nn[3], ss[3];
static const float *t0, *v0;

static int fermat(const void *a, const void *b)
/* comparison for traveltime sorting from small to large */
{
    float ta, tb;

    ta = t0[*(int *)a];
    tb = t0[*(int *)b];

    if (ta >  tb) return 1;
    if (ta == tb) return 0;
    return -1;
}

void eikods_dt (int l, float* dl1, float* ds1, float* dl2, float* ds2, float* t, float *v, float* d);
void eikods_dt1(int l, float* dl1, float* ds1, int i0, float *d);
void eikods_dt2(int l, float* dl1, float* dl2, float* ds1, float* ds2, int i0, float *d);

void eikods_init (int n3,int n2,int n1) 
/*< Initialize data dimensions >*/
{
    int maxband;
    
    maxband = 0;
    if (n1 > 1) maxband += 2*n2*n3;
    if (n2 > 1) maxband += 2*n1*n3;
    if (n3 > 1) maxband += 2*n1*n2;

    sf_pqueue_init (10*maxband);

    nn[0] = n1; nn[1] = n2; nn[2] = n3;
    ss[0] = 1;  ss[1] = n1; ss[2] = n1*n2;

    nt = n1*n2*n3;

    if (nn[2] == 1) 
	ndim = 2;
    else
	ndim = 3;
}

void eikods (float* time                /* time */, 
	     float* v                   /* slowness squared */, 
	     int* in                    /* in/front/out flag */, 
	     bool* plane                /* if plane source */,
	     int   n3,  int n2,  int n1 /* dimensions */,
	     float o3,float o2,float o1 /* origin */,
	     float d3,float d2,float d1 /* sampling */,
	     float s3,float s2,float s1 /* source */,
	     int   b3,  int b2,  int b1 /* box around the source */,
	     int order                  /* accuracy order (1,2,3) */,
	     int l                      /* direction of source perturbation */, 
	     float* dl1, float* ds1     /* first-order derivatives */,
	     float* dl2, float* ds2     /* second-order derivatives */)
/*< Run fast marching eikonal solver >*/
{
    float xs[3], d[3], *p;
    int n[3], b[3], npoints, i;
    
    n[0] = n1; xs[0] = s1-o1; b[0] = b1; d[0] = d1;
    n[1] = n2; xs[1] = s2-o2; b[1] = b2; d[1] = d2;
    n[2] = n3; xs[2] = s3-o3; b[2] = b3; d[2] = d3;

    sf_pqueue_start();
    sf_neighbors_init (in, d, n, order, time);

    for (npoints =  sf_neighbors_nearsource (xs, b, d, v, plane);
	 npoints > 0;
	 npoints -= sf_neighbours(i)) {
	/* Pick smallest value in the NarrowBand
	   mark as good, decrease points_left */

	p = sf_pqueue_extract();

	if (p == NULL) {
	    sf_warning("%s: heap exausted!",__FILE__);
	    break;
	}
	
	i = p - time;

	in[i] = SF_IN;
    }

    /* derivatives */
    eikods_dt(l,dl1,ds1,dl2,ds2,time,v,d);
}

void eikods_close (void)
/*< Free allocated storage >*/
{
    sf_pqueue_close();
}

void eikods_dt(int l, float* dl1, float* ds1, float* dl2, float* ds2,
	       float* t, float *v, float* d)
/* derivatives */
{
    int i, it;
    int *sort;

    t0 = t;
    v0 = v;

    /* sort from small to large traveltime */
    sort = sf_intalloc(nt);
    for (it = 0; it < nt; it++) {
	sort[it] = it;
    }
    qsort(sort, nt, sizeof(int), fermat);
    
    for (it = 0; it < nt; it++) {
	i = sort[it];

	eikods_dt1(l,dl1,ds1,i,d);
	if (dl2 != NULL || ds2 != NULL) 
	    eikods_dt2(l,dl1,dl2,ds1,ds2,i,d);
    }
}



void eikods_dt1(int l      /* direction */,
		float* dl1 /* derivative (relative coordinate) */,
		float* ds1 /* derivative (source coordinate) */,
		int i0     /* grid location */, 
		float* d   /* sampling*/)
/* first-order derivative */
{
    int j, ii[3], a, b, i1=-1, i2=-1;
    double rhs=0., den=0., dw=0., dt=0.;

    sf_line2cart(3,nn,i0,ii);

    for (j=0; j < ndim; j++) {
	a = i0-ss[j];
	b = i0+ss[j];

	if ((ii[j] == 0) || 
	    (ii[j] != nn[j]-1 && 1==fermat(&a,&b))) {
	    i1 = b;
	} else {
	    i1 = a;
	}

	if (t0[i1] < t0[i0]) {	    
	    rhs += (t0[i0]-t0[i1])*dl1[i1]/(d[j]*d[j]);
	    den += (t0[i0]-t0[i1])/(d[j]*d[j]);

	    if (l == j) {
		if (i1 == b) {
		    i2 = i1+ss[j];

		    if (ii[j] < nn[j]-2 && 1==fermat(&i1,&i2))
			dt = (-t0[i2]+4.*t0[i1]-3.*t0[i0])/(2.*d[j]);
		    else
			dt = (t0[i1]-t0[i0])/d[j];
		} else {
		    i2 = i1-ss[j];

		    if (ii[j] > 1 && 1==fermat(&i1,&i2))
			dt = (3.*t0[i0]-4.*t0[i1]+t0[i2])/(2.*d[j]);
		    else
			dt = (t0[i0]-t0[i1])/d[j];
		}
	    }
	}	
    }

    if (den == 0.) {
	dl1[i0] = 0.;
	ds1[i0] = 0.;
	return;
    } else {
	if (0 < ii[l] && ii[l] < nn[l]-1)
	    dw = (v0[i0+ss[l]]-v0[i0-ss[l]])/(2.*d[l]);
	else if (ii[l] == 0)
	    dw = (-v0[i0+2*ss[l]]+4.*v0[i0+ss[l]]-3.*v0[i0])/(2.*d[l]);
	else
	    dw = (3.*v0[i0]-4.*v0[i0-ss[l]]+v0[i0-2*ss[l]])/(2.*d[l]);

	dl1[i0] = (0.5*dw+rhs)/den;
	ds1[i0] = dl1[i0]-dt;
    }
}

void eikods_dt2(int l                  /* direction */,
		float* dl1, float *dl2 /* derivative (relative coordinate) */,
		float* ds1, float *ds2 /* derivative (source coordinate) */,
		int i0                 /* grid location */, 
		float* d               /* sampling */)
/* second-order derivative */
{
    int j, ii[3], a, b, i1=-1, i2=-1;
    double rhs=0., den=0., dw, dldt=0., dtdt=0.;

    sf_line2cart(3,nn,i0,ii);

    for (j=0; j < ndim; j++) {
	a = i0-ss[j];
	b = i0+ss[j];

	if ((ii[j] == 0) || 
	    (ii[j] != nn[j]-1 && 1==fermat(&a,&b))) {
	    i1 = b;
	    if (ii[j] < nn[j]-2) i2 = b+ss[j];
	} else {
	    i1 = a;
	    if (ii[j] > 2)       i2 = a-ss[j];
	}

	if (t0[i1] < t0[i0]) {
	    rhs += (t0[i0]-t0[i1])*dl2[i1]/(d[j]*d[j])
		-(dl1[i0]-dl1[i1])*(dl1[i0]-dl1[i1])/(d[j]*d[j]);
	    den += (t0[i0]-t0[i1])/(d[j]*d[j]);

	    if (l == j) {
		if (i1 == b)
		    dldt = (dl1[i1]-dl1[i0])/d[j];
		else
		    dldt = (dl1[i0]-dl1[i1])/d[j];

		if (i2 != -1 && t0[i2] < t0[i1])
		    dtdt = (t0[i0]-2.*t0[i1]+t0[i2])/(d[j]*d[j]);
		else
		    dtdt = (t0[i0]-t0[i1])/(d[j]*d[j]);
	    }
	}
    }

    if (den == 0.) {
	dl2[i0] = 0.;
	return;
    } else {
	if (0 < ii[l] && ii[l] < nn[l]-1)
	    dw   = (v0[i0+ss[l]]-2.*v0[i0]+v0[i0-ss[l]])/(d[l]*d[l]);
	else if (ii[l] == 0)
	    dw   = (2.*v0[i0]-5.*v0[i0+ss[l]]+4.*v0[i0+2*ss[l]]-3.*v0[i0+3*ss[l]])/(d[l]*d[l]);
	else
	    dw   = (2.*v0[i0]-5.*v0[i0-ss[l]]+4.*v0[i0-2*ss[l]]-3.*v0[i0-3*ss[l]])/(d[l]*d[l]);
	
	dl2[i0] = (0.5*dw+rhs)/den;
	ds2[i0] = dl2[i0]-2.*dldt+dtdt;
    }
}

/* 	$Id: fastmarch.c 7107 2011-04-10 02:04:14Z ivlad $	 */
