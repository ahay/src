/* Generate simple data: spikes, planes, constants.

Takes: [n1= n2= ... d1= d2= ... o1= o2= ... label1= label2= ... k1= k2= ... mag=1,1,...] 

k1,k2,... specify the spike position (indexing starts with 1).

Inserts label1="Time (s)" and label2=label3=...="Distance (km)"
Inserts d1=0.004 and d2=d3=...=0.1

*/
/*
Copyright (C) 2004 University of Texas at Austin

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

int main(int argc, char* argv[])
{ 
    int i, is, dim, n[SF_MAX_DIM], ii[SF_MAX_DIM];
    int nsp, **k, n1, n2, i1, i2, kk;
    char key[7], *label;
    float f, *trace, *mag;
    sf_file spike;

    sf_init (argc,argv);
    spike = sf_output("out");
    sf_setformat(spike,"native_float");

    /* dimensions */
    for (i=0; i < SF_MAX_DIM; i++) {
	snprintf(key,3,"n%d",i+1);
	if (!sf_getint(key,n+i)) break;
	sf_putint(spike,key,n[i]);
    }

    if (0==i) sf_error("Need n1=");
    dim=i;
    
    /* basic parameters */
    for (i=0; i < dim; i++) {
	snprintf(key,3,"o%d",i+1);
	if (!sf_getfloat(key,&f)) f=0.;
	sf_putfloat(spike,key,f);

	snprintf(key,3,"d%d",i+1);
	if (!sf_getfloat(key,&f)) f = (i==0)? 0.004: 0.1;
	sf_putfloat(spike,key,f);

	snprintf(key,7,"label%d",i+1);
	if (NULL == (label = sf_getstring(key)))
	    label = (i==0)? "Time (s)":"Distance (km)";
	sf_putstring(spike,key,label);
    }
	
    if (!sf_getint("nsp",&nsp)) nsp=1;
    /* Number of spikes */
    mag = sf_floatalloc (nsp);
    k = sf_intalloc2 (nsp,dim);

    for (i=0; i < dim; i++) {
	snprintf(key,3,"k%d",i+1);
	if (!sf_getints(key,k[i],nsp)) {
	    for (is=0; is < nsp; is++) {
		k[i][is]=-1;
	    }
	} else {
	    for (is=0; is < nsp; is++) {
		if (k[i][is] > n[i]) 
		    sf_error("Invalid k%d[%d]=%d > n%d=%d",
			     i+1,is+1,k[i][is],i+1,n[i]);
		k[i][is]--; /* C notation */
	    }
	}
    }

    if (!sf_getfloats("mag",mag,nsp)) {
	for (is=0; is < nsp; is++) {
	    mag[is]=1.;
	}
    }

    n1 = n[0];
    n2 = sf_leftsize(spike,1);

    trace = sf_floatalloc (n[0]);

    for (i2=0; i2 < n2; i2++) {
	sf_line2cart(dim-1, n+1, i2, ii+1);
	/* zero trace */
	for (i1=0; i1 < n1; i1++) trace[i1]=0.;
	/* put spikes in it */
	for (is=0; is < nsp; is++) {
	    for (i=1; i < dim; i++) {
		kk = k[i][is];
		if ((kk < -1) || (kk >= 0 && kk != ii[i])) break;	    
	    }
	    if (i < dim) continue;
	    kk = k[0][is];
	    if (kk >= 0) { /* one spike per trace */
		trace[kk] += mag[is];
	    } else {
		for (i1=0; i1 < n1; i1++) {
		    trace[i1] += mag[is];
		}
	    }
	}
	sf_floatwrite(trace,n1,spike);
    }

    exit (0);
}

/* 	$Id: spike.c,v 1.6 2004/06/23 18:30:00 fomels Exp $	 */
