/* Generate simple data: spikes, boxes, planes, constants. 

Spike positioning is given in samples and starts with 1.
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
    int i, j, is, ip, dim, n[SF_MAX_DIM], ii[SF_MAX_DIM];
    int nsp, **k=NULL, **l=NULL, n1, n2, i1, i2, kk, ll;
    char key[7];
    const char *label, *unit;
    float f, *trace, *mag=NULL, **p=NULL, pp;
    sf_file in, spike;

    sf_init (argc,argv);

    if (!sf_stdin()) { /* no input file in stdin */
	in = NULL;
    } else {
	in = sf_input("in");
    }

    spike = sf_output("out");
    
    if (NULL == in) {
	sf_setformat(spike,"native_float");
    } else if (SF_FLOAT != sf_gettype(in)) {
	sf_error("Need float input");
    }
    
    /* dimensions */
    for (i=0; i < SF_MAX_DIM; i++) {
	snprintf(key,3,"n%d",i+1);
	if (!sf_getint(key,n+i) && 
	    (NULL == in || !sf_histint(in,key,n+i))) break;
	/*( n# size of #-th axis )*/  
	sf_putint(spike,key,n[i]);
    }

    if (0==i) sf_error("Need n1=");
    dim=i;
    
    /* basic parameters */
    for (i=0; i < dim; i++) {
	snprintf(key,3,"o%d",i+1);
	if (!sf_getfloat(key,&f) && 
	    (NULL == in || !sf_histfloat(in,key,&f))) f=0.;
	/*( o#=[0,0,...] origin on #-th axis )*/  
	sf_putfloat(spike,key,f);

	snprintf(key,3,"d%d",i+1);
	if (!sf_getfloat(key,&f) &&
	    (NULL == in || !sf_histfloat(in,key,&f))) f = (i==0)? 0.004: 0.1;
	/*( d#=[0.004,0.1,0.1,...] sampling on #-th axis )*/  
	sf_putfloat(spike,key,f);

	snprintf(key,7,"label%d",i+1);
	if (NULL == (label = sf_getstring(key)) &&
	    (NULL == in || NULL == (label = sf_histstring(in,key))))
	    label = (i==0)? "Time":"Distance";
	/*( label#=[Time,Distance,Distance,...] label on #-th axis )*/  
	if (*label != '\0' && (*label != ' ' || *(label+1) != '\0')) 	
	    sf_putstring(spike,key,label);

	snprintf(key,6,"unit%d",i+1);
	if (NULL == (unit = sf_getstring(key)) &&
	    (NULL == in || NULL == (unit = sf_histstring(in,key))))
	    unit = (i==0)? "s":"km";
        /*( unit#=[s,km,km,...] unit on #-th axis )*/  
	if (*unit != '\0' && (*unit != ' ' || *(unit+1) != '\0')) 	
	    sf_putstring(spike,key,unit);
    }
	
    if (NULL != (label = sf_getstring("title")))
	sf_putstring(spike,"title",label);
    /* title for plots */

    if (!sf_getint("nsp",&nsp)) nsp=1;
    /* Number of spikes */

    if (nsp >= 1) { 
	mag = sf_floatalloc (nsp);
	k = sf_intalloc2 (nsp,dim);
	l = sf_intalloc2 (nsp,dim);
	p = sf_floatalloc2 (nsp,dim);
    
	for (i=0; i < dim; i++) {
	    snprintf(key,3,"k%d",i+1);
	    if ( !sf_getints(key,k[i],nsp)) {
		/*( k#=[0,...] spike starting position )*/
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
	    snprintf(key,3,"l%d",i+1);
	    if (!sf_getints(key,l[i],nsp)) {
		/*( l#=[k1,k2,...] spike ending position )*/
		for (is=0; is < nsp; is++) {
		    l[i][is]=k[i][is];
		}
	    } else {
		for (is=0; is < nsp; is++) {
		    if (l[i][is] > n[i]) 
			sf_error("Invalid l%d[%d]=%d > n%d=%d",
				 i+1,is+1,l[i][is],i+1,n[i]);
		    l[i][is]--; /* C notation */
		}
	    }	
	    snprintf(key,3,"p%d",i+1);
	    if (!sf_getfloats(key,p[i],nsp)) {
		/*( p#=[0,...] spike inclination (in samples) )*/
		for (is=0; is < nsp; is++) {
		    p[i][is]=0.;
		}
	    }
	}
	
	if (!sf_getfloats("mag",mag,nsp)) {
	    /* spike magnitudes */
	    for (is=0; is < nsp; is++) {
		mag[is]=1.;
	    }
	}
    }

    n1 = n[0];
    n2 = sf_leftsize(spike,1);

    trace = sf_floatalloc (n[0]);

    for (i2=0; i2 < n2; i2++) { /* loop over traces */
	sf_line2cart(dim-1, n+1, i2, ii+1);
	/* zero trace */
	for (i1=0; i1 < n1; i1++) trace[i1]=0.;
	/* put spikes in it */
	for (is=0; is < nsp; is++) { /* loop over spikes */
	    pp = 0.;
	    for (i=1; i < dim; i++) {
		kk = k[i][is];
		ll = l[i][is];
		if ((kk < -1 && ll < -1) || 
		    (kk >= 0 && ll >= 0 && 
		     (kk > ii[i] || ll < ii[i]))) break;
		pp += p[i][is]*(ii[i]-k[i][is]-1);
	    }
	    if (i < dim) continue; /* skip this spike */

	    /* linear interpolation */
	    ip = floorf(pp);
	    pp = 1.-(pp-ip);

	    kk = k[0][is];
	    ll = l[0][is];
	    if (kk >= 0) { /* one segment per trace */
		kk = SF_MAX(kk+ip,0);
		ll = SF_MIN(ll+ip,n1-1);
	    } else {
		kk = SF_MAX(ip,0);
		ll = SF_MIN(n1-1+ip,n1-1);
	    }

	    for (j=kk; j <= ll; j++) {
		trace[j] += pp*mag[is];
		if (j+1 < n1) trace[j+1] += (1.-pp)*mag[is];
	    }
	}
	sf_floatwrite(trace,n1,spike);
    }

    exit (0);
}

/* 	$Id: spike.c 8430 2012-05-01 15:29:37Z sfomel $	 */
