#include <rsf.h>

#include "fastmarch.h"
#include "pqueue.h"
#include "neighbors.h"

void fastmarch_init (int n3,int n2,int n1) 
{
    int maxband;
    
    maxband = 0;
    if (n1 > 1) maxband += 2*n2*n3;
    if (n2 > 1) maxband += 2*n1*n3;
    if (n3 > 1) maxband += 2*n1*n2;

    pqueue_init (maxband);
}

void fastmarch (float* time, float* v, int* in,
		int n3,int n2,int n1,
		float o3,float o2,float o1,
		float d3,float d2,float d1,
		float s3,float s2,float s1,
		int b3, int b2, int b1,
		int order)
{
    float xs[3], d[3], *p;
    int n[3], b[3], npoints, i;
    
    n[0] = n1; xs[0] = s1-o1; b[0] = b1; d[0] = d1;
    n[1] = n2; xs[1] = s2-o2; b[1] = b2; d[1] = d2;
    n[2] = n3; xs[2] = s3-o3; b[2] = b3; d[2] = d3;

    pqueue_start();
    neighbors_init (in, d, n, order, time);

    for (npoints =  nearsource (xs, b, d, v);
	 npoints > 0;
	 npoints -= neighbours(i)) {
	/* Pick smallest value in the NarrowBand
	   mark as good, decrease points_left */

	p = pqueue_extract();

	if (p == NULL) {
	    sf_warning("%s: heap exausted!",__FILE__);
	    break;
	}

	i = p - time;
	in[i] = FMM_IN;
    }
}

void fastmarch_close (void)
{
    pqueue_close();
}

/* 	$Id: fastmarch.c,v 1.2 2003/09/30 14:30:52 fomels Exp $	 */
