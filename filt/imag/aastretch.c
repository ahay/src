#include "aastretch.h"

#include <rsf.h>

struct Aamap {
    int nt, nd, **x;
    float t0,dt, **w, *a, *tmp;
    bool *m;
};

aamap aastretch_init (int n1, float o1, float d1, int nd)
{
    aamap str;

    str = (aamap) sf_alloc(1,sizeof(*str));
    str->nt = n1; 
    str->t0 = o1; 
    str->dt = d1; 
    str->nd = nd;

    str->x = sf_intalloc2(3,nd);
    str->m = sf_boolalloc(nd);
    str->w = sf_floatalloc2(3,nd);
    str->a = sf_floatalloc(nd);

    str->tmp = sf_floatalloc(n1*3);

    return str;
}

void aastretch_define (aamap str, 
		       const float *coord, 
		       const float *dt, 
		       const float *amp)
{
    int id, ix[3], j;
    float rx[3];

    for (id = 0; id < str->nd; id++) {
	str->m[id] = false;

	rx[0] = coord[id] - dt[id] - str->dt;
	rx[1] = coord[id];          
	rx[2] = coord[id] + dt[id] + str->dt;
	for (j=0; j < 3; j++) {
	    rx[j] = (rx[j] - str->t0)/str->dt; 
	    ix[j] = rx[j]; 
	    if ((ix[j] < 0) || (ix[j] > str->nt - 2)) {
		str->m[id] = true;
		break;
	    }
	    str->w[id][j] = rx[j] - ix[j];	    
	}
	if (str->m[id]) continue;

	str->x[id][0] = ix[0];
	str->x[id][1] = ix[1] + str->nt;
	str->x[id][2] = ix[2] + 2*str->nt;

	str->a[id] = amp[id]*str->dt*str->dt/
	    (dt[id]*dt[id] + str->dt*str->dt); 
    }
}

void aastretch_apply (aamap str, const float *ord, float *modl)
{
    int id, i1, i2, j, it, n;
    float w1, w2, t, w, *tmp;

    tmp = str->tmp;

    for (it=0; it < 3*str->nt; it++) {
	tmp[it]=0.;
    }

    n = str->nt;
    for (id = 0; id < str->nd; id++) {
	if (str->m[id]) continue;
	
	w = str->a[id] * ord[id];
        for (j=0; j < 3; j++) {
	    i1 = str->x[id][j]; 
	    i2 = i1 + 1;

	    w2 = str->w[id][j]; 
	    w1 = 1. - w2;

	    tmp[i1] += w1 * w;
	    tmp[i2] += w2 * w;
	}
    }

    for (it=0; it < n; it++) {
	modl[it] = 2.* tmp[it+n] - tmp[it] - tmp[it+2*n];
    }

    t = 0.;
    for (it = 0; it < n; it++) {
	t += modl[it];
	tmp[it] = t;
    }
    t = 0.;
    for (it = n-1; it >=0; it--) {
	t += tmp[it];
	modl[it] = t;
    }
}

void aastretch_close (aamap str)
{
    free (str->x[0]);
    free (str->x);
    free (str->m);
    free (str->w[0]);
    free (str->w);
    free (str->a);
    free (str->tmp);
    free (str);
}



