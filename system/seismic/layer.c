/* Layered medium */
/*
  Copyright (C) 2006 University of Texas at Austin
   
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

#include "layer.h"
#include "kirmod.h"

static int n;                   /* number of layers */
static float *y;                /* ray points */
static sf_eno *lay, *dip, *cur; /* interfaces, dips, and curvatures */
static char *type;              /* type of velocity distribution */
static ktable table;            /* traveltime attributes table */
static float *v0, **g, **x0;    /* velocity, gradient, reference */     

void layer_init(int n_in        /* number of layers */,
		char *types     /* [n] velocity model types */, 
		float *v0_in    /* [n] reference velocity  */,
		float **g_in    /* [n] velocity gradient */,
		float **x0_in   /* [n] gradient reference point */,
		int order       /* interpolation order */,
		int nx          /* horizontal sampling */,
		float **lays, 
		float **dips, 
		float **curs    /* [n+1][nx] interfaces */) 
/*< initialize >*/
{
    int k;

    n = n_in;

    y = sf_floatalloc(n-1);
    type = types;
    lay = (sf_eno*) sf_alloc(n,sizeof(sf_eno));
    dip = (sf_eno*) sf_alloc(n,sizeof(sf_eno));
    cur = (sf_eno*) sf_alloc(n,sizeof(sf_eno));
    table = (ktable) sf_alloc(1,sizeof(*table));

    for (k=0; k < n+1; k++) {
	lay[k] = sf_eno_init (order,nx); sf_eno_set (lay[k],lays[k]);
	dip[k] = sf_eno_init (order,nx); sf_eno_set (dip[k],dips[k]);
	cur[k] = sf_eno_init (order,nx); sf_eno_set (cur[k],curs[k]);
    }    

    v0 = v0_in;
    g = g_in;
    x0 = x0_in;
}

void layer_close(void)
/*< free allocated memory >*/
{
    int k;

    free(y);
    free(table);

    for (k=0; k < n+1; k++) {
	sf_eno_close(lay[k]);
	sf_eno_close(dip[k]);
	sf_eno_close(cur[k]);
    }

    free(lay);
    free(dip);
    free(cur);
}

void layer_set(float* start /* [2] start ray point */, 
	       float* end   /* [2] end ray point */)
/*< initialize intersections >*/
{
    int i;
    float dx;

    dx = (end[0]-start[0])/n;

    for (i=1; i < n; i++) {
	y[i-1] = start[0] + i*dx;
    }
}

float layer_time(float* start /* [2] start ray point */, 
		 float* end   /* [2] end ray point */)
{
    int ix, i;
    float t, x, x2, z, z2, v, v2;

    /* first intersection */
    x = y[0];
    ix = floorf(x);
    sf_eno_apply (lay[0],ix,x-ix,&z,NULL,FUNC);

    /* velocities */
    v  = v0[0]+g[0][0]*(start[0]-x0[0][0])+g[1][0]*(start[1]-x0[1][0]);
    v2 = v0[0]+g[0][0]*(x       -x0[0][0])+g[1][0]*(z       -x0[1][0]);

    /* first layer */
    kirmod_table(type[0], true,
		 z-start[1], 
		 x-start[0],0.,
		 hypotf(g[0][0],g[1][0]),
		 g[0][0],0.,g[1][0], /* gradient */
		 v,v2,v,0.,
		 0., 0., 0., 1.,
		 table);
    t = table->t;

    for (i=1; i < n-1; i++) {
	/* next intersection */
	x2 = y[i];
	ix = floorf(x2);
	sf_eno_apply (lay[i],ix,x2-ix,&z2,NULL,FUNC);
	v  = v0[i]+g[0][i]*(x -x0[0][i])+g[1][0]*(z -x0[1][i]);
	v2 = v0[i]+g[0][i]*(x2-x0[0][i])+g[1][0]*(z2-x0[1][i]);
	
	/* current layer */
	kirmod_table(type[i], true, 
		     z2-z,x2-x,0.,
		     hypotf(g[0][i],g[1][i]),
		     g[0][i],0.,g[1][i], /* gradient */
		     v,v2,v,0.,
		     0., 0., 0., 1.,
		     table);
	t += table->t;
	
	x = x2;
	z = z2;
    }
     
    /* endvelocity */
    v  = v0[n-1]+g[0][n-1]*(x     -x0[0][n-1])+g[1][n-1]*(z     -x0[1][n-1]);
    v2 = v0[n-1]+g[0][n-1]*(end[0]-x0[0][n-1])+g[1][n-1]*(end[1]-x0[1][n-1]);
    
    /* last layer */
    kirmod_table(type[n-1], true, 
		 end[1]-z,end[0]-x,0.,
		 hypotf(g[0][n-1],g[1][n-1]),
		 g[0][n-1],0.,g[1][n-1], /* gradient */
		 v,v2,v,0.,
		 0., 0., 0., 1.,
		 table);
    t += table->t;

    return t;
}
