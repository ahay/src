#include <stdlib.h>

#include <rsf.h>

#include "kirmod.h"

#ifndef _kirmod_h

typedef enum {CONST, S2LIN, VLIN} maptype;
/* type of velocity distribution */
/*^*/

#endif

typedef struct Surface {
    int is, ih;
    float x, **ta;
} *surface;

static int ny, **map;
static surface y;

static int surface_comp(const void *a, const void *b)
/* compare by the surface coordinate */
{
    float aa, bb;

    aa = ((surface) a)->x;
    bb = ((surface) b)->x;

    if (aa <  bb) return (-1);
    if (aa == bb) return 0;
    return 1;
}

void kirmod_init(int ns, float s0, float ds /* source axis */,
		 int nh, float h0, float dh /* offset axis */)
/*< Initialize surface locations >*/ 
{
    int is, ih, iy;
    float s;
    surface yi;

    ny = ns*(nh+1);
    y = (surface) sf_alloc(ny,sizeof(*y));
    map = sf_intalloc2(ns,nh+1);

    yi = y;
    for (is=0; is < ns; is++, yi++) {
	s = s0 + is*ds;
	for (ih=0; ih < nh; ih++, yi++) {
	    yi->x = s + h0 + ih*dh;
	    yi->is = is;
	    yi->ih = ih;
	}
	yi->x = s;
	yi->is = is;
	yi->ih = nh;	
    }

    qsort(y,ny,sizeof(*y),surface_comp);

    for (iy=0; iy < ny; iy++) {
	yi = y+iy;
	map[yi->is][yi->ih] = iy;
    }
}

void kirmod_close(void) 
/*< Free allocated storage >*/
{
    int iy;
    float **ta;

    for (iy=0; iy < ny; iy++) {
	ta = y[iy].ta;
	if (NULL != ta) {
	    free(*ta);
	    free(ta);
	    ta = NULL;
	}
    }
    free (y);
    free (map[0]);
    free (map);
}

void kirmod_table (maptype type               /* velocity distribution */, 
		   int nx, float x0, float dx /* reflector axis */,
		   float *curve               /* reflector */, 
		   float *veloc               /* velocity attributes */) 
/*< Compute traveltime/amplitude map >*/
{
    int ix, iy;
    
    for (iy=0; iy < ny; iy++) {
	for (ix=0; ix < nx; ix++) {
	    switch (type) {
		case CONST:
		    break;
		default:
		    sf_error("__FILE__: case %d is not implemented",type);
		    break;
	    }
	}
    }
}

float *kirmod_map(int is, int ih, int ix) 
/*< Extract from traveltime/amplitude map >*/
{
    return y[map[is][ih]].ta[ix];
}
