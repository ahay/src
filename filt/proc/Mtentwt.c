/* Tent-like weight for patching.

Takes: > wallwt.rsf windwt=windwt.rsf
*/

#include <rsf.h>

#include "tent.h"
#include "tent2.h"
#include "mkwallwt.h"

int main(int argc, char* argv[])
{
    int dim, i, j, *n, *w, *k, *a, *center, n12, w12;
    float *wind, *wall;
    char string[256];
    bool tnt;
    sf_file wallwt, windwt;

    sf_init (argc, argv);
    wallwt = sf_output("out");
    windwt = sf_output("windwt");

    sf_setformat(wallwt,"native_float");
    sf_setformat(windwt,"native_float");

    if (!sf_getint ("dim",&dim)) dim = 2;
    /* number of dimensions */

    n = sf_intalloc (dim);
    w = sf_intalloc (dim);
    k = sf_intalloc (dim);
    center = sf_intalloc (dim);
    a = sf_intalloc (dim);

    n12 = 1;
    for(j=0; j < dim; j++) {
	sprintf(string,"n%d",j+1);
	if (!sf_getint(string,&n[j])) n[j]=1;
	sf_putint(wallwt,string,n[j]);
	n12 *= n[j];
    }

    if (!sf_getints ("w",w,dim)) sf_error("Need w=");
    /* window size */

    w12 = 1;
    for(j=0; j < dim; j++) {
	if (w[j] > n[j]) w[j] = n[j];
	k[j] = 1.5 * n[j] / (w[j] - a[j] + 1.);
	a[j] = 1;
	center[j] = 1;
	w12 *= w[j];
	sprintf(string,"n%d",j+1);
	sf_putint(windwt,string,w[j]);
    }

    sf_getints ("k",k,dim);
    /* number of windows */

    sf_getints ("a",a,dim);
    /* filter size */
    
    for(j=1; j < dim; j++) {
	if (a[j] > 1) center[j-1] = a[j-1]/2;
    }

    sf_getints ("center",center,dim);

    if (!sf_getbool ("tent",&tnt)) tnt = true;
    /* if y, use tent-like weight; n, cosine weight */

    wall = sf_floatalloc(n12);
    wind = sf_floatalloc(w12);

    if (tnt) {
	tent (dim, w, center, a, wind);
    } else {
	tent2 (dim, w, wind);
    }

    sf_floatwrite (wind, w12, windwt);

    mkwallwt (dim, k, n, w, wind, wall);

    for (i=0; i < n12; i++) {
	if (wall[i] != 0.) wall[i] = 1. / wall[i];
    }

    sf_floatwrite (wall, n12, wallwt);

    exit (0);
}

