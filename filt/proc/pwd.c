#include <rsf.h>

#include "pwd.h"
#include "apfilt.h"

struct Pwd {
    int n, na;
    float **a;
};

pwd pwd_init(int n1, int nw)
{
    pwd w;

    w = (pwd) sf_alloc(1,sizeof(*w));
    w->n = n1;
    w->na = 2*nw+1;
 
    w->a = sf_floatalloc2 (n1,w->na);

    return w;
}

void pwd_close (pwd w)
{
    free (w->a[0]);
    free (w->a);
    free (w);
}

void pwd_define (bool adj, pwd w, const float* pp, float* diag, float** offd)
{
    int i;
    float b[3];
    
    for (i=0; i < w->n; i++) {
      passfilter (1, pp[i], b);
	
      if (adj) {
	  w->a[0][i] = b[2];
	  w->a[1][i] = b[1];	
	  w->a[2][i] = b[0];
      } else {
	  w->a[0][i] = b[0];
	  w->a[1][i] = b[1];	
	  w->a[2][i] = b[2];
      }
    }
    
    for (i=1; i < w->n-1; i++) {
	diag[i] += 
	    w->a[0][i-1]*w->a[0][i-1] + 
	    w->a[1][i  ]*w->a[1][i  ] + 
	    w->a[2][i+1]*w->a[2][i+1];
	offd[0][i] += 
	    w->a[0][i  ]*w->a[1][i  ] + 
	    w->a[1][i+1]*w->a[2][i+1];
	offd[1][i-1] +=
	    w->a[0][i]*w->a[2][i];
    }
    /* zero-slope */
    offd[0][0] += w->a[0][0]*(w->a[1][0]+w->a[2][0]) + w->a[1][1]*w->a[2][1];
    offd[0][w->n-1] += w->a[0][w->n-1]*w->a[2][w->n-1];
    diag[0] += 
	(w->a[1][0]+w->a[2][0])*(w->a[1][0]+w->a[2][0]) + 
	w->a[2][1]*w->a[2][1];
    diag[w->n-1] += 
	w->a[0][w->n-2]*w->a[0][w->n-2] + 
	(w->a[0][w->n-1] + w->a[1][w->n-1])*(w->a[0][w->n-1] + w->a[1][w->n-1]);
    /* zero value 
    offd[0][0] += 
	w->a[0][0]*w->a[1][0] + 
	w->a[1][1]*w->a[2][1];
    diag[0] += 
	w->a[1][0]*w->a[1][0] + 
	w->a[2][1]*w->a[2][1];
    diag[w->n-1] += 
	w->a[0][w->n-2]*w->a[0][w->n-2] + 
	w->a[1][w->n-1]*w->a[1][w->n-1];
    */
}

void pwd_set (pwd w, float* inp, float* out, float* tmp)
{
    int i;

    for (i=1; i < w->n-1; i++) {
	tmp[i] = 
	    w->a[0][i]*inp[i-1] +
	    w->a[1][i]*inp[i  ] +
	    w->a[2][i]*inp[i+1];
    }
    /* zero-slope */
    tmp[0] =
	w->a[0][0]*inp[0] +
	w->a[1][0]*inp[0] + 
	w->a[2][0]*inp[1];
    tmp[w->n-1] = 
	w->a[0][w->n-1]*inp[w->n-2] + 
	w->a[1][w->n-1]*inp[w->n-1] +
	w->a[2][w->n-1]*inp[w->n-1];
    /* zero value 
    tmp[0] = 
	w->a[1][0]*inp[0] +
	w->a[2][0]*inp[1];
    tmp[w->n-1] = 
	w->a[0][w->n-1]*inp[w->n-2] +
	w->a[1][w->n-1]*inp[w->n-1];
    */
    for (i=1; i < w->n-1; i++) {
	out[i] = 
	    w->a[0][i-1]*tmp[i-1] +
	    w->a[1][i  ]*tmp[i  ] +
	    w->a[2][i+1]*tmp[i+1];
    }
    /* zero slope */
    out[0] = 
	w->a[2][0]*tmp[0] +
	w->a[1][0]*tmp[0] +
	w->a[2][1]*tmp[1];
    out[w->n-1] =
	w->a[0][w->n-2]*tmp[w->n-2] +
	w->a[1][w->n-1]*tmp[w->n-1] +
	w->a[0][w->n-1]*tmp[w->n-1];
    
    /* zero value 
    out[0] = 
	w->a[1][0]*tmp[0] +
	w->a[2][1]*tmp[1];
    out[w->n-1] = 
	w->a[0][w->n-2]*tmp[w->n-2] +
	w->a[1][w->n-1]*tmp[w->n-1];
    */
}
    
