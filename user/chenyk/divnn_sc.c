/* Smooth division with several components and nonstationary smoothing, required by Mltftn.c*/
/*
  Copyright (C) 2010 University of Texas at Austin
  
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

static int n0, n1, n2;
static sf_triangle *tr;
static float *tmp;

void ntrianglen2_init(int nw      /* number of components */, 
		      int n       /* data size */,
		      int nf      /* first dimension */,
		      int *nbox   /* smoothing radius [nw] */)
/*< initialization >*/
{
    int i2;

    n0 = nf;
    n1 = n;
    n2 = nw;

    tmp = sf_floatalloc(n1*n2);

    tr = (sf_triangle *) sf_alloc(nw,sizeof(*tr));

    for (i2=0; i2 < n2; i2++) {
	tr[i2] = sf_triangle_init (nbox[i2],n0,false);
    }
}

// void ntrianglen2n_init(int nw      /* number of components */, 
// 		      int n       /* data size */,
// 		      int nf      /* first dimension */,
// 		      int *nbox   /* smoothing radius [nw] */)
// /*< initialization (N-dimensional smoothing)>*/
// {
//     int i2;
// 
//     n0 = nf;
//     n1 = n;
//     n2 = nw;
// 
//     tmp = sf_floatalloc(n1*n2);
// 
//     tr = (sf_triangle *) sf_alloc(nw,sizeof(*tr));
// 
//     for (i2=0; i2 < n2; i2++) {
// 	tr[i2] = sf_triangle_init (nbox[i2],n0,false);
//     }
//     
//     for (i=0; i < dim; i++) {
// 	tr[i] = (nbox[i] > 1)? ntriangle_init (nbox[i],ndat[i]): NULL;
// 	s[i] = nd;
// 	nd *= ndat[i];
//     }
//     
// }

void ntriangle2_close(void) 
/*< free allocated storage >*/
{
    int i2;

    for (i2=0; i2 < n2; i2++) {
		sf_triangle_close(tr[i2]);
    }

    free(tr);
    free(tmp);
}

void ntriangle2_lop (bool adj, bool add, int nx, int ny, float *x, float *y)
/*< combined linear operator (1-dimensional smoothing) >*/
{
    int i, j, i2;       
    
    if (nx != ny || nx != n1*n2) 
	sf_error("%s: Wrong size (nx=%d ny=%d n1=%d n2=%d)",
		 __FILE__,nx,ny,n1,n2);

    sf_adjnull (adj, add, nx, ny, x, y);

    if (adj) {
		for (i=0; i < nx; i++) {
			tmp[i] = y[i];
		}
    } else {
		for (i=0; i < nx; i++) {
			tmp[i] = x[i];
		}
    }

    for (i2=0; i2 < n2; i2++) {
		for (j=0; j < n1/n0; j++) {
			sf_smooth2 (tr[i2],j*n0,1,false,tmp+i2*n1);
		}
    }

    if (adj) {
    	for (i=0; i < nx; i++) {
			x[i] += tmp[i];
		}
    } else {
		for (i=0; i < nx; i++) {
			y[i] += tmp[i];
		}
    }
}

// void ntriangle2n_lop (bool adj, bool add, int nx, int ny, float *x, float *y)
// /*< combined linear operator (N-dimensional smoothing) >*/
// {
//     int i, j, i2;       
//     
//     if (nx != ny || nx != n1*n2) 
// 	sf_error("%s: Wrong size (nx=%d ny=%d n1=%d n2=%d)",
// 		 __FILE__,nx,ny,n1,n2);
// 
//     sf_adjnull (adj, add, nx, ny, x, y);
// 
//     if (adj) {
// 		for (i=0; i < nx; i++) {
// 			tmp[i] = y[i];
// 		}
//     } else {
// 		for (i=0; i < nx; i++) {
// 			tmp[i] = x[i];
// 		}
//     }
// 
//     for (i=0; i < dim; i++) {
// 	if (NULL != tr[i]) {
// 	    for (j=0; j < nd/n[i]; j++) {
// 		i0 = sf_first_index (i,j,dim,n,s);
// 		    nsmooth (tr[i], i0, s[i], false, tlen[i], tsft[i], tmp);
// 		    
// 		    sf_smooth2 (tr[i2],j*n0,1,false,tmp+i2*n1);
// 	    }
// 	}
//     }
// 
//     if (adj) {
//     	for (i=0; i < nx; i++) {
// 			x[i] += tmp[i];
// 		}
//     } else {
// 		for (i=0; i < nx; i++) {
// 			y[i] += tmp[i];
// 		}
//     }
// }


static float *p;

void divnn_sc_init(int nw       /* number of components */, \
		     int ndim     /* number of dimensions */,
		     int n        /* data size */,
		     int *ndat    /* data dimensions [ndim] */, 
		     int *nbox    /* smoothing radius [nw] */,
		     float* den   /* denominator [nw*nd] */,
		     bool verb    /* verbosity flag */)
/*< initialize >*/
{
    int n2;

    n2 = n*nw;
    
    ntrianglen2_init(nw,n,ndat[0],nbox);

    sf_conjgrad_init(n2, n2, n, n, 1., 1.e-6, verb, false);
    p = sf_floatalloc (n2);
    sf_weight2_init(nw,n,den);
}

void divnn_sc_init2(int nw       /* number of components */, \
		     int ndim     /* number of dimensions */,
		     int n        /* data size */,
		     int *ndat    /* data dimensions [ndim] */, 
		     int *nbox    /* smoothing radius [nw] */,
		  	 float **rct /* triangle lengths [ndim][nd] */,
          	 int **sft /* triangle shifts [ndim][nd] */,
		     float* den   /* denominator [nw*nd] */,
		     bool verb    /* verbosity flag */)
/*< initialize >*/
{
    int n2;

    n2 = n*nw;
    
//     ntrianglen2_init(nw,n,ndat[0],nbox);

	/*initialization for the non-stationary triangle smoothing operator*/
    sf_ntrianglen_init(ndim,nbox,ndat,rct,sft,1);
    
    sf_conjgrad_init(n2, n2, n, n, 1., 1.e-6, verb, false);
    p = sf_floatalloc (n2);
    sf_weight2_init(nw,n,den);
}

void divnn_sc_close (void)
/*< free allocated storage >*/
{
    ntriangle2_close();
    
    sf_conjgrad_close();
    free (p);
    sf_weight2_close();
}

void divnn_sc (float* num  /* numerator */, 
		   float* rat  /* ratio */, 
		   int niter   /* number of iterations */)
/*< smoothly divide num/rat >*/
{
    sf_conjgrad(NULL,sf_weight2_lop,ntriangle2_lop,p,rat,num,niter);
}

void divnn_sc2 (float* num  /* numerator */, 
		   float* rat  /* ratio */, 
		   int niter   /* number of iterations */)
/*< smoothly divide num/rat >*/
{
    sf_conjgrad(NULL,sf_weight2_lop,sf_ntrianglen_lop,p,rat,num,niter);
}







