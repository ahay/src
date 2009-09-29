/* Fast mean-curvature flow. */
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

static int n1, n2;
static void rhs(void* par, float* x, float* y);
static int term(void* par, float* x);

int main(int argc, char* argv[])
{
    int n12, i, i1, i2, rect, order, niter;
    float **img, **grad, **grad1, **grad2, *y, **x;
    float med, t, t2, tol[2], band;
    sf_eno2 g[2];
    sf_file in, out;
    
    sf_init(argc,argv);
    in = sf_input("in");
    out = sf_output("out");

    if (!sf_histint(in,"n1",&n1)) sf_error("No n1= in input");
    if (!sf_histint(in,"n2",&n2)) sf_error("No n2= in input");
    n12 = n1*n2;

    if (!sf_getint("rect",&rect)) rect=3;
    /* smoothing radius */
    t = (rect*rect-1)/6.0;
    sf_runge_init(2,1,t/2);

    if (!sf_getint("order",&order)) order=3;
    /* interpolation order */
    if (!sf_getfloat("tol",&tol[0])) tol[0]=0.1;
    /* error tolerance */
    tol[1]=tol[0];
    if (!sf_getint("niter",&niter)) niter=100;
    /* number of iterations */
    if (!sf_getfloat("band",&band)) band=1.;
    /* narrow band */

    img = sf_floatalloc2(n1,n2);
    grad = sf_floatalloc2(n1,n2);
    grad1 = sf_floatalloc2(n1,n2);
    grad2 = sf_floatalloc2(n1,n2);
    y = sf_floatalloc(n12);
    x = sf_floatalloc2(2,n12);

    g[0] = sf_eno2_init (order, n1, n2);
    g[1] = sf_eno2_init (order, n1, n2);

    sf_floatread(img[0],n12,in);

    /* compute gradient squared */
    sf_sobel2 (n1,n2,img,grad);

    /* normalize by median */
    for (i=0; i < n12; i++) {
	y[i] = grad[0][i];
    }

    med = sf_quantile(n12/2,n12,y);

    for (i=0; i < n12; i++) {
	grad[0][i] = -0.5*logf(1.+grad[0][i]/med);
    }

    /* compute gradient */
    sf_sobel (n1,n2,grad,grad1,grad2);
    sf_eno2_set (g[0], grad1);
    sf_eno2_set (g[1], grad2);

    /* convection */
    i=0;
    for (i2=0; i2 < n2; i2++) {
	for (i1=0; i1 < n1; i1++) {
	    if (fabsf(img[i2][i1]) < band) {
		x[i][0] = i1;
		x[i][1] = i2;
		t2 = sf_ode23 (t,tol, x[i], g, rhs, term);
		if (t2 >= t) {
		    y[i] = img[i2][i1];
		    i++;
		}
	    }
	}
    }

    /* diffusion */
    sf_int2_init (x, 0,0,1,1,n1,n2, sf_lin_int, 2, i);
    sf_triangle2_init(rect,rect, n1, n2, 1);
    sf_conjgrad_init(n12, n12, i, i, 1.0/i, 1.e-8, true, false);
    sf_conjgrad(NULL, 
		sf_int2_lop, sf_triangle2_lop, grad[0], img[0], y, niter);

    sf_floatwrite(img[0],n12,out);

    exit(0);
}

static void rhs(void* par, float* x, float* y)
{
    int i, j;
    sf_eno2 *g;

    g = (sf_eno2 *) par;
    i = x[0];
    j = x[1];

    sf_eno2_apply (g[0], i, j, x[0]-i, x[1]-j, &y[0], NULL, FUNC);
    sf_eno2_apply (g[1], i, j, x[0]-i, x[1]-j, &y[1], NULL, FUNC);
}

static int term (void* par, float* x)
{
    return (x[0] < 0 || x[0] > n1-1 || 
	    x[1] < 0 || x[1] > n2-1);
}
