/* 1-D covariance estimator. */
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

int main(int argc, char* argv[])
{
    bool verb;
    char key[6];
    int n1, n2, n12, i1, i2;
    int n[SF_MAX_DIM], rect[SF_MAX_DIM], niter, dim, i;
    float **data, **top, **bot, eps;
    const float b1=16, b2=2, t1=6, t2=3;
    sf_file inp, rct;

    sf_init(argc, argv);
    inp = sf_input("in");
    rct = sf_output("out");

    if (SF_FLOAT != sf_gettype(inp)) sf_error("Need float input");

    dim = sf_filedims (inp,n);
    n12 = 1;
    for (i=0; i < dim; i++) {
	snprintf(key,6,"rect%d",i+1);
	if (!sf_getint(key,rect+i)) rect[i]=1;
	/*( rect#=(1,1,...) smoothing radius on #-th axis )*/ 
	n12 *= n[i];
    }
    n1 = n[0];
    n2 = n12/n[0];

    if (!sf_getint("niter",&niter)) niter=100;
    /* number of iterations */

    if (!sf_getbool("verb",&verb)) verb=true;
    /* verbosity */

    if (!sf_getfloat("eps",&eps)) eps=0.0f;
    /* regularization */

    sf_divn_init(dim, n12, n, rect, niter, verb);

    data = sf_floatalloc2(n1,n2);
    top  = sf_floatalloc2(n1,n2);
    bot  = sf_floatalloc2(n1,n2);

    n12 = n1*n2;

    sf_floatread(data[0],n12,inp);

    for (i2=0; i2 < n2; i2++) {
	bot[i2][0]=top[i2][0]=0.0f;
	bot[i2][1]=top[i2][1]=0.0f;
	for (i1=2; i1 < n1-2; i1++) {
	    bot[i2][i1]= 
		b1*(data[i2][i1+1]-data[i2][i1-1]) - 
		b2*(data[i2][i1+2]-data[i2][i1-2]);
	    top[i2][i1]= 
		t1*(data[i2][i1+1]-data[i2][i1-1]) - 
		t2*(data[i2][i1+2]-data[i2][i1-2]);
	}
	bot[i2][n1-2]=top[i2][n1-2]=0.0f;
	bot[i2][n1-1]=top[i2][n2-1]=0.0f;
    }

    sf_divne (top[0], bot[0], data[0], eps);

    sf_floatwrite(data[0],n12,rct);

    exit(0);
}

