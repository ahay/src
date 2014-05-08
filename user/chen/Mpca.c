/* KL transform. */

/*
  Copyright (C) 2012 Zhonghuan Chen, UT Austin, Tsinghua University
  
  This program is free software; you can redistribute it and/or modify
  it under the terms of the GNU General Public License as published by
  the Free Software Foundation; either version 2 of the License, or
  (at your option) any later version.
  
  This program is distributed in the hope that it will be useful,
  but WITHOUT ANY WA:RRANTY; without even the implied warranty of
  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
  GNU General Public License for more details.
  
  You should have received a copy of the GNU General Public License
  along with this program; if not, write to the Free Software
  Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA
*/


#include <rsf.h>

#include "svd.h"
#include "pca.h"


int main(int argc, char* argv[])
{
	sf_file in, out;
	int n1, n2, n3, nc, nc0, n0;
	int i3;
	float eta, **u1, **u2, **u, *s, **v;
	bool verb;

    sf_init(argc, argv);

    in  = sf_input("in");	/* 2D input : column vectors */
    out = sf_output("out");	/* 2D output */

    if (SF_COMPLEX == sf_gettype(in)) 
		sf_error("Complex pca not surport yet");


    if (!sf_histint(in,"n1",&n1)) sf_error("No n1= in input");
    if (!sf_histint(in,"n2",&n2)) sf_error("No n2= in input");

    if (!sf_getbool("verb",&verb) ) verb=true;
    /* verbosity */
    if (!sf_getint("nc",&nc0) ) nc0=n1;
    /* number of components [ 0 < nc < n1 ] */
	if (!sf_getfloat("eta",&eta)) eta=0.9;
	/* energy ratio for signal subspace [ 0 eta < 1 ]*/

    u1 = sf_floatalloc2(n1, n2);
    u2 = sf_floatalloc2(n1, n2);
    u = sf_floatalloc2(n1, n1);
    v = sf_floatalloc2(n2, n2);

	n0 = n1<n2? n1 : n2;
    s = sf_floatalloc(n0);

	n3 = sf_leftsize(in, 2);
/*
#ifdef _OPENMP
#pragma omp parallel for  ordered       \
	schedule(dynamic,n3/10+1)          \
	private(i3)                  
#endif */
	for(i3=0; i3<n3; i3++)
	{
		sf_floatread(u1[0], n1*n2, in);
		svd(n1, n2, u1, true, u, s, v);
		if(nc0 == 0) nc = pca_rank(n0,s, eta);
		else nc=nc0;
		pca_klt(n1, n2, nc, s, v, u1);
		pca_iklt(n1, n2, nc, u, u1, n1, u2);
		sf_floatwrite(u2[0], n1*n2, out);
		if(verb) sf_warning("%d of %d", i3, n3);
	}

    free(u1[0]);
	free(u1);
    free(u2[0]);
	free(u2);
    free(u[0]);
    free(v[0]);
    free(u);
    free(v);
    free(s);
    return 0;
}

