/* 2-D dip estimation using analytical equation. */

/*
  Copyright (C) 2012 Zhonghuan Chen, UT Austin, Tsinghua University
  
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

#include "fastpwd.h"

int main(int argc, char* argv[])
{
    bool verb;
    int n1, n2, n12, i3, n3, rect[2], niter, n[2];
    float **wav, **dip, **num, **den, eps;
    sf_file inp, out;

    sf_init(argc,argv);
    inp = sf_input("in");
    out = sf_output("out");

    if (SF_FLOAT != sf_gettype(inp)) sf_error("Need float input");
    
    if (!sf_histint(inp,"n1",&n1)) sf_error("No n1= in input");
    if (!sf_histint(inp,"n2",&n2)) sf_error("No n2= in input");
    n3 = sf_leftsize(inp,2);
    n12 = n1*n2;
    n[0] = n1;
    n[1] = n2;

    if (!sf_getint("rect1",&rect[0])) rect[0]=1;
    if (!sf_getint("rect2",&rect[1])) rect[1]=1;
    /* smoothing radius */

    if (!sf_getint("niter",&niter)) niter=100;
    /* number of iterations */

    if (!sf_getbool("verb",&verb)) verb=true;
    /* verbosity */

    if (!sf_getfloat("eps",&eps)) eps=0.0f;
    /* regularization */

    wav = sf_floatalloc2(n1,n2);
    dip = sf_floatalloc2(n1,n2);
    num = sf_floatalloc2(n1,n2);
    den = sf_floatalloc2(n1,n2);

    fastpwd_init(n1,n2);
    sf_divn_init(2, n12, n, rect, niter, verb);

    for (i3=0; i3 < n3; i3++) {
	sf_floatread(wav[0],n12,inp);

	fastpwd(wav,num,den);

	/* smooth division */
	sf_divne (num[0], den[0], dip[0], eps);

	sf_floatwrite(dip[0],n12,out);
    }

    exit(0);
}
    
