/* Signal and noise separation using adaptive PEFs. */
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

#include "apefsignoi.h"

int main(int argc, char* argv[])
{
    int niter, n1, n2, i;
    int sf1, sf2, sf3, sf4, nf1, nf2, nf3, nf4, n3, i3, sf, nf;
    float *mm, *dd, *sfilter, *nfilter, eps;
    bool verb;
    sf_file in, out, sfilt, nfilt;

    sf_init (argc,argv);
    in = sf_input("in");
    sfilt = sf_input("sfilt");
    nfilt = sf_input("nfilt");
    out = sf_output("out");

    if (!sf_getint("niter",&niter)) niter=100;
    /* Number of iterations */

    if (!sf_getfloat("eps",&eps)) eps=0.;
    /* regularization parameter */

    if (!sf_getbool("verb",&verb)) verb = false;
    /* verbosity flag */

    /* read input file parameters */
    if (!sf_histint(in,"n1",&n1)) sf_error("No n1= in input");
    if (!sf_histint(in,"n2",&n2)) sf_error("No n2= in input");
    n3 = sf_leftsize(in,2);

    if (1==n1 && 1==n2) {
	if (!sf_histint(in,"n3",&n1)) sf_error("No n3= in input");
	if (!sf_histint(in,"n4",&n2)) sf_error("No n4= in input");
	n3 = sf_leftsize(in,4);
    }

    sf_fileflush(out,in);  /* copy data dimensions */

    /* input filter */
    if (!sf_histint(sfilt,"n1",&sf1)) sf_error("No n1= in sfilt");
    if (!sf_histint(sfilt,"n2",&sf2)) sf_error("No n2= in sfilt");
    if (!sf_histint(sfilt,"n3",&sf3)) sf_error("No n3= in sfilt");
    if (!sf_histint(sfilt,"n4",&sf4)) sf_error("No n4= in sfilt");
    
    if (!sf_histint(nfilt,"n1",&nf1)) sf_error("No n1= in nfilt");
    if (!sf_histint(nfilt,"n2",&nf2)) sf_error("No n2= in nfilt");
    if (!sf_histint(nfilt,"n3",&nf3)) sf_error("No n3= in nfilt");
    if (!sf_histint(nfilt,"n4",&nf4)) sf_error("No n4= in nfilt");

    if (nf3!=n1 || nf4!=n2 || sf3!=n1 || sf4!=n2) 
	sf_error("need n1==sf3==nf3 && n2==sf4==nf4!");

    sf = sf1*sf2*sf3*sf4;
    nf = nf1*nf2*nf3*nf4;
    sfilter = sf_floatalloc(sf);
    nfilter = sf_floatalloc(nf);
    mm = sf_floatalloc(n1*n2);
    dd = sf_floatalloc(n1*n2);

    for (i=0; i < n1*n2; i++) {
	mm[i] = 0.;
	dd[i] = 0.;
    }

    for (i3=0; i3 < n3; i3++) {
	sf_warning("slice %d of %d",i3+1,n3);
	sf_floatread(dd,n1*n2,in);
	sf_floatread(sfilter,sf,sfilt);
	sf_floatread(nfilter,nf,nfilt);
	
	apefsignoi_lop(niter, sf1, sf2, nf1, nf2, n1, n2, 
		       sfilter, nfilter, mm, dd, eps, verb);
	
	sf_floatwrite(mm,n1*n2,out);
    }

    exit (0);
}

/* 	$Id$	 */
