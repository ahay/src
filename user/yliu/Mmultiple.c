/* 2-D shot gather multiple prediction */
/*
  Copyright (C) 2009 University of Texas at Austin
   
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

int main(int argc, char *argv[])
{
    int nw, n1, n2, n4, iw, i1, i2, i4, s1, x1, s2, x2, m;
    float d1;
    bool verb;

    sf_complex *dd, *mm;
    sf_file in, out;

    sf_init(argc,argv);

    in = sf_input("in");
    out = sf_output("out");

    if (SF_COMPLEX != sf_gettype(in)) sf_error("Need complex input");

    if (!sf_histint(in,"n1",&n1)) sf_error("No n1= in input");
    if (!sf_histint(in,"n2",&n2)) sf_error("No n2= in input");
    if (!sf_histint(in,"n3",&nw)) sf_error("No n3= in input");

    if (!sf_histfloat(in,"d1",&d1)) sf_error("No d1= in input");

    n4 = sf_leftsize(in,3);

    if (!sf_getbool("verb",&verb)) verb=false;
    /* verbosity flag */

    sf_putint(out,"n1",(2*n1-1));
    sf_putfloat(out,"d1",d1);
    sf_putfloat(out,"o1",(1-n1)*d1);
    (void) sf_shiftdim(in, out, 1);

    dd = sf_complexalloc(n1*n2);            /* data space */
    mm = sf_complexalloc((2*n1-1)*n1*n2);   /* multiple space */

    /* loop over n4 */
    for (i4=0; i4 < n4; i4++) {
	for (iw=0; iw < nw; iw++) { /* loop over frequency */
	    if (verb) sf_warning("frequency %d of %d",iw+1,nw);
	    sf_complexread(dd,n1*n2,in);

	    for (i2=0; i2 < n2; i2++) {
		for (i1=0; i1 < n1; i1++) {
		    for (m=(-1*n1+1); m < n1; m++) {
			
			if (m < 0) {
			    x1 = -1*m;
			    s1 = i2 + m;
			    x2 = i1 - m;
			    s2 = m + i2;
			} else if ((i1-m) < 0) {
			    x1 = m;
			    s1 = i2;
			    x2 = m - i1;
			    s2 = i2 + i1;
			} else {
			    x1 = m;
			    s1 = i2;
			    x2 = i1 - m;
			    s2 = m + i2;
			}
			if (s1 >= 0 && s1 < n2 && x1 >= 0 && x1 < n1 && 
			    s2 >= 0 && s2 < n2 && x2 >= 0 && x2 < n1 ) {
#ifdef SF_HAS_COMPLEX_H
			    mm[i2*(2*n1-1)*n1+i1*(2*n1-1)+m+n1-1] = dd[s1*n1+x1]*dd[s2*n1+x2];
#else
			    mm[i2*(2*n1-1)*n1+i1*(2*n1-1)+m+n1-1] = sf_cmul(dd[s1*n1+x1],dd[s2*n1+x2]);
#endif
			} else {
			    mm[i2*(2*n1-1)*n1+i1*(2*n1-1)+m+n1-1] = 0.;
			}
		    }
		}
	    }
	    sf_complexwrite(mm,(2*n1-1)*n1*n2,out);
	}       
    }
		
    exit(0);
}
/* 	$Id$	 */
