/* Generate shifts with offset for 1-D regularized autoregression. */
/*
  Copyright (C) 2012 Jilin University
  
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
    int n1, is, ns, nf, esize, pad;
    off_t i2, n2;
    char *trace, *zero;
    sf_file in, shift;

    sf_init(argc,argv);
    in = sf_input("in");
    shift = sf_output("out");

    if (!sf_histint(in,"n1",&n1)) sf_error("No n1= in input");
    n2 = sf_leftsize(in,1);

    if (!sf_histint(in,"esize",&esize)) {
	esize=sf_esize(in);
    } else if (0>=esize) {
	sf_error("cannot handle esize=%d",esize);
    }

    if (!sf_getint("ns",&ns)) sf_error("Need ns=");
    /* number of shifts */

    if (!sf_getint("nf",&nf)) nf = 1;
    /* offset of first shift */

    if ((nf+ns-1) >= n1) sf_error("(nf+ns-1)=%d is too large",nf+ns-1);

    sf_putint(shift,"n2",ns);
    sf_shiftdim(in, shift, 2);

    sf_fileflush(shift,in);
    sf_setform(in,SF_NATIVE);
    sf_setform(shift,SF_NATIVE);

    n1 *= esize;

    trace = sf_charalloc(n1);
    zero = sf_charalloc(n1);
    memset(zero,0,n1);

    for (i2=0; i2 < n2; i2++) {
	sf_charread(trace,n1,in);
	for (is=nf; is <= (nf+ns-1); is++) {
	    pad = is*esize;
	    sf_charwrite(zero,pad,shift);
	    sf_charwrite(trace,n1-pad,shift);
	}
    }

    exit(0);
}
