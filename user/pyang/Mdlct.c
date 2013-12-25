/* discrete linear chirp transfrom (DLCT)
Note: In my implementation:, to make the adjoint as same as the inverse,
I normalized the forward transform of DLCT with a factor sqrt(N*L).
*/
/*
  Copyright (C) 2013  Xi'an Jiaotong University, UT Austin (Pengliang Yang)

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

#include "dlct.h"

int main(int argc, char* argv[])
{
    bool inv, verb;
    int L,N, n2;
    float C;
    float *sig;
    sf_complex *Sc;
    sf_file in, out;

    /* input and output variables */
    sf_init(argc,argv);
    in = sf_input("in");
    out = sf_output("out");

    if (!sf_getbool("inv",&inv)) inv=false;
    /* if y, do inverse transform (Here adjoint is the same as inverse!) */
    if (!sf_getbool("verb",&verb)) verb = false;/* verbosity flag */
    if (!sf_getfloat("C",&C)) C=0.005;/* C=2*Lambda/L, unit slice */
    if (!sf_histint(in,"n1",&N)) sf_error("No n1= in input"); /*length of signal */
    /* N is assumed to be 2^k */ 	
 
    if(!inv){
      	/* then: in is signal itself, out will be DCLT coefficients. */
	n2 = sf_leftsize(in,1);
      	if (!sf_getint("L",&L)) sf_error("No L");
	sf_shiftdim(in, out, 1);
      	sf_putint(out,"n1",N);
      	sf_putint(out,"n2",L);
	sf_settype(out,SF_COMPLEX);
    }else{
	/*then: in is DLCT coefficients, out will be signal.*/
	n2 = sf_leftsize(in,2);
	if (!sf_histint(in,"n1",&N)) sf_error("No n1= in input");
	if (!sf_histint(in,"n2",&L)) sf_error("No n2= in input");
	sf_unshiftdim(in, out, 2);
	sf_putint(out,"n1",N);
	sf_settype(out,SF_FLOAT);
    }
	
    sig = sf_floatalloc(N);
    Sc = sf_complexalloc(N*L);

    for (int i2=0; i2 < n2; i2++)  {
	sf_warning("slice %d of %d;",i2+1,n2);
	    if(!inv){
		sf_floatread(sig,N,in);
		forward_dlct(N, L, C, sig, Sc);
		sf_complexwrite(Sc, L*N, out);
	    }else {
		sf_complexread(Sc,L*N,in);
		inverse_dlct(N, L, C, sig, Sc);
		sf_floatwrite(sig, N, out);
	    }
    }
    sf_warning(".");

    free(sig);
    free(Sc);

    exit(0);
}
