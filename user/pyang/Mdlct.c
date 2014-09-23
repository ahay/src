/* Discrete linear chirp transfrom (DLCT)
Note: I normalized the forward transform of DLCT with a factor sqrt(N*L)
 to ensure the adjoint is the same as the inverse!
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

  Reference: Alkishriwo, Osama A., and Luis F. Chaparro. "A discrete 
	linear chirp transform (DLCT) for data compression." Information 
	Science, Signal Processing and their Applications (ISSPA), 2012 
	11th International Conference on. IEEE, 2012.
*/
#include <rsf.h>
#include <complex.h>

#include "dlct.h"

int main(int argc, char* argv[])
{
    bool inv, verb;
    int L, n1, n2, i1, i2;
    float C;
    float *sig;
    sf_complex *mm, *dd;
    sf_file in, out;

    /* input and output variables */
    sf_init(argc,argv);
    in = sf_input("in");
    out = sf_output("out");

    if (!sf_getbool("inv",&inv)) inv=false;
    /* if y, do inverse transform (Here adjoint is the same as inverse!) */
    if (!sf_getbool("verb",&verb)) verb = false;
    /* verbosity flag */
    if (!sf_getfloat("C",&C)) C=0.005;/* C=2*Lambda/L, unit slice */
    if (!sf_histint(in,"n1",&n1)) sf_error("No n1= in input"); 
    /*length of signal, n1 is assumed to be 2^k */ 	
 
    if(!inv){
      	/* then: in is signal itself, out will be DCLT coefficients. */
	n2 = sf_leftsize(in,1);
      	if (!sf_getint("L",&L)) sf_error("No L");
	/* number of discretized chirp rate*/
	sf_shiftdim(in, out, 1);
      	sf_putint(out,"n1",n1);
      	sf_putint(out,"n2",L);
	sf_settype(out,SF_COMPLEX);
    }else{
	/* then: in is DLCT coefficients, out will be signal.*/
	n2 = sf_leftsize(in,2);
	if (!sf_histint(in,"n1",&n1)) sf_error("No n1= in input");
	/* length of signal in freq domain, from -f to +f */
	if (!sf_histint(in,"n2",&L)) sf_error("No n2= in input");
	/* number of discretized chirp rate*/
	sf_unshiftdim(in, out, 2);
	sf_putint(out,"n1",n1);
	sf_settype(out,SF_FLOAT);
    }
	
    sig = sf_floatalloc(n1);
    dd = sf_complexalloc(n1);
    mm = sf_complexalloc(n1*L);

    dlct_init(n1, L, C);

    for (i2=0; i2 < n2; i2++)  {
	if (verb) sf_warning("trace %d of %d;",i2+1,n2);
	    if(inv){
		sf_complexread(mm,L*n1,in);
		dlct_lop(true, false, n1*L, n1, mm, dd);
		for(i1=0; i1<n1; i1++) sig[i1]=dd[i1];
		sf_floatwrite(sig, n1, out);
	    }else {
		sf_floatread(sig,n1,in);
		for(i1=0; i1<n1; i1++) dd[i1]=sig[i1];
		dlct_lop(false, false, n1*L, n1, mm, dd);
		sf_complexwrite(mm, L*n1, out);
	    }
    }


    dlct_close();
    free(sig);
    free(mm);
    free(dd);

    exit(0);
}
