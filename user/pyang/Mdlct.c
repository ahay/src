/* discrete linear chirp transfrom (DLCT)
*/
/*
  Copyright (C) 2013  Xi'an Jiaotong University (Pengliang Yang)

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
#include <complex.h>
#include <fftw3.h>

#ifdef _OPENMP
#include <omp.h>
#endif


int main(int argc, char* argv[])
{
  bool inv, verb;
  int N,L;
  float C;
  float *sig;
  sf_complex **Sc,**hg;
    sf_file in, out;

    /* input and output variables */
    sf_init(argc,argv);
    in = sf_input("in");
    out = sf_output("out");

    if (!sf_getbool("inv",&inv)) inv=false; /* if y, do inverse transform */
    if (!sf_getbool("verb",&verb)) verb = false;/* verbosity flag */
    if (!sf_getfloat("C",&C)) C=0.01;/* C=2*Lambda/nL, unit slice */

    if(!inv){
      	/* then: in is signal itself, out will be DCLT coefficients. */
     	if (!sf_histint(in,"n1",&N)) sf_error("No n1= in input");
      	if (!sf_getint("L",&L)) sf_error("No L");
      	/* number of freq slices each freq*/
      	sf_putint(out,"n1",N);
      	sf_putint(out,"n2",L);
    }else{
	/*then: in is DLCT coefficients, out will be signal.*/
	if (!sf_histint(in,"n1",&N)) sf_error("No n1= in input");
	if (!sf_histint(in,"n2",&L)) sf_error("No n2= in input");
	sf_putint(out,"n1",N);
    }


    fftwf_plan p;
    sig=sf_floatalloc(N);
    Sc=sf_complexalloc2(L,N);
    hg=sf_complexalloc2(L,N);

    if(!inv){
      sf_floatread(sig,N,in);

      for(int l=-L/2;l<L/2;l++){
	for(int n=0;n<N;n++){
	  hg[l+L/2][n]=sig[n]*expf(-I*2*SF_PI*C*l*n*n/N);
	}
      }
      for(int l=-L/2;l<L/2;l++){
	p = fftwf_plan_dft_1d(N, hg[l+L/2], Sc[l+L/2], FFTW_FORWARD, FFTW_ESTIMATE);
	fftwf_execute(p); 
      }

	sf_complexwrite(Sc[0],L*N,out);
    }else{
	sf_complexread(Sc[0],L*N,in);

      for(int l=-L/2;l<L/2;l++){
	p = fftwf_plan_dft_1d(N, Sc[l+L/2], hg[l+L/2], FFTW_BACKWARD, FFTW_ESTIMATE);
	fftwf_execute(p); 
      }

      for(int n=0;n<N;n++){
	sig[n]=0.0;
        for(int l=-L/2;l<L/2;l++){
	    sig[n]+=crealf(hg[l+L/2][n]*expf(I*2*SF_PI*C*l*n*n/N));
	}
      }
	
	sf_floatwrite(sig,N,out);
    }

    fftwf_destroy_plan(p);

  return 0;
}
