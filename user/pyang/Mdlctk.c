/* DLCT-K iterative shrinkage
discrete linear chirp transfrom for t-->f, fourier transform for x-->kx
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
#include "pthresh.h"

#include "dlctk.h"
int main(int argc, char* argv[])
{
    bool inv, verb;
    char *mode;
    int L, n1, n2, niter, iter, i, nthr;
    float C, thr, p, pclip, *dat, *tmp;
    sf_complex *coeffs;
    sf_file in, out;

    /* input and output variables */
    sf_init(argc,argv);
    in = sf_input("in");
    out = sf_output("out");

    if (!sf_getbool("inv",&inv)) inv=false;
    /* if y, do inverse transform (Here adjoint is the same as inverse!) */
    if (!sf_getbool("verb",&verb)) verb = false;/* verbosity flag */
    if (!sf_getint("niter",&niter)) niter=1;
    if (!sf_histint(in,"n1",&n1)) sf_error("No n1= in input"); /*trace length */
    if (!sf_histint(in,"n2",&n2)) sf_error("No n2= in input"); /*number of traces */
    /* n1,n2 are assumed to be 2^k */ 
    if (!sf_getint("L",&L)) sf_error("No L");
    if (!sf_getfloat("C",&C)) C=0.005;/* C=2*Lambda/L, unit slice */
    if (!sf_getfloat("pclip",&pclip)) 	pclip=99.;
    /* starting data clip percentile (default is 99)*/
    if ( !(mode=sf_getstring("mode")) ) mode = "exp";
    /* thresholding mode: 'hard', 'soft','pthresh','exp';
	'hard', hard thresholding;	'soft', soft thresholding; 
	'pthresh', generalized quasi-p; 'exp', exponential shrinkage */
    if (!sf_getfloat("p",&p)) 		p=0.35;
    /* norm=p, where 0<p<=1 */;
    if (strcmp(mode,"soft") == 0) 	p=1;
    else if (strcmp(mode,"hard") == 0) 	p=0;

    dat=sf_floatalloc(n1*n2);
    coeffs=sf_complexalloc(n1*n2*L);
    tmp=sf_floatalloc(n1*n2*L);// save thresholds
    sf_floatread(dat, n1*n2, in);

    for(iter=0; iter<niter; iter++)
    {
	dlctk_lop(true, false, n1*n2*L, n1*n2, coeffs, dat);

	for(i=0; i<n1*n2*L; i++) tmp[i]=cabsf(coeffs[i]);
   	nthr = 0.5+n1*n2*L*(1.-0.01*pclip); 
    	if (nthr < 0) nthr=0;
    	if (nthr >= n1*n2*L) nthr=n1*n2*L-1;
	thr=sf_quantile(nthr,n1*n2*L,tmp);
	thr*=(float)(niter-iter)/niter;
    	for(i=0; i<n1*n2*L; i++) sf_cpthresh(coeffs, n1*n2*L, thr, p, mode);

	dlctk_lop(false, false, n1*n2*L, n1*n2, coeffs, dat);

	if (verb) sf_warning("iteration %d", iter);
    }

    sf_floatwrite(dat, n1*n2, out);

    free(dat);
    free(coeffs);
    free(tmp);

    exit(0);
}
