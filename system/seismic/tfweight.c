/* time and frequency weights */
/*
  Copyright (C) 2004 University of Texas at Austin
  
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

static int nt, nw, nf, nx;
static float *w, *f, *tmp;
static kiss_fftr_cfg forw, invs;
static kiss_fft_cpx *cdata;

void tfweight_init(int n1, int nw1, int n2,
		   float* ww /* [n1*n2] time weight */,
		   float* ff /* [nw1*n2] frequency weight */)
/*< initialoze >*/
{
    nt = n1;
    nw = nw1;
    nx = n2;
    w = ww;
    f = ff;

    nf = 2*(nw-1);

    forw = kiss_fftr_alloc(nf,0,NULL,NULL);
    invs = kiss_fftr_alloc(nf,1,NULL,NULL);

    tmp = sf_floatalloc(nf);    
    cdata =  (kiss_fft_cpx*) sf_alloc(nw,sizeof(kiss_fft_cpx));
}

void tfweight_close(void)
/*< free allocated storage >*/
{
    free(tmp);
    free(cdata);
    free(forw);
    free(invs);
}

void tfweight_lop(bool adj, bool add, int nxx, int nyy, float* x, float* y) 
/*< linear operator >*/
{
    int iw, i1, i2;

    if (nxx != nt*nx || nyy != nt*nx) sf_error("%s: Wrong size",__FILE__);

    sf_adjnull(adj,add,nxx,nyy,x,y);
    
    for (i2=0; i2 < nx; i2++) {
	if (adj) {
	    for (i1=0; i1 < nt; i1++) {
		tmp[i1] = w[i1+i2*nt]*y[i1+i2*nt];
	    }
	    for (i1=nt; i1 < nf; i1++) {
		tmp[i1] = 0.0f;
	    }
	    kiss_fftr(forw, tmp, cdata);
	    for (iw=0; iw < nw; iw++) {
		cdata[iw]=sf_crmul(cdata[iw],f[iw+i2*nw]/nf);
	    }
	    kiss_fftri(invs, cdata, tmp);
	    for (i1=0; i1 < nt; i1++) {
		x[i1+i2*nt] += tmp[i1];
	    }
	} else {
	    for (i1=0; i1 < nt; i1++) {
		tmp[i1] = x[i1+i2*nt];
	    }
	    for (i1=nt; i1 < nf; i1++) {
		tmp[i1] = 0.0f;
	    }
	    kiss_fftr(forw, tmp, cdata);
	    for (iw=0; iw < nw; iw++) {
		cdata[iw]=sf_crmul(cdata[iw],f[iw+i2*nw]/nf);
	    }
	    kiss_fftri(invs, cdata, tmp);
	    for (i1=0; i1 < nt; i1++) {
		y[i1+i2*nt] += w[i1+i2*nt]*tmp[i1];
	    }
	}
    }
}
	
