/* Frequency domain Radon transform. */

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
#ifdef _OPENMP
#include <omp.h>
#endif

#include "frt.h"


int main(int argc, char* argv[])
{
    int n[SF_MAX_DIM], dim;
    int n1, n2, nn, nw;
    int iw;
    float d1, d2, o1, o2, ow, dw;
    sf_file in, out, ref;
    bool inv;
    int curv, mtd, niter;
    float x0, maxoff, mu, eta, rfreq, t1, t2;
    char buf[8];
    void *h;

    sf_complex **ibuf, **obuf, **rbuf=NULL;

    sf_init(argc, argv);

    in  = sf_input("in");	/* f-s-x domain input */
    out = sf_output("out");	/* f-s-p domain output */

    if (SF_COMPLEX != sf_gettype(in)) 
	sf_error("FX input need complex type");

    dim = sf_filedims(in,n);
    if(dim<2) sf_error("input dim should >= 2");
    nw = n[dim-1];
    if(dim == 2 ) nn=1;
    else {
	nn = n[1];
	for(n1=2; n1<dim-1; n1++) nn *= n[n1];
    }
    n1 = n[0];

    sprintf(buf,"o%d",dim);
    if (!sf_histfloat(in, buf, &ow)) sf_error("No o%d= in input", dim);
    sprintf(buf,"d%d",dim);
    if (!sf_histfloat(in, buf ,&dw)) sf_error("No d%d= in input", dim);

    if (!sf_histfloat(in,"o1",&o1)) sf_error("No o1= in input");
    if (!sf_histfloat(in,"d1",&d1)) sf_error("No d1= in input");

    if (!sf_getint("curv",&curv)) curv=0;
    /* 0: linear; 1:parabolic */
    if (!sf_getbool("inv",&inv)) inv=false;
    /* if y, perform inverse operation */
    if (!sf_getint("mtd",&mtd)) mtd=0;
    /* 0: fart; 1:firt; 2:fhrt; 3:fcrt */
    if(mtd>0)
    {
	if (!sf_getfloat("mu",&mu)) mu=0.05; 
	/* mu: for firt, fhrt, fcrt */
    }
    if(mtd>1)
    {
	if (!sf_getfloat("eta",&eta)) eta=0.05; 
	/* eta: for fhrt, fcrt */
	if (!sf_getint("niter",&niter)) niter=5; 
	/* sparse iterations: for fhrt, fcrt */
    }
    if(mtd>2)
    {
	ref = sf_input("ref");
	/* reference frequency slice for fcrt */
    } else {
	ref = NULL;
    }


    if (!sf_getint("np",&n2)) n2=0; 
    /* stepout number */
    if (!sf_getfloat("op",&o2)) o2=0; 
    /* first stepout (moveout at "ref") */
    if (!sf_getfloat("dp",&d2)) d2=0; 
    /* stepout interval */

	if(inv)
	{
		t1 = fabs(o2+d2*(n2-1));
		t2 = fabs(o2);
		maxoff=(t1>t2?t1:t2);  
	}else{
		t1 = fabs(o1+d1*(n1-1));
		t2 = fabs(o1);
		maxoff=(t1>t2?t1:t2);  
	}
//    maxoff = fabs((inv?(o2+(n2-1)*d2):(o1+d1*(n1-1))));
    if (!sf_getfloat("x0",&x0)) x0 = maxoff; 
    /* reference offset */


    if(n2 <= 0) sf_error("np should larger than 0");
    sf_putint  (out, "n1", n2);
    sf_putfloat(out, "o1", o2);
    sf_putfloat(out, "d1", d2);

    ibuf = sf_complexalloc2(n1, nn);
    obuf = sf_complexalloc2(n2, nn);

    if(inv)
    {	
	h =	sf_frt_init(mtd, curv,
			    o2/x0, d2/x0, n2, 	// offset
			    o1, d1, n1				// stepout
	    );
    }else{
	h =	sf_frt_init(mtd, curv,
			    o1/x0, d1/x0, n1, 	// offset
			    o2, d2, n2				// stepout
	    );
	if(mtd == 3)
	{
	    sprintf(buf,"o%d",dim);
	    if (!sf_histfloat(ref, buf, &rfreq)) 
		sf_error("No o%d= in input", dim);
	    sf_complexread(ibuf[0], n1*nn, ref);
	    rbuf = sf_complexalloc2(n2, nn);
	    sf_fhrt(h, rfreq, ibuf, rbuf, nn, mu, eta, niter);
	    sf_fhrt_reg(rbuf, n2, nn, mu, eta);
	}
    }

#ifdef _OPENMP
#pragma omp parallel for  ordered       \
    schedule(dynamic,5)          \
    private(iw)                  
#endif
    for(iw=0; iw<nw; iw++)
    {
	sf_complexread(ibuf[0], n1*nn, in);
	if(inv)	
	    sf_ifart(h, ow+dw*iw, ibuf, obuf, nn);
	else
	{
	    switch(mtd)
	    {
		case 1:
		    sf_firt(h, ow+dw*iw, ibuf, obuf, nn, mu);
		    break;
		case 2:
		    sf_fhrt(h, ow+dw*iw, ibuf, obuf, nn, mu, eta, niter);
		    break;
		case 3:
		    sf_fcrt(h, ow+dw*iw, ibuf, obuf, rbuf, nn);
		    break;
		default:
		    sf_fart(h, ow+dw*iw, ibuf, obuf, nn);
		    break;
	    }
	}

	sf_complexwrite(obuf[0], n2*nn, out);
    }

    sf_frt_close(h);

    free(ibuf[0]);
    free(ibuf);
    free(obuf[0]);
    free(obuf);
    return 0;
}

