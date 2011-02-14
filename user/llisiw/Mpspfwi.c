/* Full-waveform Inversion by Preconditioned Helmholtz Solver */
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

#include  "pspfwi.h"

int main(int argc, char* argv[])
{
    int n[SF_MAX_DIM], dim, i, nt, it, iter, niter, is, ns, iw, nw;
    float o[SF_MAX_DIM], d[SF_MAX_DIM], dw, ow;
    float *s, *ds, **source;
    sf_complex ***recv, *a, *b, *e;
    char key[4];
    sf_file shot, reco;

    sf_init(argc,argv);
    in  = sf_input("in");
    out = sf_output("out");

    /* read dimensions */
    dim = sf_filedims(in,n);

    nt = 1;
    for (i=0; i < dim; i++) {
	sprintf(key,"d%d",i+1);
	if (!sf_histfloat(in,key,d+i)) sf_error("No %s= in input",key);
	sprintf(key,"o%d",i+1);
	if (!sf_histfloat(in,key,o+i)) o[i]=0.;
	nt *= n[i];
    }

    if (dim < 3) {
	n[2] = 1; o[2] = o[1]; d[2] = d[1];
    }

    /* read initial slowness */
    s  = sf_floatalloc(nt);
    ds = sf_floatalloc(nt);
    sf_floatread(s,nt,in);

    /* read shot file */
    if (NULL == sf_getstring("shot")) sf_error("Need source");
    
    shot = sf_input("shot");

    if (!sf_histint(shot,"n2",&ns)) ns=1;
    
    source = sf_floatalloc2(nt,ns);
    sf_floatread(source[0],nt*ns,shot);
    sf_fileclose(shot);

    /* read record file */
    if (NULL == sf_getstring("record")) sf_error("Need record");
 
    reco = sf_input("record");

    if (SF_COMPLEX != sf_gettype(reco)) sf_error("Record must be complex");

    if (!sf_histint(reco,"n1",nw)) sf_error("No nw in record");
    if (!sf_histfloat(reco,"d1",dw)) sf_error("No dw in record");
    if (!sf_histfloat(reco,"o1",ow)) sf_error("No ow in record");

    recv = sf_complexalloc3(nw,n[1]*n[2],ns);
    sf_complexread(recv[0][0],nw*n[1]*n[2]*ns,reco);
    sf_fileclose(reco);

    /* temporary arrays */
    a = sf_complexalloc(nt);
    b = sf_complexalloc(nt);
    e = sf_complexalloc(nt);

    for (it=0; it < nt; it++) {
	e[it] = sf_cmplx(0.,0.);
    }

    if (!sf_getint("niter",niter)) niter=1;
    /* number of inversion iterations */

    /* loop over iterations */
    for (iter=0; iter < niter; iter++) {
	
	for (it=0; it < nt; it++) {
	    ds[it] = 0.;
	}

	/* loop over shots */
	for (is=0; is < ns; is++) {
	    
	    /* loop over frequencies */
	    for (iw=0; iw < nw; iw++) {
		w = w0 + iw*dw;

		b = helmholtz(w,s,source[is]);

		/* right-hand side: e = d-b */
		for (it=0; it < n[1]*[2]; it++) {
		    e[it*n[1]*n[2]] = recv[is][it][iw] - b[it*n[1]*n[2]];
		}
		
		a = helmholtz(w,s,e);
		
		/* slowness update: ds = -w*w*a*b */
		for (it=0; it < nt; it++) {
		    ds[it] += -w*w*a[it]*b[it];
		}
	    }
	}
	
	for (it=0; it < nt; it++) {
	    s[it] += ds[it];
	}

    }

    exit(0);
}
