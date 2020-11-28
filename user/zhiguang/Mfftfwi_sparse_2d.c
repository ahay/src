/* 2D frequency domain full waveform inversion with sparsity regularization. */
/*
  Copyright (C) 2014 University of Texas at Austin
  
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

#include "fwisparse2d.h"

int main(int argc, char* argv[])
{
 
    bool hermite_false, hermite_true, sparsity;
    int n1, n2, npml, pad1, pad2, ns, nw;
    float d1, d2, ds, os, dw, ow, pclip;
    double omega;
    sf_file in, out, vout=NULL, misfit, source, receiver, record, dip;
    char *order, *type = NULL;
    int uts, mts, iw, niter, par;
    float **v, **vnew, **slope=NULL, **recloc, *error;
    sf_complex ***f, ***obs;

    sf_init(argc, argv);

    in = sf_input("in");
    out = sf_output("out");

	if (!sf_getbool("sparsity", &sparsity)) sparsity=true;
	/* if true, sparsity constriant; if false, normal FWI */
    if (!sf_getint("niter", &niter)) niter=10;
	/* number of iteration */
	if (!sf_getint("npml", &npml)) npml=20;
	/* PML width */
	if (NULL == (order = sf_getstring("order"))) order="j";
	/* discretization scheme (default optimal 9-point) */
	if (sparsity && (NULL == (type = sf_getstring("type")))) type="b";
	/* [haar,linear,biorthogonal] wavelet type, the default is biorthogonal */
	if (sparsity && !sf_getint("par", &par)) par=1;
	/* seislet transform accuracy order */
	if (sparsity && !sf_getfloat("pclip", &pclip)) pclip=8.;
	/* soft thresholding parameter */

	if (!sf_getint("uts",&uts)) uts=0;
    mts = 1;
    uts = (uts < 1)? mts: uts;
	/* number of threads */

    hermite_false=false;
    hermite_true=true;
    /* Hermite operator */

    /* read input dimension */
    if (!sf_histint(in,"n1",&n1)) sf_error("No n1= in input.");
    if (!sf_histint(in,"n2",&n2)) sf_error("No n2= in input.");

    if (!sf_histfloat(in,"d1",&d1)) sf_error("No d1= in input.");
    if (!sf_histfloat(in,"d2",&d2)) sf_error("No d2= in input.");

	/* initial velocity */
    v = sf_floatalloc2(n1,n2);
    sf_floatread(v[0],n1*n2,in);
    
	/* PML padding */
	pad1 = n1+2*npml;
	pad2 = n2+2*npml;

    /* read receiver */ 
    if (NULL == sf_getstring("receiver")) sf_error("Need receiver="); 
    receiver = sf_input("receiver"); 
    recloc = sf_floatalloc2(n1,n2);
    sf_floatread(recloc[0],n1*n2,receiver);

    /* read source */
    if (NULL == sf_getstring("source")) sf_error("Need source=");
    source = sf_input("source");

    if (!sf_histint(source,"n3",&ns)) sf_error("No ns=.");
    if (!sf_histfloat(source,"d3",&ds)) ds=d2;
    if (!sf_histfloat(source,"o3",&os)) os=0.;

	/* read data */
    if (NULL == sf_getstring("record")) sf_error("Need record=");
    record = sf_input("record");

    if (!sf_histint(record,"n4",&nw)) sf_error("No nw=.");
    if (!sf_histfloat(record,"d4",&dw)) sf_error("No dw=.");
    if (!sf_histfloat(record,"o4",&ow)) sf_error("No ow=."); 

    /* input slope */
    if (sparsity) {
        if (NULL == sf_getstring("dip")) sf_error("Need dip=");
        dip = sf_input("dip");
        slope = sf_floatalloc2(n1,n2);
        sf_floatread(slope[0],n1*n2,dip);
    }

	/* output the last iteration */
	if (NULL != sf_getstring("vout")){
		vout=sf_output("vout");
	}

	/* misfit */
	if (NULL == sf_getstring("misfit")) sf_error("Need misfit=");
	misfit=sf_output("misfit");
	error=sf_floatalloc(niter*nw);

	/* set up output file dimension */
    sf_putint(out,"n1",n1);
    sf_putint(out,"n2",n2);
    sf_putint(out,"n3",niter*nw);
	sf_putint(misfit, "n1", niter*nw);
	sf_putint(misfit, "n2", 1);

	/* storage allocation */
    obs = sf_complexalloc3(n1,n2,ns);
    f = sf_complexalloc3(n1,n2,ns);
	vnew = sf_floatalloc2(n1,n2);

	/* parameter passing */
	fwisparse_init(out, v, vnew, slope, recloc, order, niter, uts,
			hermite_false, hermite_true, sparsity, error, type, 
			n1, n2, npml, pad1, pad2, ns, d1, d2, par, pclip);

    /* Loop over frequency */
    for ( iw = 0; iw < nw; iw ++ ) { 
        omega=(double) 2.*SF_PI*(ow+iw*dw); 
            
        sf_warning("Calculating frequency %d out of %d for %f HZ.\n", iw+1,nw,ow+iw*dw);

        sf_complexread(f[0][0],n1*n2*ns,source);
        sf_complexread(obs[0][0],n1*n2*ns,record);

		/* program execution for one frequency slice */
		fwisparse_exec(omega, f, obs);

		sf_warning("Ending frequency %d out of %d for %f HZ. \n\n", iw+1, nw, ow+iw*dw);
    } /* end frequency */

	if(vout != NULL) sf_floatwrite(v[0], n1*n2, vout);
	sf_floatwrite(error, niter*nw, misfit);

	fwisparse_free();
	sf_fileclose(out);
	free(**f); free(*f); free(f);
	free(**obs); free(*obs); free(obs);
       
    exit(0);
}
