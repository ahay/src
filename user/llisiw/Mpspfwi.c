/* Full-waveform Inversion by Parallel Helmholtz Solver */
/*
  Copyright (C) 2011 University of Texas at Austin
  
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

#include "helm.h"
#include "pspfwi.h"

int main(int argc, char* argv[])
{
    int n[SF_MAX_DIM], dim, i, nt, it, is, ns, iw, nw;
    int shift, **mask, water, iter, niter, cgiter, istep, nstep;
    float o[SF_MAX_DIM], d[SF_MAX_DIM], dw, ow, *s, *temps, step, misnorm, misnorm0, misnorm1;
    sf_complex **source, **data, *ds, *dp;
    char key[4], *what;
    sf_file in, out, shot, recv, reco;

    sf_init(argc,argv);
    in  = sf_input("in");
    out = sf_output("out");

    if (NULL == (what = sf_getstring("what"))) what="inversion";
    /* what to compute (default inversion) */

    /* read model dimensions [n0*n1*n2] */
    dim = sf_filedims(in,n);

    nt = 1;
    for (i=0; i < dim; i++) {
	sprintf(key,"d%d",i+1);
	if (!sf_histfloat(in,key,d+i)) sf_error("No %s= in input.",key);
	sprintf(key,"o%d",i+1);
	if (!sf_histfloat(in,key,o+i)) o[i]=0.;
	nt *= n[i];
    }

    if (dim < 3) {
	n[2] = 1; o[2] = o[1]; d[2] = d[1];
    }

    /* read initial slowness */
    s = sf_floatalloc(nt);
    sf_floatread(s,nt,in);

    /* read receiver [n1*n2][ns] */
    mask = sf_intalloc2(n[1]*n[2],ns);
    
    if (NULL == sf_getstring("recv")) {
	recv = NULL;

	/* receiver everywhere */
	for (is=0; is < ns; is++) {
	    for (it=0; it < n[1]*n[2]; it++) {
		mask[is][it] = 1;
	    }
	}
    } else {
	recv = sf_input("recv");

	if (SF_INT != sf_gettype(recv)) sf_error("Receiver must be integer.");

	sf_intread(mask[0],n[1]*n[2]*ns,recv);
	sf_fileclose(recv);
    }

    if (!sf_getint("shift",shift)) shift=0;
    /* receiver shift */

    /* read shot file [n0*n1*n2][ns][nw] */
    if (NULL == sf_getstring("shot")) sf_error("Need source.");   
    shot = sf_input("shot");

    if (SF_COMPLEX != sf_gettype(shot)) sf_error("Shot must be complex.");

    if (!sf_histint(shot,"n4",&ns)) ns=1;

    /* NOTE: shot frequency by frequency */
    source = sf_complexalloc2(nt,ns);

    /* NOTE: data frequency by frequency */
    data = sf_complexalloc2(n[1]*n[2],ns);
    
    /* temporary array */
    ds = sf_complexalloc(nt);

    switch (what[0]) {
	case 'f': /* foward modeling */

	    reco = NULL;
	    dp = NULL;

	    if (!sf_getint("nw",nw)) nw=1;
	    /* number of frequency */
	    if (!sf_getfloat("dw",dw)) dw=1.;
	    /* frequency increment */
	    if (!sf_getfloat("ow",ow)) ow=5.;
	    /* starting frequency */

	    /* write output dimensions [n1][n2][ns][nw] */
	    sf_settype(out,SF_COMPLEX);
	    
	    sf_putint(out,"n1",n[1]); sf_putfloat(out,"d1",d[1]); sf_putfloat(out,"o1",o[1]);
	    sf_putint(out,"n2",n[2]); sf_putfloat(out,"d2",d[2]); sf_putfloat(out,"o2",o[2]);
	    sf_putint(out,"n3",ns);   sf_putfloat(out,"d3",1.);   sf_putfloat(out,"o3",1.);
	    sf_putint(out,"n4",nw);   sf_putfloat(out,"d4",dw);   sf_putfloat(out,"o4",ow);

	    /* loop over frequency */
	    for (iw=0; iw < nw; iw++) {
		w = ow + iw*dw;

		sf_complexread(source[0],nt*ns,shot);

		/* helmholtz setup */
		helm_setup(s,w);

		/* loop over shot */
		for (is=0; is < ns; is++) {

		    /* helmholtz solve */
		    helm_solve(source[is],ds);

		    /* receiver */
		    for (it=0; it < n[1]*n[2]; it++) {
			if (mask[is][it])
			    data[is][it] = ds[shift+it*n[0]];
			else
			    data[is][it] = sf_cmplx(0.,0.);
		    }
		}

		sf_complexwrite(data[0],n[1]*n[2]*ns,out);
	    }

	    break;
	    
	case 'i': /* inversion */

	    /* read record file [n1*n2*ns][nw] */
	    if (NULL == sf_getstring("reco")) sf_error("Need record.");	    
	    reco = sf_input("reco");
	    
	    if (SF_COMPLEX != sf_gettype(reco)) sf_error("Record must be complex.");
	    
	    if (!sf_histint(reco,"n4",nw)) sf_error("No nw in record.");
	    if (!sf_histfloat(reco,"d4",dw)) sf_error("No dw in record.");
	    if (!sf_histfloat(reco,"o4",ow)) sf_error("No ow in record.");
	    
	    if (!sf_getint("niter",niter)) niter=1;
	    /* number of inversion iterations */

	    if (!sf_getint("cgiter",cgiter)) cgiter=20;
	    /* number of conjugate-gradient iterations */
	    
	    if (!sf_getint("nstep",nstep)) nstep=3;
	    /* number of line-search iterations */

	    if (!sf_getint("water",water)) water=0;
	    /* water layer depth */

	    /* temporary array */
	    dp = sf_complexalloc(n[1]*n[2]*ns);
	    temps = sf_floatalloc(nt);

	    /* initialize */
	    fwi_init(nt,n,ns,shift,mask,water);
	    
	    /* loop over frequency */
	    for (iw=0; iw < nw; iw++) {
		w = ow + iw*dw;

		sf_warning("frequency %g start:",w);
		
		sf_complexread(source[0],nt*ns,shot);
		sf_complexread(data[0],n[1]*n[2]*ns,reco);

		/* fwi setup */
		fwi_setup(s,w);
		
		/* forward */
		fwi_forward(source,data,dp);

		/* initial misfit norm */
		misnorm = cblas_scnrm2(n[1]*n[2]*ns,dp,1);
		misnorm0 = misnorm;

		/* iteration */
		for (iter=0; iter < niter; iter++) {

		    /* compute update */
		    sf_csolver(fwi_operator,sf_ccgstep,nt,n[1]*n[2]*ns,ds,dp,cgiter,"end");
		    sf_ccgstep_close();

		    /* line search */
		    for (istep=0, step=1.; istep < nstep; istep++, step *= 0.5) {

			/* take only real part of ds */
			for (it=0; it < nt; it++) {
			    temps[it] = s[it]+step*crealf(ds[it]);
			}

			/* compare misfit */
			fwi_setup(temps,w);
			fwi_forward(source,data,dp);

			misnorm1 = cblas_scnrm2(n[1]*n[2]*ns,dp,1);

			/* break if misfit norm decreases */
			if (misnorm1 < misnorm0) {
			    for (it=0; it < nt; it++) {
				s[it] = temps[it];
			    }

			    sf_warning("relative misfit %g after iteration %d (line search %d)",misnorm1/misnorm,iter+1,istep);

			    misnorm0 = misnorm1;
			    break;
			}
		    }

		    /* line search failure */
		    if (istep == nstep) {
			sf_warning("line search failed... jump to next frequency...");
			break;
		    }
		}
	    }

	    break;
    }

    exit(0);
}
