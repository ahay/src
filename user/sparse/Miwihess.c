/* Image-domain waveform tomography (approximate Hessian). */
/*
  Copyright (C) 2012 University of Texas at Austin
  
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

#ifdef _OPENMP
#include <omp.h>
#endif

int main(int argc, char* argv[])
{
    bool verb;
    int n1, n2, nh, ns, nw;
    int ih, is, iw, **pp, i, ip;
    int zz, xx;
    double omega;
    float dw, ow, ***wght;
    float *ahess, ***ahesss, ***ahessr;
    sf_complex *f, ***swave, ***rwave;
    sf_complex ***stemp, ***rtemp;
    sf_file in, out, list, us, ur, wvlt;
    sf_file weight;
    int uts, mts, its;
    sf_timer timer;

    sf_init(argc,argv);
    in  = sf_input("in");
    out = sf_output("out");

    if (!sf_getbool("verb",&verb)) verb=false;
    /* verbosity flag */

    if (verb)
	timer = sf_timer_init();
    else
	timer = NULL;

    if (!sf_getint("nh",&nh)) nh=0;
    /* horizontal space-lag */

    if (!sf_getint("uts",&uts)) uts=0;
    /* number of OMP threads */

#ifdef _OPENMP
    mts = omp_get_max_threads();
#else
    mts = 1;
#endif

    uts = (uts < 1)? mts: uts;

    /* read model dimensions */
    if (!sf_histint(in,"n1",&n1)) sf_error("No n1= in input.");
    if (!sf_histint(in,"n2",&n2)) sf_error("No n2= in input.");

    /* read source wavefield */
    if (NULL == sf_getstring("us"))
	sf_error("Need source wavefield us=");
    us = sf_input("us");

    if (!sf_histint(us,"n3",&ns)) sf_error("No ns=.");
    if (!sf_histint(us,"n4",&nw)) sf_error("No nw=.");
    if (!sf_histfloat(us,"d4",&dw)) sf_error("No dw=.");
    if (!sf_histfloat(us,"o4",&ow)) sf_error("No ow=.");

    /* read receiver wavefield */	
    if (NULL == sf_getstring("ur"))
	sf_error("Need receiver wavefield ur=");
    ur = sf_input("ur");
    
    /* read image weight */
    if (NULL == sf_getstring("weight")) {
	weight = NULL;
	wght = NULL;
    } else {
	weight = sf_input("weight");
	wght = sf_floatalloc3(n1,n2,2*nh+1);
	sf_floatread(wght[0][0],n1*n2*(2*nh+1),weight);
	sf_fileclose(weight);
    }

    /* read wavelet */
    if (NULL == sf_getstring("wvlt"))
	sf_error("Need wvlt=");
    wvlt = sf_input("wvlt");
    
    f = sf_complexalloc(nw);
    sf_complexread(f,nw,wvlt);
    sf_fileclose(wvlt);

    /* read list */
    if (NULL == sf_getstring("list"))
	sf_error("Need list=");
    list = sf_input("list");

    pp = sf_intalloc2(2,ns);
    sf_intread(pp[0],2*ns,list);
    sf_fileclose(list);

    /* allocate memory */
    swave = sf_complexalloc3(n1,n2,ns);
    rwave = sf_complexalloc3(n1,n2,ns);
    stemp = sf_complexalloc3(2*nh+1,ns,ns);
    rtemp = sf_complexalloc3(2*nh+1,ns,ns);

    ahess = sf_floatalloc(n1*n2);
    ahesss = sf_floatalloc3(2*nh+1,ns,n1*n2); /* NOTE: we cannot afford this! */
    ahessr = sf_floatalloc3(2*nh+1,ns,n1*n2); /* NOTE: we cannot afford this! */

    /* loop over frequency */
    for (iw=0; iw < nw; iw++) {
	omega = (double) 2.*SF_PI*(ow+iw*dw);

	if (verb) {
	    sf_warning("Frequency %d of %d.",iw+1,nw);
	    sf_timer_start(timer);
	}

	/* read wavefields */
	sf_complexread(swave[0][0],n1*n2*ns,us);
	sf_complexread(rwave[0][0],n1*n2*ns,ur);

#ifdef _OPENMP
#pragma omp parallel num_threads(uts) private(its,is,ip,ih,zz,xx,i)
#endif
	{
#ifdef _OPENMP
	    its = omp_get_thread_num();
#else
	    its = 0;
#endif

	    /* temparary wavefields */
	    /* NOTE: negletible cost */
#ifdef _OPENMP
#pragma omp for
#endif
	    for (is=0; is < ns; is++) {
		for (ip=0; ip < ns; ip++) {
		    for (ih=-nh; ih < nh+1; ih++) {
			zz = pp[ip][0];

			/* source */
			xx = pp[ip][1]+ih;
			if (xx+ih >= 0 && xx+ih < n2)
			    stemp[is][ip][ih+nh] = -omega*omega/conjf(f[iw])
				*(wght==NULL? 1.: wght[ih+nh][xx][zz])*rwave[is][xx+ih][zz];
			
			/* receiver */
			xx = pp[ip][1]-ih;
			if (xx-ih >= 0 && xx-ih < n2)
			    rtemp[is][ip][ih+nh] = -omega*omega/conjf(f[iw])
				*(wght==NULL? 1.: wght[ih+nh][xx][zz])*conjf(swave[is][xx-ih][zz]);
		    }
		}
	    }

	    /* loop over model */
	    /* NOTE: unacceptable cost */
#ifdef _OPENMP
#pragma omp for
#endif
	    for (i=0; i < n1*n2; i++) {
		for (ih=-nh; ih < nh+1; ih++) {
		    for (is=0; is < ns; is++) {
			for (ip=0; ip < ns; ip++) {
			    ahesss[i][ip][ih+nh] += 
				crealf(conjf(swave[ip][0][i]*swave[is][0][i])*stemp[is][ip][ih+nh]);
			    ahessr[i][ip][ih+nh] += 
				crealf(conjf(swave[ip][0][i])*rwave[is][0][i]*rtemp[is][ip][ih+nh]);
			}
		    }
		}
	    }
	}

	if (verb) {
	    sf_timer_stop (timer);
	    sf_warning("Finished in %g seconds.",sf_timer_get_diff_time(timer)/1.e3);
	}
    }

    /* assemble approximate Hessian */
    /* NOTE: to be finished */

    exit(0);
}
