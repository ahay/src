/* Image-domain waveform tomography (approximate Hessian). */
/*
  Copyright (C) 2013 University of Texas at Austin
  
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
    int n1, n2, ns, nw;
    int is, iw, **pp, i, ip;
    double omega;
    float dw, ow;
    float *ahess, **ahesss, **ahessr;
    sf_complex **f, ***swave, ***rwave;
    sf_complex **stemp, **rtemp;
    sf_file in, out, list, us, ur, wvlt;
    int uts, mts;

    sf_init(argc,argv);
    in  = sf_input("in");
    out = sf_output("out");

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
    
    /* read wavelet */
    if (NULL == sf_getstring("wvlt"))
	sf_error("Need wvlt=");
    wvlt = sf_input("wvlt");
    
    f = sf_complexalloc2(nw,ns);
    sf_complexread(f[0],nw*ns,wvlt);
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
    stemp = sf_complexalloc2(ns,ns);
    rtemp = sf_complexalloc2(ns,ns);

    ahesss = sf_floatalloc2(n1*n2,ns);
    ahessr = sf_floatalloc2(n1*n2,ns);

    ahess = sf_floatalloc(n1*n2);

    /* loop over frequency */
    for (iw=0; iw < nw; iw++) {
	omega = (double) 2.*SF_PI*(ow+iw*dw);	

	/* read wavefields */
	sf_complexread(swave[0][0],n1*n2*ns,us);
	sf_complexread(rwave[0][0],n1*n2*ns,ur);

#ifdef _OPENMP
#pragma omp parallel num_threads(uts) private(is,ip,i)
#endif
	{
#ifdef _OPENMP
#pragma omp for
#endif
	    for (is=0; is < ns; is++) {
		for (ip=0; ip < ns; ip++) {
		    /* temps */
		    stemp[is][ip] = -omega*omega/conjf(f[ip][iw])
			*rwave[is][pp[ip][1]][pp[ip][0]];
		    
		    /* tempr */
		    rtemp[is][ip] = -omega*omega/conjf(f[ip][iw])
			*conjf(swave[is][pp[ip][1]][pp[ip][0]]);
		}
	    }

	    /* loop over model */
#ifdef _OPENMP
#pragma omp for
#endif
	    for (i=0; i < n1*n2; i++) {
		for (is=0; is < ns; is++) {
		    for (ip=0; ip < ns; ip++) {
			ahesss[ip][i] += crealf(
			    conjf(swave[ip][0][i]*swave[is][0][i])*stemp[is][ip]);
			ahessr[ip][i] += crealf(
			    conjf(swave[ip][0][i])*rwave[is][0][i]*rtemp[is][ip]);
		    }
		}
	    }
	}	
    }

    /* assemble */
#ifdef _OPENMP
#pragma omp parallel for num_threads(uts) private(i,ip)
#endif
    for (i=0; i < n1*n2; i++) {
	for (ip=0; ip < ns; ip++) {
	    ahess[i] += powf(ahesss[ip][i]+ahessr[ip][i],2.);
	}
    }

    /* output hessian */
    sf_floatwrite(ahess,n1*n2,out);

    exit(0);
}
