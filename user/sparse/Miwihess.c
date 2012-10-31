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
    int n1, n2, nh, ns, nw;
    int ih, is, iw, **pos, i, j, ipos;
    double omega;
    float dw, ow, ***wght, ***ahess;
    sf_complex ****swave, ****rwave, *f;
    sf_complex **green, temps, tempr;
    sf_file in, out, list, us, ur, wvlt;
    sf_file weight;
    int uts, mts, its;

    sf_init(argc,argv);
    in  = sf_input("in");
    out = sf_output("out");

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
    
    swave = sf_complexalloc4(n1,n2,ns,nw);
    sf_complexread(swave[0][0][0],n1*n2*ns*nw,us);
    sf_fileclose(us);

    /* read receiver wavefield */	
    if (NULL == sf_getstring("ur"))
	sf_error("Need receiver wavefield ur=");
    ur = sf_input("ur");

    rwave = sf_complexalloc4(n1,n2,ns,nw);
    sf_complexread(rwave[0][0][0],n1*n2*ns*nw,ur);
    sf_fileclose(ur);

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

    pos = sf_intalloc2(2,ns);
    sf_intread(pos[0],2*ns,list);
    sf_fileclose(list);

    /* allocate memory */
    green = sf_complexalloc2(ns,uts);
    ahess = sf_floatalloc3(n1,n2,uts);

    /* loop over frequency */
#ifdef _OPENMP
#pragma omp parallel num_threads(uts) private(its,omega,j,i,ipos,ih,is,temps,tempr)
#endif
    {
#ifdef _OPENMP
	its = omp_get_thread_num();
#else
	its = 0;
#endif

#ifdef _OPENMP
#pragma omp for
#endif
	for (iw=0; iw < nw; iw++) {
	    omega = (double) 2.*SF_PI*(ow+iw*dw);	    
	    
	    /* loop over model */
	    for (j=0; j < n2; j++) {
		for (i=0; i < n1; i++) {
		    
		    /* green's function */
		    for (ipos=0; ipos < ns; ipos++) {
			green[its][ipos] = swave[iw][ipos][j][i]/f[iw];
		    }

		    /* loop over shots */
		    for (is=0; is < ns; is++) {
			/* loop over space-lags */
			for (ih=-nh; ih < nh+1; ih++) {
			    /* loop over green's functions */
			    for (ipos=0; ipos < ns; ipos++) {

				if (pos[ipos][1]+2*ih >= 0 && pos[ipos][1]+2*ih < n2) {
				    temps = conjf(swave[iw][is][j][i]*green[its][ipos])
					*rwave[iw][is][pos[ipos][1]+2*ih][pos[ipos][0]];
				} else {
				    temps = sf_cmplx(0.,0.);
				}

				if (pos[ipos][1]-2*ih >= 0 && pos[ipos][1]-2*ih < n2) {
				    tempr = conjf(swave[iw][is][pos[ipos][1]-2*ih][pos[ipos][0]])
					*rwave[iw][is][j][i]*conjf(green[its][ipos]);
				} else {
				    tempr = sf_cmplx(0.,0.);
				}

				if (wght != NULL) {
				    temps *= wght[ih+nh][pos[ipos][1]+ih][pos[ipos][0]];
				    tempr *= wght[ih+nh][pos[ipos][1]-ih][pos[ipos][0]];
				}
				
				ahess[its][j][i] += pow(-omega*omega*creal(temps+tempr),2.);
			    }
			}
		    }
		}
	    }
	}

#ifdef _OPENMP
#pragma omp for
#endif
	for (j=0; j < n2; j++) {
	    for (i=0; i < n1; i++) {
		for (its=1; its < uts; its++) {
		    ahess[0][j][i] += ahess[its][j][i];
		}
	    }
	}
    }
    
    /* write output */
    sf_floatwrite(ahess[0][0],n1*n2,out);

    exit(0);
}
