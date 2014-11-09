/* 2D Frequency Domain Reverse Time Migration. */
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
#include <umfpack.h>

#ifdef _OPENMP
#include <omp.h>
#endif

#include "fdprep.h"
#include "sparsesolver.h"
#include "optimization.h" 

int main(int argc, char* argv[])
{
 
    bool hermite_false, hermite_true;
    int n1, n2, npml, pad1, pad2, ns, nw, nh;
    float d1, d2, **v, ds, os, dw, ow;
    double omega;
    sf_complex ***f, ***srcw, ***recw, ***obs, ***obs_cut;
    sf_file in, out, source, receiver, record;
    int uts, mts;
    char *order;
    int is, i, j, iw, ih;
    float ***image, **recloc;

    sf_init(argc, argv);

    in = sf_input("in");
    out = sf_output("out");

    if (!sf_getint("nh",&nh)) nh=0;

    if (!sf_getint("uts",&uts)) uts=0;
//#ifdef _OPENMP
//    mts = omp_get_max_threads();
//#else
    mts = 1;
//#endif
    uts = (uts < 1)? mts: uts;

    hermite_false=false;
    hermite_true=true;
    /* Hermite operator */
    
    if (!sf_getint("npml",&npml)) npml=20;
    /* PML width */

    if (NULL == (order = sf_getstring("order"))) order="c";
    /* discretization scheme (default optimal 9-point) */

    fdprep_order(order);

    /* read input dimension */
    if (!sf_histint(in,"n1",&n1)) sf_error("No n1= in input.");
    if (!sf_histint(in,"n2",&n2)) sf_error("No n2= in input.");

    if (!sf_histfloat(in,"d1",&d1)) sf_error("No d1= in input.");
    if (!sf_histfloat(in,"d2",&d2)) sf_error("No d2= in input.");

    v = sf_floatalloc2(n1,n2);
    sf_floatread(v[0],n1*n2,in);
    
	/* PML padding */
	pad1 = n1+2*npml;
	pad2 = n2+2*npml;

    /* read receiver */ 
    if (NULL == sf_getstring("receiver")) sf_error("Need receiver="); 
    receiver = sf_input("receiver"); 
    recloc=sf_floatalloc2(n1,n2);
    sf_floatread(recloc[0],n1*n2,receiver);

    /* read source */
    if (NULL == sf_getstring("source")) sf_error("Need source=");
    source = sf_input("source");

    if (!sf_histint(source,"n3",&ns)) sf_error("No ns=.");
    if (!sf_histfloat(source,"d3",&ds)) ds=d2;
    if (!sf_histfloat(source,"o3",&os)) os=0.;

    f = sf_complexalloc3(n1,n2,ns);

    /* read observed data */
    if (NULL == sf_getstring("record")) sf_error("Need record=");
    record = sf_input("record");

    if (!sf_histint(record,"n4",&nw)) sf_error("No nw=.");
    if (!sf_histfloat(record,"d4",&dw)) sf_error("No dw=.");
    if (!sf_histfloat(record,"o4",&ow)) sf_error("No ow=."); 

    obs = sf_complexalloc3(n1,n2,ns);
    obs_cut = sf_complexalloc3(n1,n2,ns);

    srcw = sf_complexalloc3(n1,n2,ns);
    recw = sf_complexalloc3(n1,n2,ns);
    image = sf_floatalloc3(n1,n2,2*nh+1);

    /* Loop over frequency */
    for (iw=0; iw<nw; iw++ ) { 
        omega=(double) 2.*SF_PI*(ow+iw*dw); 
            
        sf_warning("Calculating frequency %d out of %d for %f HZ.",iw+1,nw,ow+iw*dw);

    	sf_complexread(f[0][0],n1*n2*ns,source);
        sf_complexread(obs[0][0],n1*n2*ns,record);
        
        /* generate adjoint source for reverse time migration */
        genadjsrc_rtm(obs, obs_cut, recloc, n1, n2, ns); 

        /* initialize sparse solver */ 
        sparse_init(uts, pad1, pad2);

        /* factorize matrix, change according to different frequencies and models */
        sparse_factor(omega,n1,n2,d1,d2,v,npml,pad1,pad2);

        for (is=0; is < ns; is++ ) { 
        for (j=0; j < n2; j++ ) { 
        for (i=0; i < n1; i++ ) { 
            srcw[is][j][i]=f[is][j][i];
            recw[is][j][i]=obs_cut[is][j][i];
        }
        }
        }

        /* sparse solver for source wavefield */ 
        sparse_solve(npml, pad1, pad2, srcw, hermite_false, ns, uts);

        /* sparse solver for receiver wavefield */
        sparse_solve(npml, pad1, pad2, recw, hermite_true, ns, uts);

        /* imaging condition */
        for (ih=-nh; ih < nh+1; ih++ ) {
        for (j=0; j<n2; j++ ) {
        for (i=0; i< n1; i++ ) {
        for (is=0; is < ns; is++ ) {
        if (j-abs(ih) >= 0 && j+abs(ih) < n2) {
            image[ih+nh][j][i] += crealf(omega*omega*conjf(srcw[is][j-ih][i])*recw[is][j+ih][i]/(v[j][i]*v[j][i]));
        }
        }
        }
        }
        }

        /* free memory */
        sparse_free(uts);
    } /* end frequency */

    sf_putint(out,"n1",n1);
    sf_putint(out,"n2",n2);
    sf_putint(out,"n3",2*nh+1); 
    sf_putfloat(out,"d3",d2);
    sf_putfloat(out,"o3", (float) -nh*d2);

    sf_floatwrite(image[0][0],n1*n2*(2*nh+1),out);

    exit(0);
}
