/* 2D Helmholtz forward solver by LU factorization. */
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


int main(int argc, char* argv[])
{
 
    bool hermite;
    int n1, n2, npml, pad1, pad2, ns, nw, iw;
    float d1, d2, **v, ds, os, ow, dw;
    double omega;
    sf_complex ***f;
    sf_file in, out, source;
    int uts, mts;
    char *order;
    
    sf_init(argc, argv);

    in = sf_input("in");
    out = sf_output("out");

    if (!sf_getint("uts",&uts)) uts=0;
//#ifdef _OPENMP
//    mts = omp_get_max_threads();
//#else
    mts = 1;
//#endif
    uts = (uts < 1)? mts: uts;

    if (!sf_getbool("hermite",&hermite)) hermite=false;
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

    /* read source */
    if (NULL == sf_getstring("source"))
	sf_error("Need source=");
    source = sf_input("source");

    if (!sf_histint(source,"n3",&ns)) sf_error("No ns=.");
    if (!sf_histfloat(source,"d3",&ds)) ds=d2;
    if (!sf_histfloat(source,"o3",&os)) os=0.;
    if (!sf_histint(source,"n4",&nw)) sf_error("No nw=.");
    if (!sf_histfloat(source,"d4",&dw)) sf_error("No dw=.");
    if (!sf_histfloat(source,"o4",&ow)) sf_error("No ow=.");

    f = sf_complexalloc3(n1,n2,ns);

    /* write out forward simulation */
    sf_settype(out,SF_COMPLEX);

    sf_putint(out,"n3",ns);
    sf_putfloat(out,"d3",ds);
    sf_putfloat(out,"o3",os);
    sf_putstring(out,"label3","Shot");
    sf_putstring(out,"unit3","");
    sf_putint(out,"n4",nw);
    sf_putfloat(out,"d4",dw);
    sf_putfloat(out,"o4",ow);
    sf_putstring(out,"label4","Frequency");
    sf_putstring(out,"unit4","Hz");

    /* Loop over frequency */
    for (iw=0; iw<nw; iw++ ) { 
        omega=(double) 2.*SF_PI*(ow+iw*dw); 

        sf_warning("Calculating frequency %d out of %d for %f HZ.",iw+1,nw,ow+iw*dw);

        /* read in source */
        sf_complexread(f[0][0],n1*n2*ns,source);

        /* initialize sparse solver */ 
        sparse_init(uts, pad1, pad2);

        /* factorize matrix, change according to different frequencies and models */
        sparse_factor(omega, n1, n2, d1, d2, v, npml, pad1, pad2);
    
        /* sparse solver */
        sparse_solve(npml, pad1, pad2, f, hermite, ns, uts);

        /* write out wavefield */
        sf_complexwrite(f[0][0],n1*n2*ns,out);

        /* free memory */
        sparse_free(uts);
    }

    exit(0);
}
