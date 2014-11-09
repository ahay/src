/* 2D Frequency Domain Least Squares Reverse Time Migration. */
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
#include "waveoperator.h"

int main(int argc, char* argv[])
{
 
    bool hermite_false, hermite_true;
    int n1, n2, npml, pad1, pad2, ns, nw;
    float d1, d2, **v, ds, os, dw, ow;
    sf_complex ****f,  ****obs; 
    sf_file in, out, misfit, source, receiver, record;
    char *order;
    int uts, mts, is, i, j, iw, iter, niter;
    float **recloc;
    float **m_old, **m_new;
    float **d_new, **d_old;
    float **g_old, **g_new; 
    sf_complex ****r_new, ****r_old;
    sf_complex ****Fg; 
    float alpha, beta, gnorm, rnorm; 
    float *datamisfit;

    sf_init(argc, argv);

    in = sf_input("in");
    out = sf_output("out");
    misfit = sf_output("misfit");

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
    
    if (!sf_getint("niter",&niter)) niter=0; 
    /* Number of iterations */

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

    /* read observed data */
    if (NULL == sf_getstring("record")) sf_error("Need record=");
    record = sf_input("record");

    if (!sf_histint(record,"n4",&nw)) sf_error("No nw=.");
    if (!sf_histfloat(record,"d4",&dw)) sf_error("No dw=.");
    if (!sf_histfloat(record,"o4",&ow)) sf_error("No ow=."); 

    f = sf_complexalloc4(n1,n2,ns,nw);
    obs = sf_complexalloc4(n1,n2,ns,nw);

    sf_complexread(f[0][0][0],n1*n2*ns*nw,source);
    sf_complexread(obs[0][0][0],n1*n2*ns*nw,record);
   
    /* allocate variables */
    m_old = sf_floatalloc2(n1,n2);
    m_new = sf_floatalloc2(n1,n2);
    d_old = sf_floatalloc2(n1,n2);
    d_new = sf_floatalloc2(n1,n2);
    g_old = sf_floatalloc2(n1,n2);
    g_new = sf_floatalloc2(n1,n2);
    
    r_old = sf_complexalloc4(n1,n2,ns,nw);
    r_new = sf_complexalloc4(n1,n2,ns,nw);
    Fg = sf_complexalloc4(n1,n2,ns,nw);

    /* set output dimension */
    sf_putint(out,"n1",n1);
    sf_putint(out,"n2",n2);
    sf_putint(out,"n3",niter);
    
    sf_putint(misfit,"n1",niter);
    sf_putfloat(misfit,"d1",1);
    sf_putint(misfit,"n2",1);
    datamisfit = sf_floatalloc(niter);

    rnorm = 0.0; 
    for ( iw = 0; iw < nw; iw ++ ) { 
    for ( is = 0; is < ns; is ++) { 
    for ( j = 0; j < n2 ; j++ ) { 
    for ( i = 0; i < n1; i++ ) { 
        r_old[iw][is][j][i] = obs[iw][is][j][i];
        if ( recloc[j][i] > 0.0 ) { 
            rnorm += cabsf( r_old[iw][is][j][i] * r_old[iw][is][j][i] ); 
        }
    }
    }
    }
    }
    sf_warning("rnorm = %g.",rnorm);


    sf_warning("Adjoint calculation for the first iteration.");
    /* adjoint calculation */
    adjlsm_operator(nw, ow, dw, ns, n1, n2, 
                     uts, pad1, pad2, npml, d1, d2,
                     hermite_false, hermite_true, v, f, recloc, 
                     r_old, g_old);

    /* set starting valuables */
    for (j = 0; j < n2; j++) { 
    for (i = 0; i < n2; i++ ) { 
        d_old[j][i] = g_old[j][i];
        m_old[j][i] = 0.0;
    }
    }
    
    for ( iter = 0; iter < niter; iter ++ ) { 

        sf_warning("Calculating iteration %d out of %d.",iter,niter);

        /* born forward operator */
        bornsyn_operator(nw, ow, dw, ns, n1, n2, 
                         uts, pad1, pad2, npml, d1, d2,
                         hermite_false, v, f, recloc, 
                         d_old, Fg); 

        /* calculate alpha value */
        alpha = calc_alpha(g_old, Fg, recloc, n1, n2, ns, nw);
        
        /* update model */ 
        update_model_lsm(m_old, m_new, d_old, alpha, n1, n2); 

        /* update residual */ 
        update_residual(r_old, r_new, Fg, recloc, alpha, n1, n2, ns, nw);  

        /* adjoint operator */ 
        adjlsm_operator(nw, ow, dw, ns, n1, n2, 
                         uts, pad1, pad2, npml, d1, d2,
                         hermite_false, hermite_true, v, f, recloc, 
                         r_new, g_new);

        /* update direction */
        beta = direction_cg_fletcher(g_old, d_old, g_new, d_new, n1, n2); 

        sf_warning("alpha = %g, beta = %g.",alpha, beta);

        /* update vectors */
        gnorm = 0.0 ;  
        for (j = 0; j < n2; j++ ) { 
        for (i = 0; i < n1; i++ ) { 
            d_old[j][i] = d_new[j][i]; 
            g_old[j][i] = g_new[j][i];
            m_old[j][i] = m_new[j][i]; 
            gnorm += g_old[j][i] * g_old[j][i]; 
        }
        }

        rnorm = 0.0 ; 
        for (iw = 0; iw < nw; iw ++ ) { 
        for (is = 0; is < ns; is ++ ) { 
        for (j = 0; j < n2; j ++ ) { 
        for (i = 0; i < n1; i ++ ) { 
            r_old[iw][is][j][i] = r_new[iw][is][j][i]; 
            if ( recloc[j][i] > 0.0 ) { 
                rnorm += cabsf( r_old[iw][is][j][i] * r_old[iw][is][j][i] ); 
            }
        }
        }
        }
        }
        sf_warning("gnorm = %g; rnorm = %g.",gnorm, rnorm);

        datamisfit[iter] = rnorm; 

        sf_floatwrite(m_old[0],n1*n2,out);

    }  /* end iteration */

    sf_floatwrite(datamisfit, niter, misfit); 

    exit(0);
}
