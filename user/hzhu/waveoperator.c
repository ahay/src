/* Subroutines for waveoperators */
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

#include "sparsesolver.h" 
#include "optimization.h"

/* subroutines for fwi */ 
float forward_operator(const int uts, 
                      const int pad1,
                      const int pad2, 
                      const double omega, 
                      const int n1,
                      const int n2, 
                      const float d1, 
                      const float d2, 
                      const int npml,
                      const int ns,
                      sf_complex ***f, 
                      sf_complex ***obs,
                      bool hermite_false, 
                      float **recloc,
                      float **input) 
/*< forward operator for fwi for a single frequency >*/ 
{

    int is, i, j;
    sf_complex ***syn;
    float misfit; 

    syn = sf_complexalloc3(n1, n2, ns);

    /* initialize sparse solver */ 
    sparse_init(uts, pad1, pad2);
            
    /* factorize matrix, change according to different frequencies and models */
    sparse_factor(omega, n1, n2, d1, d2, input, npml, pad1, pad2);
    
    for ( is = 0; is < ns; is++ ) { 
    for ( j = 0; j < n2; j++ ) { 
    for ( i = 0; i < n1; i++ ) { 
        syn[is][j][i]=f[is][j][i];
    }
    }
    }
    /* sparse solver for source wavefield */ 
    sparse_solve(npml, pad1, pad2, syn, hermite_false, ns, uts);

    misfit = calc_misfit(obs, syn, recloc, n1, n2, ns);

    /* free memory */
    sparse_free(uts);

    return misfit; 
}

float adjfwi_operator(const int uts, 
                      const int pad1, 
                      const int pad2,
                      const double omega, 
                      const int n1,
                      const int n2,
                      const float d1,
                      const float d2,
                      const int npml,
                      const int ns, 
                      sf_complex ***f,
                      sf_complex ***obs,
                      bool hermite_false, 
                      bool hermite_true,
                      float **recloc, 
                      float **input,
                      float **output) 
/*< adjoint operator for fwi for a single frequency >*/
{

    int is, i, j;
    sf_complex ***syn, ***adj; 
    float misfit, weight; 

    syn = sf_complexalloc3(n1, n2, ns);
    adj = sf_complexalloc3(n1, n2, ns);

    /* initialize sparse solver */ 
    sparse_init(uts, pad1, pad2);
            
    /* factorize matrix, change according to different frequencies and models */
    sparse_factor(omega, n1, n2, d1, d2, input, npml, pad1, pad2);

    for (is=0; is < ns; is++ ) { 
    for (j=0; j < n2; j++ ) { 
    for (i=0; i < n1; i++ ) { 
        syn[is][j][i]=f[is][j][i];
    }
    }
    }
    
    /* sparse solver for source wavefield */ 
    sparse_solve(npml, pad1, pad2, syn, hermite_false, ns, uts);
            
    misfit = genadjsrc_fwi(obs, syn, adj, recloc, n1, n2, ns);
            
    /* sparse solver for receiver wavefield */
    sparse_solve(npml, pad1, pad2, adj, hermite_true, ns, uts);

    /* calculate gradient */
    for (j=0; j<n2; j++ ) {
    for (i=0; i< n1; i++ ) {
        output[j][i]=0.0;
        weight=-1.0;
        for (is=0; is < ns; is++ ) {
            output[j][i] += weight*crealf(conjf(syn[is][j][i])*adj[is][j][i]);
        }
    }
    }
 
    /* free memory */
    sparse_free(uts);

    return misfit; 
}


/* subroutines for lsm */ 
void bornsyn_operator(const int nw,
                      const float ow, 
                      const float dw, 
                      const int ns, 
                      const int n1,
                      const int n2, 
                      const int uts, 
                      const int pad1, 
                      const int pad2, 
                      const int npml,
                      const float d1, 
                      const float d2,
                      bool hermite_false, 
                      float **v, 
                      sf_complex ****f, 
                      float **recloc, 
                      float **input, 
                      sf_complex ****output)
/*< born forward operator >*/
{        

    int iw, is, j, i; 
    double omega; 
    sf_complex ***fin, ***fborn; 

    fin = sf_complexalloc3(n1,n2,ns);
    fborn = sf_complexalloc3(n1,n2,ns);

    /* Loop over frequency */
    for ( iw = 0; iw < nw; iw++ ) { 
        omega=(double) 2. * SF_PI * ( ow + iw * dw ); 

        /* initialize sparse solver */ 
        sparse_init(uts, pad1, pad2);

        /* factorize matrix, change according to different frequencies and models */
        sparse_factor(omega, n1, n2, d1, d2, v, npml, pad1, pad2);

        for ( is = 0; is < ns; is ++ ) { 
        for ( j = 0; j < n2; j ++) { 
        for ( i = 0; i < n1; i++ ) { 
            fin[is][j][i]=f[iw][is][j][i];
        }
        }
        }

        /* solver the Helmholtz equation */
        sparse_solve(npml, pad1, pad2, fin, hermite_false, ns, uts);

        /* Born sources */
        for ( is = 0; is < ns; is++ ) { 
        for ( j = 0; j < n2; j++ ) {
        for ( i = 0; i < n1; i++ ) {
            fborn[is][j][i] = omega*omega*fin[is][j][i]*input[j][i]/(v[j][i]*v[j][i]);
        }
        }
        }
       
        /* calculate Born synthetic */
        sparse_solve(npml, pad1, pad2, fborn, hermite_false, ns, uts);
            
        for ( is = 0; is < ns; is++ ) { 
        for ( j = 0; j < n2; j++ ) {
        for ( i = 0; i < n1; i++ ) {
            output[iw][is][j][i] = fborn[is][j][i];
        }
        }
        }

        /* free memory */
        sparse_free(uts);
    }
}

void adjlsm_operator(const int nw, 
                     const float ow,
                     const float dw,
                     const int ns,
                     const int n1,
                     const int n2,
                     const int uts,
                     const int pad1,
                     const int pad2,
                     const int npml,
                     const float d1,
                     const float d2,
                     bool hermite_false,
                     bool hermite_true,
                     float **v,
                     sf_complex ****f,
                     float **recloc,
                     sf_complex ****input,
                     float **output) 
/*< adjoint operator for lsm >*/ 
{
    int iw, is, j, i; 
    double omega; 
    sf_complex ***obs_in, ***obs_cut, ***srcw, ***recw; 

    obs_in = sf_complexalloc3(n1,n2,ns);
    obs_cut = sf_complexalloc3(n1,n2,ns);
    srcw = sf_complexalloc3(n1,n2,ns);
    recw = sf_complexalloc3(n1,n2,ns);

    for (iw = 0; iw < nw; iw ++ ) { 
        omega=(double) 2.* SF_PI * ( ow + iw * dw ); 
     
        for ( is = 0; is < ns; is ++ ) { 
        for ( j = 0; j < n2; j ++) { 
        for ( i = 0; i < n1; i++ ) { 
            obs_in[is][j][i] = input[iw][is][j][i];
        }
        }
        }

        /* generate adjoint source for reverse time migration */
        genadjsrc_rtm(obs_in, obs_cut, recloc, n1, n2, ns); 

        /* initialize sparse solver */ 
        sparse_init(uts, pad1, pad2);

        /* factorize matrix, change according to different frequencies and models */
        sparse_factor(omega,n1,n2,d1,d2,v,npml,pad1,pad2);

        for (is = 0; is < ns; is++ ) { 
        for (j = 0; j < n2; j++ ) { 
        for (i = 0; i < n1; i++ ) { 
            srcw[is][j][i] = f[iw][is][j][i];
            recw[is][j][i] = obs_cut[is][j][i];
        }
        }
        }

        /* sparse solver for source wavefield */ 
        sparse_solve(npml, pad1, pad2, srcw, hermite_false, ns, uts);

        /* sparse solver for receiver wavefield */
        sparse_solve(npml, pad1, pad2, recw, hermite_true, ns, uts);

        /* imaging condition */
        for ( j = 0; j < n2; j++ ) {
        for ( i = 0; i < n1; i++ ) {
            output[j][i] = 0.0; 
            for ( is = 0; is < ns; is++ ) {
                output[j][i] += crealf(omega*omega*conjf(srcw[is][j][i])*recw[is][j][i]/(v[j][i]*v[j][i]));
            }
        }
        }

        /* free memory */
        sparse_free(uts);
    } /* end frequency */
}


