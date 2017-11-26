/* Subroutines used for optimization */
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

/* subroutines for rtm */
void genadjsrc_rtm(sf_complex ***obs,
                   sf_complex ***adj, 
                   float **recloc,
                   const int n1, 
                   const int n2, 
                   const int ns)
/*< generate adjoint source for rtm >*/
{
    int i, j, k;
    for ( k = 0; k < ns; k++ ) {
    for ( j = 0; j < n2; j++ ) { 
    for ( i = 0; i < n1; i++ ) { 
        adj[k][j][i] = sf_cmplx(0.0,0.0); 
            
        /* pick out receiver location */
        if (recloc[j][i] > 0.0 ) { 
            adj[k][j][i] = obs[k][j][i];
        }
    }
    }
    }
}

/* subroutines for fwi */
float genadjsrc_fwi(sf_complex ***obs,
                    sf_complex ***syn,
                    sf_complex ***adj,
                    float **recloc,
                    const int n1, 
                    const int n2,
                    const int ns)
/*< generate adjoint source for fwi >*/
{
    int i, j, k;
    sf_complex temp;
    float misfit;

    misfit = 0.0; 
    for ( k = 0; k < ns; k++ ){ 
    for ( j = 0; j < n2; j++) {
    for ( i = 0; i < n1; i++ ){ 
        adj[k][j][i] = sf_cmplx(0.0,0.0); 

        if ( recloc[j][i] > 0.0 ) { 
            temp = obs[k][j][i] - syn[k][j][i];
            misfit += cabsf(temp) * cabsf(temp); 
            adj[k][j][i] = temp;
        }
    }  /* for i  */ 
    }  /* for j  */
    }  /* for ns */

    return misfit;
}

float calc_misfit(sf_complex ***obs,
                  sf_complex ***syn,
                  float **recloc,
                  const int n1, 
                  const int n2, 
                  const int ns)
/*< calculate misfit value >*/
{
    int i, j, k;
    float misfit; 
    sf_complex temp;

    misfit = 0.0; 
    for ( k = 0; k < ns; k++ ){ 
    for ( j = 0; j < n2; j++) {
    for ( i = 0; i < n1; i++ ){ 
        if ( recloc[j][i] > 0.0 ) { 
            temp = obs[k][j][i] - syn[k][j][i];
            misfit += cabsf(temp) * cabsf(temp); 
        }
    }  /* for i */
    }  /* for j */
    }  /* for ns */
    return misfit;
}

void direction_sd(float **g,
                  float **d,
                  const int n1, 
                  const int n2)
/*< calculate steepest descent direction >*/
{
   int i, j;
   for ( j = 0; j < n2; j++ ){
   for ( i = 0; i < n1; i++ ){
        d[j][i] = -g[j][i];
   }
   }
}

float direction_cg_polak(float **g0,
                         float **d0,
                         float **g1,
                         float **d1,
                         const int n1, 
                         const int n2)
/*< calculate conjugate gradient direction >*/ 
{
    int i,j; 
    float diff, beta_top, beta_bot, beta;

    beta_top=0.0;
    beta_bot=0.0;

    for ( j = 0; j < n2; j++ ) {
    for ( i = 0; i < n1; i++ ) {
        diff = g1[j][i] - g0[j][i];
        beta_top += g1[j][i] * diff;
        beta_bot += g0[j][i] * g0[j][i];
    }
    }

    beta=beta_top / beta_bot;
    
    if (beta < 0.0 ) { 
        beta =0.0; 
    }

    for ( j = 0; j < n2; j++ ){
    for ( i = 0; i < n1; i++ ) { 
        d1[j][i] = -g1[j][i] + beta * d0[j][i];
    }
    }
    return beta; 
}


void update_model_fwi(float **m0, 
                      float **m1,
                      float **d, 
                      const float alpha,
					  const float max,
                      const int n1, 
                      const int n2)
/*< update model for fwi, nonlinear >*/ 
{
    int i, j;
    float dmax;

    dmax=0.0;
    for ( j = 0; j < n2; j++ ) { 
    for ( i = 0; i < n1; i++ ) { 
    if ( fabs( d[j][i] ) > dmax ) {
        dmax = fabs( d[j][i] );
    }
    }
    }

    for ( j = 0; j < n2; j++ ){ 
    for ( i = 0; i < n1; i++ ) { 
        m1[j][i] = m0[j][i] + alpha * d[j][i]*max / dmax;
    }
    }
}


/* subroutines for lsm */ 
float calc_alpha(float **g, 
                 sf_complex ****Fg,
                 float **recloc,
                 const int n1,
                 const int n2,
                 const int ns,
                 const int nw) 
/*< calculate alpha value for linear inversion >*/
{
    float alpha_top, alpha_bot, alpha; 
    int i, j, iw, is; 

    alpha_top = 0.0 ; 
    alpha_bot = 0.0 ; 

    for ( j = 0; j < n2; j++) { 
    for ( i = 0; i < n1; i++) { 
        alpha_top += g[j][i] * g[j][i];
    }
    }

    for ( iw = 0; iw < nw; iw++ ) {
    for ( is = 0; is < ns; is++ ) { 
    for ( j = 0; j < n2; j++ ) { 
    for ( i = 0; i < n1; i++ ) { 
        if ( recloc[j][i] > 0.0 ) { 
            alpha_bot += cabsf( conjf( Fg[iw][is][j][i]) * Fg[iw][is][j][i] );
        }
    }
    }
    }
    }
    
    alpha=alpha_top / alpha_bot; 

    return alpha; 
}

void update_model_lsm(float **m0, 
                      float **m1,
                      float **d, 
                      const float alpha,
                      const int n1, 
                      const int n2)
/*< update model for lsm, linear >*/ 
{
    int i, j;

    for ( j = 0; j < n2; j++ ){ 
    for ( i = 0; i < n1; i++ ) { 
        m1[j][i] = m0[j][i] + alpha * d[j][i];
    }
    }
}

void update_residual(sf_complex ****r0,
                     sf_complex ****r1,
                     sf_complex ****Fg,
                     float **recloc, 
                     const float alpha, 
                     const int n1,
                     const int n2,
                     const int ns, 
                     const int nw)
/*< update data residual >*/
{
    int iw, is, j, i; 

    for ( iw = 0; iw < nw; iw ++ ) { 
    for ( is = 0; is < ns; is ++ ) {
    for ( j = 0; j < n2; j++ ) { 
    for ( i = 0; i < n1; i++ ) {
        if (recloc[j][i] > 0.0 ) { 
            r1[iw][is][j][i] = r0[iw][is][j][i] - alpha * Fg[iw][is][j][i];
        }
    }
    }
    }
    }
}

float direction_cg_fletcher(float **g0,
                            float **d0,
                            float **g1,
                            float **d1,
                            const int n1, 
                            const int n2)
/*< calculate conjugate gradient direction >*/ 
{
    int i,j; 
    float beta_top, beta_bot, beta;

    beta_top=0.0;
    beta_bot=0.0;
    for ( j = 0; j < n2; j++ ) {
    for ( i = 0; i < n1; i++ ) {
        beta_top += g1[j][i] * g1[j][i];
        beta_bot += g0[j][i] * g0[j][i];
    }
    }

    beta=beta_top / beta_bot;
    
    if (beta < 0.0 ) { 
        beta =0.0; 
    }

    for ( j = 0; j< n2; j++ ){
    for ( i = 0; i < n1; i++ ) { 
        d1[j][i] = g1[j][i] + beta * d0[j][i];
    }
    }
    return beta; 
}

