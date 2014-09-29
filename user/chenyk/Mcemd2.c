/* Bivariate empirical mode decomposition using second algorithm. */

/*
  Copyright (C) 2013 the University of Texas at Austin

  This program is the Madagascar version of cemdc2 program written by G. Rilling (2007).

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
#include "cemdutil2.h"
#include <stdlib.h>
#include <stdio.h>

int main(int argc, char* argv[])
{
    /* declarations */
    int i,n,nb_imfs,max_imfs,max_iterations,iteration_counter,stop_status,allocated_x,stop_EMD,nbphases;
    float dt;
    double *x,*a;

    float threshold, tolerance;
    extrema_t ex;
    envelope_t env;
    stop_t stop_params;
    COMPLEX_T *y,*m,*z;
    imf_list_t list;
    kiss_fft_cpx *ind, *imf;
    sf_file inp, outp; 

    sf_init(argc,argv);
    inp  = sf_input("in");
    outp = sf_output("out");

    if (SF_COMPLEX != sf_gettype(inp)) sf_error("Need complex input");
    if (!sf_histint(inp,"n1",&n)) sf_error("No n1= in input");
    if (!sf_histfloat(inp,"d1",&dt)) sf_error("No d1= in input");

    ind=(kiss_fft_cpx*) sf_complexalloc(n);
    imf=(kiss_fft_cpx*) sf_complexalloc(n*20);
    sf_floatread((float*)ind,n*2,inp);

    x = (double *)malloc(n*sizeof(double));
    y=(COMPLEX_T *)malloc(n*sizeof(COMPLEX_T));

    for(i=0;i<n;i++)
	{
	y[i].r=ind[i].r;y[i].i=ind[i].i;
	x[i]=i*dt;
	}

    if(!sf_getfloat("threshold",&threshold)) threshold=DEFAULT_THRESHOLD;
    /* Sifting stoping parameter: threshold, the default is 0.05. */
    stop_params.threshold=threshold;

    if(!sf_getfloat("tolerance",&tolerance)) tolerance=DEFAULT_TOLERANCE;
    /* Sifting stoping parameter: tolerance, the default is 0.05. */
    stop_params.tolerance=tolerance;

    /* input checking */
    if (stop_params.threshold <= 0 || stop_params.threshold >=1)
        sf_warning("threshold must be a real number in [O,1]");
    if (stop_params.tolerance < 0 || stop_params.tolerance >= 1)
        sf_warning("tolerance must be a real number in [O,1]");
    
    if(!sf_getint("miter", &max_iterations)) max_iterations = MAX_ITERATIONS;
    /* Maximum number of iterations during sifting, the default is 1000. */
    
    if(!sf_getint("mimf", &max_imfs)) max_imfs = 0;
    /* Maximum number of IMFs, the default is as many as possible. */

    if(!sf_getint("nbdir", &nbphases)) nbphases = DEFAULT_NBPHASES;
    /* Number of directions used to compute the local mean, the default is 4. */

    allocated_x=0;

    /* initialisations */
    list=init_imf_list(n);
    z=(COMPLEX_T *)malloc(n*sizeof(COMPLEX_T));
    m=(COMPLEX_T *)malloc(n*sizeof(COMPLEX_T));
    a=(double *)malloc(n*sizeof(double));
    ex=init_extr(n+2*NBSYM);
    env=init_local_mean(n+2*NBSYM);
  
    /* MAIN LOOP */
    nb_imfs=0;
    stop_EMD=0;
    while ((!max_imfs || (nb_imfs < max_imfs)) && !stop_EMD) {
    /* initialisation */
    for (i=0;i<n;i++) z[i]=y[i];
    for (i=0;i<n;i++) m[i]=y[i];
    iteration_counter=0;
    
    stop_status = mean_and_amplitude(x,z,m,a,n,nbphases,&ex,&env);
   
    /* SIFTING LOOP */    
    while (!stop_status && !stop_sifting(m,a,&ex,&stop_params,n,iteration_counter,max_iterations)) {
      /* subtract the local mean */
      #ifdef C99_OK
      for (i=0;i<n;i++) z[i]=z[i]-m[i];
      #else
      for (i=0;i<n;i++) {
        z[i].r=z[i].r-m[i].r;
        z[i].i=z[i].i-m[i].i;
      }
      #endif
      iteration_counter++;
      stop_status = mean_and_amplitude(x,z,m,a,n,nbphases,&ex,&env);    
    }
    
    /* save current IMF into list if at least     */
    /* one sifting iteration has been performed */
    if (iteration_counter) {
      add_imf(&list,z,iteration_counter);
      nb_imfs++;
      #ifdef C99_OK
      for (i=0;i<n;i++) y[i]=y[i]-z[i];
      #else
      for (i=0;i<n;i++) {
        y[i].r=y[i].r-z[i].r;
        y[i].i=y[i].i-z[i].i;
	imf[i].r=z[i].r;
	imf[i].i=z[i].i;
      }
      #endif
    }
    else
      stop_EMD = 1;    
    }

    /* save the residual into list */
    add_imf(&list,y,0);

    for (i=0;i<n;i++) {imf[i+nb_imfs*n].r=y[i].r; imf[i+nb_imfs*n].i=y[i].i;};	

    sf_settype(outp, SF_COMPLEX);
    sf_putint(outp,"n1",n);
    sf_putint(outp,"n2",nb_imfs+1);
    sf_putfloat(outp,"d1",dt);
    sf_putfloat(outp,"d2",1);
    sf_putfloat(outp,"o1",0);
    sf_putfloat(outp,"o2",0);  

    /* output  */
    sf_complexwrite((sf_complex*) imf,n*(nb_imfs+1),outp);
  
    /* free allocated memory */
    if (allocated_x)
       free(x);
    free(y);
    free(m);
    free(a);
    free_local_mean(env);
    free(z);
    free_imf_list(list);
    free_extr(ex);

    exit (0);
}
