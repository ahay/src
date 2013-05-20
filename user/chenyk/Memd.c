/* Empirical Mode Decomposition */
/*
  Copyright (C) 2013 the University of Texas at Austin

  This program is the Madagascar version of emd program written by G. Rilling (2007).

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
#include "emdutil.h"
#include <stdlib.h>
#include <stdio.h>

int main(int argc, char* argv[])
{
    /* declarations */
    int i,n,nb_imfs,max_imfs,iteration_counter,max_iterations,stop_status,allocated_x,stop_EMD;
    extrema_t ex;
    envelope_t env;
    stop_t stop_params;
    float threshold, tolerance;
    double *x,*y,*z,*m,*a;
    float *dat,*imf,dt;
    imf_list_t list;
    sf_file inp, outp; 

    sf_init(argc,argv);
    inp  = sf_input("in");
    outp = sf_output("out");

    if (SF_FLOAT != sf_gettype(inp)) sf_error("Need float input");
    if (!sf_histint(inp,"n1",&n)) sf_error("No n1= in input");
    if (!sf_histfloat(inp,"d1",&dt)) sf_error("No d1= in input");
    x=(double*)sf_alloc(n,sizeof(double));
    y=(double*)sf_alloc(n,sizeof(double)); 

    dat=sf_floatalloc(n);
    imf=sf_floatalloc(n*20);

    sf_floatread(dat,n,inp);
    for(i=0;i<n;i++){y[i]=dat[i];x[i]=i*dt;}

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

    allocated_x=0;

    /* initialisations */
    ex=init_extr(n+2*NBSYM);
    list=init_imf_list(n);
    z=(double *)malloc(n*sizeof(double));
    m=(double *)malloc(n*sizeof(double));
    a=(double *)malloc(n*sizeof(double));
    env=init_local_mean(n+2*NBSYM);
  
    /* MAIN LOOP */
    nb_imfs=0;
    stop_EMD=0;
    while ((!max_imfs || (nb_imfs < max_imfs)) && !stop_EMD) {   
        /* initialisation */
	for (i=0;i<n;i++) z[i]=y[i];
	for (i=0;i<n;i++) m[i]=y[i];
	iteration_counter=0;
	stop_status = mean_and_amplitude(x,z,m,a,n,&ex,&env);
        /* SIFTING LOOP */
	while (!stop_status && !stop_sifting(m,a,&ex,&stop_params,n,iteration_counter, max_iterations)) {
            /* subtract the local mean */
	    for (i=0;i<n;i++) z[i]=z[i]-m[i];
	    iteration_counter++;
	    stop_status = mean_and_amplitude(x,z,m,a,n,&ex,&env);      
	}  
        /* save current IMF into list if at least   */
        /* one sifting iteration has been performed */
	if (iteration_counter) {
	    add_imf(&list,z,iteration_counter);

	    nb_imfs++;
	    for (i=0;i<n;i++){
		y[i]=y[i]-z[i];
		imf[i+(nb_imfs-1)*n]=z[i];	
	    }      
	}
	else
	    stop_EMD = 1;    
    }
  
    /* save the residual into list */
    add_imf(&list,y,0);
    for (i=0;i<n;i++) imf[i+nb_imfs*n]=y[i];	

    sf_putint(outp,"n1",n);
    sf_putint(outp,"n2",nb_imfs+1);
    sf_putfloat(outp,"d1",dt);
    sf_putfloat(outp,"d2",1);
    sf_putfloat(outp,"o1",0);
    sf_putfloat(outp,"o2",0);  

    /* output  */
    sf_floatwrite(imf,n*(nb_imfs+1),outp);
  
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
