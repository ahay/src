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

//#define STOP_DEFAULT {.threshold = 0.05, .tolerance = 0.05}
#define DEFAULT_THRESHOLD 0.05
#define DEFAULT_TOLERANCE 0.05
#define MAX_ITERATIONS 1000
#define LIM_GMP 30000
#define NBSYM 2

int main(int argc, char* argv[])
{
    /* declarations */
    int i,n,nb_imfs,max_imfs,iteration_counter,stop_status,allocated_x,stop_EMD;
    extrema_t ex;
    input_t input;
    envelope_t env;
    stop_t stop_params;
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

    /* get input data */
    input.stop_params.threshold=DEFAULT_THRESHOLD;
    input.stop_params.tolerance=DEFAULT_TOLERANCE;
    input.allocated_x=0;
#ifdef _ALT_MEXERRMSGTXT
    input.error_flag=0;
#endif  
    input.max_imfs=0;
    input.n=n;  
    input.x=x;
    input.y=y;

    max_imfs=input.max_imfs;
    stop_params=input.stop_params;
    allocated_x=input.allocated_x;
    x=input.x;
    y=input.y;
  
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
    
	while (!stop_status && !stop_sifting(m,a,&ex,&stop_params,n,iteration_counter)) {
      
            /* subtract the local mean */
	    for (i=0;i<n;i++) z[i]=z[i]-m[i];
	    iteration_counter++;
      
	    stop_status = mean_and_amplitude(x,z,m,a,n,&ex,&env);
      
      
	}
    
        /* save current IMF into list if at least     */
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
