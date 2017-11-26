/* Fast sweeping TTI eikonal solver (2D) */
/*
  Copyright (C) 2016 Princeton University
  
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

#include <math.h>
#include <rsf.h>

#include "ttieikonal.h"
#include "derivative.h"

void sf_derivative_2D (float *,float *, float *,
                       double , double ,
                       int , int , int );

int main (int argc,char* argv[]) {
    int n1, n2, n3 = 1, i, j, loop, nshot, is, n12, ndim, niter, nfpi;
    float o1, o2, o3 , d1, d2, d3;
    float **s, *t, *tn, *tn1, *vz, *vx, *eps, *del, *theta;
    char *sfile = NULL, *file = NULL;
    int sy, sz; 
    float *py, *pz, pydash, pzdash, *rhs, sum;
    sf_file vzf, time, shots, epsilonf, deltaf, thetaf;



    sf_init (argc, argv);
    vzf = sf_input ("in");
    time = sf_output ("out");


    if (SF_FLOAT != sf_gettype (vzf))
        sf_error("Need float input");
    if (!sf_histint (vzf, "n1", &n1)) sf_error ("No n1= in input");
    if (!sf_histint (vzf, "n2", &n2)) sf_error ("No n2= in input");
    if (!sf_histint (vzf, "n3", &n3)) n3 = 1;

    if (!sf_histfloat (vzf, "d1", &d1)) sf_error ("No d1= in input");
    if (!sf_histfloat (vzf, "d2", &d2)) sf_error ("No d2= in input");
    if (!sf_histfloat (vzf, "d3", &d3)) d3 = d2;

    if (!sf_histfloat (vzf, "o1", &o1)) o1 = 0.;
    if (!sf_histfloat (vzf, "o2", &o2)) o2 = 0.;
    if (!sf_histfloat (vzf, "o3", &o3)) o3 = 0.;


    if (!sf_getint ("niter", &niter)) niter = 4;
    /* number of sweeping iterations */

    if (!sf_getint ("nfpi", &nfpi)) nfpi = 3;
    /* number of fixed-point iterations */


    /* Assigning file pointers for input parameters */

    if(NULL!=(file = sf_getstring("epsilon"))){
	epsilonf = sf_input(file);
	free(file);
    }
    else{
	epsilonf = NULL;
    }


    if(NULL!=(file = sf_getstring("delta"))){
	deltaf = sf_input(file);
	free(file);
    }
    else{
	deltaf = NULL;
    }

    if(NULL!=(file = sf_getstring("theta"))){
	thetaf = sf_input(file);
	free(file);
    }
    else{
	thetaf = NULL;
    }



    sfile = sf_getstring ("shotfile");
    /* File with shot locations (n2=number of shots, n1=3) */

    if (NULL != sfile) {
        shots = sf_input (sfile);

        if (SF_FLOAT != sf_gettype (shots)) 
            sf_error ("Need float shotfile");
        if (!sf_histint (shots, "n1", &ndim) || ndim != 3)
            sf_error ("Need n1=3 in shotfile");
        nshot = sf_leftsize (shots, 1);

        s = sf_floatalloc2 (ndim, nshot);
        sf_floatread (s[0], nshot * ndim, shots);
        sf_fileclose (shots);

        sf_putint (time, 1 == n3 ? "n3" : "n4", nshot);
        free (sfile); sfile = NULL;
    } else {
        nshot = 1;
        ndim = 2;

        s = sf_floatalloc2 (ndim, nshot);

        if (!sf_getfloat ("zshot", &s[0][0])) s[0][0] = 0.; 
        /* Shot location (used if no shotfile) */
        if (!sf_getfloat ("yshot", &s[0][1])) s[0][1] = o2 + 0.5*(n2-1)*d2;

        sf_warning ("Shooting from zshot=%g yshot=%g",
                    s[0][0], s[0][1]);
    }


    n12 = n1*n2;

    t  = sf_floatalloc (n12);
    vz  = sf_floatalloc (n12);
    vx  = sf_floatalloc (n12);
    eps  = sf_floatalloc (n12);
    del  = sf_floatalloc (n12);
    theta  = sf_floatalloc (n12);
    py  = sf_floatalloc (n12);
    pz  = sf_floatalloc (n12);
    rhs = sf_floatalloc (n12);
    tn  = sf_floatalloc (n12);
    tn1  = sf_floatalloc (n12);



    /* Reading input parameters */

    sf_floatread (vz, n12, vzf);


    if( epsilonf != NULL){
	    sf_floatread(eps,n12,epsilonf);
	    sf_fileclose(epsilonf);
	}
    else{
        for(i=0;i<n12;i++){
            eps[i] = 0.0;
	    }
    }


    if( deltaf != NULL){
	    sf_floatread(del,n12,deltaf);
	    sf_fileclose(deltaf);
	}
    else{
        for(i=0;i<n12;i++){
            del[i] = eps[i];
	    }
    }

    if( thetaf != NULL){
	    sf_floatread(theta,n12,thetaf);
	    sf_fileclose(thetaf);
	}
    else{
        for(i=0;i<n12;i++){
            theta[i] = 0.0;
	    }
    }

    /* Convert angles from degrees to radians */

    for(i=0;i<n12;i++){
        
        vx[i] = vz[i]*sqrtf(1+2*eps[i]); /* Compute horizontal velocity */
        theta[i] = theta[i]*SF_PI/180.0; /* Convert angle from degrees to radians */
        rhs[i] = 1.0;     /* Inititialize rhs to 1 */
        tn[i] = 0.;       /* tn is the current traveltime estimate, and tn1 is the estimate from previous iteration */

    }


    if (sfile) {
        free (sfile); sfile = NULL;
    }

    sf_warning ("Performing 2-D sweeps");

    /* loop over shots */
    for (is = 0; is < nshot; is++) {

        /* Converting source location into grid points */
        sy = (int)((s[is][1] - o2) / d2 + 0.5f);
        sz = (int)((s[is][0] - o1) / d1 + 0.5f);


        for(loop=0;loop<nfpi;loop++){
    
            sf_warning("Fixed-point iteration %d of %d", loop+1, nfpi);


            if (nshot > 1)
                sf_warning ("Calculating shot %d of %d", is + 1, nshot);
            if (false == sf_init_fast_sweep (t, n2, n1,
                                             o2, o1, d2, 
                                             d1, sy, sz))
                sf_error ("Incorrect shot location");

            sf_run_fast_sweep (t, vz, vx, 
                               theta, rhs,
                               niter, n2, n1,
                               o2, o1, d2, d1,
                               sy, sz);

            sf_derivative_2D(t,pz,py,d2,d1,n2,n1,2);

            /* Calculating L1 norm of the update made by fixed point iteration */

            sum = 0.;

            for(i=0;i<n12;i++) { 
         
                tn1[i] = tn[i];   
                tn[i] = t[i]; 
                sum = sum + fabs(tn[i]-tn1[i]);
            } 

            sf_warning("==================================");
            sf_warning("L1 norm of update is %f",sum/n12);
            sf_warning("==================================");   


            /* Evaluating the rhs term */

            for(i=0;i<n1;i++){
                for(j=0;j<n2;j++){
 
                    pydash = cos(theta[j*n1+i])*py[j*n1+i] + sin(theta[j*n1+i])*pz[j*n1+i];
                    pzdash = cos(theta[j*n1+i])*pz[j*n1+i] - sin(theta[j*n1+i])*py[j*n1+i];
                    rhs[j*n1+i] = 1 + 2*((eps[j*n1+i]-del[j*n1+i])/(1+2*eps[j*n1+i]))
                                *vx[j*n1+i]*vx[j*n1+i]*vz[j*n1+i]*vz[j*n1+i]*pydash*pydash*pzdash*pzdash;

                }
            }

        } // for(loop=0;loop<nfpi;loop++)


        sf_floatwrite (t, n12, time); /* Writing the traveltime into the file */

    } //for (is = 0; is < nshot; is++)


    free(t); free(vz); free(vx); free(eps);
    free(del); free(theta);  
    free(py); free(pz); free(rhs);
    free(tn); free(tn1);
    free(*s); free(s);


    return 0;
}
