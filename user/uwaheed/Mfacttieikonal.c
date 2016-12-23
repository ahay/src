/* Fast sweeping factored TTI eikonal solver (2D) */
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

#include "facttieikonal.h"
#include "derivative.h"

void sf_derivative_2D (float *,float *, float *,
                       double , double ,
                       int , int , int );

int main (int argc,char* argv[]) {
    int n1, n2, n3, i, j, loop, nshot, is, n12, ndim, niter, nfpi, fac;
    float o1, o2, o3, d1, d2, d3;
    float **s, *t, *tau, *vz, *vx, *eps, *del, *theta, *T0, *py0, *pz0, *tn, *tn1, a0, b0, c0, vx0, vz0, st0, ct0;
    char *sfile = NULL, *file = NULL;
    int sy, sz, loc=0; 
    float *py, *pz, pydash, pzdash, *rhs, sum;
    bool optloc;
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

    if (!sf_getint ("fac", &fac)) fac = 1;
    /* Type of factorization: (0)Additive, (1)Multiplicative */
    /* Multiplicative factorization is more stable */

    if (!sf_getbool ("optloc", &optloc)) optloc = false;
    /* Selects optimal location for homogeneous medium parameter */
    /* Useful for stability of additive factorization when the highest velocity in the medium is 
       much larger than the velocity at the source point */



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
    tau  = sf_floatalloc (n12);
    vz  = sf_floatalloc (n12);
    vx  = sf_floatalloc (n12);
    eps  = sf_floatalloc (n12);
    del  = sf_floatalloc (n12);
    theta  = sf_floatalloc (n12);
    T0  = sf_floatalloc (n12);
    py0  = sf_floatalloc (n12);
    pz0  = sf_floatalloc (n12);
    py  = sf_floatalloc (n12);
    pz  = sf_floatalloc (n12);
    rhs = sf_floatalloc (n12);
    tn = sf_floatalloc (n12);
    tn1 = sf_floatalloc (n12);


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
        theta[i] = theta[i]*SF_PI/180.0;       /* Convert angle from degrees to radians */
        rhs[i] = 1.0;     /* Inititialize rhs to 1 */
        tn[i] = 0.;       /* tn is the current time, and tn1 is the time from previous iteration */

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


			if(optloc) loc = (int)(n1/2.);

            if(optloc){
                vx0 = vx[sy*n1+loc]; vz0 = vz[sy*n1+loc]; 
                ct0 = cos(theta[sy*n1+loc]); st0 = sin(theta[sy*n1+loc]);

                a0 = (vx0*vx0*ct0*ct0 + vz0*vz0*st0*st0)/rhs[sy*n1+loc];
                b0 = (vx0*vx0*st0*st0 + vz0*vz0*ct0*ct0)/rhs[sy*n1+loc];
                c0 = ((vx0*vx0 - vz0*vz0)*st0*ct0)/rhs[sy*n1+loc];

                sf_warning("Homogeneous model parameters taken at (%f,%f)"
                           ,o2+(sy)*d2,o1+(loc)*d1);
            }
            else{
                vx0 = vx[sy*n1+sz]; vz0 = vz[sy*n1+sz]; 
                ct0 = cos(theta[sy*n1+sz]); st0 = sin(theta[sy*n1+sz]);

                a0 = (vx0*vx0*ct0*ct0 + vz0*vz0*st0*st0)/rhs[sy*n1+sz];
                b0 = (vx0*vx0*st0*st0 + vz0*vz0*ct0*ct0)/rhs[sy*n1+sz];
                c0 = ((vx0*vx0 - vz0*vz0)*st0*ct0)/rhs[sy*n1+sz];

                sf_warning("Homogeneous model parameters taken at (%f,%f)"
                           ,o2+(sy)*d2,o1+(sz)*d1);
            }

            




            /* Traveltime and derivative computation for homogeneous TEA medium */
            for(i=0;i<n1;i++){
                for(j=0;j<n2;j++){

                    if(i==sz && j==sy){
                        T0[j*n1+i] = 0.; 
                        py0[j*n1+i] =0;  pz0[j*n1+i]=0;
                        continue;
                    } 


                    T0[j*n1+i] = sqrtf((b0*(j-sy)*(j-sy)*d2*d2 - 2*c0*(j-sy)*(i-sz)*d1*d2 
                                 + a0*(i-sz)*(i-sz)*d1*d1)/(a0*b0-c0*c0));


                    py0[j*n1+i] = (b0*(j-sy)*d2 - c0*(i-sz)*d1)/(sqrtf((b0*(j-sy)*(j-sy)*d2*d2
                                   - 2*c0*(j-sy)*(i-sz)*d1*d2 + a0*(i-sz)*(i-sz)*d1*d1)
				    			   *(a0*b0-c0*c0)));

                    pz0[j*n1+i] = (a0*(i-sz)*d1 - c0*(j-sy)*d2)/(sqrtf((b0*(j-sy)*(j-sy)*d2*d2
                                   - 2*c0*(j-sy)*(i-sz)*d1*d2 + a0*(i-sz)*(i-sz)*d1*d1)
						    	   *(a0*b0-c0*c0)));

                }
            }


            if (nshot > 1)
                sf_warning ("Calculating shot %d of %d", is + 1, nshot);
            if (false == sf_init_fast_sweep (tau,
                                             n2, n1,
                                             o2, o1,
                                             d2, d1,
                                             sy, sz,fac))
                sf_error ("Incorrect shot location");

            sf_run_fast_sweep (tau, T0, py0, pz0, 
                               vz, vx, theta,
                               rhs, niter,
                               n2, n1,
                               o2, o1,
                               d2, d1,
                               sy, sz, fac);

            for(i=0;i<n12;i++){

                if(fac==0) {t[i] = T0[i]+tau[i]; tn1[i] = tn[i]; tn[i] = t[i];}
                else if(fac==1){ t[i] = T0[i]*tau[i]; tn1[i] = tn[i]; tn[i] = t[i];}
                else sf_error("Choose fac=0(Additive) or fac=1(Multiplicative) factorization");
            
            }

            sf_derivative_2D(tau,pz,py,d2,d1,n2,n1,2);

            sum = 0.;

            /* dT/dx and dT/dz computation */
            for(i=0;i<n12;i++){
                if(fac==0){
                    py[i] = py[i]+py0[i]; pz[i] = pz[i]+pz0[i];
                }
                else if(fac==1){
                    py[i] = T0[i]*py[i] + tau[i]*py0[i]; pz[i] = T0[i]*pz[i] + tau[i]*pz0[i];
                }
                
                sum = sum + fabs(tn[i]-tn1[i]);

            }

            sf_warning("==================================");
            sf_warning("L1 norm of update is %f",sum/n12);
            sf_warning("==================================");

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


    free(t); free(tau); free(vz); free(vx); free(eps);
    free(del); free(theta); free(T0); free(py0);
    free(pz0); free(py); free(pz); free(rhs);
    free(tn); free(tn1);
    free(*s); free(s);


    return 0;
}
