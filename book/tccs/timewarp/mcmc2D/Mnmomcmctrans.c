/* 2D NMO GMA MCMC transdimensional inversion with Metropolis rule (Mosegaard and Tarantola, 1995) */
/*
  Copyright (C) 2009 University of Texas at Austin
  
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
#include <math.h>
#include <time.h>

/* Global variables */
int i,j;
float testval,x,t0,L,w1,A1,B1,C1,sigmapc,A,B,C;
static const double two_pi = 2.0*3.14159265358979323846;

/* Function declaration */
static float likelihood(float *offsetx, int noffset, float lmx /* maximum offset to compute likelihood */,
				float datarms /* RMS average of data*/,
				float *dattime, /*true time (data) at the time (depth) slice*/
				float *t0sq, /*t0 sq at the time (depth) slice*/
				float *esttime, /*est time at the time (depth) slice*/				
				float *newcoef /*new coeff drawn from uniform sampling*/);

/* Main */
int main(int argc, char* argv[])
{
  int noffset, nt, npara, nmodel, nc, seed, reject, saveiter, index;
  float dt, tini, currentL, newL;
  float *t0sq, *dattime, *esttime, *rangecoef, *newcoef, *currentcoef, *drangecoef, *midcoef, *offx, **final, *limitx, *drms; 
  int k1,k2,k3,k4,l;
  bool prior,wgauss,eta;

  sf_file inp, t0sqf, out ,rangecoeff, ofx, limx, datrms;

  sf_init(argc, argv);
    inp = sf_input("in");
    t0sqf = sf_input("t0sq");
    rangecoeff = sf_input("rangecoef");
    ofx = sf_input("offsetx");
    limx = sf_input("limitx");
    datrms = sf_input("datrms");
    out = sf_output("out");

    /* input vector data t^2' */
    if (!sf_histint(inp,"n1",&noffset)) sf_error("No n1=");
    if (!sf_histint(inp,"n2",&nt)) sf_error("No n2=");
    if (!sf_histfloat(inp,"d2",&dt)) sf_error("No d2=");
    if (!sf_histfloat(inp,"o2",&tini)) sf_error("No o2=");

	
	if (!sf_getint("seed",&seed)) seed = time(NULL);
    /* random seed */
    init_genrand((unsigned long) seed);
    
    if (!sf_getint("saveiter",&saveiter)) saveiter = 20;
    /* save state every iter  */
    
    if (!sf_getbool("prior",&prior)) prior = false;
    /* generate prior or posterior */
    
    if (!sf_getbool("wgauss",&wgauss)) wgauss = false;
    /* use gaussian distribution for W */
    
    if (!sf_getbool("eta",&eta)) eta = false;
    /* use eta to constrain A B and C */

    /*Number of fitting parameters*/
    npara=5; /* W A B C sigma*/
    
    if (!sf_getint("nmodel",&nmodel)) nmodel=1000;
    /*Get the number of MC models*/
        
    /*Memory allocation*/
    rangecoef = sf_floatalloc(2*npara);
	drangecoef = sf_floatalloc(npara);
	midcoef = sf_floatalloc(npara);
    currentcoef = sf_floatalloc(npara);
    newcoef = sf_floatalloc(npara);
    dattime = sf_floatalloc(noffset);
    t0sq = sf_floatalloc(noffset);
    esttime = sf_floatalloc(noffset);
    offx = sf_floatalloc(noffset);
    limitx = sf_floatalloc(nt);
    drms = sf_floatalloc(nt);
    final = sf_floatalloc2(npara,nmodel);

    if (!sf_histint(rangecoeff,"n1",&nc) || nc != 2*npara) 
	sf_error("Need n1=%d in rangecoeff",2*npara);
    
     /*Output dimension*/
     sf_putint(out,"n1",npara);
     sf_putint(out,"d1",1);
     sf_putint(out,"o1",0);
     sf_putint(out,"n2",nmodel);
     sf_putint(out,"d2",1);
     sf_putint(out,"o2",0);
     sf_putint(out,"n3",nt);
     sf_putint(out,"d3",dt);
     sf_putint(out,"o3",tini);

     /*Read offx */
     sf_floatread(offx,noffset,ofx);
     sf_floatread(limitx,nt,limx);
     sf_floatread(drms,nt,datrms);
     
     /*Read and define the range of coefficients*/
	sf_floatread(rangecoef,2*npara,rangecoeff);
    
     /* Loop over time slices */
    for(k1=0; k1 < nt; k1++) {

		/* Change mu and SD for different events */
        if(wgauss && k1!=0) sf_floatread(rangecoef,2*npara,rangecoeff);
		 
		for(l=0; l < npara; l++) {
			if (l==0 && wgauss){ // Case of Gaussian W
				drangecoef[0] = rangecoef[1]; // sigma (SD)
				midcoef[0] = rangecoef[0]; // Mu (mean)
			} else {
				drangecoef[l] = rangecoef[2*l+1] - rangecoef[2*l];
				midcoef[l] = (rangecoef[2*l+1] + rangecoef[2*l])/2;
			}
		
			/* Step 1: Initial model */
			if (l==0 && wgauss) currentcoef[0] = drangecoef[0]*(sqrt(-2.0 * log(genrand_real1())) * cos(two_pi * genrand_real1())) + midcoef[0];
			else currentcoef[l] = drangecoef[l]*(genrand_real1()-0.5) + midcoef[l];
		
		 }
	 
		 if (eta) {
				currentcoef[2] = -currentcoef[1]/currentcoef[0] + currentcoef[0] + currentcoef[0]*currentcoef[1]/(currentcoef[1]-2*currentcoef[0]*currentcoef[0]);
				currentcoef[3] = 4*powf(currentcoef[0],6)/powf(currentcoef[1]-2*currentcoef[0]*currentcoef[0],2);
		 }
     	
    
        sf_floatread(dattime,noffset,inp);
        sf_floatread(t0sq,noffset,t0sqf);
     	
        /* Initial likelihood */
        currentL = likelihood(offx,noffset,limitx[k1],drms[k1],dattime,t0sq,esttime,currentcoef);
        reject = 0;
        
//         sf_warning("%g %g",currentL,limitx[k1]);
        
		/* Loop over randomly generated models */
		for(k2=0; k2 < nmodel; k2++) {

			/* Step 4: Save state every saveiter loops to increase sampling independence */
			for(k3=0; k3 < saveiter; k3++){

				/* Step 2: Find new model */
				/* random which parameter to update */
// 				if (eta) index = (int)round((genrand_real1())); /* either W or A -- B and C follow will follow */
// 				index = (int)round((npara-1)*(genrand_real1()));
// // 				
// 				if (index==0 && wgauss) {
// 					newcoef[0] = drangecoef[0]*(sqrt(-2.0 * log(genrand_real1())) * cos(two_pi * genrand_real1())) + midcoef[0];
// 				} else {
// 					newcoef[index] = drangecoef[index]*(genrand_real1()-0.5) + midcoef[index];	
// 				}
				for(k4=0; k4 < npara; k4++) {
					if (k4==0 && wgauss) {
						newcoef[0] = drangecoef[0]*(sqrt(-2.0 * log(genrand_real1())) * cos(two_pi * genrand_real1())) + midcoef[0];
					} else {
						newcoef[k4] = drangecoef[k4]*(genrand_real1()-0.5) + midcoef[k4];		
					}
				}
				
				if (eta) {
					newcoef[2] = -newcoef[1]/newcoef[0] + newcoef[0] + newcoef[0]*newcoef[1]/(newcoef[1]-2*newcoef[0]*newcoef[0]);
     				newcoef[3] = 4*powf(newcoef[0],6)/powf(newcoef[1]-2*newcoef[0]*newcoef[0],2);
     			}

// powf(currentcoef[4]/newcoef[4],2)*
				/* Step 3: Metropolis rule */
				if(!prior) {
				
					newL = likelihood(offx,noffset,limitx[k1],drms[k1],dattime,t0sq,esttime,newcoef);
// if( log(genrand_real1()) < (noffset/2)*(log(currentcoef[4])-log(newcoef[4])) + (currentL-newL))
					if( log(genrand_real1()) < (noffset/2)*(log(currentcoef[4]/newcoef[4])) + (currentL-newL)) { /* Do it in log form to make sure the value is manageable */
						for(l=0; l < npara; l++) currentcoef[l] = newcoef[l];
						currentL = newL;
					} else {
						reject++;
					}
				} else {
					for(l=0; l < npara; l++) currentcoef[l] = newcoef[l];
				}
			
				/* Record state */
				if (k3==saveiter-1) for(l=0; l < npara; l++) final[k2][l] = currentcoef[l];
			} // k3
			
			sf_warning(" Time step: %d of %d Model : %d of %d Reject : %d of %d;",k1+1,nt,k2+1,nmodel,reject,saveiter*nmodel);
		
		} // k2
        sf_floatwrite(final[0],npara*nmodel,out);
        sf_warning("\n");
    } // k1



  /* free memory */
	free(t0sq);
	free(dattime);
	free(esttime);
	free(rangecoef); 
	free(newcoef);
	free(currentcoef);
	free(drangecoef); 
	free(midcoef); 
	free(offx); 
	free(final); 
	free(limitx);  
	free(drms); 
	 
  exit(0);
}


/* Compute the likelihood function of the new proposed coeff */
static float likelihood(float *offsetx, int noffset, float lmx /* max offset to compute likelihood */,
				float datarms /* RMS average of data*/,
				float *dattime, /*true time (data) at the time (depth) slice*/
				float *t0sq, /*t0 sq at the time (depth) slice*/
				float *esttime, /*est time at the time (depth) slice*/				
				float *newcoef /*new coeff drawn from uniform sampling*/)
{
	float sigm;
	
	L = 0.0; /* Initial likelihood */
	
	w1 = newcoef[0];A1 = newcoef[1];B1 = newcoef[2];C1 = newcoef[3]; sigmapc=newcoef[4];

	
	/* Loops for each offset (x,y) */		
	for(i=0;i<noffset;i++){
		x = offsetx[i];
		if (fabsf(x)<lmx) {
		    
		    t0 = sqrt(t0sq[i]);

			A = A1*pow(x,4);
			B = B1*pow(x,2);
			C = C1*pow(x,4);

			testval = pow(t0,2) + B + sqrt(pow(t0,4) + C + 2*pow(t0,2)*B);

			/* Avoid dividing by zero*/
			if (fabsf(testval)> 0.001) {

			/* Compute traveltime (t^2)*/
			esttime[i] = pow(t0,2) + w1*pow(x,2) + A/(pow(t0,2) + B + sqrt(pow(t0,4) + C + 2*pow(t0,2)*B));
		
			} else {
	
			esttime[i] = 0.0; /* Generate large misfit */
	
			}
		// 		sf_warning("true %g predict %g",dattime[i],esttime[i]);
			/* Compute the likelihood */
			sigm = sigmapc/100*datarms;
			L = L + 0.5*powf(fabs(dattime[i]-esttime[i])/sigm,2);
		}
	}
	return L;
}








