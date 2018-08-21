/* 3D NMO GMA MCMC inversion for spiral sorted gather with Metropolis rule (Mosegaard and Tarantola, 1995) */
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
float testval,x,y,t0,L,w1,w2,w3,A1,A2,A3,A4,A5,B1,B2,B3,C1,C2,C3,C4,C5,A,B,C;

/* Function declaration */
static float likelihood(float *offsetx, float *offsety, int noffset, float sigma, /*Noise variance */
				float *dattime, /*true time (data) at the time (depth) slice*/
				float *t0sq, /*t0 sq at the time (depth) slice*/
				float *esttime, /*est time at the time (depth) slice*/				
				float *newcoef /*new coeff drawn from uniform sampling*/);


/* Main */
int main(int argc, char* argv[])
{
  int noffset, nt, npara, nmodel, nc, seed, reject, getin, saveiter;
  float dt, tini, currentL, newL, sigma;
  float *t0sq, *dattime, *esttime, *rangecoef, *newcoef, *currentcoef, *drangecoef, *midcoef, *offx, *offy, **final; 
  int k1,k2,k3,k4,l;
  bool prior;

  sf_file inp, t0sqf, out ,rangecoeff, ofx, ofy;

  sf_init(argc, argv);
    inp = sf_input("in");
    t0sqf = sf_input("t0sq");
    rangecoeff = sf_input("rangecoef");
    ofx = sf_input("offsetx");
    ofy = sf_input("offsety");
    out = sf_output("out");

    /* input vector data t^2' */
    if (!sf_histint(inp,"n1",&noffset)) sf_error("No n1=");
    if (!sf_histint(inp,"n2",&nt)) sf_error("No n2=");
    if (!sf_histfloat(inp,"d2",&dt)) sf_error("No d2=");
    if (!sf_histfloat(inp,"o2",&tini)) sf_error("No o2=");

	
	if (!sf_getint("seed",&seed)) seed = time(NULL);
    /* random seed */
    init_genrand((unsigned long) seed);
    
    if (!sf_getfloat("sigma",&sigma)) sigma = 1.0;
    /* noise variance */
    
    if (!sf_getint("saveiter",&saveiter)) saveiter = 20;
    /* save state every iter  */
    
    if (!sf_getbool("prior",&prior)) prior = false;
    /* generate prior or posterior */

    /*Number of fitting parameters*/
    npara=16;
    
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
    offy = sf_floatalloc(noffset);
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

     /*Read offx and offy */
     sf_floatread(offx,noffset,ofx);
      sf_floatread(offy,noffset,ofy);
     
     /*Read and define the range of coefficients*/
     sf_floatread(rangecoef,2*npara,rangecoeff);
     for(l=0; l < npara; l++) {
     	drangecoef[l] = rangecoef[2*l+1] - rangecoef[2*l];
     	midcoef[l] = (rangecoef[2*l+1] + rangecoef[2*l])/2;
     	
     	/* Step 1: Initial model */
     	currentcoef[l] = drangecoef[l]*(genrand_real1()-0.5) + midcoef[l]; 
     }

    /* Loop over time slices */
    for(k1=0; k1 < nt; k1++) {
        sf_floatread(dattime,noffset,inp);
        sf_floatread(t0sq,noffset,t0sqf);
        
        /* Initial likelihood */
        currentL = likelihood(offx,offy,noffset,sigma,dattime,t0sq,esttime,currentcoef);
        reject = 0;
        
		/* Loop over randomly generated models */
		for(k2=0; k2 < nmodel; k2++) {

			/* Step 4: Save state every saveiter loops to increase sampling independence */
			for(k4=0; k4 < saveiter; k4++){

				/* Step 2: Find new model */
				/* Loop over all parameters */
				for(k3=0; k3 < npara; k3++) {
						newcoef[k3] = drangecoef[k3]*(genrand_real1()-0.5) + midcoef[k3];
// 						getin = 0; /* first time getting in */
// 						while (newcoef[k3] > rangecoef[2*k3+1] || newcoef[k3] < rangecoef[2*k3] || getin == 0) {
// 							newcoef[k3] = drangecoef[k3]*(genrand_real1()-0.5) + currentcoef[k3];
// 							getin = 1;
// 						}
						
				}
			
				/* Step 3: Metropolis rule */
				if(!prior) {
				
					newL = likelihood(offx,offy,noffset,sigma,dattime,t0sq,esttime,newcoef);
					
					if( genrand_real1() < newL/currentL) {
						for(l=0; l < npara; l++) currentcoef[l] = newcoef[l];
						currentL = newL;
					} else {
						reject++;
					}
				} 
			
				/* Record state */
				if (k4==saveiter-1) for(l=0; l < npara; l++) final[k2][l] = currentcoef[l];
			} // k4
			
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
	free(offy);  
	free(final); 
 
  exit(0);
}


/* Compute the likelihood function of the new proposed coeff */
static float likelihood(float *offsetx, float *offsety, int noffset, float sigma, /*Noise variance */
				float *dattime, /*true time (data) at the time (depth) slice*/
				float *t0sq, /*t0 sq at the time (depth) slice*/
				float *esttime, /*est time at the time (depth) slice*/				
				float *newcoef /*new coeff drawn from uniform sampling*/)
{
	L = 0.0; /* Initial likelihood */
	
	w1 = newcoef[0];w2 = newcoef[1];w3 = newcoef[2];
	A1 = newcoef[3];A2 = newcoef[4];A3 = newcoef[5];A4 = newcoef[6];A5 = newcoef[7];
	B1 = newcoef[8];B2 = newcoef[9];B3 = newcoef[10];
	C1 = newcoef[11];C2 = newcoef[12];C3 = newcoef[13];C4 = newcoef[14];C5 = newcoef[15];

	
	/* Loops for each offset (x,y) */		
	for(i=0;i<noffset;i++){
		x = offsetx[i];
		y = offsety[i];
	
		t0 = sqrt(t0sq[i]);
	
		A = A1*pow(x,4) + A2*pow(x,3)*y + A3*pow(x,2)*pow(y,2) + A4*x*pow(y,3) + A5*pow(y,4);
		B = B1*pow(x,2) + B2*x*y + B3*pow(y,2);
		C = C1*pow(x,4) + C2*pow(x,3)*y + C3*pow(x,2)*pow(y,2) + C4*x*pow(y,3) + C5*pow(y,4);
	
		testval = pow(t0,2) + B + sqrt(pow(t0,4) + C + 2*pow(t0,2)*B);
	
		/* Avoid dividing by zero*/
		if (fabsf(testval)> 0.001) {
	
		/* Compute traveltime (t^2)*/
		esttime[i] = pow(t0,2) + w1*pow(x,2) + w2*x*y + w3*pow(y,2) + A/(pow(t0,2) + B + sqrt(pow(t0,4) + C + 2*pow(t0,2)*B));
			
		} else {
		
		esttime[i] = 0.0; /* Generate large misfit */
		
		}
		
		/* Compute the likelihood */
		L = L + 0.5*powf(fabs(dattime[i]-esttime[i]),2);
		
	}
	
	L = exp(-L/(sigma*sigma));
	return L;
}










