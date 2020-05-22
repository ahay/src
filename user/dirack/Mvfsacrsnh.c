/* Version 1.0 - Zero offset CRS parameter inversion (RN, RNIP, BETA) with Very Fast Simulated Aneeling (VFSA) Global Optimization

This program the Non-Hyperbolic CRS approximation to fit data cube and get the parameters (Fomel, 2013).

Programer: Rodolfo A. C. Neves (Dirack) 19/09/2019

Email:  rodolfo_profissional@hotmail.com

License: GPL-3.0 <https://www.gnu.org/licenses/gpl-3.0.txt>.

 */

#include "vfsacrsnh_lib.h"

int main(int argc, char* argv[])
{

	float m0; // central CMP
	float om; // CMP axis origin
	float dm; // CMP sampling
	int nm; // Number of CMP's
	float oh; // Offset axis origin
	float dh; // Offset sampling
	int nh; // Number of Offsets
	int nt; // Number of time samples
	float ot; // Time axis origin
	float dt; // Time sampling
	bool verb; // Key to turn On/Off active mode
	float v0; // Near surface velocity
	float t0; // Normal ray time travel
	float cnew[3]; // Temporary parameters vector - actual iteration
	float c[3]; // Temporary parameters vector - last iteration
	float *otm; // Optimazed parameters
	float otrn, otrnip, otbeta, otsemb; // Optimazed parameters - actual iteration
	float deltaE, PM; // Metrópolis criteria
	float Em0=0; // Major semblance
	float u; // Random number
	float ***t; // Data cube A(m,h,t)
	int q, i; // loop counter
	float semb; // Semblance - actual iteration
	float RN, RNIP, BETA; // CRS parameters
	float semb0; // Inicial semblance value
	float c0; // VFSA damping factor
	float temp0; // inicial VFSA temperature
	float temp; // VFSA temperature
	int repeat; // perform VFSA optimization more than once

	/* RSF files I/O */  
	sf_file in, out;

	/* RSF files axis */
	sf_axis ax,ay,az;

	sf_init(argc,argv); 

	in = sf_input("in");
	out = sf_output("out");

	if (!sf_getfloat("m0",&m0)) m0=0;
	/* central CMP of the approximation (Km) */

	if (!sf_getfloat("v0",&v0)) v0=1.5;
	/* Near surface velocity (Km/s) */

	if (!sf_getfloat("t0",&t0)) t0=1.5;
	/* Normal ray traveltime (s) */

	if (!sf_getfloat("c0",&c0)) c0=0.5;
	/* damping factor of VFSA */

	if (!sf_getfloat("temp0",&temp0)) temp0=10;
	/* inicial VFSA temperature */

	if(!sf_getint("repeat",&repeat)) repeat=1;
	/* How many times to perform VFSA global optimization */

	if (!sf_histint(in,"n1",&nt)) sf_error("No n1= in input");
	if (!sf_histfloat(in,"d1",&dt)) sf_error("No d1= in input");
	if (!sf_histfloat(in,"o1",&ot)) sf_error("No o1= in input");
	if (!sf_histint(in,"n2",&nh)) sf_error("No n2= in input");
	if (!sf_histfloat(in,"d2",&dh)) sf_error("No d2= in input");
	if (!sf_histfloat(in,"o2",&oh)) sf_error("No o2= in input");
	if (!sf_histint(in,"n3",&nm)) sf_error("No n3= in input");
	if (!sf_histfloat(in,"d3",&dm)) sf_error("No d3= in input");
	if (!sf_histfloat(in,"o3",&om)) sf_error("No o3= in input");

	if(! sf_getbool("verb",&verb)) verb=0;
	/* 1: active mode; 0: quiet mode */

	if (verb) {

		sf_warning("Active mode on!!!");
		sf_warning("Command line parameters: "); 
		sf_warning("m0=%f v0=%f t0=%f c0=%f temp0=%f repeat=%i",m0,v0,t0,c0,temp0,repeat);
		sf_warning("Input file parameters: ");
		sf_warning("n1=%i d1=%f o1=%f",nt,dt,ot);
		sf_warning("n2=%i d2=%f o2=%f",nh,dh,oh);
		sf_warning("n3=%i d3=%f o3=%f",nm,dm,om);
	}
	
	c[0] = 0;
	c[1] = 0;
	c[2] = 0;
	cnew[0] = 0;
	cnew[1] = 0;
	cnew[2] = 0;

	srand(time(NULL));

	/* Read seismic data cube */
	t=sf_floatalloc3(nt,nh,nm);
	sf_floatread(t[0][0],nm*nh*nt,in);

	sf_fileclose(in);

	semb0=0;

	for(i=0;i<repeat;i++){

		for (q=0; q <ITMAX; q++){
				
			/* calculate VFSA temperature for this iteration */
			temp=getVfsaIterationTemperature(q,c0,temp0);
							
			/* parameter disturbance */
			disturbParameters(temp,cnew,c);
																	
			RN = cnew[0];
			RNIP = cnew[1];
			BETA = cnew[2];
		
			semb=0;
			
			/* Calculate semblance: Non-hyperbolic CRS approximation with data */		
			semb=semblance(m0,dm,om,oh,dh,dt,nt,t0,v0,RN,RNIP,BETA,t);
			
			
			/* VFSA parameters convergence condition */		
			if(fabs(semb) > fabs(semb0) ){
				otsemb = semb;
				otrn = RN;
				otrnip = RNIP;
				otbeta = BETA;
				semb0 = semb;			
			}

			/* VFSA parameters update condition */
			deltaE = -semb - Em0;
			
			/* Metrópolis criteria */
			PM = expf(-deltaE/temp);
			
			if (deltaE<=0){
				c[0] = cnew[0];
				c[1] = cnew[1];
				c[2] = cnew[2];
				Em0 = -semb;
			} else {
				u=getRandomNumberBetween0and1();
				if (PM > u){
					c[0] = cnew[0];
					c[1] = cnew[1];
					c[2] = cnew[2];
					Em0 = -semb;
				}	
			}	
			
		} /* loop over iterations */

	} /* repeat VFSA global optimization */

	free(t);

	/* Save optimized parameters in 'param' file */
	otm=sf_floatalloc(8);
	otm[0] = otrn;
	otm[1] = otrnip;
	otm[2] = otbeta;
	otm[3] = otsemb;
	otm[4] = c0;
	otm[5] = temp0;
	otm[6] = t0;
	otm[7] = m0;

	/* Show optimized parameters on screen before save them */
	sf_warning("Optimized parameters:\n RN=%f, RNIP=%f, BETA=%f, SEMB=%f",otrn,otrnip,otbeta,otsemb);

	/* axis = sf_maxa(n,o,d)*/
	ax = sf_maxa(8, 0, 1);
	ay = sf_maxa(1, 0, 1);
	az = sf_maxa(1, 0, 1);

	/* sf_oaxa(file, axis, axis index) */
	sf_oaxa(out,ax,1);
	sf_oaxa(out,ay,2);
	sf_oaxa(out,az,3);
	sf_floatwrite(otm,8,out);

	sf_close();
	exit(0);
}
