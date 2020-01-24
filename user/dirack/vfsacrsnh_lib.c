/*
	 vfsacrsnh_lib.c (c)
	 
	 Purpose: 'Mvfsacrsnh.c' library.
	 	 
	 Version 1.0
	 
	 Site: http://www.dirackslounge.online
	 
	 Programer: Rodolfo A. C. Neves (Dirack) 19/09/2019

	 Email:  rodolfo_profissional@hotmail.com

	 License: GPL-3.0 <https://www.gnu.org/licenses/gpl-3.0.txt>.

*/

#define Beta_MAX 1
#define Beta_MIN -1
#define Rnip_MAX 5
#define Rnip_MIN 0
#define Rn_MAX 5
#define Rn_MIN 0
#define hMAX 50
#define mMAX 50
#define ITMAX 5000
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <rsf.h>
/*^*/

float signal(float s) { 
/*< Signal function >*/

	if(s >= 0){
		
		return s = 1;		
	}
	
	return s=-1;
		
}

float getRandomNumberBetween0and1(){
/*< Function to get a random number between 0 and 1 >*/

	float randomNumber;
	int u;

	u = rand()%1000;
			
	randomNumber = (float)u/1000;

	return randomNumber;

}

float getVfsaIterationTemperature(int iteration,float dampingFactor,float inicialTemperature){
/*< Temperature function for VFSA algorithm >*/

	float temperature;

	temperature=inicialTemperature*expf(-dampingFactor*pow(iteration,0.25));

	return temperature;
}

void disturbParameters(float temperature, float* disturbedParameter, float* parameter){
/*< Perturbar os parâmetros da iteração anterior >*/

	float u;
	float disturbance;
	float aperture;

	u=getRandomNumberBetween0and1();
			
	disturbance = signal(u - 0.5) * temperature * (pow( (1+temperature),fabs(2*u-1) )-1);

	/* RN */
	aperture = Rn_MAX - Rn_MIN;

	disturbedParameter[0] = parameter[0] + disturbance * (aperture);
				
	if (disturbedParameter[0] >= Rn_MAX || disturbedParameter[0] <= Rn_MIN) {

		disturbedParameter[0] = (aperture) * getRandomNumberBetween0and1() + Rn_MIN;
		
	}

	/* RNIP */
	aperture = Rnip_MAX - Rnip_MIN;

	disturbedParameter[1] = parameter[1] + disturbance * (aperture);
				
	if (disturbedParameter[1] >= Rnip_MAX || disturbedParameter[1] <= Rnip_MIN) {

		disturbedParameter[1] = (aperture) * getRandomNumberBetween0and1() + Rnip_MIN;
		
	}

	/* BETA */
	aperture = Beta_MAX - Beta_MIN;

	disturbedParameter[2] = parameter[2] + disturbance * (aperture);

	if (disturbedParameter[2] >= Beta_MAX || disturbedParameter[2] <= Beta_MIN) {

		disturbedParameter[2] = (aperture) * getRandomNumberBetween0and1() + Beta_MIN;

	}		

}

float nonHyperbolicCRSapp(float m, float h, float t0, float v0, float RN, float RNIP, float BETA){
/*< Non hyperbolic CRS approximation (FOMEL; KAZINNIK, 2013) >*/
	float t;
	float a1, a2, b2, c1, Fd, Fd1, Fd2;
			
	a1=(2*sin(BETA))/(v0);		
	a2=(2*cos(BETA)*cos(BETA)*t0)/(v0*RN);
	b2=(2*cos(BETA)*cos(BETA)*t0)/(v0*RNIP);
	c1=2*b2+a1*a1-a2;
												
	Fd=(t0+a1*m)*(t0+a1*m)+a2*m*m;				
	Fd2=(t0+a1*(m-h))*(t0+a1*(m-h))+a2*(m-h)*(m-h);
	Fd1=(t0+a1*(m+h))*(t0+a1*(m+h))+a2*(m+h)*(m+h);					
	return t=sqrt((Fd+c1*h*h+sqrt(Fd2*Fd1))*0.5); 

}

float semblance(float m0, float dm, float om, float oh, float dh, float dt, int nt,float t0, float v0,float RN, float RNIP, float BETA, float*** t){
/*< Semblance: Non Hyperbolic CRS approximation with data >*/

	int im, ih, numSamples=0;
	float m, h;
	float amplitude=0.;
	float amplitudeSampleSum=0.;
	float amplitudeSquaredSampleSum=0.;
	float semblance=0;
	int tetai;
	float teta;
	int m0_index;

	m0_index = (int)(m0/dm);

	for (im=m0_index-mMAX; im < m0_index+mMAX; im++){
			
		for(ih=0;ih<hMAX;ih++){

			m=im*dm+om;
	
			m=m-m0;
			
			h=ih*dh+oh;

			teta = nonHyperbolicCRSapp(m,h,t0,v0,RN,RNIP,BETA);

			tetai=teta/dt;

			if(tetai>=0 && tetai < nt){
				amplitude = t[im][ih][tetai];
			}else{
				amplitude=0.;
			}
	
			amplitudeSampleSum=amplitudeSampleSum+amplitude;
					
			amplitudeSquaredSampleSum=amplitudeSquaredSampleSum+(amplitude*amplitude);
				
			numSamples++;
		}
		
	}		

	if(amplitudeSquaredSampleSum==0 || amplitudeSampleSum==0)		
	return semblance=0;
	else
	return semblance=(amplitudeSampleSum*amplitudeSampleSum)/(numSamples*amplitudeSquaredSampleSum);

}

