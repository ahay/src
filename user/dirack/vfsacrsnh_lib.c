/*
	 vfsacrsnh_lib.c (c)
	 
	 Purpose: 'Mvfsacrsnh.c' library.
	 	 
	 Version 2.0.1
	 
	 Site: https://dirack.github.io
	 
	 Programmer: Rodolfo A. C. Neves (Dirack) 13/08/2021

	 Email:  rodolfo_profissional@hotmail.com

	 License: GPL-3.0 <https://www.gnu.org/licenses/gpl-3.0.txt>.

*/

#define Beta_MAX 1 // Beta maximum value
#define Beta_MIN -1 // Beta minimun value
#define BETA_APERTURE Beta_MAX-Beta_MIN 
#define Rnip_MAX 4 // RNIP maximum value
#define Rnip_MIN 0 // RNIP minimum value
#define RNIP_APERTURE Rnip_MAX-Rnip_MIN
#define Rn_MAX 5 // RN maximum value
#define Rn_MIN 0 // RN minimum value
#define RN_APERTURE Rn_MAX-Rn_MIN
#define hMAX 50 // Max of samples to stack in half-offset
#define mMAX 50 // Max of samples to stack in CMP
#define ITMAX 5000 // Maximum number of iterations in VFSA
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <rsf.h>
/*^*/

#define signal(s) ((s<0)?(-1.):(1.))
/*< Signal function >*/
/*^*/

float getRandomNumberBetween0and1(){
/*< Function to get a random number between 0 and 1 >*/

	return (float)(rand()%1000)/1000;
}

float getVfsaIterationTemperature(	int iteration /* Number of the current iteration */,
					float dampingFactor /* C0, the VFSA damping factor */,
					float initialTemperature /* T0, the VFSA Initial temperature */)
/*< Temperature function for VFSA algorithm >*/
{

	return initialTemperature*expf(-dampingFactor*pow(iteration,0.25));

}

void disturbParameters(  float temperature /* Temperature of the current iteration */,
			 float* disturbedParameter /* disturbed parameters vector */,
			 float* parameter /* parameters vector */)
/*< Disturb parameters from the previous VFSA iteration

Note: It receives a parameter vector and distubs it accordingly to 
VFSA algorithm disturb parameters step.
 >*/
{
	float u;
	float disturbance;

	u=getRandomNumberBetween0and1();
			
	disturbance = signal(u - 0.5) * temperature * (pow( (1+temperature),fabs(2*u-1) )-1);

	/* Disturb RN */

	disturbedParameter[0] = parameter[0] + disturbance * (RN_APERTURE);
				
	if (disturbedParameter[0] >= Rn_MAX || disturbedParameter[0] <= Rn_MIN) {

		disturbedParameter[0] = (RN_APERTURE) * getRandomNumberBetween0and1() + Rn_MIN;
		
	}

	/* Disturb RNIP */

	disturbedParameter[1] = parameter[1] + disturbance * (RNIP_APERTURE);
				
	if (disturbedParameter[1] >= Rnip_MAX || disturbedParameter[1] <= Rnip_MIN) {

		disturbedParameter[1] = (RNIP_APERTURE) * getRandomNumberBetween0and1() + Rnip_MIN;
		
	}

	/* Disturb BETA */

	disturbedParameter[2] = parameter[2] + (disturbance/10.) * (BETA_APERTURE);

	if (disturbedParameter[2] >= Beta_MAX || disturbedParameter[2] <= Beta_MIN) {

		disturbedParameter[2] = (BETA_APERTURE) * getRandomNumberBetween0and1() + Beta_MIN;

	}		

}

void nonHyperbolicCRSapp(	float t[2*mMAX+1][hMAX] /* non-hyperbolic CRS traveltime surface */,
				float m0 /* Central CMP of the approximation */,
				float dm, /* CMP sampling */
				float om /* CMP axis origin */,
				float dh /* half-offset sampling */,
				float oh /* half-offset axis origin */,
				float t0 /* Normal ray traveltime */,
				float v0 /* Near surface velocity */,
				float RN /* RN, CRS parameter */,
				float RNIP /* RNIP, CRS parameter */,
				float BETA /* BETA, CRS parameter */)
/*< Returns the Non hyperbolic CRS approximation surface (FOMEL; KAZINNIK, 2013) >*/
{
	float m0_index=(int)(m0/dm);
	float a1, a2, b2, c1, Fd, Fd1, Fd2;
	int im, ih;
	float m, h, mmh, mph;
	float sinB=sin(BETA),cosB=cos(BETA);
	
	om = om+(m0_index-mMAX)*dm;

	a1=(2*sinB)/(v0);		
	a2=(2*cosB*cosB*t0)/(v0*RN);
	b2=(2*cosB*cosB*t0)/(v0*RNIP);
	c1=2*b2+a1*a1-a2;

	for (im=0; im < 2*mMAX+1; im++){
			
		m=(im*dm+om)-m0;

		for(ih=0;ih<hMAX;ih++){
			
			h=ih*dh+oh;
			mmh=m-h;
			mph=m+h;

			Fd=(t0+a1*m)*(t0+a1*m)+a2*m*m;				
			Fd2=(t0+a1*(mmh))*(t0+a1*(mmh))+a2*(mmh)*(mmh);
			Fd1=(t0+a1*(mph))*(t0+a1*(mph))+a2*(mph)*(mph);					
			t[im][ih]=sqrt((Fd+c1*h*h+sqrt(Fd2*Fd1))*0.5); 
		}
	}

}

float semblance(float m0 /* Central CMP of the approximation */,
		float dm /* CMP sampling */,
		float om /* CMP axis origin */,
		float oh /* half-offset axis origin */,
		float dh /* half-offset axis sampling */,
		float dt /* time sampling */,
		int nt /* number of time samples */,
		float t0 /* Normal ray traveltime */,
		float v0 /* Near surface velocity */,
		float RN /* RN, CRS parameter */,
		float RNIP /* RNIP, CRS parameter */,
		float BETA /* BETA, CRS parameter */,
		float*** t /* reflection data cube A(m,h,t) */)
/*< Calculate semblance between the Non Hyperbolic CRS approximation surface and reflection data >*/
{

	int im, ih, numSamples=0;
	float amplitude=0.;
	float amplitudeSampleSum=0.;
	float amplitudeSquaredSampleSum=0.;
	float semblance=0;
	int tetai;
	float teta[2*mMAX+1][hMAX];
	int m0_index_init, m0_index_end;

	m0_index_init = (int)(m0/dm)-mMAX;
	m0_index_end = (int)(m0/dm)+mMAX;

	nonHyperbolicCRSapp(teta,m0,dm,om,dh,oh,t0,v0,RN,RNIP,BETA);

	for (im=m0_index_init; im < m0_index_end; im++){
			
		for(ih=0;ih<hMAX;ih++){

			tetai=teta[im-m0_index_init][ih]/dt;

			if(tetai>=0 && tetai < nt){
				amplitude = t[im][ih][tetai];
				
				amplitudeSampleSum += amplitude;
					
				amplitudeSquaredSampleSum += (amplitude*amplitude);
			}
				
			numSamples++;
		}
		
	}		

	if(amplitudeSampleSum==0)		
		return semblance=0;
	else
		return semblance=(amplitudeSampleSum*amplitudeSampleSum)/(numSamples*amplitudeSquaredSampleSum);

}

