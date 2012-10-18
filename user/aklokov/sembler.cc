#include "sembler.hh"
#include <string.h>

Sembler::Sembler () {
}

Sembler::~Sembler () {
}

void Sembler::getSemblanceForTrace (int tracesNum, float* data, float* dataSq, int zNum, int sembWindow, float* semb, int k) {
   
	const int halfWindow = sembWindow / 2;

	const int zNumFull = zNum + 2 * halfWindow;
	const int numInsideWindow = 1 + 2 * halfWindow; 	
	
	float* traceSumOutput = new float [zNumFull];
	float* traceSumInput  = new float [zNumFull];
    memset (traceSumOutput, 0, zNumFull * sizeof (float));   
    memset (traceSumInput,  0, zNumFull * sizeof (float));	

	for (int iz = 0; iz < zNum; ++iz) {
		const int ind = iz + halfWindow; 
		traceSumOutput[ind] += data [iz] * data [iz];
		traceSumInput [ind] += dataSq [iz];
	}

    for (int iz = 0; iz < zNum; ++iz) {
        double sumOutput (0.f);
        double sumInput  (0.f);
		const int temp = iz + numInsideWindow;
		for (int j = iz; j < temp; ++j){
		    sumOutput += traceSumOutput[j];
		    sumInput  += traceSumInput[j];
		}
		semb[iz] = sumInput ? sumOutput / (tracesNum * k * sumInput) : 0.f;
    }

    delete traceSumOutput;
    delete traceSumInput;

    return;
}
