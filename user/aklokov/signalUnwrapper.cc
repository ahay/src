#include <rsf.hh>
#include "signalUnwrapper.hh"

SignalUnwrapper::SignalUnwrapper () {
}

SignalUnwrapper::~SignalUnwrapper () {
}

void SignalUnwrapper::unwrap (float* ptrZO) {

	float vel = 2000;

	int pNum = 201;
	float pStart = 0.f;
	float pStep = 10.0;


	const int   tNum = dp_->zNum;
	const float tStart = dp_->zStart;
	const float tStep = dp_->zStep;

	const int   hNum = dp_->hNum;
	const float hStart = dp_->hStart;
	const float hStep = dp_->hStep;

	const int   xNum = dp_->xNum;
	const float xStart = dp_->xStart;
	const float xStep = dp_->xStep;

	for (int ix = 0; ix < xNum; ++ix) {		
		for (int ih = 0; ih < hNum; ++ih) {
			const float curOffset = hStart + ih * hStep;
			const float halfOffset = curOffset / 2.f;
			for (int it = 0; it < tNum; ++it) {
				const float curTime = tStart + it * tStep;

				const int dataShift = ix * hNum * tNum + ih * tNum + it;

				const float sample = *(ptrToData_ + dataShift);		
				// unwrapping
				for (int ip = 0; ip < pNum; ++ip) {
					const float curPos = pStart + ip * pStep;
					const float a = pow (curTime, 2) - pow (curOffset / vel, 2);
					const float b = pow (curPos - halfOffset, 2);
					const float c = pow (halfOffset, 2);

					const float t0 = c ? sqrt ( a * (1 - b / c) ) : curTime;
					const int tInd = t0 / tStep;
		
					if (tInd < 0 || tInd >= tNum) continue; 
	
					const int indZO   = ip * tNum + tInd;
					*(ptrZO + indZO) += sample;
		
//					sf_warning ("%f ", sample);
				}
			}
		}
	}	

	return;
}
