#ifndef CURVE_DEFINER_DIP_OFFSET_H
#define CURVE_DEFINER_DIP_OFFSET_H

#include "curveDefinerBase.hh"

class CurveDefinerDipOffset : public CurveDefinerBase {

public:

    CurveDefinerDipOffset ();
   ~CurveDefinerDipOffset ();

	void  getEscapePoint   (const float curY, const float curX, const float curZeroTime, 
	  				        const float curDip, const float curAz, const float migVel,
							float& time, float& xRayExit, float& yRayExit);
								   
    float getEscapeTime    (const float dy, const float dx, const float t0,
						    const float migVel, float& p);
		    
};
#endif
