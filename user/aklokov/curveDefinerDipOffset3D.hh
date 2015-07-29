#ifndef CURVE_DEFINER_DIP_OFFSET_3D_H
#define CURVE_DEFINER_DIP_OFFSET_3D_H

#include "curveDefinerBase.hh"

class CurveDefinerDipOffset3D : public CurveDefinerBase {

public:

    CurveDefinerDipOffset3D ();
   ~CurveDefinerDipOffset3D ();

	void  getEscapePoint   (const float curY, const float curX, const float curZeroTime, 
	  				        const float curDip, const float curAz, const float migVel,
							float& time, float& xRayExit, float& yRayExit);
								   
    float getEscapeTime    (const float dy, const float dx, const float t0,
						    const float migVel, float& p);
		    
};
#endif
