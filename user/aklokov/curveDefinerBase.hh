#ifndef CURVE_DEFINER_BASE_H
#define CURVE_DEFINER_BASE_H

class CurveDefinerBase {

public:

                  CurveDefinerBase ();
                 virtual ~CurveDefinerBase ();

	virtual void  getEscapePoint   (const float curY, const float curX, const float curZeroTime, 
							       const float curDip, const float curAz, const float migVel,
								   float& time, float& xRayExit, float& yRayExit) = 0;
								   
	virtual float getEscapeTime    (const float dy, const float dx, const float t0,
									const float migVel, float& p) = 0;
									
	float curOffset_;						    	
};
#endif
