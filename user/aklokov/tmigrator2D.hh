#ifndef TIME_MIGRATOR_2D_H
#define TIME_MIGRATOR_2D_H

#include "tmigratorBase.hh"

class TimeMigrator2D : public TimeMigratorBase {

public:

	    TimeMigrator2D ();
	   ~TimeMigrator2D ();

	void  processGather (Point2D& curGatherCoords, float curOffset,  const float* const velTrace, const bool isAzDip,
								  float* curoffsetGather, float* curoffsetImage, float* curoffsetImageSq);

	void getStackTaper (const float edgeTaper, const bool dummy);

private:

	float getSampleFromData (const float geoY, const float geoX, const float ti, const float p);
	int   getSampleByBeam   (const float yCIG, const float xCIG, const float curZeroTime, 
				 const float curDip, const float curAz, const float migVel, const bool isAzDip,
				 float& sample);
	float getSampleByRay (const float yCIG, const float xCIG, const float curZeroTime, 
					      const float curDip, const float curAz, const float migVel, const bool isAzDip,
					      const float yEmAngle, const float xEmAngle);
};
#endif
