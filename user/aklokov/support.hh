#ifndef SUPPORT_H
#define SUPPORT_H

const double CONVRATIO = 3.14159265 / 180.f;

class Point2D {
public:
    Point2D ();
    Point2D (float x, float y);
   ~Point2D ();

    Point2D& operator= (const Point2D& point);
   	bool operator < (const Point2D& p) const
    {
        return (x_ < p.x_);
    }

    void  setX (float x) {x_ = x;}
    float getX ()        {return x_;}
    void  setY (float y) {y_ = y;}
    float getY ()        {return y_;}
 
//private:
    float x_;
    float y_;
};

class EscapePoint {
public:
    EscapePoint ();
    EscapePoint (float x, float z, float t, float p, float offset, float startDir, bool isSurf);
   ~EscapePoint ();

    EscapePoint& operator= (const EscapePoint& point);

	float x;	
	float z;
	float t;
	float p;
	float offset;    // offset
	float startDir;  // scattering direction of the ray
	bool  isSurf;    // if ray reached the daylight surface
};


struct RunParamsTmigda {


    bool              isDipAz;   // 0 - "crossline/inline" mode
                             // 1 - "azimuth/dip" mode

    bool             isDag;
    bool             isCig;
    bool             isSemb;
	bool             isMCig;  // super (dip-angle and scattering-angle) CIG
	bool             isTT;

    
    bool              is3D;   // 0 - 2D mode (by default)
    						  // 1 - 3D mode

	bool              isAA;   // 0 - no anti-aliasing filter
							  // 1 - anti-aliasing filter after Lumley-Claerbout-Bevc (by default)

	bool              isCMP;  // 0 - if traces have coordinates of a receiver 
							  // 1 - if traces have coordinates of CMP
	bool              isVelMS; 

	int               hMigNum;
	int               sembWindow;

	float             edgeTaper;
	
};

struct VolumeParams {
  
    int   zNum;
    float zStart;
    float zStep; 

    int   hNum;
    float hStart;
    float hStep; 

    int   xNum;
    float xStart;
    float xStep; 

    int   yNum;
    float yStart;
    float yStep; 
};

struct GatherParams {

    int   zNum;
    float zStart;
    float zStep; 

    int   dipNum;
    float dipStart;
    float dipStep; 

    int   scatNum;
    float scatStart;
    float scatStep; 

    int   sdipNum;
    float sdipStart;
    float sdipStep; 
};

#endif
