#include <string.h>
#include "support.hh"

// -- class Point2D ---
Point2D::Point2D () : x_ (0),
		      y_ (0) {
}

Point2D::Point2D (float x, float y) : x_ (x),
			 	      y_ (y) {
}

Point2D::~Point2D () {

}

Point2D& Point2D::operator= (const Point2D& point) {
    if (this == &point)
        return *this;
    x_ = point.x_;
    y_ = point.y_;

    return *this;
}

// -- class EscapePoint ---

EscapePoint::EscapePoint () : x (0.f),
						      z (0.f),
						      t (0.f), 
						      p (0.f),
						      offset   (0.f), 
							  startDir (0.f),
						      isSurf (false) {
}

EscapePoint::EscapePoint (float x1, float z1, float t1, float p1, 
						  float offset1, float startDir1, bool isSurf1) : x (x1),
							 									          z (z1),
																 	      t (t1), 	
																	      p (p1),
																	      offset   (offset1), 
																	      startDir (startDir1),
																	      isSurf   (isSurf1) {
}

EscapePoint::~EscapePoint () {
}

EscapePoint& EscapePoint::operator= (const EscapePoint& point) {
    if (this == &point)
        return *this;
    x = point.x;
    z = point.z;
    t = point.t;
	p = point.p;    
	offset   = point.offset;
	startDir = point.startDir;
	isSurf   = point.isSurf;

    return *this;
}
